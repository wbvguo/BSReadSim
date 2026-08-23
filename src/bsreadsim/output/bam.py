"""Details-alignment projection and SAM 1.6 record encoding.

The Python layer derives exact linear alignments from the Full Details reference
projection. This module owns alignment semantics and SAM representation only;
the HTSlib process lifecycle and BAM artifact belong to ``output.session``.
"""

from __future__ import annotations
from contextlib import suppress

from dataclasses import dataclass
import hashlib
import os
from pathlib import Path
import subprocess


from ..process.batch import (
    BaseState,
    MethylationContext,
    MethylationSource,
    NO_VARIANT_INDEX,
    ProcessedFragment,
    ProcessedMate,
    ReadNameError,
    VariantKind,
    _CompactBaseStates,
    format_fragment_identifier,
)

from ..htsim.protocol import Header
from .errors import OutputError

from .._cext import format_sam_batch as _cext_format_sam_batch
from .._cext import format_sam_columns as _cext_format_sam_columns


BAM_CONTRACT = "bsreadsim-bam-v3"
BAM_MAPQ = 60
ANNOTATION_STATE_ALPHABET = (
    "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"
)
ANNOTATION_STATE_SCHEMA = "state64-v1"
ANNOTATION_READ_SUMMARY_SCHEMA = "u16x12-v1"
ANNOTATION_FRAGMENT_SUMMARY_SCHEMA = "u16x12-v1"
ANNOTATION_FRAGMENT_REALIZATION_SCHEMA = "packed-b64url-v1"
_MAX_REFERENCE_LENGTH = (1 << 31) - 1
_MAX_CIGAR_OPERATION_LENGTH = (1 << 28) - 1
_MAX_TEMPLATE_LENGTH = (1 << 31) - 1
_REFERENCE_FIRST = frozenset(
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$%&+./:;?@^_|~-"
)
_REFERENCE_REST = _REFERENCE_FIRST | frozenset("*=")
_SEQUENCE_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class BamError(ValueError):
    """Details projection cannot be represented as one interoperable BAM record."""


def validate_bam_header(header: Header) -> None:
    """Fail unless a protocol header can form a SAM 1.6 reference dictionary."""

    if not isinstance(header, Header):
        raise BamError("BAM requires a protocol Header")
    if not header.has_details:
        raise BamError("BAM requires Full Details protocol columns")
    if not header.contigs:
        raise BamError("BAM requires at least one reference sequence")
    observed = set()
    for contig in header.contigs:
        _validate_reference_name(contig.name)
        if contig.name in observed:
            raise BamError("BAM reference names must be unique")
        observed.add(contig.name)
        if (
            isinstance(contig.length, bool)
            or not isinstance(contig.length, int)
            or not 1 <= contig.length <= _MAX_REFERENCE_LENGTH
        ):
            raise BamError(
                "BAM reference lengths must be in [1, 2^31-1]"
            )


def build_sam_header(
    header: Header,
    *,
    sample_name: str,
    program_version: str,
    fragment_summary: bool = False,
    fragment_realization: bool = False,
) -> bytes:
    """Build the complete SAM header streamed to the HTSlib BAM writer."""

    validate_bam_header(header)
    _validate_header_value("sample name", sample_name)
    _validate_header_value("run identifier", header.run_id)
    _validate_header_value("program version", program_version)
    if not isinstance(fragment_summary, bool):
        raise BamError("BAM fragment_summary must be a boolean")
    if not isinstance(fragment_realization, bool):
        raise BamError("BAM fragment_realization must be a boolean")
    if fragment_realization and not fragment_summary:
        raise BamError("BAM fragment_realization requires fragment_summary")
    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    lines.extend(
        "@SQ\tSN:{}\tLN:{}".format(contig.name, contig.length)
        for contig in header.contigs
    )
    lines.append("@RG\tID:{}\tSM:{}".format(header.run_id, sample_name))
    lines.append(
        "@PG\tID:bsreadsim\tPN:bsreadsim\tVN:{}".format(program_version)
    )
    lines.append(
        "@CO\tDetails alignments; MAPQ 60 denotes simulated origin, "
        "not calibrated mapping confidence"
    )
    lines.extend(
        (
            "@CO\tBSREADSIM_BAM_CONTRACT={}".format(BAM_CONTRACT),
            "@CO\tBSREADSIM_CONFIG_SHA256={}".format(
                header.normalized_config_sha256.hex()
            ),
            "@CO\tAS_SCHEME=details-max-v1;AS_MAX=query_length",
            "@CO\tBSREADSIM_ZT={};ALPHABET={}".format(
                ANNOTATION_STATE_SCHEMA, ANNOTATION_STATE_ALPHABET
            ),
            "@CO\tBSREADSIM_ZR={};REQUIRED=1".format(
                ANNOTATION_READ_SUMMARY_SCHEMA
            ),
            "@CO\tBSREADSIM_ZF={};ENABLED={}".format(
                ANNOTATION_FRAGMENT_SUMMARY_SCHEMA,
                1 if fragment_summary else 0,
            ),
            "@CO\tBSREADSIM_ZX={};ENABLED={};BIT_ORDER=LSB0".format(
                ANNOTATION_FRAGMENT_REALIZATION_SCHEMA,
                1 if fragment_realization else 0,
            ),
        )
    )
    return ("\n".join(lines) + "\n").encode("ascii")


def format_sam_fragment(
    fragment: ProcessedFragment,
    *,
    paired_end: bool,
    read_group_id: str,
    contig_length: int,
    include_fragment_summary: bool = False,
    include_fragment_realization: bool = False,
) -> tuple[bytes, ...]:
    """Format one SE/PE details alignment, keeping paired records adjacent."""

    if not isinstance(fragment, ProcessedFragment):
        raise BamError("BAM fragment must be a ProcessedFragment")
    _validate_header_value("read group identifier", read_group_id)
    _validate_reference_name(fragment.contig_name)
    if not isinstance(include_fragment_summary, bool):
        raise BamError("BAM fragment-summary policy must be boolean")
    if not isinstance(include_fragment_realization, bool):
        raise BamError("BAM fragment-realization policy must be boolean")
    if include_fragment_realization and not include_fragment_summary:
        raise BamError("BAM fragment realization requires fragment summary")
    if (
        isinstance(contig_length, bool)
        or not isinstance(contig_length, int)
        or not 1 <= contig_length <= _MAX_REFERENCE_LENGTH
    ):
        raise BamError("BAM contig length is invalid")
    expected_indices = (0, 1) if paired_end else (0,)
    mates = tuple(sorted(fragment.mates, key=lambda value: value.mate_index))
    if tuple(mate.mate_index for mate in mates) != expected_indices:
        raise BamError("BAM mate cardinality disagrees with output mode")
    try:
        qname = format_fragment_identifier(
            fragment.contig_name,
            fragment.reference_start,
            fragment.reference_end,
            fragment.fragment_ordinal,
        )
    except ReadNameError as error:
        raise BamError(str(error)) from error
    _validate_query_name(qname)

    alignments = tuple(
        _alignment_for_mate(mate, contig_length=contig_length) for mate in mates
    )
    template_lengths = _template_lengths(alignments) if paired_end else (0,)
    base_state_tag_suffixes = _site_state_suffixes(
        fragment,
        mates,
        include_fragment_summary=include_fragment_summary,
        include_fragment_realization=include_fragment_realization,
    )
    records = []
    for index, (mate, alignment) in enumerate(zip(mates, alignments)):
        flag = 0
        mate_fields = ("*", 0)
        tags = [
            "RG:Z:{}".format(read_group_id),
            "AS:i:{}".format(len(alignment.sequence)),
        ]
        if paired_end:
            other = alignments[1 - index]
            flag |= 0x1 | 0x2 | (0x40 if mate.mate_index == 0 else 0x80)
            if mate.reverse_complement:
                flag |= 0x10
            if mates[1 - index].reverse_complement:
                flag |= 0x20
            mate_fields = ("=", other.position)
            tags.append("MQ:i:{}".format(BAM_MAPQ))
            tags.append("MC:Z:{}".format(other.cigar))
        elif mate.reverse_complement:
            flag |= 0x10
        fields = (
            qname,
            str(flag),
            fragment.contig_name,
            str(alignment.position),
            str(BAM_MAPQ),
            alignment.cigar,
            mate_fields[0],
            str(mate_fields[1]),
            str(template_lengths[index]),
            alignment.sequence,
            alignment.quality,
            *tags,
        )
        records.append(
            "\t".join(fields).encode("ascii")
            + base_state_tag_suffixes[index]
            + b"\n"
        )
    return tuple(records)


def format_sam_batch(
    fragments: tuple[ProcessedFragment, ...],
    *,
    paired_end: bool,
    read_group_id: str,
    contig_lengths: tuple[int, ...],
    include_fragment_summary: bool = False,
    include_fragment_realization: bool = False,
) -> tuple[bytes, tuple[int, ...]]:
    """Format one trusted worker batch with the required C-extension encoder."""

    if not isinstance(fragments, tuple) or not fragments:
        raise BamError("details SAM batch must be a non-empty tuple")
    if not isinstance(paired_end, bool):
        raise BamError("details SAM paired_end must be a boolean")
    if not isinstance(contig_lengths, tuple) or len(contig_lengths) != len(fragments):
        raise BamError("details SAM contig lengths disagree with fragment count")
    if not isinstance(include_fragment_summary, bool):
        raise BamError("BAM fragment-summary policy must be boolean")
    if not isinstance(include_fragment_realization, bool):
        raise BamError("BAM fragment-realization policy must be boolean")
    if include_fragment_realization and not include_fragment_summary:
        raise BamError("BAM fragment realization requires fragment summary")
    tag_suffixes = tuple(
        suffix
        for fragment in fragments
        for suffix in _site_state_suffixes(
            fragment,
            tuple(sorted(fragment.mates, key=lambda value: value.mate_index)),
            include_fragment_summary=include_fragment_summary,
            include_fragment_realization=include_fragment_realization,
        )
    )
    return _cext_format_sam_batch(
        fragments,
        paired_end,
        read_group_id,
        contig_lengths,
        tag_suffixes,
    )


def format_sam_columns(
    batch,
    header: Header,
    formatted,
    *,
    paired_end: bool,
    include_fragment_summary: bool = False,
    include_fragment_realization: bool = False,
):
    """Encode one validated Full-Details batch without typed fragment objects."""

    if not isinstance(header, Header) or not header.has_details:
        raise BamError("columnar details SAM requires a Full Details header")
    if (
        batch.details is None
        or formatted.mate_sequences is None
        or formatted.base_state_codes is None
        or formatted.read_summaries is None
    ):
        raise BamError("columnar details SAM inputs are incomplete")
    if not isinstance(paired_end, bool):
        raise BamError("columnar details SAM paired_end must be a boolean")
    if not isinstance(include_fragment_summary, bool):
        raise BamError("BAM fragment-summary policy must be boolean")
    if not isinstance(include_fragment_realization, bool):
        raise BamError("BAM fragment-realization policy must be boolean")
    if include_fragment_realization and not include_fragment_summary:
        raise BamError("BAM fragment realization requires fragment summary")
    fragment_summaries = getattr(formatted, "fragment_summaries", None)
    fragment_realizations = getattr(formatted, "fragment_realizations", None)
    if include_fragment_summary != (fragment_summaries is not None):
        raise BamError("columnar fragment-summary data disagrees with policy")
    if include_fragment_realization != (fragment_realizations is not None):
        raise BamError("columnar fragment-realization data disagrees with policy")
    details = batch.details
    columns = (
        batch.reference_starts.raw,
        batch.reference_ends.raw,
        batch.mate_offsets.raw,
        batch.mate_template_starts.raw,
        batch.mate_template_ends.raw,
        batch.mate_indices.raw,
        batch.mate_reverse_complements.raw,
        details.projection_offsets.raw,
        details.projection_template_starts.raw,
        details.projection_template_ends.raw,
        details.projection_reference_starts.raw,
        details.variant_offsets.raw,
        details.variant_reference_starts.raw,
        details.variant_template_starts.raw,
        details.variant_template_ends.raw,
    )
    return _cext_format_sam_columns(
        tuple(header.contigs[index].name for index in batch.contig_indices),
        columns,
        formatted.mate_sequences,
        formatted.base_state_codes,
        formatted.read_summaries,
        fragment_summaries,
        fragment_realizations,
        formatted.read_length,
        formatted.quality_byte,
        int(paired_end),
        batch.first_fragment_ordinal,
        header.run_id,
        tuple(header.contigs[index].length for index in batch.contig_indices),
    )


_CONTEXT_CODES = {
    MethylationContext.CG_C: 1,
    MethylationContext.CG_G: 1,
    MethylationContext.CHG_C: 2,
    MethylationContext.CHG_G: 2,
    MethylationContext.CHH_C: 3,
    MethylationContext.CHH_G: 3,
}
_FLAG_HAS_ASM = 1 << 7
_FLAG_HAS_SNV = 1 << 8
_FLAG_HAS_INS = 1 << 9
_FLAG_HAS_DEL = 1 << 10
_FLAG_HAS_METH_CONV = 1 << 11
_FLAG_COUNT_OVERFLOW = 1 << 13
_SUMMARY_MAX = 0xFFFF


def _site_state_suffixes(
    fragment: ProcessedFragment,
    mates: tuple[ProcessedMate, ...],
    *,
    include_fragment_summary: bool,
    include_fragment_realization: bool,
) -> tuple[bytes, ...]:
    variant_by_index = {variant.index: variant for variant in fragment.variants}
    site_by_index = {site.site_index: site for site in fragment.site_states}
    mate_values = tuple(
        _read_summary_values(fragment, mate, variant_by_index, site_by_index)
        for mate in mates
    )
    fragment_summary = (
        _finalize_summary(_fragment_summary_values(fragment, mate_values))
        if include_fragment_summary
        else None
    )
    fragment_realization = (
        fragment.fragment_realization if include_fragment_realization else None
    )
    if include_fragment_realization and fragment_realization is None:
        raise BamError("fragment realization is missing")
    suffixes = []
    for mate, (state_text, raw_summary) in zip(mates, mate_values):
        if mate.reverse_complement:
            state_text = state_text[::-1]
        fields = [
            "zt:Z:{}".format(state_text),
            _summary_tag("zr", _finalize_summary(raw_summary)),
        ]
        if fragment_summary is not None:
            fields.append(_summary_tag("zf", fragment_summary))
        if fragment_realization is not None:
            try:
                fields.append("zx:Z:" + fragment_realization.decode("ascii"))
            except UnicodeDecodeError as error:
                raise BamError("fragment realization must be ASCII") from error
        suffixes.append(("\t" + "\t".join(fields)).encode("ascii"))
    return tuple(suffixes)


def _read_summary_values(
    fragment: ProcessedFragment,
    mate: ProcessedMate,
    variant_by_index: dict[int, object],
    site_by_index: dict[int, object],
) -> tuple[str, list[int]]:
    base_states = mate.base_states
    state_characters = []
    context_counts = [0, 0, 0]
    methylated_counts = [0, 0, 0]
    n_converted = 0
    n_conversion_failed = 0
    n_sequencing_errors = 0
    n_false_methylation = 0
    variant_ids = set()
    has_asm = False
    if isinstance(base_states, _CompactBaseStates):
        rows = (
            (
                base_states.site_indices[index],
                base_states.methylated[index],
                base_states.succeeded[index],
                base_states.attempted[index],
                base_states.variant_indices[index],
                base_states.error_flags[index],
                base_states.post_conversion_bases[index],
                base_states.final_bases[index],
            )
            for index in range(len(base_states))
        )
    else:
        rows = (
            (
                annotation.site_index,
                annotation.methylated,
                annotation.conversion_succeeded,
                annotation.conversion_attempted,
                annotation.variant_index,
                annotation.sequencing_error,
                annotation.post_conversion_base,
                annotation.final_base,
            )
            for annotation in base_states
        )
    for (
        site_index,
        methylated,
        converted,
        attempted,
        variant_index,
        sequencing_error,
        post_conversion_base,
        final_base,
    ) in rows:
        context_code = 0
        if site_index is not None:
            site = site_by_index.get(site_index)
            if site is None:
                raise BamError("BAM annotation refers to an unknown site")
            context_code = _CONTEXT_CODES[site.context]
            has_asm = has_asm or site.methylation_source is MethylationSource.ASM
            context_counts[context_code - 1] += 1
            if methylated is True:
                methylated_counts[context_code - 1] += 1
        if variant_index != NO_VARIANT_INDEX:
            if variant_index not in variant_by_index:
                raise BamError("BAM annotation refers to an unknown variant")
            variant_ids.add(variant_index)
        n_converted += int(bool(converted))
        n_conversion_failed += int(bool(attempted) and not bool(converted))
        n_sequencing_errors += int(bool(sequencing_error))
        false_methylation = bool(sequencing_error) and (
            (int(mate.conversion_mode) == 0 and post_conversion_base == 3 and final_base == 1)
            or (
                int(mate.conversion_mode) == 1
                and post_conversion_base == 0
                and final_base == 2
            )
        )
        n_false_methylation += int(false_methylation)
        state = (
            context_code
            | (int(methylated is True) << 2)
            | (int(bool(converted)) << 3)
            | (int(variant_index != NO_VARIANT_INDEX) << 4)
            | (int(bool(sequencing_error)) << 5)
        )
        state_characters.append(ANNOTATION_STATE_ALPHABET[state])

    flags = (
        (int(fragment.haplotype) & 0x3)
        | ((int(fragment.capture_strand) & 0x3) << 2)
        | ((int(mate.conversion_mode) & 0x7) << 4)
    )
    if has_asm:
        flags |= _FLAG_HAS_ASM
    for variant_index in variant_ids:
        flags |= _variant_flag(variant_by_index[variant_index].kind)
    if n_converted:
        flags |= _FLAG_HAS_METH_CONV
    summary = [
        flags,
        context_counts[0],
        methylated_counts[0],
        context_counts[1],
        methylated_counts[1],
        context_counts[2],
        methylated_counts[2],
        n_converted,
        n_conversion_failed,
        len(variant_ids),
        n_sequencing_errors,
        n_false_methylation,
    ]
    return "".join(state_characters), summary


def _fragment_summary_values(
    fragment: ProcessedFragment,
    mate_values: tuple[tuple[str, list[int]], ...],
) -> list[int]:
    context_counts = [0, 0, 0]
    methylated_counts = [0, 0, 0]
    has_asm = False
    for site in fragment.site_states:
        context_code = _CONTEXT_CODES[site.context]
        context_counts[context_code - 1] += 1
        methylated_counts[context_code - 1] += int(bool(site.methylated))
        has_asm = has_asm or site.methylation_source is MethylationSource.ASM
    flags = (
        (int(fragment.haplotype) & 0x3)
        | ((int(fragment.capture_strand) & 0x3) << 2)
        | ((int(fragment.fragment_conversion_mode) & 0x7) << 4)
    )
    if has_asm:
        flags |= _FLAG_HAS_ASM
    for event in fragment.variants:
        flags |= _variant_flag(event.kind)
    if fragment.fragment_conversion_success_count:
        flags |= _FLAG_HAS_METH_CONV
    return [
        flags,
        context_counts[0],
        methylated_counts[0],
        context_counts[1],
        methylated_counts[1],
        context_counts[2],
        methylated_counts[2],
        fragment.fragment_conversion_success_count,
        fragment.fragment_conversion_failure_count,
        len(fragment.variants),
        sum(value[1][10] for value in mate_values),
        sum(value[1][11] for value in mate_values),
    ]


def _variant_flag(kind) -> int:
    if kind is VariantKind.SNV:
        return _FLAG_HAS_SNV
    if kind is VariantKind.INSERTION:
        return _FLAG_HAS_INS
    if kind is VariantKind.DELETION:
        return _FLAG_HAS_DEL
    raise BamError("BAM variant kind is unsupported")


def _finalize_summary(values: list[int]) -> tuple[int, ...]:
    if len(values) != 12 or any(value < 0 for value in values):
        raise BamError("BAM summary has invalid values")
    result = list(values)
    if any(value > _SUMMARY_MAX for value in result[1:]):
        result[0] |= _FLAG_COUNT_OVERFLOW
    result[1:] = [min(value, _SUMMARY_MAX) for value in result[1:]]
    if result[0] > _SUMMARY_MAX:
        raise BamError("BAM summary flags exceed uint16")
    return tuple(result)


def _summary_tag(name: str, values: tuple[int, ...]) -> str:
    return "{}:B:S,{}".format(name, ",".join(str(value) for value in values))


@dataclass(frozen=True)
class _MateAlignment:
    cigar: str
    position: int
    reference_start: int
    reference_end: int
    sequence: str
    quality: str
    mate_index: int


def _alignment_for_mate(
    mate: ProcessedMate,
    *,
    contig_length: int,
) -> _MateAlignment:
    if not isinstance(mate, ProcessedMate):
        raise BamError("BAM mates must be ProcessedMate values")
    if not mate.sequence or len(mate.sequence) != len(mate.quality):
        raise BamError("BAM sequence and quality lengths disagree")
    positions = _reference_projection_positions(mate)
    if len(positions) != len(mate.sequence):
        raise BamError("BAM annotations must cover every read base")
    if mate.reverse_complement:
        positions = tuple(reversed(positions))
        sequence = mate.sequence.translate(_SEQUENCE_COMPLEMENT)[::-1]
        quality = mate.quality[::-1]
    else:
        sequence = mate.sequence
        quality = mate.quality

    operations = []  # type: list[tuple[str, int]]
    previous_mapped = None  # type: int | None
    mapped_positions = []
    for position in positions:
        if isinstance(position, bool) or not isinstance(position, int):
            raise BamError("BAM reference positions must be integers")
        if position == -1:
            _append_cigar_operation(operations, "I", 1)
            continue
        if not 0 <= position < contig_length:
            raise BamError("BAM reference position is outside its contig")
        if previous_mapped is not None:
            gap = position - previous_mapped - 1
            if gap < 0:
                raise BamError(
                    "BAM reference projection is not strictly increasing"
                )
            if gap:
                _append_cigar_operation(operations, "D", gap)
        _append_cigar_operation(operations, "M", 1)
        previous_mapped = position
        mapped_positions.append(position)

    if mapped_positions:
        reference_start = mapped_positions[0]
        reference_end = mapped_positions[-1] + 1
        if (
            mate.reference_start != reference_start
            or mate.reference_end != reference_end
        ):
            raise BamError(
                "BAM mate bounds disagree with its reference projection"
            )
    else:
        if mate.reference_start != mate.reference_end:
            raise BamError("insertion-only BAM mate has non-empty bounds")
        if not 0 <= mate.reference_start <= contig_length:
            raise BamError("insertion-only BAM anchor is outside its contig")
        # A pure-I CIGAR consumes no reference. Keep POS inside [1, LN],
        # including for an insertion anchored immediately after the final base.
        reference_start = min(mate.reference_start, contig_length - 1)
        reference_end = reference_start + 1
    query_length = sum(
        length for operation, length in operations if operation in ("M", "I")
    )
    if query_length != len(sequence):
        raise BamError("BAM CIGAR does not consume the complete read")
    cigar = "".join(
        "{}{}".format(length, operation) for operation, length in operations
    )
    return _MateAlignment(
        cigar=cigar,
        position=reference_start + 1,
        reference_start=reference_start,
        reference_end=reference_end,
        sequence=sequence,
        quality=quality,
        mate_index=mate.mate_index,
    )


def _reference_projection_positions(mate: ProcessedMate) -> tuple[int, ...]:
    base_states = mate.base_states
    if isinstance(base_states, _CompactBaseStates):
        return tuple(base_states.reference_positions)
    if not isinstance(base_states, tuple):
        raise BamError("BAM base states must be immutable")
    positions = []
    for offset, annotation in enumerate(base_states):
        if (
            not isinstance(annotation, BaseState)
            or annotation.read_offset != offset
        ):
            raise BamError("BAM base states are not in read order")
        positions.append(annotation.reference_pos)
    return tuple(positions)


def _append_cigar_operation(
    operations: list[tuple[str, int]], operation: str, length: int
) -> None:
    while length:
        if operations and operations[-1][0] == operation:
            available = _MAX_CIGAR_OPERATION_LENGTH - operations[-1][1]
            if available:
                consumed = min(available, length)
                operations[-1] = (operation, operations[-1][1] + consumed)
                length -= consumed
                continue
        consumed = min(_MAX_CIGAR_OPERATION_LENGTH, length)
        operations.append((operation, consumed))
        length -= consumed


def _template_lengths(alignments: tuple[_MateAlignment, ...]) -> tuple[int, int]:
    left_index = min(
        range(2),
        key=lambda index: (
            alignments[index].reference_start,
            alignments[index].mate_index,
        ),
    )
    span = max(item.reference_end for item in alignments) - min(
        item.reference_start for item in alignments
    )
    if not 0 < span <= _MAX_TEMPLATE_LENGTH:
        raise BamError("BAM template length exceeds signed 32-bit range")
    values = [0, 0]
    values[left_index] = span
    values[1 - left_index] = -span
    return values[0], values[1]


def _validate_reference_name(value: str) -> None:
    if (
        not isinstance(value, str)
        or not value
        or value[0] not in _REFERENCE_FIRST
        or any(character not in _REFERENCE_REST for character in value[1:])
    ):
        raise BamError(
            "reference name is outside the SAM 1.6 character contract: {!r}".format(
                value
            )
        )


def _validate_query_name(value: str) -> None:
    if (
        not isinstance(value, str)
        or not 1 <= len(value) <= 254
        or any(
            not 33 <= ord(character) <= 126 or character == "@"
            for character in value
        )
    ):
        raise BamError("BAM query name violates the SAM 1.6 contract")


def _validate_header_value(name: str, value: str) -> None:
    if (
        not isinstance(value, str)
        or not value
        or any(not 32 <= ord(character) <= 126 for character in value)
        or "\t" in value
    ):
        raise BamError("BAM {} is not printable ASCII".format(name))


BytesLike = bytes | bytearray | memoryview


@dataclass(frozen=True)
class BamConfig:
    """Validated streaming contract for one HTSlib-backed BAM."""

    writer_argv: tuple[str, ...]
    sam_header: bytes
    references: tuple[tuple[str, int], ...]
    read_group_id: str
    fragment_summary: bool = False
    fragment_realization: bool = False

    def __post_init__(self) -> None:
        if (
            not isinstance(self.writer_argv, tuple)
            or not self.writer_argv
            or any(
                not isinstance(value, str) or not value or "\x00" in value
                for value in self.writer_argv
            )
        ):
            raise OutputError("BAM writer argv must be non-empty text")
        if type(self.sam_header) is not bytes or not self.sam_header.startswith(
            b"@HD\t"
        ):
            raise OutputError("BAM SAM header is invalid")
        if not self.sam_header.endswith(b"\n"):
            raise OutputError("BAM SAM header must end with a newline")
        if not isinstance(self.references, tuple) or not self.references:
            raise OutputError("BAM references must be non-empty")
        observed = set()
        for reference in self.references:
            if (
                not isinstance(reference, tuple)
                or len(reference) != 2
                or not isinstance(reference[0], str)
                or not reference[0]
                or isinstance(reference[1], bool)
                or not isinstance(reference[1], int)
                or reference[1] <= 0
                or reference[0] in observed
            ):
                raise OutputError("BAM references are invalid")
            observed.add(reference[0])
        if not isinstance(self.read_group_id, str) or not self.read_group_id:
            raise OutputError("BAM read group identifier is invalid")
        if not isinstance(self.fragment_summary, bool):
            raise OutputError("BAM fragment_summary must be a boolean")
        if not isinstance(self.fragment_realization, bool):
            raise OutputError("BAM fragment_realization must be a boolean")
        if self.fragment_realization and not self.fragment_summary:
            raise OutputError("BAM fragment_realization requires fragment_summary")

    def reference_length(self, name: str) -> int:
        for reference_name, length in self.references:
            if reference_name == name:
                return length
        raise OutputError("BAM fragment refers to an unknown contig")


class BamOutput:
    """Stream SAM into the bundled HTSlib writer and identify its BAM bytes."""

    _WAIT_SECONDS = 30

    def __init__(self, path: Path, config: BamConfig) -> None:
        self.path = path
        self.closed = False
        self._completed = False
        self._size_bytes = 0
        self._sha256 = ""
        self._stderr_path = path.with_name(path.name + ".writer-stderr")
        self.raw = None  # type: BinaryIO | None
        self.stderr = None  # type: BinaryIO | None
        self.process = None  # type: subprocess.Popen | None
        try:
            self.raw = path.open("xb")
            self.stderr = self._stderr_path.open("xb")
            self.process = subprocess.Popen(
                config.writer_argv,
                stdin=subprocess.PIPE,
                stdout=self.raw,
                stderr=self.stderr,
                close_fds=True,
            )
            if self.process.stdin is None:
                raise OutputError("BAM writer has no SAM input stream")
            self.write_bytes(config.sam_header)
        except Exception:
            self._stop_process()
            self._close_files()
            self._unlink_stderr()
            raise

    def write_bytes(self, value: BytesLike) -> None:
        if self.closed:
            raise OutputError("cannot write to a closed BAM")
        if not isinstance(value, (bytes, bytearray, memoryview)):
            raise OutputError("BAM accepts bytes-like SAM records only")
        if self.process is None or self.process.stdin is None:
            raise OutputError("BAM writer is unavailable")
        try:
            self.process.stdin.write(value)
        except (BrokenPipeError, OSError) as error:
            raise OutputError(
                "BAM writer closed its SAM input early{}".format(
                    self._stderr_suffix()
                )
            ) from error

    def close(self) -> None:
        if self.closed:
            return
        error = None  # type: BaseException | None
        status = None
        try:
            if self.process is None or self.process.stdin is None:
                raise OutputError("BAM writer is unavailable")
            with suppress((BrokenPipeError, OSError)):
                self.process.stdin.close()
            try:
                status = self.process.wait(timeout=self._WAIT_SECONDS)
            except subprocess.TimeoutExpired as timeout_error:
                self._stop_process()
                raise OutputError(
                    "BAM writer did not terminate"
                ) from timeout_error
            if self.raw is None:
                raise OutputError("BAM staged file is unavailable")
            self.raw.flush()
            os.fsync(self.raw.fileno())
            if status != 0:
                raise OutputError(
                    "BAM writer exited with status {}{}".format(
                        status, self._stderr_suffix()
                    )
                )
        except BaseException as observed:
            error = observed
        finally:
            self._close_files()
            self.closed = True

        if error is None:
            try:
                digest = hashlib.sha256()
                size = 0
                with self.path.open("rb") as staged:
                    while True:
                        block = staged.read(1024 * 1024)
                        if not block:
                            break
                        size += len(block)
                        digest.update(block)
                self._size_bytes = size
                self._sha256 = digest.hexdigest()
                self._completed = True
            except BaseException as observed:
                error = observed
        self._unlink_stderr()
        if error is not None:
            if isinstance(error, OutputError):
                raise error
            raise OutputError("cannot finalize BAM: {}".format(error)) from error

    def _stop_process(self) -> None:
        process = self.process
        if process is None:
            return
        if process.stdin is not None:
            with suppreii(OSError):
                process.stdin.close()
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)

    def _close_files(self) -> None:
        for stream in (self.raw, self.stderr):
            if stream is not None and not stream.closed:
                with suppress(OSError):
                    stream.close()

    def _stderr_suffix(self) -> str:
        try:
            if self.stderr is not None and not self.stderr.closed:
                self.stderr.flush()
            with self._stderr_path.open("rb") as stream:
                value = stream.read(8192).decode("utf-8", errors="replace").strip()
        except OSError:
            return ""
        return ": {}".format(value) if value else ""

    def _unlink_stderr(self) -> None:
        with suppress(FileNotFoundError):
            self._stderr_path.unlink()

    @property
    def size_bytes(self) -> int:
        if not self._completed:
            raise OutputError("staged BAM identity requires a closed stream")
        return self._size_bytes

    @property
    def sha256(self) -> str:
        if not self._completed:
            raise OutputError("staged BAM identity requires a closed stream")
        return self._sha256


__all__ = [
    "BAM_CONTRACT",
    "BAM_MAPQ",
    "BamConfig",
    "BamError",
    "BamOutput",
    "build_sam_header",
    "format_sam_batch",
    "format_sam_columns",
    "format_sam_fragment",
    "validate_bam_header",
]
