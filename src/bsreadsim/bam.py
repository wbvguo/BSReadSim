"""SAM 1.6 formatting for BSReadSim truth-aligned BAM output.

The Python layer derives exact linear alignments from the Full Truth reference
projection.  HTSlib owns only the final SAM-to-BAM serialization boundary.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

from .postprocess import (
    BaseAnnotation,
    ProcessedFragment,
    ProcessedMate,
    _CompactAnnotations,
)
from .protocol import Header, TruthMode
from .read_names import ReadNameError, format_fragment_identifier


TRUTH_BAM_CONTRACT = "bsreadsim-truth-bam-v1"
TRUTH_BAM_MAPQ = 60
_MAX_REFERENCE_LENGTH = (1 << 31) - 1
_MAX_CIGAR_OPERATION_LENGTH = (1 << 28) - 1
_MAX_TEMPLATE_LENGTH = (1 << 31) - 1
_REFERENCE_FIRST = frozenset(
    "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz!#$%&+./:;?@^_|~-"
)
_REFERENCE_REST = _REFERENCE_FIRST | frozenset("*=")
_SEQUENCE_COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


class TruthBamError(ValueError):
    """Truth projection cannot be represented as one interoperable BAM record."""


def validate_truth_bam_header(header: Header) -> None:
    """Fail unless a protocol header can form a SAM 1.6 reference dictionary."""

    if not isinstance(header, Header):
        raise TruthBamError("truth BAM requires a protocol Header")
    if header.truth_columns is not TruthMode.FULL:
        raise TruthBamError("truth BAM requires Full Truth protocol columns")
    if not header.contigs:
        raise TruthBamError("truth BAM requires at least one reference sequence")
    observed = set()
    for contig in header.contigs:
        _validate_reference_name(contig.name)
        if contig.name in observed:
            raise TruthBamError("truth BAM reference names must be unique")
        observed.add(contig.name)
        if (
            isinstance(contig.length, bool)
            or not isinstance(contig.length, int)
            or not 1 <= contig.length <= _MAX_REFERENCE_LENGTH
        ):
            raise TruthBamError(
                "truth BAM reference lengths must be in [1, 2^31-1]"
            )


def build_truth_sam_header(
    header: Header,
    *,
    sample_name: str,
    program_version: str,
) -> bytes:
    """Build the complete SAM header streamed to the HTSlib BAM writer."""

    validate_truth_bam_header(header)
    _validate_header_value("sample name", sample_name)
    _validate_header_value("run identifier", header.run_id)
    _validate_header_value("program version", program_version)
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
        "@CO\tBSReadSim truth alignments; MAPQ 60 denotes simulated origin, "
        "not calibrated mapping confidence"
    )
    return ("\n".join(lines) + "\n").encode("ascii")


def format_truth_sam_fragment(
    fragment: ProcessedFragment,
    *,
    paired_end: bool,
    read_group_id: str,
    contig_length: int,
) -> Tuple[bytes, ...]:
    """Format one SE/PE truth alignment, keeping paired records adjacent."""

    if not isinstance(fragment, ProcessedFragment):
        raise TruthBamError("truth BAM fragment must be a ProcessedFragment")
    _validate_header_value("read group identifier", read_group_id)
    _validate_reference_name(fragment.contig_name)
    if (
        isinstance(contig_length, bool)
        or not isinstance(contig_length, int)
        or not 1 <= contig_length <= _MAX_REFERENCE_LENGTH
    ):
        raise TruthBamError("truth BAM contig length is invalid")
    expected_indices = (0, 1) if paired_end else (0,)
    mates = tuple(sorted(fragment.mates, key=lambda value: value.mate_index))
    if tuple(mate.mate_index for mate in mates) != expected_indices:
        raise TruthBamError("truth BAM mate cardinality disagrees with output mode")
    try:
        qname = format_fragment_identifier(
            fragment.contig_name,
            fragment.reference_start,
            fragment.reference_end,
            fragment.fragment_ordinal,
        )
    except ReadNameError as error:
        raise TruthBamError(str(error)) from error
    _validate_query_name(qname)

    alignments = tuple(
        _alignment_for_mate(mate, contig_length=contig_length) for mate in mates
    )
    template_lengths = _template_lengths(alignments) if paired_end else (0,)
    records = []
    for index, (mate, alignment) in enumerate(zip(mates, alignments)):
        flag = 0
        mate_fields = ("*", 0)
        tags = ["RG:Z:{}".format(read_group_id), "PG:Z:bsreadsim"]
        if paired_end:
            other = alignments[1 - index]
            flag |= 0x1 | 0x2 | (0x40 if mate.mate_index == 0 else 0x80)
            if mate.reverse_complement:
                flag |= 0x10
            if mates[1 - index].reverse_complement:
                flag |= 0x20
            mate_fields = ("=", other.position)
            tags.append("MC:Z:{}".format(other.cigar))
        elif mate.reverse_complement:
            flag |= 0x10
        fields = (
            qname,
            str(flag),
            fragment.contig_name,
            str(alignment.position),
            str(TRUTH_BAM_MAPQ),
            alignment.cigar,
            mate_fields[0],
            str(mate_fields[1]),
            str(template_lengths[index]),
            alignment.sequence,
            alignment.quality,
            *tags,
        )
        records.append(("\t".join(fields) + "\n").encode("ascii"))
    return tuple(records)


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
        raise TruthBamError("truth BAM mates must be ProcessedMate values")
    if not mate.sequence or len(mate.sequence) != len(mate.quality):
        raise TruthBamError("truth BAM sequence and quality lengths disagree")
    positions = _annotation_reference_positions(mate)
    if len(positions) != len(mate.sequence):
        raise TruthBamError("truth BAM annotations must cover every read base")
    if mate.reverse_complement:
        positions = tuple(reversed(positions))
        sequence = mate.sequence.translate(_SEQUENCE_COMPLEMENT)[::-1]
        quality = mate.quality[::-1]
    else:
        sequence = mate.sequence
        quality = mate.quality

    operations = []  # type: List[Tuple[str, int]]
    previous_mapped = None  # type: Optional[int]
    mapped_positions = []
    for position in positions:
        if isinstance(position, bool) or not isinstance(position, int):
            raise TruthBamError("truth BAM reference positions must be integers")
        if position == -1:
            _append_cigar_operation(operations, "I", 1)
            continue
        if not 0 <= position < contig_length:
            raise TruthBamError("truth BAM reference position is outside its contig")
        if previous_mapped is not None:
            gap = position - previous_mapped - 1
            if gap < 0:
                raise TruthBamError(
                    "truth BAM reference projection is not strictly increasing"
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
            raise TruthBamError(
                "truth BAM mate bounds disagree with its reference projection"
            )
    else:
        if mate.reference_start != mate.reference_end:
            raise TruthBamError("insertion-only truth BAM mate has non-empty bounds")
        if not 0 <= mate.reference_start <= contig_length:
            raise TruthBamError("insertion-only truth BAM anchor is outside its contig")
        # A pure-I CIGAR consumes no reference. Keep POS inside [1, LN],
        # including for an insertion anchored immediately after the final base.
        reference_start = min(mate.reference_start, contig_length - 1)
        reference_end = reference_start + 1
    query_length = sum(
        length for operation, length in operations if operation in ("M", "I")
    )
    if query_length != len(sequence):
        raise TruthBamError("truth BAM CIGAR does not consume the complete read")
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


def _annotation_reference_positions(mate: ProcessedMate) -> Tuple[int, ...]:
    annotations = mate.annotations
    if isinstance(annotations, _CompactAnnotations):
        return tuple(annotations.reference_positions)
    if not isinstance(annotations, tuple):
        raise TruthBamError("truth BAM annotations must be immutable")
    positions = []
    for offset, annotation in enumerate(annotations):
        if (
            not isinstance(annotation, BaseAnnotation)
            or annotation.read_offset != offset
        ):
            raise TruthBamError("truth BAM annotations are not in read order")
        positions.append(annotation.reference_pos)
    return tuple(positions)


def _append_cigar_operation(
    operations: List[Tuple[str, int]], operation: str, length: int
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


def _template_lengths(alignments: Tuple[_MateAlignment, ...]) -> Tuple[int, int]:
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
        raise TruthBamError("truth BAM template length exceeds signed 32-bit range")
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
        raise TruthBamError(
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
        raise TruthBamError("truth BAM query name violates the SAM 1.6 contract")


def _validate_header_value(name: str, value: str) -> None:
    if (
        not isinstance(value, str)
        or not value
        or any(not 32 <= ord(character) <= 126 for character in value)
        or "\t" in value
    ):
        raise TruthBamError("truth BAM {} is not printable ASCII".format(name))


__all__ = [
    "TRUTH_BAM_CONTRACT",
    "TRUTH_BAM_MAPQ",
    "TruthBamError",
    "build_truth_sam_header",
    "format_truth_sam_fragment",
    "validate_truth_bam_header",
]
