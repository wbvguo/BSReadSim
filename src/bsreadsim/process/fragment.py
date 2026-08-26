"""Acquire decoded fragments and orchestrate the ordered processing flow."""

from __future__ import annotations

import base64

from enum import Enum


import numpy as np

from .batch import (
    CaptureStrand,
    ConversionMode,
    Fragment,
    MethylationModelBatch,
    FastqFragment,
    FastqMate,
    EncodedFastqBatch,
    ColumnarFragmentBatch,
    ColumnarReadBatch,
    ProcessError,
    UniformQuality,
    ProcessedFragment,
    ProcessedMate,
    UniformError,
    _ConvertedFragment,
    _U1,
    _U2,
    _U4,
    _U8,
    _F4,
    _bernoulli_from_pair,
    _contig_fragment_groups,
    _owners_from_offsets,
    _philox_pairs,
)
from .config import ProcessConfig
from .methylation import (
    BernoulliStateModel,
    _materialize_site_states,
    _sample_methylation_batch_values,
    _sample_site_states,
)
from .bisulfite import (
    _convert_fragment,
    _draw_fragment_base_bernoulli,
)
from .sequencing import (
    _apply_errors,
    _derive_mates,
    _draw_uniform_errors,
    _generate_quality,
    _finalize_read,
    _mate_template_offsets,
)
from ..rng import RNGStage, _bernoulli_unchecked as bernoulli, derive_key
from ..htsim.protocol import (
    DecodedBatchView,
    Header,
)
from .._cext import (
    decode_protocol_fragments as _cext_decode_fragments,
    format_fastq_batch as _cext_format_fastq_batch,
    pack_protocol_common_columns as _cext_pack_common_columns,
)

_BASE_ASCII = np.frombuffer(b"ACGTN", dtype=np.uint8)
_SITE_STATE_ASCII = np.frombuffer(
    b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_",
    dtype=np.uint8,
)
_CONTEXT_STATE = np.asarray(
    (0, 1, 0, 2, 0, 0, 0, 3, 0, 1, 0, 2, 0, 0, 0, 3),
    dtype=np.uint8,
)
_COMPLEMENT = np.asarray((3, 2, 1, 0, 4), dtype=np.uint8)
_ALTERNATIVE_BASES = np.asarray(
    ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)),
    dtype=np.uint8,
)


def _resolve_fragment_mode(
    fragment: Fragment,
    config: ProcessConfig,
) -> ConversionMode:
    """Resolve molecule orientation before bisulfite chemistry is applied."""
    if not config.bisulfite:
        return ConversionMode.NONE
    if fragment.capture_strand is CaptureStrand.FORWARD:
        return ConversionMode.C2T
    if fragment.capture_strand is CaptureStrand.REVERSE:
        return ConversionMode.G2A
    if fragment.capture_strand is not CaptureStrand.UNKNOWN:
        raise ProcessError("fragment has an unsupported capture strand")
    if config.directional:
        return ConversionMode.C2T

    key = derive_key(
        config.master_seed,
        RNGStage.LIBRARY_ORIENTATION,
        fragment.contig_index,
    )
    reverse = bernoulli(key, fragment.fragment_ordinal, 0, 0.5)
    return ConversionMode.G2A if reverse else ConversionMode.C2T


def _resolve_fragment_modes(
    capture_strands: np.ndarray,
    ordinals: np.ndarray,
    contig_indices: np.ndarray,
    config: ProcessConfig,
) -> np.ndarray:
    """Resolve vectorized molecule orientations at the fragment boundary."""
    if np.any(capture_strands > 2):
        raise ProcessError("NumPy capture strand contains an invalid enum")
    if not config.bisulfite:
        return np.full(len(capture_strands), int(ConversionMode.NONE), dtype=np.uint8)
    result = np.zeros(len(capture_strands), dtype=np.uint8)
    result[capture_strands == 2] = 1
    unknown = np.flatnonzero(capture_strands == 0)
    if unknown.size and not config.directional:
        groups = _contig_fragment_groups(contig_indices)
        if len(groups) == 1:
            contig_index = next(iter(groups))
            key = derive_key(
                config.master_seed,
                RNGStage.LIBRARY_ORIENTATION,
                contig_index,
            )
            pair0, _ = _philox_pairs(
                key,
                ordinals[unknown],
                np.zeros(unknown.size, dtype=np.uint64),
            )
            result[unknown] = _bernoulli_from_pair(pair0, 0.5)
            return result
        for contig_index, fragments in groups.items():
            positions = np.intersect1d(unknown, fragments, assume_unique=True)
            if positions.size == 0:
                continue
            key = derive_key(
                config.master_seed,
                RNGStage.LIBRARY_ORIENTATION,
                contig_index,
            )
            pair0, _ = _philox_pairs(
                key,
                ordinals[positions],
                np.zeros(positions.size, dtype=np.uint64),
            )
            result[positions] = _bernoulli_from_pair(pair0, 0.5)
    return result


def decode_fragments(
    batch: DecodedBatchView,
    header: Header,
) -> tuple[Fragment, ...]:
    """Reconstruct typed processing values through the required C extension."""

    if not isinstance(batch, DecodedBatchView):
        raise ProcessError("batch must be a decoded protocol batch")
    if not isinstance(header, Header):
        raise ProcessError("header must be a protocol Header")
    has_details = batch.details is not None
    if header.has_details is not has_details:
        raise ProcessError("batch details columns disagree with the header")

    common_columns = (
        batch.contig_indices.raw,
        batch.reference_starts.raw,
        batch.reference_ends.raw,
        batch.template_offsets.raw,
        batch.mate_offsets.raw,
        batch.site_offsets.raw,
        batch.mate_template_starts.raw,
        batch.mate_template_ends.raw,
        batch.site_template_offsets.raw,
        batch.site_probabilities.raw,
        batch.haplotypes.raw,
        batch.capture_strands.raw,
        batch.mate_indices.raw,
        batch.mate_reverse_complements.raw,
        batch.site_contexts.raw,
        batch.methylation_sources.raw,
        batch.site_alleles.raw,
        batch.template_bases.raw,
    )
    has_details = None
    if batch.details is not None:
        details = batch.details
        has_details = (
            details.projection_offsets.raw,
            details.variant_offsets.raw,
            details.original_n_offsets.raw,
            details.projection_template_starts.raw,
            details.projection_template_ends.raw,
            details.projection_reference_starts.raw,
            details.variant_indices.raw,
            details.variant_id_offsets.raw,
            details.variant_reference_starts.raw,
            details.variant_reference_ends.raw,
            details.variant_template_starts.raw,
            details.variant_template_ends.raw,
            details.variant_ref_offsets.raw,
            details.variant_alt_offsets.raw,
            details.site_reference_positions.raw,
            details.original_n_template_offsets.raw,
            details.variant_sources.raw,
            details.variant_kinds.raw,
            details.variant_phased_haplotypes.raw,
            details.variant_ids.raw,
            details.variant_ref_bases.raw,
            details.variant_alt_bases.raw,
        )
    try:
        return _cext_decode_fragments(
            common_columns,
            has_details,
            batch.first_fragment_ordinal,
            tuple(contig.length for contig in header.contigs),
            header.mates_per_fragment,
            header.read_length_r1,
            header.read_length_r2,
        )
    except (BufferError, OverflowError, TypeError, ValueError) as error:
        raise ProcessError(str(error)) from error


def decode_common_numpy_batch(
    batch: DecodedBatchView,
    header: Header,
) -> ColumnarFragmentBatch:
    """Build exact common-column inputs with a Python reference implementation."""

    if not isinstance(batch, DecodedBatchView):
        raise ProcessError("batch must be a decoded protocol batch")
    if not isinstance(header, Header):
        raise ProcessError("header must be a protocol Header")
    if header.has_details is not (batch.details is not None):
        raise ProcessError("common-column Details policy is inconsistent")

    common_columns = (
        batch.contig_indices.raw,
        batch.reference_starts.raw,
        batch.reference_ends.raw,
        batch.template_offsets.raw,
        batch.mate_offsets.raw,
        batch.site_offsets.raw,
        batch.mate_template_starts.raw,
        batch.mate_template_ends.raw,
        batch.site_template_offsets.raw,
        batch.site_probabilities.raw,
        batch.haplotypes.raw,
        batch.capture_strands.raw,
        batch.mate_indices.raw,
        batch.mate_reverse_complements.raw,
        batch.site_contexts.raw,
        batch.methylation_sources.raw,
        batch.site_alleles.raw,
        batch.template_bases.raw,
    )
    try:
        ordinals, contig_names, model_columns, fragment_columns = (
            _cext_pack_common_columns(
                common_columns,
                batch.first_fragment_ordinal,
                tuple(contig.name for contig in header.contigs),
                tuple(contig.length for contig in header.contigs),
                header.mates_per_fragment,
                header.read_length_r1,
                header.read_length_r2,
            )
        )
    except (BufferError, OverflowError, TypeError, ValueError) as error:
        raise ProcessError(str(error)) from error
    model = MethylationModelBatch(ordinals, contig_names, *model_columns)
    details = batch.details
    return ColumnarFragmentBatch(
        model,
        *fragment_columns,
        variant_offsets=(
            _as_bytes(details.variant_offsets, _U4) if details is not None else b""
        ),
        variant_template_starts=(
            _as_bytes(details.variant_template_starts, _U4)
            if details is not None
            else b""
        ),
        variant_template_ends=(
            _as_bytes(details.variant_template_ends, _U4)
            if details is not None
            else b""
        ),
        variant_kinds=(
            _as_bytes(details.variant_kinds, _U1) if details is not None else b""
        ),
    )


def _as_bytes(values, dtype: np.dtype) -> bytes:
    raw = getattr(values, "raw", None)
    if (
        raw is not None
        and raw.format == "B"
        and raw.nbytes == len(values) * dtype.itemsize
    ):
        if dtype in (_U1, _U4, _F4):
            return raw.tobytes()
    return np.asarray(tuple(values), dtype=dtype).tobytes()


def process_fragment_batch(
    fragments: tuple[Fragment, ...],
    contig_names: tuple[str, ...],
    config: ProcessConfig,
    *,
    compact_base_states: bool = False,
    include_details: bool = True,
    include_fragment_realization: bool = False,
) -> tuple[ProcessedFragment, ...]:
    """Apply Python stages to one batch with exactly one model batch call.

    Model results are validated against the exact input ordinal vector and
    site cardinalities before any row is assigned to a fragment.  Subsequent
    conversion, quality, error, and format-ready stages retain their existing
    fragment-local deterministic RNG addresses.
    """

    _validate_fragment_batch_request(
        fragments,
        contig_names,
        config,
        compact_base_states=compact_base_states,
        include_details=include_details,
        include_fragment_realization=include_fragment_realization,
    )
    sampled_batch = (
        _sample_methylation_batch_values(fragments, config)
        if config.bisulfite
        else tuple(() for _ in fragments)
    )
    return tuple(
        _process_fragment_with_states(
            fragment,
            contig_name,
            config,
            sampled_methylation,
            compact_base_states=compact_base_states,
            include_details=include_details,
            include_fragment_realization=include_fragment_realization,
        )
        for fragment, contig_name, sampled_methylation in zip(
            fragments,
            contig_names,
            sampled_batch,
            strict=True,
        )
    )


def _process_fragment_with_states(
    fragment: Fragment,
    contig_name: str,
    config: ProcessConfig,
    sampled_methylation: tuple[bool, ...],
    *,
    compact_base_states: bool,
    include_details: bool,
    include_fragment_realization: bool = False,
) -> ProcessedFragment:
    fragment_mode = _resolve_fragment_mode(fragment, config)
    if config.bisulfite:
        site_states = (
            _materialize_site_states(fragment, sampled_methylation)
            if include_details
            else ()
        )
        converted_fragment = _convert_fragment(
            fragment,
            fragment_mode,
            sampled_methylation,
            config,
        )
    else:
        site_states = ()
        template_length = len(fragment.template_bases)
        converted_fragment = _ConvertedFragment(
            fragment_mode,
            fragment.template_bases,
            (None,) * template_length,
            (None,) * template_length,
            (False,) * template_length,
            (False,) * template_length,
        )
    converted_mates = _derive_mates(fragment, converted_fragment)
    quality_mates = tuple(
        _generate_quality(fragment, mate, config)
        for mate in converted_mates
    )
    errored_mates = tuple(
        _apply_errors(fragment, mate, config)
        for mate in quality_mates
    )
    processed_mates = tuple(
        _finalize_read(
            mate,
            compact_base_states=compact_base_states,
            include_details=include_details,
        )
        for mate in errored_mates
    )

    return ProcessedFragment(
        fragment.fragment_ordinal,
        contig_name,
        fragment.reference_start,
        fragment.reference_end,
        fragment.haplotype,
        fragment_mode,
        fragment.variants,
        site_states,
        processed_mates,
        fragment.capture_strand,
        sum(bool(value) for value in converted_fragment.succeeded),
        sum(
            bool(attempted) and not bool(succeeded)
            for attempted, succeeded in zip(
                converted_fragment.attempted,
                converted_fragment.succeeded,
                strict=True,
            )
        ),
        _encode_typed_fragment_realization(sampled_methylation, converted_fragment)
        if include_fragment_realization and config.bisulfite
        else None,
    )


def _validate_fragment_batch_request(
    fragments: tuple[Fragment, ...],
    contig_names: tuple[str, ...],
    config: ProcessConfig,
    *,
    compact_base_states: bool,
    include_details: bool,
    include_fragment_realization: bool = False,
) -> None:
    if not isinstance(fragments, tuple) or not fragments:
        raise ProcessError("fragments must be a non-empty immutable tuple")
    if not isinstance(contig_names, tuple):
        raise ProcessError("contig_names must be an immutable tuple")
    if len(contig_names) != len(fragments):
        raise ProcessError("contig_names count disagrees with fragments")
    if not isinstance(config, ProcessConfig):
        raise ProcessError("config must be ProcessConfig")
    if not isinstance(compact_base_states, bool):
        raise ProcessError("compact_base_states must be a boolean")
    if not isinstance(include_details, bool):
        raise ProcessError("include_details must be a boolean")
    if not isinstance(include_fragment_realization, bool):
        raise ProcessError("include_fragment_realization must be a boolean")
    if include_fragment_realization and not include_details:
        raise ProcessError("fragment realization requires Full Details")
    if compact_base_states and not include_details:
        raise ProcessError("compact_base_states requires include_details")
    for fragment, contig_name in zip(fragments, contig_names, strict=True):
        if not isinstance(fragment, Fragment):
            raise ProcessError("fragment must be a decoded protocol Fragment")
        if (
            not isinstance(contig_name, str)
            or not contig_name
            or "\x00" in contig_name
        ):
            raise ProcessError("contig_name must be non-empty text without NUL")


class _ResultMode(Enum):
    FRAGMENTS = "fragments"
    FORMATTED_FASTQ = "formatted-fastq"
    COLUMNAR_READS = "columnar-reads"


def supports_common_processing(config: ProcessConfig) -> bool:
    """Whether common columns support the uniform vectorized process path."""

    return (
        isinstance(config.methylation_model, BernoulliStateModel)
        and isinstance(config.quality, UniformQuality)
        and isinstance(config.error, UniformError)
    )


def generate_columnar_reads(
    batch: ColumnarFragmentBatch,
    config: ProcessConfig,
    *,
    include_fragment_summary: bool = False,
    include_fragment_realization: bool = False,
) -> ColumnarReadBatch:
    """Apply the common path without constructing unpublished FASTQ bytes."""

    if not supports_common_processing(config):
        raise ProcessError(
            "configuration has no exact common-column process path"
        )
    if not isinstance(include_fragment_summary, bool):
        raise ProcessError("include_fragment_summary must be a boolean")
    if not isinstance(include_fragment_realization, bool):
        raise ProcessError("include_fragment_realization must be a boolean")
    if include_fragment_realization and not include_fragment_summary:
        raise ProcessError("fragment realization requires fragment summary")
    result = _generate_columnar_read_batch(
        batch,
        config,
        result_mode=_ResultMode.COLUMNAR_READS,
        include_fragment_summary=include_fragment_summary,
        include_fragment_realization=include_fragment_realization,
    )
    if not isinstance(result, ColumnarReadBatch):
        raise ProcessError("common read processing returned an invalid batch")
    return result


def _pack_fragment_summaries(
    batch: ColumnarFragmentBatch,
    states: np.ndarray,
    fragment_modes: np.ndarray,
    mate_fragment: np.ndarray,
    read_summary: np.ndarray,
    config: ProcessConfig,
    *,
    include_fragment_realization: bool = False,
) -> tuple[bytes, tuple[bytes, ...] | None]:
    """Build fixed-width zf values without materializing typed fragments."""

    model = batch.model
    fragment_count = batch.fragment_count
    summary = np.zeros((fragment_count, 12), dtype=np.uint64)
    haplotypes = batch.array(batch.haplotypes, _U1)
    capture_strands = batch.array(batch.capture_strands, _U1)
    summary[:, 0] = (
        haplotypes.astype(np.uint64)
        | (capture_strands.astype(np.uint64) << np.uint64(2))
        | (fragment_modes.astype(np.uint64) << np.uint64(4))
    )

    site_offsets = batch.array(model.site_offsets, _U8)
    site_owner = _owners_from_offsets(site_offsets)
    raw_contexts = batch.array(model.site_contexts, _U1)
    if np.any(raw_contexts >= len(_CONTEXT_STATE)):
        raise ProcessError("NumPy methylation context is invalid")
    contexts = _CONTEXT_STATE[raw_contexts]
    if np.any(contexts == 0):
        raise ProcessError("NumPy methylation context is unsupported")
    methylated = states == 1
    for context_code, total_column, methylated_column in (
        (1, 1, 2),
        (2, 3, 4),
        (3, 5, 6),
    ):
        selected = contexts == context_code
        summary[:, total_column] = np.bincount(
            site_owner[selected], minlength=fragment_count
        )
        summary[:, methylated_column] = np.bincount(
            site_owner[np.logical_and(selected, methylated)],
            minlength=fragment_count,
        )

    has_asm = np.zeros(fragment_count, dtype=np.bool_)
    np.logical_or.at(
        has_asm,
        site_owner,
        batch.array(model.methylation_sources, _U1) == 2,
    )
    summary[:, 0] |= has_asm.astype(np.uint64) << np.uint64(7)

    if batch.variant_offsets:
        variant_offsets = batch.array(batch.variant_offsets, _U4)
        event_owner = _owners_from_offsets(variant_offsets)
        variant_kinds = batch.array(batch.variant_kinds, _U1)
        if np.any(np.logical_or(variant_kinds < 1, variant_kinds > 3)):
            raise ProcessError("NumPy variant kind is invalid")
        kind_bits = np.asarray((0, 1 << 8, 1 << 9, 1 << 10), dtype=np.uint64)
        fragment_event_flags = np.zeros(fragment_count, dtype=np.uint64)
        np.bitwise_or.at(
            fragment_event_flags,
            event_owner,
            kind_bits[variant_kinds],
        )
        summary[:, 0] |= fragment_event_flags
        summary[:, 9] = np.diff(variant_offsets)

    template_offsets = batch.array(model.template_offsets, _U8)
    template_bases = batch.array(model.template_bases, _U1)
    template_owner = _owners_from_offsets(template_offsets)
    local_template_offsets = (
        np.arange(template_bases.size, dtype=np.uint64)
        - template_offsets[template_owner]
    )
    template_methylation = np.full(template_bases.size, -1, dtype=np.int8)
    site_template_offsets = batch.array(model.site_template_offsets, _U4)
    template_methylation[
        template_offsets[site_owner] + site_template_offsets
    ] = states
    owner_modes = fragment_modes[template_owner]
    conversion_targets = np.logical_or(
        np.logical_and(owner_modes == 0, template_bases == 1),
        np.logical_and(owner_modes == 1, template_bases == 2),
    )
    attempts = np.flatnonzero(
        np.logical_and(conversion_targets, template_methylation != 1)
    )
    attempt_owner = template_owner[attempts]
    attempt_counts = np.bincount(attempt_owner, minlength=fragment_count)
    conversion_success = np.zeros(template_bases.size, dtype=np.bool_)
    if attempts.size:
        converted = _draw_fragment_base_bernoulli(
            config.master_seed,
            RNGStage.CONVERSION,
            batch.array(model.contig_indices, _U4),
            attempt_owner,
            batch.array(model.fragment_ordinal_bytes, _U8)[attempt_owner],
            local_template_offsets[attempts],
            config.conversion_rate,
        )
        converted_counts = np.bincount(
            attempt_owner[converted], minlength=fragment_count
        )
        conversion_success[attempts[converted]] = True
    else:
        converted_counts = np.zeros(fragment_count, dtype=np.int64)
    summary[:, 7] = converted_counts
    summary[:, 8] = attempt_counts - converted_counts
    summary[:, 0] |= (
        (converted_counts > 0).astype(np.uint64) << np.uint64(11)
    )

    np.add.at(summary[:, 10], mate_fragment, read_summary[:, 10])
    np.add.at(summary[:, 11], mate_fragment, read_summary[:, 11])
    overflow = np.any(summary[:, 1:] > np.uint64(0xFFFF), axis=1)
    summary[:, 0] |= overflow.astype(np.uint64) << np.uint64(13)
    summary[:, 1:] = np.minimum(summary[:, 1:], np.uint64(0xFFFF))
    if np.any(summary[:, 0] > np.uint64(0xFFFF)):
        raise ProcessError("NumPy fragment summary flags exceed uint16")
    realizations = None
    if include_fragment_realization:
        values = []
        for fragment_index in range(fragment_count):
            site_begin = int(site_offsets[fragment_index])
            site_end = int(site_offsets[fragment_index + 1])
            template_begin = int(template_offsets[fragment_index])
            template_end = int(template_offsets[fragment_index + 1])
            target_mask = conversion_targets[template_begin:template_end]
            values.append(
                _encode_fragment_realization(
                    states[site_begin:site_end],
                    conversion_success[template_begin:template_end][target_mask],
                )
            )
        realizations = tuple(values)
    return np.ascontiguousarray(summary.astype(_U2)).tobytes(), realizations


def _encode_bits_lsb0(values) -> bytes:
    bits = np.asarray(values, dtype=np.uint8)
    if bits.size == 0:
        return b""
    packed = np.packbits(bits, bitorder="little")
    return base64.urlsafe_b64encode(packed.tobytes()).rstrip(b"=")


def _encode_fragment_realization(states, conversions) -> bytes:
    return b".".join(
        (
            format(len(states), "x").encode("ascii"),
            format(len(conversions), "x").encode("ascii"),
            _encode_bits_lsb0(states),
            _encode_bits_lsb0(conversions),
        )
    )


def _encode_typed_fragment_realization(
    sampled_methylation: tuple[bool, ...],
    converted: _ConvertedFragment,
) -> bytes:
    conversion_values = tuple(
        bool(succeeded)
        for base, methylated, attempted, succeeded in zip(
            converted.bases,
            converted.methylated,
            converted.attempted,
            converted.succeeded,
            strict=True,
        )
        if attempted
        or (
            methylated is True
            and (
                (
                    int(converted.conversion_mode) == 0
                    and base in (1, ord("C"))
                )
                or (
                    int(converted.conversion_mode) == 1
                    and base in (2, ord("G"))
                )
            )
        )
    )
    return _encode_fragment_realization(sampled_methylation, conversion_values)


def encode_fastq_batch(
    batch: ColumnarFragmentBatch,
    config: ProcessConfig,
    paired_end: bool,
    *,
    retain_sequences: bool = False,
) -> EncodedFastqBatch:
    """Apply the exact common path and format one consecutive FASTQ batch."""

    if not supports_common_processing(config):
        raise ProcessError(
            "configuration has no exact common-column process path"
        )
    if not isinstance(paired_end, bool):
        raise ProcessError("paired_end must be a boolean")
    if not isinstance(retain_sequences, bool):
        raise ProcessError("retain_sequences must be a boolean")
    result = _generate_columnar_read_batch(
        batch,
        config,
        result_mode=_ResultMode.FORMATTED_FASTQ,
        paired_end=paired_end,
        retain_sequences=retain_sequences,
    )
    if not isinstance(result, EncodedFastqBatch):
        raise ProcessError("C-extension FASTQ formatter returned an invalid batch")
    if paired_end != (result.read2 is not None):
        raise ProcessError("C-extension FASTQ formatter changed mate cardinality")
    return result


def _generate_columnar_read_batch(
    batch: ColumnarFragmentBatch,
    config: ProcessConfig,
    *,
    result_mode: _ResultMode,
    paired_end: bool | None = None,
    retain_sequences: bool = False,
    include_fragment_summary: bool = False,
    include_fragment_realization: bool = False,
) -> (
    tuple[FastqFragment, ...]
    | EncodedFastqBatch
    | ColumnarReadBatch
):
    if not isinstance(batch, ColumnarFragmentBatch):
        raise ProcessError("batch must be ColumnarFragmentBatch")
    if not supports_common_processing(config):
        raise ProcessError("configuration has no exact NumPy process path")
    if not isinstance(result_mode, _ResultMode):
        raise ProcessError("result_mode must be a _ResultMode")
    if not isinstance(retain_sequences, bool):
        raise ProcessError("retain_sequences must be a boolean")
    if not isinstance(include_fragment_summary, bool):
        raise ProcessError("include_fragment_summary must be a boolean")
    if not isinstance(include_fragment_realization, bool):
        raise ProcessError("include_fragment_realization must be a boolean")
    if include_fragment_realization and not include_fragment_summary:
        raise ProcessError("fragment realization requires fragment summary")
    if include_fragment_summary and result_mode is not _ResultMode.COLUMNAR_READS:
        raise ProcessError("fragment summaries require columnar details output")

    format_fastq = result_mode is _ResultMode.FORMATTED_FASTQ
    columnar_result = result_mode is _ResultMode.COLUMNAR_READS
    fastq_ready = not columnar_result
    if format_fastq and not isinstance(paired_end, bool):
        raise ProcessError("formatted FASTQ paired_end must be a boolean")
    if columnar_result:
        retain_sequences = True

    states = (
        _sample_site_states(batch, config)
        if config.bisulfite
        else np.zeros(0, dtype=np.uint8)
    )
    model = batch.model
    ordinals = batch.array(model.fragment_ordinal_bytes, _U8)
    contig_indices = batch.array(model.contig_indices, _U4)
    template_offsets = batch.array(model.template_offsets, _U8)
    template_bases = batch.array(model.template_bases, _U1)
    site_offsets = batch.array(model.site_offsets, _U8)
    haplotypes = batch.array(batch.haplotypes, _U1)
    capture_strands = batch.array(batch.capture_strands, _U1)
    mate_offsets = batch.array(batch.mate_offsets, _U8)
    mate_indices = batch.array(batch.mate_indices, _U1)
    reverse = batch.array(batch.mate_reverse_complements, _U1).astype(
        np.bool_, copy=False
    )
    template_starts = batch.array(batch.mate_template_starts, _U4)
    template_ends = batch.array(batch.mate_template_ends, _U4)
    reference_starts = batch.array(batch.mate_reference_starts, _U8)
    reference_ends = batch.array(batch.mate_reference_ends, _U8)
    site_ref_offsets = batch.array(batch.site_ref_offsets, _U8)
    site_ref_read_offsets = batch.array(batch.site_ref_read_offsets, _U4)
    site_ref_site_indices = batch.array(batch.site_ref_site_indices, _U4)

    mate_fragment = _owners_from_offsets(mate_offsets)
    read_lengths = template_ends.astype(np.int64) - template_starts.astype(np.int64)
    if read_lengths.size == 0 or np.any(read_lengths <= 0):
        raise ProcessError("NumPy mate slices must be non-empty")
    read_length = int(read_lengths[0])
    if np.any(read_lengths != read_length):
        raise ProcessError("NumPy path requires one configured read length")

    cycle = np.arange(read_length, dtype=np.uint64)
    starts = template_offsets[mate_fragment] + template_starts.astype(np.uint64)
    gather = starts[:, None] + cycle[None, :]
    if gather.size and int(gather.max()) >= len(template_bases):
        raise ProcessError("NumPy mate slice exceeds template bases")
    oriented = template_bases[gather].copy()
    if np.any(oriented > 4):
        raise ProcessError("NumPy oriented bases contain an invalid base code")
    if reverse.any():
        oriented[reverse] = _COMPLEMENT[oriented[reverse, ::-1]]

    fragment_modes = _resolve_fragment_modes(
        capture_strands,
        ordinals,
        contig_indices,
        config,
    )
    mate_modes = (
        np.bitwise_xor(fragment_modes[mate_fragment], reverse.astype(np.uint8))
        if config.bisulfite
        else np.full(len(mate_fragment), int(ConversionMode.NONE), dtype=np.uint8)
    )

    methylated = np.full(oriented.shape, -1, dtype=np.int8)
    context_codes = np.zeros(oriented.shape, dtype=np.uint8)
    reference_mate = _owners_from_offsets(site_ref_offsets)
    if config.bisulfite and reference_mate.size:
        if np.any(site_ref_read_offsets >= read_length):
            raise ProcessError("NumPy site reference exceeds read length")
        reference_fragment = mate_fragment[reference_mate]
        local_site_counts = np.diff(site_offsets)[reference_fragment]
        if np.any(site_ref_site_indices >= local_site_counts):
            raise ProcessError("NumPy site reference exceeds fragment sites")
        global_site = site_offsets[reference_fragment] + site_ref_site_indices
        methylated[reference_mate, site_ref_read_offsets] = states[global_site]
        raw_contexts = batch.array(model.site_contexts, _U1)[global_site]
        if np.any(raw_contexts >= len(_CONTEXT_STATE)):
            raise ProcessError("NumPy methylation context is invalid")
        context_codes[reference_mate, site_ref_read_offsets] = _CONTEXT_STATE[
            raw_contexts
        ]

    targets = np.logical_or(
        np.logical_and(mate_modes[:, None] == 0, oriented == 1),
        np.logical_and(mate_modes[:, None] == 1, oriented == 2),
    )
    implicit = np.logical_and(targets, methylated < 0)
    methylated[implicit] = 0
    attempted = np.logical_and(targets, methylated != 1)
    converted = oriented.copy()
    succeeded = np.zeros(oriented.shape, dtype=np.bool_)
    attempt_rows, attempt_cycles = np.nonzero(attempted)
    if attempt_rows.size:
        attempt_template_offsets = _mate_template_offsets(
            template_starts,
            template_ends,
            reverse,
            attempt_rows,
            attempt_cycles,
        )
        success = _draw_fragment_base_bernoulli(
            config.master_seed,
            RNGStage.CONVERSION,
            contig_indices,
            mate_fragment[attempt_rows],
            ordinals[mate_fragment[attempt_rows]],
            attempt_template_offsets,
            config.conversion_rate,
        )
        success_rows = attempt_rows[success]
        success_cycles = attempt_cycles[success]
        succeeded[success_rows, success_cycles] = True
        converted[success_rows, success_cycles] = np.where(
            mate_modes[success_rows] == 0,
            np.uint8(3),
            np.uint8(0),
        )

    final_bases = converted.copy()
    error_flags = np.zeros(oriented.shape, dtype=np.bool_)
    eligible_rows, eligible_cycles = np.nonzero(converted != 4)
    if eligible_rows.size:
        errors, alternatives = _draw_uniform_errors(
            config.master_seed,
            contig_indices,
            mate_fragment[eligible_rows],
            ordinals[mate_fragment[eligible_rows]],
            mate_indices[eligible_rows],
            eligible_cycles,
            config.error.rate,
        )
        error_rows = eligible_rows[errors]
        error_cycles = eligible_cycles[errors]
        if error_rows.size:
            error_flags[error_rows, error_cycles] = True
            original = converted[error_rows, error_cycles]
            final_bases[error_rows, error_cycles] = _ALTERNATIVE_BASES[
                original,
                alternatives[errors],
            ]

    sequence_ascii = _BASE_ASCII[final_bases]
    base_state_codes = None
    read_summaries = None
    fragment_summaries = None
    fragment_realizations = None
    if retain_sequences:
        variant_bases, variant_counts, variant_flags = _variant_annotations_for_mates(
            batch,
            mate_offsets,
            template_starts,
            template_ends,
            reverse,
            read_length,
        )
        state_values = (
            context_codes
            | ((methylated == 1).astype(np.uint8) << np.uint8(2))
            | (succeeded.astype(np.uint8) << np.uint8(3))
            | (variant_bases.astype(np.uint8) << np.uint8(4))
            | (error_flags.astype(np.uint8) << np.uint8(5))
        )
        base_state_codes = _SITE_STATE_ASCII[state_values].tobytes()

        summary = np.zeros((batch.mate_count, 12), dtype=np.uint64)
        summary[:, 0] = (
            haplotypes[mate_fragment].astype(np.uint64)
            | (capture_strands[mate_fragment].astype(np.uint64) << np.uint64(2))
            | (mate_modes.astype(np.uint64) << np.uint64(4))
            | variant_flags.astype(np.uint64)
        )
        if reference_mate.size:
            methylation_sources = batch.array(model.methylation_sources, _U1)[global_site]
            asm = np.zeros(batch.mate_count, dtype=np.bool_)
            np.logical_or.at(asm, reference_mate, methylation_sources == 2)
            summary[:, 0] |= asm.astype(np.uint64) << np.uint64(7)
        summary[:, 0] |= (
            np.any(succeeded, axis=1).astype(np.uint64) << np.uint64(11)
        )
        for context_code, total_column, methylated_column in (
            (1, 1, 2),
            (2, 3, 4),
            (3, 5, 6),
        ):
            selected = context_codes == context_code
            summary[:, total_column] = np.count_nonzero(selected, axis=1)
            summary[:, methylated_column] = np.count_nonzero(
                np.logical_and(selected, methylated == 1), axis=1
            )
        summary[:, 7] = np.count_nonzero(succeeded, axis=1)
        summary[:, 8] = np.count_nonzero(
            np.logical_and(attempted, np.logical_not(succeeded)), axis=1
        )
        summary[:, 9] = variant_counts
        summary[:, 10] = np.count_nonzero(error_flags, axis=1)
        false_methylation = np.logical_and(
            error_flags,
            np.logical_or(
                np.logical_and(
                    mate_modes[:, None] == 0,
                    np.logical_and(converted == 3, final_bases == 1),
                ),
                np.logical_and(
                    mate_modes[:, None] == 1,
                    np.logical_and(converted == 0, final_bases == 2),
                ),
            ),
        )
        summary[:, 11] = np.count_nonzero(false_methylation, axis=1)
        if include_fragment_summary:
            fragment_summaries, fragment_realizations = _pack_fragment_summaries(
                batch,
                states,
                fragment_modes,
                mate_fragment,
                summary,
                config,
                include_fragment_realization=include_fragment_realization,
            )
        overflow = np.any(summary[:, 1:] > np.uint64(0xFFFF), axis=1)
        summary[:, 0] |= overflow.astype(np.uint64) << np.uint64(13)
        summary[:, 1:] = np.minimum(summary[:, 1:], np.uint64(0xFFFF))
        read_summaries = np.ascontiguousarray(summary.astype(_U2)).tobytes()
    quality_byte = config.quality.phred + 33
    sequence_buffer = np.ascontiguousarray(sequence_ascii)
    if columnar_result:
        if base_state_codes is None or read_summaries is None:
            raise ProcessError("columnar read details state is incomplete")
        return ColumnarReadBatch(
            batch.fragment_count,
            batch.mate_count,
            sequence_buffer.tobytes(),
            read_length,
            quality_byte,
            base_state_codes,
            read_summaries,
            fragment_summaries,
            fragment_realizations,
        )
    if format_fastq:
        try:
            read1, read2, record_lengths = _cext_format_fastq_batch(
                model.contig_names,
                model.fragment_ordinal_bytes,
                batch.mate_offsets,
                batch.mate_indices,
                batch.mate_reference_starts,
                batch.mate_reference_ends,
                sequence_buffer,
                read_length,
                quality_byte,
                int(bool(paired_end)),
            )
        except (BufferError, OverflowError, TypeError, ValueError) as error:
            raise ProcessError(str(error)) from error
        return EncodedFastqBatch(
            read1,
            read2,
            record_lengths,
            sequence_buffer.tobytes() if retain_sequences else None,
            read_length if retain_sequences else 0,
            quality_byte if retain_sequences else 0,
            base_state_codes,
            read_summaries,
        )
    quality_text = (bytes((quality_byte,)) * read_length).decode("ascii")
    processed = []
    for fragment_index in range(batch.fragment_count):
        mates = []
        first_mate = int(mate_offsets[fragment_index])
        final_mate = int(mate_offsets[fragment_index + 1])
        mate_order = sorted(
            range(first_mate, final_mate),
            key=lambda value: int(mate_indices[value]),
        )
        for mate_index in mate_order:
            sequence = sequence_ascii[mate_index].tobytes().decode("ascii")
            if fastq_ready:
                mates.append(
                    FastqMate(
                        int(mate_indices[mate_index]),
                        sequence,
                        quality_text,
                    )
                )
            else:
                mates.append(
                    ProcessedMate(
                        int(mate_indices[mate_index]),
                        bool(reverse[mate_index]),
                        ConversionMode(int(mate_modes[mate_index])),
                        int(reference_starts[mate_index]),
                        int(reference_ends[mate_index]),
                        sequence,
                        quality_text,
                        (),
                    )
                )
        if fastq_ready:
            processed.append(
                FastqFragment(
                    int(ordinals[fragment_index]),
                    model.contig_names[fragment_index],
                    int(reference_starts[first_mate]),
                    int(reference_ends[first_mate]),
                    tuple(mates),
                )
            )
        else:
            processed.append(
                ProcessedFragment(
                    int(ordinals[fragment_index]),
                    model.contig_names[fragment_index],
                    int(reference_starts[first_mate]),
                    int(reference_ends[first_mate]),
                    int(haplotypes[fragment_index]),
                    ConversionMode(int(fragment_modes[fragment_index])),
                    (),
                    (),
                    tuple(mates),
                    CaptureStrand(int(capture_strands[fragment_index])),
                )
            )
    return tuple(processed)


def _variant_annotations_for_mates(
    batch: ColumnarFragmentBatch,
    mate_offsets: np.ndarray,
    template_starts: np.ndarray,
    template_ends: np.ndarray,
    reverse: np.ndarray,
    read_length: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Project selected variant events onto read bases and read summaries."""

    variant_bases = np.zeros((batch.mate_count, read_length), dtype=np.bool_)
    variant_counts = np.zeros(batch.mate_count, dtype=np.uint64)
    variant_flags = np.zeros(batch.mate_count, dtype=np.uint16)
    if not batch.variant_offsets:
        return variant_bases, variant_counts, variant_flags
    variant_offsets = batch.array(batch.variant_offsets, _U4)
    event_begins = batch.array(batch.variant_template_starts, _U4)
    event_ends = batch.array(batch.variant_template_ends, _U4)
    variant_kinds = batch.array(batch.variant_kinds, _U1)
    kind_bits = {1: 1 << 8, 2: 1 << 9, 3: 1 << 10}
    for fragment_index in range(batch.fragment_count):
        for variant_index in range(
            int(variant_offsets[fragment_index]),
            int(variant_offsets[fragment_index + 1]),
        ):
            event_begin = int(event_begins[variant_index])
            event_end = int(event_ends[variant_index])
            kind = int(variant_kinds[variant_index])
            if kind not in kind_bits:
                raise ProcessError("NumPy variant kind is invalid")
            for mate_row in range(
                int(mate_offsets[fragment_index]),
                int(mate_offsets[fragment_index + 1]),
            ):
                mate_begin = int(template_starts[mate_row])
                mate_end = int(template_ends[mate_row])
                if event_begin == event_end:
                    observed = mate_begin < event_begin < mate_end
                else:
                    observed = max(event_begin, mate_begin) < min(event_end, mate_end)
                if not observed:
                    continue
                variant_counts[mate_row] += np.uint64(1)
                variant_flags[mate_row] |= np.uint16(kind_bits[kind])
                if event_begin == event_end:
                    continue
                overlap_begin = max(event_begin, mate_begin)
                overlap_end = min(event_end, mate_end)
                offsets = np.arange(overlap_begin, overlap_end, dtype=np.int64)
                cycles = (
                    mate_end - 1 - offsets
                    if reverse[mate_row]
                    else offsets - mate_begin
                )
                variant_bases[mate_row, cycles] = True
    return variant_bases, variant_counts, variant_flags


__all__ = [
    "decode_common_numpy_batch",
    "decode_fragments",
    "encode_fastq_batch",
    "generate_columnar_reads",
    "process_fragment_batch",
    "supports_common_processing",
]
