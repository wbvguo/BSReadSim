"""NumPy-owned batch post-processing over packed fragment columns.

The native extension stops at the representation boundary: it validates one
CRC-verified protocol payload block and decodes it directly into immutable
little-endian buffers.  Model inference and all subsequent orientation,
conversion, quality, and sequencing-error decisions remain Python-owned and
are expressed as NumPy array operations rather than per-base Python loops.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np

from .postprocess import (
    ConversionMode,
    MethylationModelBatch,
    PostprocessConfig,
    PostprocessError,
    ProcessedFragment,
    ProcessedMate,
    UniformError,
    UniformQuality,
)
from .read_names import ReadNameError, format_fragment_identifier
from .rng import RNGStage, derive_key

try:
    from ._native import philox_pairs as _native_philox_pairs
except ImportError:
    _native_philox_pairs = None

try:
    from ._native import format_fastq_batch as _native_format_fastq_batch
except ImportError:
    _native_format_fastq_batch = None


_U1 = np.dtype("u1")
_U4 = np.dtype("<u4")
_U8 = np.dtype("<u8")
_I8 = np.dtype("<i8")
_F4 = np.dtype("<f4")
_UINT64 = np.dtype(np.uint64)
_MASK32 = np.uint64(0xFFFFFFFF)
_PHILOX_M0 = np.uint64(0xD2511F53)
_PHILOX_M1 = np.uint64(0xCD9E8D57)
_PHILOX_W0 = 0x9E3779B9
_PHILOX_W1 = 0xBB67AE85
_TWO_NEGATIVE_53 = 2.0 ** -53
_THREE_WAY_THRESHOLD_1 = np.uint64(6148914691236517206)
_THREE_WAY_THRESHOLD_2 = np.uint64(12297829382473034411)
_COMPLEMENT = np.asarray((3, 2, 1, 0, 4), dtype=np.uint8)
_BASE_ASCII = np.frombuffer(b"ACGTN", dtype=np.uint8)
_ALTERNATIVE_BASES = np.asarray(
    ((1, 2, 3), (0, 2, 3), (0, 1, 3), (0, 1, 2)),
    dtype=np.uint8,
)


@dataclass(frozen=True)
class NumpyFragmentBatch:
    """Complete FASTQ-only fragment batch backed by immutable byte columns."""

    model: MethylationModelBatch
    haplotypes: bytes
    capture_strands: bytes
    mate_offsets: bytes
    mate_indices: bytes
    mate_reverse_complements: bytes
    mate_template_starts: bytes
    mate_template_ends: bytes
    mate_reference_starts: bytes
    mate_reference_ends: bytes
    site_ref_offsets: bytes
    site_ref_read_offsets: bytes
    site_ref_site_indices: bytes

    def __post_init__(self) -> None:
        values = (
            self.haplotypes,
            self.capture_strands,
            self.mate_offsets,
            self.mate_indices,
            self.mate_reverse_complements,
            self.mate_template_starts,
            self.mate_template_ends,
            self.mate_reference_starts,
            self.mate_reference_ends,
            self.site_ref_offsets,
            self.site_ref_read_offsets,
            self.site_ref_site_indices,
        )
        if any(type(value) is not bytes for value in values):
            raise PostprocessError("NumPy fragment batch buffers must be bytes")
        fragment_count = self.model.fragment_count
        if (
            len(self.haplotypes) != fragment_count
            or len(self.capture_strands) != fragment_count
            or len(self.mate_offsets) != (fragment_count + 1) * 8
        ):
            raise PostprocessError("NumPy fragment columns disagree on fragment count")
        mate_offsets = self.array(self.mate_offsets, _U8)
        _validate_offsets(mate_offsets, "mate")
        mate_count = int(mate_offsets[-1])
        if (
            len(self.mate_indices) != mate_count
            or len(self.mate_reverse_complements) != mate_count
            or len(self.mate_template_starts) != mate_count * 4
            or len(self.mate_template_ends) != mate_count * 4
            or len(self.mate_reference_starts) != mate_count * 8
            or len(self.mate_reference_ends) != mate_count * 8
            or len(self.site_ref_offsets) != (mate_count + 1) * 8
        ):
            raise PostprocessError("NumPy fragment columns disagree on mate count")
        site_ref_offsets = self.array(self.site_ref_offsets, _U8)
        _validate_offsets(site_ref_offsets, "site-reference")
        site_ref_count = int(site_ref_offsets[-1])
        if (
            len(self.site_ref_read_offsets) != site_ref_count * 4
            or len(self.site_ref_site_indices) != site_ref_count * 4
        ):
            raise PostprocessError(
                "NumPy fragment columns disagree on site-reference count"
            )

    @staticmethod
    def array(buffer: bytes, dtype: np.dtype) -> np.ndarray:
        """Return a read-only zero-copy one-dimensional NumPy view."""

        result = np.frombuffer(buffer, dtype=dtype)
        result.flags.writeable = False
        return result

    @property
    def fragment_count(self) -> int:
        return self.model.fragment_count

    @property
    def mate_count(self) -> int:
        return len(self.mate_indices)


@dataclass(frozen=True)
class NumpyFastqMate:
    """FASTQ-ready mate without unavailable no-Truth coordinate claims."""

    mate_index: int
    sequence: str
    quality: str


@dataclass(frozen=True)
class NumpyFastqFragment:
    """FASTQ-ready fragment produced by the common-column lane."""

    fragment_ordinal: int
    contig_name: str
    reference_start: int
    reference_end: int
    mates: Tuple[NumpyFastqMate, ...]


@dataclass(frozen=True)
class NumpyFormattedFastqBatch:
    """Consecutive FASTQ records formatted without per-mate Python objects."""

    read1: bytes
    read2: Optional[bytes]
    record_lengths: Tuple[Tuple[int, int, int], ...]

    def __post_init__(self) -> None:
        if type(self.read1) is not bytes or (
            self.read2 is not None and type(self.read2) is not bytes
        ):
            raise PostprocessError("formatted FASTQ batches must contain bytes")
        if not isinstance(self.record_lengths, tuple) or not self.record_lengths:
            raise PostprocessError("formatted FASTQ record lengths must be non-empty")
        for lengths in self.record_lengths:
            if (
                not isinstance(lengths, tuple)
                or len(lengths) != 3
                or any(
                    isinstance(value, bool)
                    or not isinstance(value, int)
                    or value < 0
                    for value in lengths
                )
                or lengths[0] == 0
                or lengths[2] != 0
            ):
                raise PostprocessError("formatted FASTQ record lengths are invalid")
        if sum(value[0] for value in self.record_lengths) != len(self.read1):
            raise PostprocessError("formatted read1 bytes disagree with record lengths")
        expected_read2 = sum(value[1] for value in self.record_lengths)
        if (self.read2 is None) != (expected_read2 == 0) or (
            self.read2 is not None and len(self.read2) != expected_read2
        ):
            raise PostprocessError("formatted read2 bytes disagree with record lengths")


def _supports_numpy_policies(
    config: PostprocessConfig,
    *,
    include_truth: bool,
) -> bool:
    return (
        not include_truth
        and isinstance(config.quality, UniformQuality)
        and isinstance(config.error, UniformError)
    )


def supports_common_postprocess(
    config: PostprocessConfig,
    *,
    include_truth: bool,
) -> bool:
    """Whether common columns alone preserve the selected semantics."""

    return _supports_numpy_policies(config, include_truth=include_truth)


def process_common_fastq_batch(
    batch: NumpyFragmentBatch,
    config: PostprocessConfig,
) -> Tuple[NumpyFastqFragment, ...]:
    """Apply the exact common path without materializing false provenance."""

    if not supports_common_postprocess(config, include_truth=False):
        raise PostprocessError(
            "configuration has no exact common-column post-process path"
        )
    return _process_numpy_fragment_batch(batch, config, fastq_ready=True)


def format_common_fastq_batch(
    batch: NumpyFragmentBatch,
    config: PostprocessConfig,
    paired_end: bool,
) -> NumpyFormattedFastqBatch:
    """Apply the exact common path and format one consecutive FASTQ batch."""

    if not supports_common_postprocess(config, include_truth=False):
        raise PostprocessError(
            "configuration has no exact common-column post-process path"
        )
    if not isinstance(paired_end, bool):
        raise PostprocessError("paired_end must be a boolean")
    if _native_format_fastq_batch is None:
        fragments = _process_numpy_fragment_batch(
            batch,
            config,
            fastq_ready=True,
        )
        read1_parts = []
        read2_parts = []
        record_lengths = []
        for fragment in fragments:
            records = {}
            for mate in fragment.mates:
                try:
                    identifier = format_fragment_identifier(
                        fragment.contig_name,
                        fragment.reference_start,
                        fragment.reference_end,
                        fragment.fragment_ordinal,
                        pair_number=mate.mate_index + 1,
                    )
                except ReadNameError as error:
                    raise PostprocessError(str(error)) from error
                record = "@{}\n{}\n+\n{}\n".format(
                    identifier,
                    mate.sequence,
                    mate.quality,
                ).encode("utf-8")
                records[mate.mate_index] = record
            expected_indices = {0, 1} if paired_end else {0}
            if set(records) != expected_indices:
                raise PostprocessError("FASTQ formatter changed mate cardinality")
            read1_parts.append(records[0])
            if paired_end:
                read2_parts.append(records[1])
            record_lengths.append(
                (len(records[0]), len(records.get(1, b"")), 0)
            )
        return NumpyFormattedFastqBatch(
            b"".join(read1_parts),
            b"".join(read2_parts) if paired_end else None,
            tuple(record_lengths),
        )
    result = _process_numpy_fragment_batch(
        batch,
        config,
        fastq_ready=True,
        format_fastq=True,
        paired_end=paired_end,
    )
    if not isinstance(result, NumpyFormattedFastqBatch):
        raise PostprocessError("native FASTQ formatter returned an invalid batch")
    if paired_end != (result.read2 is not None):
        raise PostprocessError("native FASTQ formatter changed mate cardinality")
    return result


def _process_numpy_fragment_batch(
    batch: NumpyFragmentBatch,
    config: PostprocessConfig,
    *,
    fastq_ready: bool,
    format_fastq: bool = False,
    paired_end: Optional[bool] = None,
):

    if not isinstance(batch, NumpyFragmentBatch):
        raise PostprocessError("batch must be NumpyFragmentBatch")
    if not _supports_numpy_policies(config, include_truth=False):
        raise PostprocessError("configuration has no exact NumPy post-process path")
    if not isinstance(format_fastq, bool):
        raise PostprocessError("format_fastq must be a boolean")
    if format_fastq and not isinstance(paired_end, bool):
        raise PostprocessError("formatted FASTQ paired_end must be a boolean")

    states = _sample_site_states(batch, config)
    model = batch.model
    ordinals = batch.array(model.fragment_ordinal_data, _U8)
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
        raise PostprocessError("NumPy mate slices must be non-empty")
    read_length = int(read_lengths[0])
    if np.any(read_lengths != read_length):
        raise PostprocessError("NumPy path requires one configured read length")

    cycle = np.arange(read_length, dtype=np.uint64)
    starts = template_offsets[mate_fragment] + template_starts.astype(np.uint64)
    gather = starts[:, None] + cycle[None, :]
    if gather.size and int(gather.max()) >= len(template_bases):
        raise PostprocessError("NumPy mate slice exceeds template bases")
    oriented = template_bases[gather].copy()
    if np.any(oriented > 4):
        raise PostprocessError("NumPy oriented bases contain an invalid base code")
    if reverse.any():
        oriented[reverse] = _COMPLEMENT[oriented[reverse, ::-1]]

    fragment_modes = _fragment_conversion_modes(
        capture_strands,
        ordinals,
        contig_indices,
        config,
    )
    mate_modes = np.bitwise_xor(fragment_modes[mate_fragment], reverse.astype(np.uint8))

    methylated = np.full(oriented.shape, -1, dtype=np.int8)
    reference_mate = _owners_from_offsets(site_ref_offsets)
    if reference_mate.size:
        if np.any(site_ref_read_offsets >= read_length):
            raise PostprocessError("NumPy site reference exceeds read length")
        reference_fragment = mate_fragment[reference_mate]
        local_site_counts = np.diff(site_offsets)[reference_fragment]
        if np.any(site_ref_site_indices >= local_site_counts):
            raise PostprocessError("NumPy site reference exceeds fragment sites")
        global_site = site_offsets[reference_fragment] + site_ref_site_indices
        methylated[reference_mate, site_ref_read_offsets] = states[global_site]

    targets = np.logical_or(
        np.logical_and(mate_modes[:, None] == 0, oriented == 1),
        np.logical_and(mate_modes[:, None] == 1, oriented == 2),
    )
    implicit = np.logical_and(targets, methylated < 0)
    methylated[implicit] = 0
    attempted = np.logical_and(targets, methylated != 1)
    converted = oriented.copy()
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
        converted[success_rows, success_cycles] = np.where(
            mate_modes[success_rows] == 0,
            np.uint8(3),
            np.uint8(0),
        )

    final_bases = converted.copy()
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
            original = converted[error_rows, error_cycles]
            final_bases[error_rows, error_cycles] = _ALTERNATIVE_BASES[
                original,
                alternatives[errors],
            ]

    sequence_ascii = _BASE_ASCII[final_bases]
    quality_byte = config.quality.phred + 33
    if format_fastq:
        if _native_format_fastq_batch is None:
            raise PostprocessError("native FASTQ batch formatter is unavailable")
        try:
            read1, read2, record_lengths = _native_format_fastq_batch(
                model.contig_names,
                model.fragment_ordinal_data,
                batch.mate_offsets,
                batch.mate_indices,
                batch.mate_reference_starts,
                batch.mate_reference_ends,
                np.ascontiguousarray(sequence_ascii),
                read_length,
                quality_byte,
                int(bool(paired_end)),
            )
        except (BufferError, OverflowError, TypeError, ValueError) as error:
            raise PostprocessError(str(error)) from error
        return NumpyFormattedFastqBatch(read1, read2, record_lengths)
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
                    NumpyFastqMate(
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
                NumpyFastqFragment(
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
                )
            )
    return tuple(processed)


def _sample_site_states(
    batch: NumpyFragmentBatch,
    config: PostprocessConfig,
) -> np.ndarray:
    model = batch.model
    probabilities = batch.array(model.site_probabilities, _F4)
    site_indices = batch.array(model.site_indices, _U4).astype(np.uint64, copy=False)
    site_offsets = batch.array(model.site_offsets, _U8)
    site_fragment = _owners_from_offsets(site_offsets)
    result = np.empty(model.site_count, dtype=np.uint8)
    contig_indices = batch.array(model.contig_indices, _U4)
    fragment_groups = _contig_fragment_groups(contig_indices)
    if len(fragment_groups) == 1:
        contig_index = next(iter(fragment_groups))
        key = derive_key(config.master_seed, RNGStage.SITE_STATE, contig_index)
        pair0, _ = _philox_pairs(
            key,
            batch.array(model.fragment_ordinal_data, _U8)[site_fragment],
            site_indices,
        )
        result[:] = _bernoulli_from_pair(pair0, probabilities)
        result.flags.writeable = False
        return result
    for contig_index, fragments in fragment_groups.items():
        selected = np.isin(site_fragment, fragments, assume_unique=False)
        positions = np.flatnonzero(selected)
        if positions.size == 0:
            continue
        key = derive_key(config.master_seed, RNGStage.SITE_STATE, contig_index)
        pair0, _ = _philox_pairs(
            key,
            batch.array(model.fragment_ordinal_data, _U8)[site_fragment[positions]],
            site_indices[positions],
        )
        result[positions] = _bernoulli_from_pair(
            pair0,
            probabilities[positions],
        )
    result.flags.writeable = False
    return result


def _fragment_conversion_modes(
    capture_strands: np.ndarray,
    ordinals: np.ndarray,
    contig_indices: np.ndarray,
    config: PostprocessConfig,
) -> np.ndarray:
    if np.any(capture_strands > 2):
        raise PostprocessError("NumPy capture strand contains an invalid enum")
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


def _mate_template_offsets(
    template_starts: np.ndarray,
    template_ends: np.ndarray,
    reverse: np.ndarray,
    mate_rows: np.ndarray,
    cycles: np.ndarray,
) -> np.ndarray:
    cycles_u64 = cycles.astype(np.uint64, copy=False)
    starts_u64 = template_starts[mate_rows].astype(np.uint64, copy=False)
    ends_u64 = template_ends[mate_rows].astype(np.uint64, copy=False)
    return np.where(
        reverse[mate_rows],
        ends_u64 - np.uint64(1) - cycles_u64,
        starts_u64 + cycles_u64,
    ).astype(np.uint64, copy=False)


def _draw_fragment_base_bernoulli(
    master_seed: int,
    stage: RNGStage,
    contig_indices: np.ndarray,
    fragment_indices: np.ndarray,
    ordinals: np.ndarray,
    template_offsets: np.ndarray,
    probability: float,
) -> np.ndarray:
    result = np.empty(len(ordinals), dtype=np.bool_)
    if probability == 0.0:
        result.fill(False)
        return result
    if probability == 1.0:
        result.fill(True)
        return result
    local_indices = template_offsets.astype(np.uint64, copy=False)
    groups = _contig_fragment_groups(contig_indices)
    if len(groups) == 1:
        contig_index = next(iter(groups))
        key = derive_key(master_seed, stage, contig_index)
        pair0, _ = _philox_pairs(key, ordinals, local_indices)
        return _bernoulli_from_pair(pair0, probability).astype(
            np.bool_, copy=False
        )
    for contig_index, fragments in groups.items():
        positions = np.flatnonzero(np.isin(fragment_indices, fragments))
        if positions.size == 0:
            continue
        key = derive_key(master_seed, stage, contig_index)
        pair0, _ = _philox_pairs(
            key,
            ordinals[positions],
            local_indices[positions],
        )
        result[positions] = _bernoulli_from_pair(pair0, probability).astype(
            np.bool_, copy=False
        )
    return result


def _draw_uniform_errors(
    master_seed: int,
    contig_indices: np.ndarray,
    fragment_indices: np.ndarray,
    ordinals: np.ndarray,
    mate_indices: np.ndarray,
    cycles: np.ndarray,
    probability: float,
) -> Tuple[np.ndarray, np.ndarray]:
    errors = np.zeros(len(ordinals), dtype=np.bool_)
    alternatives = np.zeros(len(ordinals), dtype=np.uint8)
    if probability == 0.0:
        return errors, alternatives
    local_indices = np.left_shift(
        mate_indices.astype(np.uint64),
        np.uint64(32),
    ) | cycles.astype(np.uint64)
    groups = _contig_fragment_groups(contig_indices)
    if len(groups) == 1:
        contig_index = next(iter(groups))
        key = derive_key(
            master_seed,
            RNGStage.SEQUENCING_ERROR,
            contig_index,
        )
        pair0, pair1 = _philox_pairs(key, ordinals, local_indices)
        errors[:] = _bernoulli_from_pair(pair0, probability).astype(
            np.bool_, copy=False
        )
        alternatives[:] = (
            (pair1 >= _THREE_WAY_THRESHOLD_1).astype(np.uint8)
            + (pair1 >= _THREE_WAY_THRESHOLD_2).astype(np.uint8)
        )
        return errors, alternatives
    for contig_index, fragments in groups.items():
        positions = np.flatnonzero(np.isin(fragment_indices, fragments))
        if positions.size == 0:
            continue
        key = derive_key(
            master_seed,
            RNGStage.SEQUENCING_ERROR,
            contig_index,
        )
        pair0, pair1 = _philox_pairs(
            key,
            ordinals[positions],
            local_indices[positions],
        )
        errors[positions] = _bernoulli_from_pair(pair0, probability).astype(
            np.bool_, copy=False
        )
        alternatives[positions] = (
            (pair1 >= _THREE_WAY_THRESHOLD_1).astype(np.uint8)
            + (pair1 >= _THREE_WAY_THRESHOLD_2).astype(np.uint8)
        )
    return errors, alternatives


def _philox_pairs(
    key: int,
    entity_ordinals: np.ndarray,
    local_indices: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray]:
    entity = np.asarray(entity_ordinals, dtype=_UINT64)
    local = np.asarray(local_indices, dtype=_UINT64)
    if entity.shape != local.shape:
        raise PostprocessError("NumPy Philox counter columns must have one shape")
    if _native_philox_pairs is not None:
        entity = np.ascontiguousarray(entity)
        local = np.ascontiguousarray(local)
        pair_0_data, pair_1_data = _native_philox_pairs(key, entity, local)
        pair_0 = np.frombuffer(pair_0_data, dtype=_UINT64).reshape(entity.shape)
        pair_1 = np.frombuffer(pair_1_data, dtype=_UINT64).reshape(entity.shape)
        pair_0.flags.writeable = False
        pair_1.flags.writeable = False
        return pair_0, pair_1
    counter_0 = entity & _MASK32
    counter_1 = entity >> np.uint64(32)
    counter_2 = local & _MASK32
    counter_3 = local >> np.uint64(32)
    key_0 = key & 0xFFFFFFFF
    key_1 = key >> 32
    for _ in range(10):
        product_0 = _PHILOX_M0 * counter_0
        product_1 = _PHILOX_M1 * counter_2
        counter_0, counter_1, counter_2, counter_3 = (
            ((product_1 >> np.uint64(32)) ^ counter_1 ^ np.uint64(key_0))
            & _MASK32,
            product_1 & _MASK32,
            ((product_0 >> np.uint64(32)) ^ counter_3 ^ np.uint64(key_1))
            & _MASK32,
            product_0 & _MASK32,
        )
        key_0 = (key_0 + _PHILOX_W0) & 0xFFFFFFFF
        key_1 = (key_1 + _PHILOX_W1) & 0xFFFFFFFF
    return (
        counter_0 | (counter_1 << np.uint64(32)),
        counter_2 | (counter_3 << np.uint64(32)),
    )


def _bernoulli_from_pair(pair: np.ndarray, probability) -> np.ndarray:
    probabilities = np.asarray(probability, dtype=np.float64)
    uniform = (pair >> np.uint64(11)).astype(np.float64) * _TWO_NEGATIVE_53
    return np.where(
        probabilities <= 0.0,
        False,
        np.where(probabilities >= 1.0, True, uniform < probabilities),
    ).astype(np.uint8, copy=False)


def _owners_from_offsets(offsets: np.ndarray) -> np.ndarray:
    _validate_offsets(offsets, "ragged")
    return np.repeat(
        np.arange(len(offsets) - 1, dtype=np.intp),
        np.diff(offsets).astype(np.intp),
    )


def _validate_offsets(offsets: np.ndarray, name: str) -> None:
    if offsets.ndim != 1 or offsets.size == 0 or int(offsets[0]) != 0:
        raise PostprocessError("{} offsets must start at zero".format(name))
    if np.any(offsets[1:] < offsets[:-1]):
        raise PostprocessError("{} offsets must be nondecreasing".format(name))


def _contig_fragment_groups(contig_indices: np.ndarray) -> Dict[int, np.ndarray]:
    groups = {}  # type: Dict[int, list]
    for index, contig_index in enumerate(contig_indices):
        groups.setdefault(int(contig_index), []).append(index)
    return {
        contig_index: np.asarray(indices, dtype=np.intp)
        for contig_index, indices in groups.items()
    }


__all__ = [
    "NumpyFastqFragment",
    "NumpyFastqMate",
    "NumpyFormattedFastqBatch",
    "NumpyFragmentBatch",
    "format_common_fastq_batch",
    "process_common_fastq_batch",
    "supports_common_postprocess",
]
