"""Shared immutable records and columnar batches for the processing flow."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, IntEnum
import math
import struct

import numpy as np

from .._cext import philox_pairs as _cext_philox_pairs

NO_VARIANT_INDEX = 0xFFFFFFFF
_COLUMN_U64 = struct.Struct("<Q")

_U1 = np.dtype("u1")
_U2 = np.dtype("<u2")
_U4 = np.dtype("<u4")
_U8 = np.dtype("<u8")
_I8 = np.dtype("<i8")
_F4 = np.dtype("<f4")
_UINT64 = np.dtype(np.uint64)
_TWO_NEGATIVE_53 = 2.0 ** -53

class CaptureStrand(IntEnum):
    UNKNOWN = 0
    FORWARD = 1
    REVERSE = 2


class VariantKind(IntEnum):
    SNV = 1
    INSERTION = 2
    DELETION = 3


class VariantSource(str, Enum):
    VCF = "vcf"
    DE_NOVO = "de_novo"

    # Match enum.StrEnum while retaining Python 3.10 compatibility.
    __str__ = str.__str__


class MethylationContext(IntEnum):
    CG_C = 1
    CHG_C = 3
    CHH_C = 7
    CG_G = 9
    CHG_G = 11
    CHH_G = 15


class MethylationSource(IntEnum):
    CGMAP = 1
    ASM = 2
    BETA = 3
    POOLED_CGMAP = 4


class MethylationAllele(IntEnum):
    SHARED = 0
    REFERENCE_HAPLOTYPE = 1
    ALTERNATE_HAPLOTYPE = 2


@dataclass(frozen=True)
class Variant:
    index: int
    id: str
    source: VariantSource
    kind: VariantKind
    phased_haplotype: int
    reference_start: int
    reference_end: int
    ref_bases: bytes
    alt_bases: bytes


@dataclass(frozen=True)
class MethylationSite:
    site_index: int
    template_offset: int
    reference_pos: int
    context: MethylationContext
    methylation_source: MethylationSource
    allele: MethylationAllele
    methylation_probability: float


@dataclass(frozen=True)
class SiteReference:
    read_offset: int
    site_index: int


@dataclass(frozen=True)
class Mate:
    mate_index: int
    reverse_complement: bool
    template_start: int
    template_end: int
    reference_start: int
    reference_end: int
    site_refs: tuple[SiteReference, ...]


@dataclass(frozen=True)
class Fragment:
    fragment_ordinal: int
    contig_index: int
    haplotype: int
    capture_strand: CaptureStrand
    reference_start: int
    reference_end: int
    template_bases: bytes
    reference_positions: tuple[int, ...]
    base_variant_indices: tuple[int, ...]
    variants: tuple[Variant, ...]
    methylation_sites: tuple[MethylationSite, ...]
    mates: tuple[Mate, ...]


@dataclass(frozen=True)
class FragmentSummary:
    """Fragment aggregates observed independently by Python workers."""

    fragment_count: int
    mate_count: int
    template_base_count: int
    methylation_site_count: int
    per_contig_fragment_counts: tuple[int, ...]


class ProcessError(ValueError):
    """A decoded fragment or baseline configuration cannot be processed."""


@dataclass(frozen=True)
class MethylationModelBatch:
    """Contiguous model inputs for one variable-length fragment block.

    Multi-byte buffers use explicit little-endian values. ``template_offsets``
    and ``site_offsets`` are uint64 prefix sums with ``fragment_count + 1``
    entries, so models can build padding masks or use ragged tensors without
    allocating one Python object per site.  Site-local arrays are flat and
    aligned by their shared site index in the batch.
    """

    fragment_ordinals: tuple[int, ...]
    contig_names: tuple[str, ...]
    fragment_ordinal_bytes: bytes
    contig_indices: bytes
    template_offsets: bytes
    template_bases: bytes
    site_offsets: bytes
    site_indices: bytes
    site_template_offsets: bytes
    site_reference_positions: bytes
    site_contexts: bytes
    methylation_sources: bytes
    site_alleles: bytes
    site_probabilities: bytes

    def __post_init__(self) -> None:
        if not isinstance(self.fragment_ordinals, tuple) or not self.fragment_ordinals:
            raise ProcessError(
                "columnar batch ordinals must be a non-empty immutable tuple"
            )
        if not isinstance(self.contig_names, tuple):
            raise ProcessError(
                "columnar batch contig names must be an immutable tuple"
            )
        fragment_count = len(self.fragment_ordinals)
        if len(self.contig_names) != fragment_count:
            raise ProcessError(
                "columnar batch contig-name count disagrees with ordinals"
            )
        buffers = (
            self.fragment_ordinal_bytes,
            self.contig_indices,
            self.template_offsets,
            self.template_bases,
            self.site_offsets,
            self.site_indices,
            self.site_template_offsets,
            self.site_reference_positions,
            self.site_contexts,
            self.methylation_sources,
            self.site_alleles,
            self.site_probabilities,
        )
        if any(type(value) is not bytes for value in buffers):
            raise ProcessError("columnar batch buffers must be immutable bytes")
        if len(self.fragment_ordinal_bytes) != fragment_count * 8:
            raise ProcessError("columnar fragment ordinal buffer has wrong length")
        if len(self.contig_indices) != fragment_count * 4:
            raise ProcessError("columnar contig-index buffer has wrong length")
        if len(self.template_offsets) != (fragment_count + 1) * 8:
            raise ProcessError("columnar template offsets have wrong length")
        if len(self.site_offsets) != (fragment_count + 1) * 8:
            raise ProcessError("columnar site offsets have wrong length")
        _validate_column_offsets(
            self.template_offsets,
            len(self.template_bases),
            "template",
        )
        site_count = len(self.site_contexts)
        _validate_column_offsets(self.site_offsets, site_count, "site")
        if (
            len(self.site_indices) != site_count * 4
            or len(self.site_template_offsets) != site_count * 4
            or len(self.site_reference_positions) != site_count * 8
            or len(self.methylation_sources) != site_count
            or len(self.site_alleles) != site_count
            or len(self.site_probabilities) != site_count * 4
        ):
            raise ProcessError("columnar site buffers disagree on site count")
        encoded_ordinals = tuple(
            _COLUMN_U64.unpack_from(self.fragment_ordinal_bytes, index * 8)[0]
            for index in range(fragment_count)
        )
        if encoded_ordinals != self.fragment_ordinals:
            raise ProcessError(
                "columnar ordinal buffer disagrees with fragment ordinals"
            )

    @property
    def fragment_count(self) -> int:
        return len(self.fragment_ordinals)

    @property
    def site_count(self) -> int:
        return len(self.site_contexts)


def _validate_column_offsets(data: bytes, final_value: int, name: str) -> None:
    previous = None
    for offset in range(0, len(data), _COLUMN_U64.size):
        value = _COLUMN_U64.unpack_from(data, offset)[0]
        if previous is None:
            if value != 0:
                raise ProcessError(
                    "columnar {} offsets must start at zero".format(name)
                )
        elif value < previous:
            raise ProcessError(
                "columnar {} offsets must be nondecreasing".format(name)
            )
        previous = value
    if previous != final_value:
        raise ProcessError(
            "columnar {} offsets disagree with flat data".format(name)
        )


class ConversionMode(IntEnum):
    """Bisulfite substitution visible in the current oriented mate."""

    C2T = 0
    G2A = 1
    NONE = 2


@dataclass(frozen=True)
class UniformQuality:
    """One constant Phred value for every cycle."""

    phred: int = 40

    def __post_init__(self) -> None:
        if (
            isinstance(self.phred, bool)
            or not isinstance(self.phred, int)
            or not 0 <= self.phred <= 93
        ):
            raise ProcessError("quality phred must be an integer in [0, 93]")


@dataclass(frozen=True)
class UniformError:
    """One constant substitution probability for every A/C/G/T call."""

    rate: float = 0.005

    def __post_init__(self) -> None:
        _require_probability("error rate", self.rate)


@dataclass(frozen=True)
class SiteState:
    """The single latent draw for one fragment-level methylation site."""

    site_index: int
    template_offset: int
    reference_pos: int
    context: MethylationContext
    methylation_source: MethylationSource
    allele: MethylationAllele
    methylated: bool
    probability: float


@dataclass(frozen=True)
class BaseState:
    """Traceable state for one base after every processing stage.

    ``site_index`` is ``None`` when the core did not declare a methylation site.
    For an undeclared base that is nevertheless a target of the mate's
    conversion mode, ``methylated`` is ``False`` to record its implicit
    unmethylated state.  Bases that are neither declared sites nor conversion
    targets retain ``None``.
    """

    read_offset: int
    reference_pos: int
    variant_index: int
    site_index: int | None
    methylated: bool | None
    oriented_base: int
    post_conversion_base: int
    final_base: int
    conversion_attempted: bool
    conversion_succeeded: bool
    sequencing_error: bool
    quality_phred: int


@dataclass(frozen=True)
class _CompactBaseStates:
    """Worker-only columnar base state for C-extension detail formatting.

    Keeping one immutable column per field avoids allocating one 12-field
    Python dataclass for every read base. The public ``process_fragment`` path
    continues to return regular ``BaseState`` tuples by default.
    """

    reference_positions: tuple[int, ...]
    variant_indices: tuple[int, ...]
    site_indices: tuple[int | None, ...]
    methylated: tuple[bool | None, ...]
    oriented_bases: bytes
    post_conversion_bases: bytes
    final_bases: bytes
    attempted: tuple[bool, ...]
    succeeded: tuple[bool, ...]
    error_flags: tuple[bool, ...]
    quality_phreds: tuple[int, ...]

    def __len__(self) -> int:
        return len(self.final_bases)


@dataclass(frozen=True)
class ProcessedMate:
    """A mate ready for a separate annotation/FASTQ formatter."""

    mate_index: int
    reverse_complement: bool
    conversion_mode: ConversionMode
    reference_start: int
    reference_end: int
    sequence: str
    quality: str
    base_states: tuple[BaseState, ...] | _CompactBaseStates


@dataclass(frozen=True)
class ProcessedFragment:
    """Atomic processing result for one SE or PE fragment."""

    fragment_ordinal: int
    contig_name: str
    reference_start: int
    reference_end: int
    haplotype: int
    fragment_conversion_mode: ConversionMode
    variants: tuple[Variant, ...]
    site_states: tuple[SiteState, ...]
    mates: tuple[ProcessedMate, ...]
    capture_strand: CaptureStrand = CaptureStrand.UNKNOWN
    fragment_conversion_success_count: int = 0
    fragment_conversion_failure_count: int = 0
    fragment_realization: bytes | None = None


@dataclass(frozen=True)
class _OrientedMate:
    mate: Mate
    conversion_mode: ConversionMode
    bases: bytes
    reference_positions: tuple[int, ...]
    variant_indices: tuple[int, ...]


@dataclass(frozen=True)
class _ConvertedFragment:
    conversion_mode: ConversionMode
    bases: bytes
    site_indices: tuple[int | None, ...]
    methylated: tuple[bool | None, ...]
    attempted: tuple[bool, ...]
    succeeded: tuple[bool, ...]


@dataclass(frozen=True)
class _ConvertedMate:
    oriented: _OrientedMate
    bases: bytes
    site_indices: tuple[int | None, ...]
    methylated: tuple[bool | None, ...]
    attempted: tuple[bool, ...]
    succeeded: tuple[bool, ...]


@dataclass(frozen=True)
class _QualityMate:
    converted: _ConvertedMate
    qualities: tuple[int, ...]


@dataclass(frozen=True)
class _ErroredMate:
    quality: _QualityMate
    bases: bytes
    error_flags: tuple[bool, ...]


@dataclass(frozen=True)
class ColumnarFragmentBatch:
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
    variant_offsets: bytes = b""
    variant_template_starts: bytes = b""
    variant_template_ends: bytes = b""
    variant_kinds: bytes = b""

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
            raise ProcessError("NumPy fragment batch buffers must be bytes")
        fragment_count = self.model.fragment_count
        if (
            len(self.haplotypes) != fragment_count
            or len(self.capture_strands) != fragment_count
            or len(self.mate_offsets) != (fragment_count + 1) * 8
        ):
            raise ProcessError("NumPy fragment columns disagree on fragment count")
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
            raise ProcessError("NumPy fragment columns disagree on mate count")
        site_ref_offsets = self.array(self.site_ref_offsets, _U8)
        _validate_offsets(site_ref_offsets, "site-reference")
        site_ref_count = int(site_ref_offsets[-1])
        if (
            len(self.site_ref_read_offsets) != site_ref_count * 4
            or len(self.site_ref_site_indices) != site_ref_count * 4
        ):
            raise ProcessError(
                "NumPy fragment columns disagree on site-reference count"
            )
        event_values = (
            self.variant_offsets,
            self.variant_template_starts,
            self.variant_template_ends,
            self.variant_kinds,
        )
        if any(type(value) is not bytes for value in event_values):
            raise ProcessError("NumPy event columns must be bytes")
        if self.variant_offsets:
            if len(self.variant_offsets) != (fragment_count + 1) * 4:
                raise ProcessError("NumPy event offsets disagree on fragment count")
            variant_offsets = self.array(self.variant_offsets, _U4)
            _validate_offsets(variant_offsets, "event")
            variant_count = int(variant_offsets[-1])
            if (
                len(self.variant_template_starts) != variant_count * 4
                or len(self.variant_template_ends) != variant_count * 4
                or len(self.variant_kinds) != variant_count
            ):
                raise ProcessError("NumPy event columns disagree on event count")
        elif any(event_values[1:]):
            raise ProcessError("NumPy event columns require event offsets")

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
class FastqMate:
    """FASTQ-ready mate without unavailable no-Details coordinate claims."""

    mate_index: int
    sequence: str
    quality: str


@dataclass(frozen=True)
class FastqFragment:
    """FASTQ-ready fragment produced by the common-column lane."""

    fragment_ordinal: int
    contig_name: str
    reference_start: int
    reference_end: int
    mates: tuple[FastqMate, ...]


@dataclass(frozen=True)
class ColumnarReadBatch:
    """Columnar read and details state retained for BAM serialization only."""

    fragment_count: int
    mate_count: int
    mate_sequences: bytes
    read_length: int
    quality_byte: int
    base_state_codes: bytes
    read_summaries: bytes
    fragment_summaries: bytes | None = None
    fragment_realizations: tuple[bytes, ...] | None = None

    def __post_init__(self) -> None:
        if (
            isinstance(self.fragment_count, bool)
            or not isinstance(self.fragment_count, int)
            or self.fragment_count <= 0
            or isinstance(self.mate_count, bool)
            or not isinstance(self.mate_count, int)
            or self.mate_count not in (self.fragment_count, self.fragment_count * 2)
            or isinstance(self.read_length, bool)
            or not isinstance(self.read_length, int)
            or self.read_length <= 0
            or isinstance(self.quality_byte, bool)
            or not isinstance(self.quality_byte, int)
            or not 33 <= self.quality_byte <= 126
        ):
            raise ProcessError("processed read-batch geometry is invalid")
        if (
            type(self.mate_sequences) is not bytes
            or len(self.mate_sequences) != self.mate_count * self.read_length
            or type(self.base_state_codes) is not bytes
            or len(self.base_state_codes) != self.mate_count * self.read_length
            or type(self.read_summaries) is not bytes
            or len(self.read_summaries) != self.mate_count * 12 * 2
            or (
                self.fragment_summaries is not None
                and (
                    type(self.fragment_summaries) is not bytes
                    or len(self.fragment_summaries)
                    != self.fragment_count * 12 * 2
                )
            )
            or (
                self.fragment_realizations is not None
                and (
                    not isinstance(self.fragment_realizations, tuple)
                    or len(self.fragment_realizations) != self.fragment_count
                    or any(
                        type(value) is not bytes
                        or not value
                        or any(byte < 32 or byte > 126 or byte == 9 for byte in value)
                        for value in self.fragment_realizations
                    )
                )
            )
        ):
            raise ProcessError("processed read-batch columns are invalid")


@dataclass(frozen=True)
class EncodedFastqBatch:
    """Consecutive FASTQ records formatted without per-mate Python objects."""

    read1: bytes
    read2: bytes | None
    record_lengths: tuple[tuple[int, int], ...]
    mate_sequences: bytes | None = None
    read_length: int = 0
    quality_byte: int = 0
    base_state_codes: bytes | None = None
    read_summaries: bytes | None = None

    def __post_init__(self) -> None:
        if type(self.read1) is not bytes or (
            self.read2 is not None and type(self.read2) is not bytes
        ):
            raise ProcessError("formatted FASTQ batches must contain bytes")
        if not isinstance(self.record_lengths, tuple) or not self.record_lengths:
            raise ProcessError("formatted FASTQ record lengths must be non-empty")
        for lengths in self.record_lengths:
            if (
                not isinstance(lengths, tuple)
                or len(lengths) != 2
                or any(
                    isinstance(value, bool)
                    or not isinstance(value, int)
                    or value < 0
                    for value in lengths
                )
                or lengths[0] == 0
            ):
                raise ProcessError("formatted FASTQ record lengths are invalid")
        if sum(value[0] for value in self.record_lengths) != len(self.read1):
            raise ProcessError("formatted read1 bytes disagree with record lengths")
        expected_read2 = sum(value[1] for value in self.record_lengths)
        if (self.read2 is None) != (expected_read2 == 0) or (
            self.read2 is not None and len(self.read2) != expected_read2
        ):
            raise ProcessError("formatted read2 bytes disagree with record lengths")
        if self.mate_sequences is None:
            if (
                self.read_length != 0
                or self.quality_byte != 0
                or self.base_state_codes is not None
                or self.read_summaries is not None
            ):
                raise ProcessError("discarded mate sequences retain geometry")
        else:
            mate_count = len(self.record_lengths) * (2 if self.read2 is not None else 1)
            if (
                type(self.mate_sequences) is not bytes
                or isinstance(self.read_length, bool)
                or not isinstance(self.read_length, int)
                or self.read_length <= 0
                or len(self.mate_sequences) != mate_count * self.read_length
                or isinstance(self.quality_byte, bool)
                or not isinstance(self.quality_byte, int)
                or not 33 <= self.quality_byte <= 126
            ):
                raise ProcessError("retained mate sequence geometry is invalid")
            if (self.base_state_codes is None) != (self.read_summaries is None):
                raise ProcessError("retained details columns must be paired")
            if self.base_state_codes is not None and (
                type(self.base_state_codes) is not bytes
                or len(self.base_state_codes) != mate_count * self.read_length
                or type(self.read_summaries) is not bytes
                or len(self.read_summaries) != mate_count * 12 * 2
            ):
                raise ProcessError("retained details-column geometry is invalid")


def _philox_pairs(
    key: int,
    entity_ordinals: np.ndarray,
    local_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    entity = np.asarray(entity_ordinals, dtype=_UINT64)
    local = np.asarray(local_indices, dtype=_UINT64)
    if entity.shape != local.shape:
        raise ProcessError("NumPy Philox counter columns must have one shape")
    entity = np.ascontiguousarray(entity)
    local = np.ascontiguousarray(local)
    pair_0_data, pair_1_data = _cext_philox_pairs(key, entity, local)
    pair_0 = np.frombuffer(pair_0_data, dtype=_UINT64).reshape(entity.shape)
    pair_1 = np.frombuffer(pair_1_data, dtype=_UINT64).reshape(entity.shape)
    pair_0.flags.writeable = False
    pair_1.flags.writeable = False
    return pair_0, pair_1


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
        raise ProcessError("{} offsets must start at zero".format(name))
    if np.any(offsets[1:] < offsets[:-1]):
        raise ProcessError("{} offsets must be nondecreasing".format(name))


def _contig_fragment_groups(contig_indices: np.ndarray) -> dict[int, np.ndarray]:
    groups = {}  # type: dict[int, list]
    for index, contig_index in enumerate(contig_indices):
        groups.setdefault(int(contig_index), []).append(index)
    return {
        contig_index: np.asarray(indices, dtype=np.intp)
        for contig_index, indices in groups.items()
    }


READ_NAME_CONTRACT = "bsreadsim-read-name"
_MAX_U64 = (1 << 64) - 1

# Cross-stage batches retain only numeric identity columns.  FASTQ and BAM
# serializers materialize this canonical text lazily at the output boundary;
# no per-fragment QNAME object or eagerly allocated string is retained.

class ReadNameError(ValueError):
    """A fragment cannot be represented by the read-name contract."""


def fragment_identifier_coordinates(
    reference_start: int,
    reference_end: int,
) -> tuple[int, int]:
    """Convert a zero-based half-open envelope to display coordinates.

    A physical insertion-only fragment has a zero-width reference envelope.
    It is represented by its one-based insertion anchor on both sides so the
    identifier remains a valid inclusive interval.
    """

    for name, value in (
        ("fragment reference start", reference_start),
        ("fragment reference end", reference_end),
    ):
        if isinstance(value, bool) or not isinstance(value, int):
            raise ReadNameError("{} must be an integer".format(name))
        if not 0 <= value <= _MAX_U64:
            raise ReadNameError("{} is outside uint64".format(name))
    if reference_start > reference_end:
        raise ReadNameError("fragment reference envelope is reversed")
    if reference_start == _MAX_U64:
        raise ReadNameError("fragment reference start cannot be converted to 1-based")

    left = reference_start + 1
    right = reference_end if reference_end > reference_start else left
    return left, right


def format_fragment_identifier(
    contig_name: str,
    reference_start: int,
    reference_end: int,
    fragment_ordinal: int,
    *,
    pair_number: int | None = None,
) -> str:
    """Return ``contig:start-end:hex`` with an optional mate suffix.

    Coordinates are one-based inclusive.  The fragment ordinal is canonical
    variable-width lowercase hexadecimal without ``0x`` or leading zeroes.
    """

    if (
        not isinstance(contig_name, str)
        or not contig_name
        or any(
            character.isspace() or not character.isprintable()
            for character in contig_name
        )
    ):
        raise ReadNameError(
            "contig name must be non-empty printable text without whitespace"
        )
    if (
        isinstance(fragment_ordinal, bool)
        or not isinstance(fragment_ordinal, int)
        or not 0 <= fragment_ordinal <= _MAX_U64
    ):
        raise ReadNameError("fragment ordinal is outside uint64")
    if pair_number is not None and (
        isinstance(pair_number, bool)
        or not isinstance(pair_number, int)
        or pair_number not in (1, 2)
    ):
        raise ReadNameError("pair number must be integer 1 or 2")

    left, right = fragment_identifier_coordinates(reference_start, reference_end)
    identifier = "{}:{}-{}:{}".format(
        contig_name,
        left,
        right,
        format(fragment_ordinal, "x"),
    )
    if pair_number is not None:
        identifier += "/{}".format(pair_number)
    return identifier


def _require_probability(name: str, value: float) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ProcessError("{} must be a finite number in [0, 1]".format(name))
    converted = float(value)
    if not math.isfinite(converted) or not 0.0 <= converted <= 1.0:
        raise ProcessError("{} must be a finite number in [0, 1]".format(name))


__all__ = [
    "BaseState",
    "CaptureStrand",
    "ConversionMode",
    "Fragment",
    "FragmentSummary",
    "Mate",
    "MethylationAllele",
    "MethylationContext",
    "MethylationModelBatch",
    "MethylationSite",
    "MethylationSource",
    "NO_VARIANT_INDEX",
    "FastqFragment",
    "FastqMate",
    "EncodedFastqBatch",
    "ColumnarFragmentBatch",
    "ColumnarReadBatch",
    "ProcessError",
    "ProcessedFragment",
    "ProcessedMate",
    "READ_NAME_CONTRACT",
    "ReadNameError",
    "SiteReference",
    "SiteState",
    "UniformError",
    "UniformQuality",
    "Variant",
    "VariantKind",
    "format_fragment_identifier",
    "fragment_identifier_coordinates",
]
