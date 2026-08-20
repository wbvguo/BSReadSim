"""Pure Python post-processing for decoded scientific fragments.

This module owns no file I/O, concurrency, process lifecycle, or FASTQ
publication. A decoded :class:`~bsreadsim.model.Fragment` is transformed
by pure functions using explicit stateless RNG counters, making uniform and
versioned empirical models independent of worker count, chunking, and
completion order.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
import math
import struct
from typing import Optional, Sequence, Tuple, Union

from .model import (
    CaptureStrand,
    Fragment,
    Mate,
    MethylationAllele,
    MethylationContext,
    MethylationSource,
    VariantEvent,
)
from .rng import (
    RNGStage,
    _bernoulli_unchecked as bernoulli,
    _u64_unchecked as u64,
    derive_key,
)
from .sequencing_models import QualityConfusionModel, QualityMarkovModel

try:
    from ._native import (
        apply_uniform_errors as _native_apply_uniform_errors,
        sample_bernoulli_sites as _native_sample_bernoulli_sites,
    )
except ImportError:
    _native_apply_uniform_errors = None
    _native_sample_bernoulli_sites = None


POSTPROCESS_STAGE_ORDER = (
    "methylation",
    "fragment-orientation",
    "fragment-conversion",
    "mate-derivation",
    "quality",
    "sequencing-error",
    "format-ready",
)

_UINT64_MAX = (1 << 64) - 1
_COMPLEMENT = (3, 2, 1, 0, 4)
_BASE_TEXT_TRANSLATION = bytes.maketrans(bytes(range(5)), b"ACGTN")
_QUALITY_TEXT_TRANSLATION = bytes.maketrans(
    bytes(range(94)),
    bytes(range(33, 127)),
)
_ALTERNATIVE_BASES = (
    (1, 2, 3),
    (0, 2, 3),
    (0, 1, 3),
    (0, 1, 2),
)
_COLUMN_U64 = struct.Struct("<Q")


class PostprocessError(ValueError):
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

    fragment_ordinals: Tuple[int, ...]
    contig_names: Tuple[str, ...]
    fragment_ordinal_data: bytes
    contig_indices: bytes
    template_offsets: bytes
    template_bases: bytes
    site_offsets: bytes
    site_indices: bytes
    site_template_offsets: bytes
    site_reference_positions: bytes
    site_contexts: bytes
    site_sources: bytes
    site_alleles: bytes
    site_probabilities: bytes

    def __post_init__(self) -> None:
        if not isinstance(self.fragment_ordinals, tuple) or not self.fragment_ordinals:
            raise PostprocessError(
                "columnar batch ordinals must be a non-empty immutable tuple"
            )
        if not isinstance(self.contig_names, tuple):
            raise PostprocessError(
                "columnar batch contig names must be an immutable tuple"
            )
        fragment_count = len(self.fragment_ordinals)
        if len(self.contig_names) != fragment_count:
            raise PostprocessError(
                "columnar batch contig-name count disagrees with ordinals"
            )
        buffers = (
            self.fragment_ordinal_data,
            self.contig_indices,
            self.template_offsets,
            self.template_bases,
            self.site_offsets,
            self.site_indices,
            self.site_template_offsets,
            self.site_reference_positions,
            self.site_contexts,
            self.site_sources,
            self.site_alleles,
            self.site_probabilities,
        )
        if any(type(value) is not bytes for value in buffers):
            raise PostprocessError("columnar batch buffers must be immutable bytes")
        if len(self.fragment_ordinal_data) != fragment_count * 8:
            raise PostprocessError("columnar fragment ordinal buffer has wrong length")
        if len(self.contig_indices) != fragment_count * 4:
            raise PostprocessError("columnar contig-index buffer has wrong length")
        if len(self.template_offsets) != (fragment_count + 1) * 8:
            raise PostprocessError("columnar template offsets have wrong length")
        if len(self.site_offsets) != (fragment_count + 1) * 8:
            raise PostprocessError("columnar site offsets have wrong length")
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
            or len(self.site_sources) != site_count
            or len(self.site_alleles) != site_count
            or len(self.site_probabilities) != site_count * 4
        ):
            raise PostprocessError("columnar site buffers disagree on site count")
        encoded_ordinals = tuple(
            _COLUMN_U64.unpack_from(self.fragment_ordinal_data, index * 8)[0]
            for index in range(fragment_count)
        )
        if encoded_ordinals != self.fragment_ordinals:
            raise PostprocessError(
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
                raise PostprocessError(
                    "columnar {} offsets must start at zero".format(name)
                )
        elif value < previous:
            raise PostprocessError(
                "columnar {} offsets must be nondecreasing".format(name)
            )
        previous = value
    if previous != final_value:
        raise PostprocessError(
            "columnar {} offsets disagree with flat data".format(name)
        )


class ConversionMode(IntEnum):
    """Bisulfite substitution visible in the current oriented mate."""

    C_TO_T = 0
    G_TO_A = 1


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
            raise PostprocessError("quality phred must be an integer in [0, 93]")


@dataclass(frozen=True)
class UniformError:
    """One constant substitution probability for every A/C/G/T call."""

    rate: float = 0.005

    def __post_init__(self) -> None:
        _require_probability("error rate", self.rate)


QualityPolicy = Union[UniformQuality, QualityMarkovModel]
ErrorPolicy = Union[UniformError, QualityConfusionModel]


@dataclass(frozen=True)
class PostprocessConfig:
    """Validated immutable settings for all released Python read stages."""

    master_seed: int
    directional: bool
    conversion_rate: float
    quality: QualityPolicy
    error: ErrorPolicy

    def __post_init__(self) -> None:
        _require_u64("master_seed", self.master_seed)
        if not isinstance(self.directional, bool):
            raise PostprocessError("directional must be a boolean")
        _require_probability("conversion_rate", self.conversion_rate)
        if not isinstance(self.quality, (UniformQuality, QualityMarkovModel)):
            raise PostprocessError("quality policy is outside the released contract")
        if not isinstance(self.error, (UniformError, QualityConfusionModel)):
            raise PostprocessError("error policy is outside the released contract")


class UniformPostprocessConfig(PostprocessConfig):
    """Convenience constructor for uniform quality and error policies."""

    def __init__(
        self,
        master_seed: int,
        directional: bool = True,
        conversion_rate: float = 0.998,
        quality_phred: int = 40,
        error_rate: float = 0.005,
    ) -> None:
        super().__init__(
            master_seed=master_seed,
            directional=directional,
            conversion_rate=conversion_rate,
            quality=UniformQuality(quality_phred),
            error=UniformError(error_rate),
        )

    @property
    def quality_phred(self) -> int:
        return self.quality.phred

    @property
    def error_rate(self) -> float:
        return self.error.rate


@dataclass(frozen=True)
class SiteState:
    """The single latent draw for one fragment-level methylation site."""

    site_index: int
    template_offset: int
    reference_pos: int
    context: MethylationContext
    source: MethylationSource
    allele: MethylationAllele
    methylated: bool
    probability: float


@dataclass(frozen=True)
class BaseAnnotation:
    """Traceable state for one base after every post-processing stage.

    ``site_index`` is ``None`` when the core did not declare a methylation site.
    For an undeclared base that is nevertheless a target of the mate's
    conversion mode, ``methylated`` is ``False`` to record its implicit
    unmethylated state.  Bases that are neither declared sites nor conversion
    targets retain ``None``.
    """

    read_offset: int
    reference_pos: int
    variant_event_id: int
    site_index: Optional[int]
    methylated: Optional[bool]
    oriented_base: int
    post_conversion_base: int
    final_base: int
    conversion_attempted: bool
    conversion_succeeded: bool
    sequencing_error: bool
    quality_phred: int


@dataclass(frozen=True)
class _CompactAnnotations:
    """Worker-only columnar annotation state for native truth formatting.

    Keeping one immutable column per field avoids allocating one 12-field
    Python dataclass for every read base.  The public ``process_fragment`` path
    continues to return regular ``BaseAnnotation`` tuples by default.
    """

    reference_positions: Tuple[int, ...]
    variant_event_ids: Tuple[int, ...]
    site_indices: Tuple[Optional[int], ...]
    methylated: Tuple[Optional[bool], ...]
    oriented_bases: bytes
    post_conversion_bases: bytes
    final_bases: bytes
    attempted: Tuple[bool, ...]
    succeeded: Tuple[bool, ...]
    error_flags: Tuple[bool, ...]
    quality_phreds: Tuple[int, ...]

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
    annotations: Union[Tuple[BaseAnnotation, ...], _CompactAnnotations]


@dataclass(frozen=True)
class ProcessedFragment:
    """Atomic post-processing result for one SE or PE fragment."""

    fragment_ordinal: int
    contig_name: str
    reference_start: int
    reference_end: int
    haplotype: int
    fragment_conversion_mode: ConversionMode
    variant_events: Tuple[VariantEvent, ...]
    site_states: Tuple[SiteState, ...]
    mates: Tuple[ProcessedMate, ...]


@dataclass(frozen=True)
class _OrientedMate:
    mate: Mate
    conversion_mode: ConversionMode
    bases: bytes
    reference_positions: Tuple[int, ...]
    variant_event_ids: Tuple[int, ...]


@dataclass(frozen=True)
class _ConvertedFragment:
    conversion_mode: ConversionMode
    bases: bytes
    site_indices: Tuple[Optional[int], ...]
    methylated: Tuple[Optional[bool], ...]
    attempted: Tuple[bool, ...]
    succeeded: Tuple[bool, ...]


@dataclass(frozen=True)
class _ConvertedMate:
    oriented: _OrientedMate
    bases: bytes
    site_indices: Tuple[Optional[int], ...]
    methylated: Tuple[Optional[bool], ...]
    attempted: Tuple[bool, ...]
    succeeded: Tuple[bool, ...]


@dataclass(frozen=True)
class _QualityMate:
    converted: _ConvertedMate
    qualities: Tuple[int, ...]


@dataclass(frozen=True)
class _ErroredMate:
    quality: _QualityMate
    bases: bytes
    error_flags: Tuple[bool, ...]


def process_fragment(
    fragment: Fragment,
    contig_name: str,
    config: PostprocessConfig,
    *,
    compact_annotations: bool = False,
    include_truth: bool = True,
) -> ProcessedFragment:
    """Apply the ordered Python-owned stages to one decoded fragment.

    The function is stateless.  ``fragment.fragment_ordinal`` is the RNG entity
    ordinal for every Python stage; local draw numbers are stable site indices
    or collision-free packed mate/base offsets.
    """
    if not isinstance(fragment, Fragment):
        raise PostprocessError("fragment must be a decoded protocol Fragment")
    if not isinstance(config, PostprocessConfig):
        raise PostprocessError("config must be PostprocessConfig")
    if not isinstance(compact_annotations, bool):
        raise PostprocessError("compact_annotations must be a boolean")
    if not isinstance(include_truth, bool):
        raise PostprocessError("include_truth must be a boolean")
    if compact_annotations and not include_truth:
        raise PostprocessError("compact_annotations requires include_truth")
    if not isinstance(contig_name, str) or not contig_name or "\x00" in contig_name:
        raise PostprocessError("contig_name must be non-empty text without NUL")

    sampled_methylation = _sample_methylation_values(
        fragment,
        config,
    )
    return _process_fragment_with_states(
        fragment,
        contig_name,
        config,
        sampled_methylation,
        compact_annotations=compact_annotations,
        include_truth=include_truth,
    )


def process_fragment_batch(
    fragments: Tuple[Fragment, ...],
    contig_names: Tuple[str, ...],
    config: PostprocessConfig,
    *,
    compact_annotations: bool = False,
    include_truth: bool = True,
) -> Tuple[ProcessedFragment, ...]:
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
        compact_annotations=compact_annotations,
        include_truth=include_truth,
    )
    sampled_batch = _sample_methylation_batch_values(
        fragments,
        config,
    )
    return tuple(
        _process_fragment_with_states(
            fragment,
            contig_name,
            config,
            sampled_methylation,
            compact_annotations=compact_annotations,
            include_truth=include_truth,
        )
        for fragment, contig_name, sampled_methylation in zip(
            fragments,
            contig_names,
            sampled_batch,
        )
    )


def _process_fragment_with_states(
    fragment: Fragment,
    contig_name: str,
    config: PostprocessConfig,
    sampled_methylation: Tuple[bool, ...],
    *,
    compact_annotations: bool,
    include_truth: bool,
) -> ProcessedFragment:
    site_states = (
        _materialize_site_states(fragment, sampled_methylation)
        if include_truth
        else ()
    )
    fragment_mode = _select_fragment_conversion_mode(fragment, config)
    converted_fragment = _convert_fragment(
        fragment,
        fragment_mode,
        sampled_methylation,
        config,
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
        _make_format_ready(
            mate,
            compact_annotations=compact_annotations,
            include_annotations=include_truth,
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
        fragment.variant_events,
        site_states,
        processed_mates,
    )


def _validate_fragment_batch_request(
    fragments: Tuple[Fragment, ...],
    contig_names: Tuple[str, ...],
    config: PostprocessConfig,
    *,
    compact_annotations: bool,
    include_truth: bool,
) -> None:
    if not isinstance(fragments, tuple) or not fragments:
        raise PostprocessError("fragments must be a non-empty immutable tuple")
    if not isinstance(contig_names, tuple):
        raise PostprocessError("contig_names must be an immutable tuple")
    if len(contig_names) != len(fragments):
        raise PostprocessError("contig_names count disagrees with fragments")
    if not isinstance(config, PostprocessConfig):
        raise PostprocessError("config must be PostprocessConfig")
    if not isinstance(compact_annotations, bool):
        raise PostprocessError("compact_annotations must be a boolean")
    if not isinstance(include_truth, bool):
        raise PostprocessError("include_truth must be a boolean")
    if compact_annotations and not include_truth:
        raise PostprocessError("compact_annotations requires include_truth")
    for fragment, contig_name in zip(fragments, contig_names):
        if not isinstance(fragment, Fragment):
            raise PostprocessError("fragment must be a decoded protocol Fragment")
        if (
            not isinstance(contig_name, str)
            or not contig_name
            or "\x00" in contig_name
        ):
            raise PostprocessError("contig_name must be non-empty text without NUL")


def _select_fragment_conversion_mode(
    fragment: Fragment,
    config: PostprocessConfig,
) -> ConversionMode:
    # Target capture metadata is authoritative.  An uncaptured directional
    # library keeps the forward/C-to-T strand.  An undirectional library draws
    # one fragment-level strand, shared by both mates.
    if fragment.capture_strand is CaptureStrand.FORWARD:
        return ConversionMode.C_TO_T
    if fragment.capture_strand is CaptureStrand.REVERSE:
        return ConversionMode.G_TO_A
    if fragment.capture_strand is not CaptureStrand.UNKNOWN:
        raise PostprocessError("fragment has an unsupported capture strand")
    if config.directional:
        return ConversionMode.C_TO_T

    key = derive_key(
        config.master_seed,
        RNGStage.LIBRARY_ORIENTATION,
        fragment.contig_index,
    )
    reverse = bernoulli(key, fragment.fragment_ordinal, 0, 0.5)
    return ConversionMode.G_TO_A if reverse else ConversionMode.C_TO_T


def _sample_methylation_batch_values(
    fragments: Tuple[Fragment, ...],
    config: PostprocessConfig,
) -> Tuple[Tuple[bool, ...], ...]:
    return tuple(
        _sample_methylation_values(fragment, config)
        for fragment in fragments
    )


def _sample_methylation_values(
    fragment: Fragment,
    config: PostprocessConfig,
) -> Tuple[bool, ...]:
    key = derive_key(
        config.master_seed,
        RNGStage.SITE_STATE,
        fragment.contig_index,
    )
    if _native_sample_bernoulli_sites is not None:
        try:
            return _native_sample_bernoulli_sites(
                fragment.methylation_sites,
                key,
                fragment.fragment_ordinal,
            )
        except (TypeError, ValueError, OverflowError) as error:
            raise PostprocessError(
                "native Bernoulli site sampling failed: {}".format(error)
            ) from error
    return tuple(
        bernoulli(
            key,
            fragment.fragment_ordinal,
            site.site_index,
            site.methylation_probability,
        )
        for site in fragment.methylation_sites
    )


def _materialize_site_states(
    fragment: Fragment,
    sampled: Sequence[Union[bool, int]],
) -> Tuple[SiteState, ...]:
    return tuple(
        SiteState(
            site.site_index,
            site.template_offset,
            site.reference_pos,
            site.context,
            site.source,
            site.allele,
            bool(state),
            float(site.methylation_probability),
        )
        for site, state in zip(fragment.methylation_sites, sampled)
    )


def _convert_fragment(
    fragment: Fragment,
    conversion_mode: ConversionMode,
    sampled_methylation: Sequence[Union[bool, int]],
    config: PostprocessConfig,
) -> _ConvertedFragment:
    key = derive_key(
        config.master_seed,
        RNGStage.CONVERSION,
        fragment.contig_index,
    )
    bases = bytearray(fragment.template_bases)
    length = len(bases)
    site_indices = [None] * length
    methylated = [None] * length
    attempted = [False] * length
    succeeded = [False] * length
    if len(sampled_methylation) != len(fragment.methylation_sites):
        raise PostprocessError("methylation-state count disagrees with fragment sites")
    for site_index, site in enumerate(fragment.methylation_sites):
        if site.site_index != site_index:
            raise PostprocessError("fragment methylation sites are not index ordered")
        offset = site.template_offset
        if not 0 <= offset < length:
            raise PostprocessError("fragment site points outside the template")
        if site_indices[offset] is not None:
            raise PostprocessError("fragment sites contain a duplicate template offset")
        site_indices[offset] = site_index
        methylated[offset] = bool(sampled_methylation[site_index])

    # Conversion is one physical fragment event. A convertible C/G without a
    # declared protocol site is implicitly unmethylated, and every template
    # offset receives at most one conversion draw even when both mates observe it.
    for offset, template_base in enumerate(fragment.template_bases):
        if not _is_conversion_target(template_base, conversion_mode):
            continue
        if site_indices[offset] is None:
            methylated[offset] = False
        elif methylated[offset]:
            continue
        attempted[offset] = True
        converted = bernoulli(
            key,
            fragment.fragment_ordinal,
            offset,
            config.conversion_rate,
        )
        if converted:
            bases[offset] = (
                3 if conversion_mode is ConversionMode.C_TO_T else 0
            )
            succeeded[offset] = True

    return _ConvertedFragment(
        conversion_mode,
        bytes(bases),
        tuple(site_indices),
        tuple(methylated),
        tuple(attempted),
        tuple(succeeded),
    )


def _derive_mates(
    fragment: Fragment,
    converted: _ConvertedFragment,
) -> Tuple[_ConvertedMate, ...]:
    derived = []
    for mate in sorted(fragment.mates, key=lambda value: value.mate_index):
        template_slice = slice(mate.template_start, mate.template_end)
        original_bases = fragment.template_bases[template_slice]
        converted_bases = converted.bases[template_slice]
        positions = tuple(fragment.reference_positions[template_slice])
        event_ids = tuple(fragment.base_event_ids[template_slice])
        site_indices = tuple(converted.site_indices[template_slice])
        methylated = tuple(converted.methylated[template_slice])
        attempted = tuple(converted.attempted[template_slice])
        succeeded = tuple(converted.succeeded[template_slice])
        mate_mode = converted.conversion_mode

        if mate.reverse_complement:
            original_bases = bytes(
                _COMPLEMENT[base] for base in reversed(original_bases)
            )
            converted_bases = bytes(
                _COMPLEMENT[base] for base in reversed(converted_bases)
            )
            positions = tuple(reversed(positions))
            event_ids = tuple(reversed(event_ids))
            site_indices = tuple(reversed(site_indices))
            methylated = tuple(reversed(methylated))
            attempted = tuple(reversed(attempted))
            succeeded = tuple(reversed(succeeded))
            mate_mode = _opposite_mode(mate_mode)

        expected_refs = {
            (read_offset, site_index)
            for read_offset, site_index in enumerate(site_indices)
            if site_index is not None
        }
        actual_refs = {
            (site_ref.read_offset, site_ref.site_index)
            for site_ref in mate.site_refs
        }
        if len(actual_refs) != len(mate.site_refs) or actual_refs != expected_refs:
            raise PostprocessError(
                "mate site references disagree with fragment site projection"
            )

        oriented = _OrientedMate(
            mate,
            mate_mode,
            original_bases,
            positions,
            event_ids,
        )
        derived.append(
            _ConvertedMate(
                oriented,
                converted_bases,
                site_indices,
                methylated,
                attempted,
                succeeded,
            )
        )
    return tuple(derived)


def _generate_quality(
    fragment: Fragment,
    converted: _ConvertedMate,
    config: PostprocessConfig,
) -> _QualityMate:
    if isinstance(config.quality, UniformQuality):
        qualities = (config.quality.phred,) * len(converted.bases)
    else:
        qualities = config.quality.sample(
            config.master_seed,
            fragment.contig_index,
            fragment.fragment_ordinal,
            converted.oriented.mate.mate_index,
            len(converted.bases),
        )
    return _QualityMate(
        converted,
        qualities,
    )


def _apply_errors(
    fragment: Fragment,
    quality: _QualityMate,
    config: PostprocessConfig,
) -> _ErroredMate:
    if isinstance(config.error, QualityConfusionModel):
        sampled = config.error.sample(
            config.master_seed,
            fragment.contig_index,
            fragment.fragment_ordinal,
            quality.converted.oriented.mate.mate_index,
            quality.converted.bases,
            quality.qualities,
        )
        return _ErroredMate(
            quality,
            bytes(sampled),
            tuple(
                before != after
                for before, after in zip(quality.converted.bases, sampled)
            ),
        )

    key = derive_key(
        config.master_seed,
        RNGStage.SEQUENCING_ERROR,
        fragment.contig_index,
    )
    if _native_apply_uniform_errors is not None:
        try:
            sampled_bases, sampled_flags = _native_apply_uniform_errors(
                quality.converted.bases,
                key,
                fragment.fragment_ordinal,
                quality.converted.oriented.mate.mate_index,
                config.error.rate,
            )
        except (TypeError, ValueError, OverflowError) as error:
            raise PostprocessError(
                "native uniform sequencing-error sampling failed: {}".format(
                    error
                )
            ) from error
        return _ErroredMate(
            quality,
            sampled_bases,
            tuple(value != 0 for value in sampled_flags),
        )
    bases = bytearray(quality.converted.bases)
    error_flags = [False] * len(bases)
    mate_index = quality.converted.oriented.mate.mate_index

    for offset, base in enumerate(bases):
        if base == 4:  # Ambiguous bases have no defined A/C/G/T confusion row.
            continue
        local_index = _pack_mate_base(mate_index, offset)
        if bernoulli(
            key,
            fragment.fragment_ordinal,
            local_index,
            config.error.rate,
        ):
            bases[offset] = _alternative_base(
                key, fragment.fragment_ordinal, local_index, base
            )
            error_flags[offset] = True

    return _ErroredMate(
        quality,
        bytes(bases),
        tuple(error_flags),
    )


def _make_format_ready(
    errored: _ErroredMate,
    *,
    compact_annotations: bool,
    include_annotations: bool,
) -> ProcessedMate:
    quality = errored.quality
    converted = quality.converted
    oriented = converted.oriented
    reference_positions = oriented.reference_positions
    variant_event_ids = oriented.variant_event_ids
    site_indices = converted.site_indices
    methylated = converted.methylated
    oriented_bases = oriented.bases
    converted_bases = converted.bases
    attempted = converted.attempted
    succeeded = converted.succeeded
    error_flags = errored.error_flags
    qualities = quality.qualities
    if not include_annotations:
        annotations = ()
    elif compact_annotations:
        annotations = _CompactAnnotations(
            reference_positions,
            variant_event_ids,
            site_indices,
            methylated,
            oriented_bases,
            converted_bases,
            errored.bases,
            attempted,
            succeeded,
            error_flags,
            qualities,
        )  # type: Union[Tuple[BaseAnnotation, ...], _CompactAnnotations]
    else:
        annotation_values = []
        for offset, final_base in enumerate(errored.bases):
            annotation_values.append(
                BaseAnnotation(
                    offset,
                    reference_positions[offset],
                    variant_event_ids[offset],
                    site_indices[offset],
                    methylated[offset],
                    oriented_bases[offset],
                    converted_bases[offset],
                    final_base,
                    attempted[offset],
                    succeeded[offset],
                    error_flags[offset],
                    qualities[offset],
                )
            )
        annotations = tuple(annotation_values)

    return ProcessedMate(
        oriented.mate.mate_index,
        oriented.mate.reverse_complement,
        oriented.conversion_mode,
        oriented.mate.reference_start,
        oriented.mate.reference_end,
        errored.bases.translate(_BASE_TEXT_TRANSLATION).decode("ascii"),
        bytes(qualities).translate(_QUALITY_TEXT_TRANSLATION).decode("ascii"),
        annotations,
    )


def _is_conversion_target(base: int, mode: ConversionMode) -> bool:
    return (mode is ConversionMode.C_TO_T and base == 1) or (
        mode is ConversionMode.G_TO_A and base == 2
    )


def _opposite_mode(mode: ConversionMode) -> ConversionMode:
    return (
        ConversionMode.G_TO_A
        if mode is ConversionMode.C_TO_T
        else ConversionMode.C_TO_T
    )


def _pack_mate_base(mate_index: int, read_offset: int) -> int:
    # The protocol bounds both values to u32. The packed value is therefore a
    # collision-free u64 local_index independent of iteration order.
    return (mate_index << 32) | read_offset


def _alternative_base(
    key: int,
    fragment_ordinal: int,
    local_index: int,
    original_base: int,
) -> int:
    # Pair zero of the Philox block decides whether an error occurs.  Pair one
    # independently selects one of the three non-original A/C/G/T bases.
    choice = (
        u64(key, fragment_ordinal, local_index, pair=1) * 3
    ) >> 64
    return _ALTERNATIVE_BASES[original_base][choice]


def _require_u64(name: str, value: int) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 0 <= value <= _UINT64_MAX
    ):
        raise PostprocessError("{} must be an unsigned 64-bit integer".format(name))


def _require_probability(name: str, value: float) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise PostprocessError("{} must be a finite number in [0, 1]".format(name))
    converted = float(value)
    if not math.isfinite(converted) or not 0.0 <= converted <= 1.0:
        raise PostprocessError("{} must be a finite number in [0, 1]".format(name))
