"""Apply fragment-level bisulfite orientation and conversion."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from .batch import (
    ConversionMode,
    Fragment,
    ProcessError,
    _ConvertedFragment,
    _bernoulli_from_pair,
    _contig_fragment_groups,
    _philox_pairs,
)
from .config import ProcessConfig
from ..rng import RNGStage, _bernoulli_unchecked as bernoulli, derive_key


def _convert_fragment(
    fragment: Fragment,
    conversion_mode: ConversionMode,
    sampled_methylation: Sequence[bool | int],
    config: ProcessConfig,
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
        raise ProcessError("methylation-state count disagrees with fragment sites")
    for site_index, site in enumerate(fragment.methylation_sites):
        if site.site_index != site_index:
            raise ProcessError("fragment methylation sites are not index ordered")
        offset = site.template_offset
        if not 0 <= offset < length:
            raise ProcessError("fragment site points outside the template")
        if site_indices[offset] is not None:
            raise ProcessError("fragment sites contain a duplicate template offset")
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
                3 if conversion_mode is ConversionMode.C2T else 0
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


def _is_conversion_target(base: int, mode: ConversionMode) -> bool:
    return (mode is ConversionMode.C2T and base == 1) or (
        mode is ConversionMode.G2A and base == 2
    )


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
