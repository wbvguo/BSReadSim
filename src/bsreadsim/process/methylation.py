"""Assign fragment methylation states with Bernoulli or future BiLSTM models."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol, runtime_checkable

import numpy as np

from .batch import (
    Fragment,
    ColumnarFragmentBatch,
    ProcessError,
    SiteState,
    _F4,
    _U4,
    _U8,
    _bernoulli_from_pair,
    _contig_fragment_groups,
    _owners_from_offsets,
    _philox_pairs,
)
from ..rng import RNGStage, derive_key

from .._native import sample_bernoulli_sites as _native_sample_bernoulli_sites

if TYPE_CHECKING:
    from .config import ProcessConfig


@runtime_checkable
class MethylationStateModel(Protocol):
    """Stable fragment-state model boundary for built-ins and future plugins."""

    contract: str

    def sample_fragment(
        self, fragment: Fragment, config: ProcessConfig
    ) -> tuple[bool, ...]: ...

    def sample_batch(
        self, batch: ColumnarFragmentBatch, config: ProcessConfig
    ) -> np.ndarray: ...


@dataclass(frozen=True)
class BernoulliStateModel:
    """Independent per-site sampling through the frozen Philox domain."""

    contract: str = "bernoulli-site-v1"

    def sample_fragment(
        self, fragment: Fragment, config: ProcessConfig
    ) -> tuple[bool, ...]:
        return _sample_methylation_values_bernoulli(fragment, config)

    def sample_batch(
        self, batch: ColumnarFragmentBatch, config: ProcessConfig
    ) -> np.ndarray:
        return _sample_site_states_bernoulli(batch, config)

def _sample_methylation_batch_values(
    fragments: tuple[Fragment, ...],
    config: ProcessConfig,
) -> tuple[tuple[bool, ...], ...]:
    return tuple(
        config.methylation_model.sample_fragment(fragment, config)
        for fragment in fragments
    )


def _sample_methylation_values(
    fragment: Fragment,
    config: ProcessConfig,
) -> tuple[bool, ...]:
    return config.methylation_model.sample_fragment(fragment, config)


def _sample_methylation_values_bernoulli(
    fragment: Fragment,
    config: ProcessConfig,
) -> tuple[bool, ...]:
    key = derive_key(
        config.master_seed,
        RNGStage.SITE_STATE,
        fragment.contig_index,
    )
    try:
        return _native_sample_bernoulli_sites(
            fragment.methylation_sites,
            key,
            fragment.fragment_ordinal,
        )
    except (TypeError, ValueError, OverflowError) as error:
        raise ProcessError(
            "native Bernoulli site sampling failed: {}".format(error)
        ) from error


def _materialize_site_states(
    fragment: Fragment,
    sampled: Sequence[bool | int],
) -> tuple[SiteState, ...]:
    return tuple(
        SiteState(
            site.site_index,
            site.template_offset,
            site.reference_pos,
            site.context,
            site.methylation_source,
            site.allele,
            bool(state),
            float(site.methylation_probability),
        )
        for site, state in zip(fragment.methylation_sites, sampled)
    )


def _sample_site_states(
    batch: ColumnarFragmentBatch,
    config: ProcessConfig,
) -> np.ndarray:
    return config.methylation_model.sample_batch(batch, config)


def _sample_site_states_bernoulli(
    batch: ColumnarFragmentBatch,
    config: ProcessConfig,
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
            batch.array(model.fragment_ordinal_bytes, _U8)[site_fragment],
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
            batch.array(model.fragment_ordinal_bytes, _U8)[site_fragment[positions]],
            site_indices[positions],
        )
        result[positions] = _bernoulli_from_pair(
            pair0,
            probabilities[positions],
        )
    result.flags.writeable = False
    return result


assign_methylation = _sample_methylation_values
assign_methylation_batch = _sample_methylation_batch_values

__all__ = [
    "BernoulliStateModel",
    "MethylationStateModel",
    "assign_methylation",
    "assign_methylation_batch",
]
