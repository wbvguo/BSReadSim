"""Test-only scalar fragment processing and uniform configuration helpers."""

from __future__ import annotations
from bsreadsim.rng import _bernoulli_unchecked as bernoulli

import bsreadsim.process.fragment as _fragment
from bsreadsim.process import (
    Fragment,
    ProcessConfig,
    ProcessedFragment,
    UniformError,
    UniformQuality,
)
from bsreadsim.rng import RNGStage, derive_key


class UniformProcessConfig(ProcessConfig):
    """Concise uniform policy constructor for tests."""

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


def process_fragment(
    fragment: Fragment,
    contig_name: str,
    config: ProcessConfig,
    *,
    compact_base_states: bool = False,
    include_details: bool = True,
) -> ProcessedFragment:
    """Scalar reference path retained only for stage-level tests."""

    _fragment._validate_fragment_batch_request(
        (fragment,),
        (contig_name,),
        config,
        compact_base_states=compact_base_states,
        include_details=include_details,
    )
    key = derive_key(
        config.master_seed,
        RNGStage.SITE_STATE,
        fragment.contig_index,
    )
    sampled = tuple(
        bernoulli(
            key,
            fragment.fragment_ordinal,
            site.site_index,
            site.methylation_probability,
        )
        for site in fragment.methylation_sites
    )
    return _fragment._process_fragment_with_states(
        fragment,
        contig_name,
        config,
        sampled,
        compact_base_states=compact_base_states,
        include_details=include_details,
    )


__all__ = ["UniformProcessConfig", "process_fragment"]
