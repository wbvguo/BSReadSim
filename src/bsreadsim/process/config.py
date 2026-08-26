"""Validated configuration shared by all fragment-processing stages."""

from __future__ import annotations

from dataclasses import dataclass, field

from .batch import (
    ProcessError,
    UniformError,
    UniformQuality,
    _require_probability,
)
from .sequencing import QualityConfusionModel, QualityMarkovModel
from .methylation import BernoulliStateModel, MethylationStateModel


_UINT64_MAX = (1 << 64) - 1

QualityPolicy = UniformQuality | QualityMarkovModel
ErrorPolicy = UniformError | QualityConfusionModel


@dataclass(frozen=True)
class ProcessConfig:
    """Validated immutable settings for all released Python read stages."""

    master_seed: int
    directional: bool
    conversion_rate: float
    quality: QualityPolicy
    error: ErrorPolicy
    methylation_model: MethylationStateModel = field(
        default_factory=BernoulliStateModel
    )
    bisulfite: bool = True

    def __post_init__(self) -> None:
        _require_u64("master_seed", self.master_seed)
        if not isinstance(self.directional, bool):
            raise ProcessError("directional must be a boolean")
        if not isinstance(self.bisulfite, bool):
            raise ProcessError("bisulfite must be a boolean")
        _require_probability("conversion_rate", self.conversion_rate)
        if not isinstance(self.quality, (UniformQuality, QualityMarkovModel)):
            raise ProcessError("quality policy is outside the released contract")
        if not isinstance(self.error, (UniformError, QualityConfusionModel)):
            raise ProcessError("error policy is outside the released contract")
        if not isinstance(self.methylation_model, MethylationStateModel):
            raise ProcessError("methylation model is outside the released contract")


def _require_u64(name: str, value: int) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 0 <= value <= _UINT64_MAX
    ):
        raise ProcessError("{} must be an unsigned 64-bit integer".format(name))


__all__ = ["ProcessConfig"]
