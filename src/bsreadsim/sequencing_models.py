"""Versioned quality and sequencing-error model components.

The parsers in this module consume already snapshotted bytes.  They perform no
filesystem I/O and reject executable serialization formats such as pickle.
Sampling is counter-addressed through the shared Philox contract, so worker
count, batch size, and completion order cannot alter a result.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
import json
from typing import Dict, Mapping, Sequence, Tuple

from .rng import RNGStage, bounded_integer, derive_key


QUALITY_MARKOV_FORMAT = "json"
QUALITY_MARKOV_VERSION = "quality-markov-v1"
QUALITY_MARKOV_SCHEMA = "bsreadsim-quality-markov-v1"
QUALITY_CONFUSION_FORMAT = "json"
QUALITY_CONFUSION_VERSION = "quality-confusion-v1"
QUALITY_CONFUSION_SCHEMA = "bsreadsim-quality-confusion-v1"
QUALITY_INITIAL_CYCLES = 5
MAX_MODEL_BYTES = 8 * 1024 * 1024

_UINT32_MAX = (1 << 32) - 1
_UINT64_MAX = (1 << 64) - 1
_BASE_COUNT = 4


class SequencingModelError(ValueError):
    """A model artifact or sampling request violates the v1 contract."""


@dataclass(frozen=True)
class WeightedRow:
    """One exact integer-weight categorical distribution."""

    cumulative_counts: Tuple[int, ...]
    total_count: int

    def __post_init__(self) -> None:
        if not isinstance(self.cumulative_counts, tuple) or not self.cumulative_counts:
            raise SequencingModelError(
                "weighted row must contain an immutable nonempty count tuple"
            )
        previous = 0
        for cumulative in self.cumulative_counts:
            if (
                isinstance(cumulative, bool)
                or not isinstance(cumulative, int)
                or cumulative < previous
                or cumulative - previous > _UINT32_MAX
                or cumulative > _UINT64_MAX
            ):
                raise SequencingModelError("weighted row cumulative counts are invalid")
            previous = cumulative
        if (
            isinstance(self.total_count, bool)
            or not isinstance(self.total_count, int)
            or self.total_count <= 0
            or self.total_count != previous
        ):
            raise SequencingModelError("weighted row total is invalid")


@dataclass(frozen=True)
class QualityMateModel:
    """The five initial distributions and one transition matrix for a mate."""

    initial: Tuple[WeightedRow, ...]
    transitions: Tuple[WeightedRow, ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.initial, tuple)
            or len(self.initial) != QUALITY_INITIAL_CYCLES
            or not all(isinstance(row, WeightedRow) for row in self.initial)
            or not isinstance(self.transitions, tuple)
            or not self.transitions
            or not all(isinstance(row, WeightedRow) for row in self.transitions)
        ):
            raise SequencingModelError("quality mate model shape is invalid")


@dataclass(frozen=True)
class QualityMarkovModel:
    """R1/R2 quality Markov model decoded from a strict JSON artifact."""

    quality_scores: Tuple[int, ...]
    mates: Tuple[QualityMateModel, QualityMateModel]

    def __post_init__(self) -> None:
        _validate_frozen_quality_scores(self.quality_scores)
        if (
            not isinstance(self.mates, tuple)
            or len(self.mates) != 2
            or not all(isinstance(mate, QualityMateModel) for mate in self.mates)
        ):
            raise SequencingModelError("quality model must contain R1 and R2")
        state_count = len(self.quality_scores)
        for mate in self.mates:
            if len(mate.transitions) != state_count or any(
                len(row.cumulative_counts) != state_count
                for row in mate.initial + mate.transitions
            ):
                raise SequencingModelError("quality model matrix shape is invalid")

    def sample(
        self,
        master_seed: int,
        contig_index: int,
        fragment_ordinal: int,
        mate_index: int,
        read_length: int,
    ) -> Tuple[int, ...]:
        """Sample one quality sequence with stable per-cycle RNG addresses."""
        _require_mate_index(mate_index)
        _require_read_length(read_length)
        _require_u64("fragment_ordinal", fragment_ordinal)
        key = _sampling_key(master_seed, RNGStage.QUALITY, contig_index)
        mate = self.mates[mate_index]
        state_index = 0
        qualities = []
        for read_offset in range(read_length):
            row = (
                mate.initial[read_offset]
                if read_offset < QUALITY_INITIAL_CYCLES
                else mate.transitions[state_index]
            )
            state_index = _sample_row(
                row,
                key,
                fragment_ordinal,
                _pack_mate_offset(mate_index, read_offset),
            )
            qualities.append(self.quality_scores[state_index])
        return tuple(qualities)


@dataclass(frozen=True)
class ConfusionMateModel:
    """Quality-indexed 4x4 base-call distributions for one mate."""

    matrices: Tuple[Tuple[WeightedRow, ...], ...]

    def __post_init__(self) -> None:
        if (
            not isinstance(self.matrices, tuple)
            or not self.matrices
            or any(
                not isinstance(matrix, tuple)
                or len(matrix) != _BASE_COUNT
                or any(
                    not isinstance(row, WeightedRow)
                    or len(row.cumulative_counts) != _BASE_COUNT
                    for row in matrix
                )
                for matrix in self.matrices
            )
        ):
            raise SequencingModelError("confusion mate model shape is invalid")


@dataclass(frozen=True)
class QualityConfusionModel:
    """R1/R2 quality-specific base transition model."""

    quality_scores: Tuple[int, ...]
    mates: Tuple[ConfusionMateModel, ConfusionMateModel]

    def __post_init__(self) -> None:
        _validate_frozen_quality_scores(self.quality_scores)
        if (
            not isinstance(self.mates, tuple)
            or len(self.mates) != 2
            or not all(isinstance(mate, ConfusionMateModel) for mate in self.mates)
            or any(
                len(mate.matrices) != len(self.quality_scores)
                for mate in self.mates
            )
        ):
            raise SequencingModelError("confusion model matrix shape is invalid")

    def sample(
        self,
        master_seed: int,
        contig_index: int,
        fragment_ordinal: int,
        mate_index: int,
        bases: Sequence[int],
        qualities: Sequence[int],
    ) -> Tuple[int, ...]:
        """Sample final A/C/G/T calls; ambiguous N bases remain unchanged."""
        _require_mate_index(mate_index)
        _require_u64("fragment_ordinal", fragment_ordinal)
        if len(bases) != len(qualities):
            raise SequencingModelError(
                "base and quality sequences must have equal length"
            )
        if not 1 <= len(bases) <= _UINT32_MAX:
            raise SequencingModelError("read length must be in [1, 2**32)")
        key = _sampling_key(
            master_seed,
            RNGStage.SEQUENCING_ERROR,
            contig_index,
        )
        score_indices = {
            score: index for index, score in enumerate(self.quality_scores)
        }
        result = []
        for read_offset, (base, quality) in enumerate(zip(bases, qualities)):
            _require_base(base)
            _require_quality_value(quality)
            try:
                quality_index = score_indices[quality]
            except KeyError as error:
                raise SequencingModelError(
                    "quality score is absent from the confusion model"
                ) from error
            if base == 4:
                result.append(4)
                continue
            row = self.mates[mate_index].matrices[quality_index][base]
            result.append(
                _sample_row(
                    row,
                    key,
                    fragment_ordinal,
                    _pack_mate_offset(mate_index, read_offset),
                )
            )
        return tuple(result)


def parse_quality_markov(payload: bytes) -> QualityMarkovModel:
    """Decode one strict ``bsreadsim-quality-markov-v1`` JSON document."""
    document = _parse_document(payload)
    _require_keys(
        document,
        {"schema", "quality_scores", "mates"},
        "$",
    )
    if document["schema"] != QUALITY_MARKOV_SCHEMA:
        raise SequencingModelError("quality model schema is unsupported")
    scores = _parse_quality_scores(document["quality_scores"])
    mate_values = _require_list(document["mates"], "$.mates")
    if len(mate_values) != 2:
        raise SequencingModelError("quality model must contain exactly R1 and R2")
    mates = tuple(
        _parse_quality_mate(value, len(scores), "$.mates[{}]".format(index))
        for index, value in enumerate(mate_values)
    )
    return QualityMarkovModel(scores, (mates[0], mates[1]))


def parse_quality_confusion(payload: bytes) -> QualityConfusionModel:
    """Decode one strict ``bsreadsim-quality-confusion-v1`` JSON document."""
    document = _parse_document(payload)
    _require_keys(
        document,
        {"schema", "quality_scores", "mates"},
        "$",
    )
    if document["schema"] != QUALITY_CONFUSION_SCHEMA:
        raise SequencingModelError("confusion model schema is unsupported")
    scores = _parse_quality_scores(document["quality_scores"])
    mate_values = _require_list(document["mates"], "$.mates")
    if len(mate_values) != 2:
        raise SequencingModelError(
            "confusion model must contain exactly R1 and R2"
        )
    mates = tuple(
        _parse_confusion_mate(value, len(scores), "$.mates[{}]".format(index))
        for index, value in enumerate(mate_values)
    )
    return QualityConfusionModel(scores, (mates[0], mates[1]))


def _parse_quality_mate(
    value: object,
    state_count: int,
    path: str,
) -> QualityMateModel:
    mapping = _require_mapping(value, path)
    _require_keys(mapping, {"initial_counts", "transition_counts"}, path)
    initial_values = _require_list(
        mapping["initial_counts"], path + ".initial_counts"
    )
    if len(initial_values) != QUALITY_INITIAL_CYCLES:
        raise SequencingModelError(
            "{} must contain exactly five cycle rows".format(
                path + ".initial_counts"
            )
        )
    initial = tuple(
        _parse_weight_row(
            row,
            state_count,
            "{}.initial_counts[{}]".format(path, index),
        )
        for index, row in enumerate(initial_values)
    )
    transition_values = _require_list(
        mapping["transition_counts"], path + ".transition_counts"
    )
    if len(transition_values) != state_count:
        raise SequencingModelError(
            "{} must contain one row per quality score".format(
                path + ".transition_counts"
            )
        )
    transitions = tuple(
        _parse_weight_row(
            row,
            state_count,
            "{}.transition_counts[{}]".format(path, index),
        )
        for index, row in enumerate(transition_values)
    )
    return QualityMateModel(initial=initial, transitions=transitions)


def _parse_confusion_mate(
    value: object,
    state_count: int,
    path: str,
) -> ConfusionMateModel:
    mapping = _require_mapping(value, path)
    _require_keys(mapping, {"base_transition_counts"}, path)
    matrix_values = _require_list(
        mapping["base_transition_counts"], path + ".base_transition_counts"
    )
    if len(matrix_values) != state_count:
        raise SequencingModelError(
            "{} must contain one matrix per quality score".format(
                path + ".base_transition_counts"
            )
        )
    matrices = []
    for quality_index, matrix_value in enumerate(matrix_values):
        rows = _require_list(
            matrix_value,
            "{}.base_transition_counts[{}]".format(path, quality_index),
        )
        if len(rows) != _BASE_COUNT:
            raise SequencingModelError(
                "every confusion matrix must contain four source-base rows"
            )
        matrices.append(
            tuple(
                _parse_weight_row(
                    row,
                    _BASE_COUNT,
                    "{}.base_transition_counts[{}][{}]".format(
                        path, quality_index, base
                    ),
                )
                for base, row in enumerate(rows)
            )
        )
    return ConfusionMateModel(matrices=tuple(matrices))


def _parse_quality_scores(value: object) -> Tuple[int, ...]:
    values = _require_list(value, "$.quality_scores")
    if not 1 <= len(values) <= 94:
        raise SequencingModelError(
            "quality_scores must contain between 1 and 94 states"
        )
    scores = []
    previous = -1
    for index, score in enumerate(values):
        _require_quality_value(score, "$.quality_scores[{}]".format(index))
        if score <= previous:
            raise SequencingModelError(
                "quality_scores must be strictly increasing and unique"
            )
        scores.append(score)
        previous = score
    return tuple(scores)


def _validate_frozen_quality_scores(value: object) -> None:
    if not isinstance(value, tuple) or not 1 <= len(value) <= 94:
        raise SequencingModelError(
            "quality_scores must be an immutable tuple with 1 to 94 states"
        )
    previous = -1
    for score in value:
        _require_quality_value(score)
        if score <= previous:
            raise SequencingModelError(
                "quality_scores must be strictly increasing and unique"
            )
        previous = score


def _parse_weight_row(value: object, width: int, path: str) -> WeightedRow:
    values = _require_list(value, path)
    if len(values) != width:
        raise SequencingModelError(
            "{} must contain exactly {} counts".format(path, width)
        )
    cumulative = []
    total = 0
    for index, count in enumerate(values):
        if (
            isinstance(count, bool)
            or not isinstance(count, int)
            or not 0 <= count <= _UINT32_MAX
        ):
            raise SequencingModelError(
                "{}[{}] must be an unsigned 32-bit integer".format(path, index)
            )
        if count > _UINT64_MAX - total:
            raise SequencingModelError("{} total exceeds uint64".format(path))
        total += count
        cumulative.append(total)
    if total == 0:
        raise SequencingModelError("{} must have positive total count".format(path))
    return WeightedRow(tuple(cumulative), total)


def _sample_row(
    row: WeightedRow,
    key: int,
    fragment_ordinal: int,
    local_index: int,
) -> int:
    draw = bounded_integer(
        key,
        fragment_ordinal,
        local_index,
        0,
        row.total_count,
    )
    index = bisect_right(row.cumulative_counts, draw)
    if index >= len(row.cumulative_counts):
        raise SequencingModelError("categorical draw escaped its row")
    return index


def _parse_document(payload: bytes) -> Mapping[str, object]:
    if not isinstance(payload, bytes):
        raise SequencingModelError("model payload must be bytes")
    if len(payload) == 0 or len(payload) > MAX_MODEL_BYTES:
        raise SequencingModelError("model payload size is outside the v1 limit")
    try:
        text = payload.decode("utf-8")
    except UnicodeDecodeError as error:
        raise SequencingModelError("model payload must be UTF-8 JSON") from error

    def reject_constant(value: str) -> object:
        raise SequencingModelError(
            "non-finite JSON number is forbidden: {}".format(value)
        )

    def unique_pairs(pairs: Sequence[Tuple[str, object]]) -> Dict[str, object]:
        result = {}  # type: Dict[str, object]
        for key, value in pairs:
            if key in result:
                raise SequencingModelError(
                    "duplicate JSON object key: {}".format(key)
                )
            result[key] = value
        return result

    try:
        value = json.loads(
            text,
            object_pairs_hook=unique_pairs,
            parse_constant=reject_constant,
        )
    except SequencingModelError:
        raise
    except (TypeError, ValueError, json.JSONDecodeError) as error:
        raise SequencingModelError("model payload is not strict JSON") from error
    return _require_mapping(value, "$")


def _require_mapping(value: object, path: str) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise SequencingModelError("{} must be an object".format(path))
    return value


def _require_list(value: object, path: str) -> Sequence[object]:
    if not isinstance(value, list):
        raise SequencingModelError("{} must be an array".format(path))
    return value


def _require_keys(
    value: Mapping[str, object], expected: set, path: str
) -> None:
    observed = set(value)
    if observed != expected:
        missing = sorted(expected - observed)
        unexpected = sorted(observed - expected)
        details = []
        if missing:
            details.append("missing {}".format(", ".join(missing)))
        if unexpected:
            details.append("unexpected {}".format(", ".join(unexpected)))
        raise SequencingModelError(
            "{} has invalid fields: {}".format(path, "; ".join(details))
        )


def _sampling_key(
    master_seed: int,
    stage: RNGStage,
    contig_index: int,
) -> int:
    _require_u32("contig_index", contig_index)
    try:
        return derive_key(master_seed, stage, contig_index)
    except ValueError as error:
        raise SequencingModelError(str(error)) from error


def _require_mate_index(value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or value not in (0, 1):
        raise SequencingModelError("mate_index must be 0 or 1")


def _require_read_length(value: int) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 1 <= value <= _UINT32_MAX
    ):
        raise SequencingModelError("read_length must be in [1, 2**32)")


def _require_u32(name: str, value: object) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 0 <= value <= _UINT32_MAX
    ):
        raise SequencingModelError(
            "{} must be an unsigned 32-bit integer".format(name)
        )


def _require_u64(name: str, value: object) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 0 <= value <= _UINT64_MAX
    ):
        raise SequencingModelError(
            "{} must be an unsigned 64-bit integer".format(name)
        )


def _require_quality_value(value: object, path: str = "quality") -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 0 <= value <= 93
    ):
        raise SequencingModelError(
            "{} must be an integer in [0, 93]".format(path)
        )


def _require_base(value: object) -> None:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or not 0 <= value <= 4
    ):
        raise SequencingModelError("base must use protocol encoding [0, 4]")


def _pack_mate_offset(mate_index: int, read_offset: int) -> int:
    return (mate_index << 32) | read_offset


__all__ = [
    "MAX_MODEL_BYTES",
    "QUALITY_CONFUSION_FORMAT",
    "QUALITY_CONFUSION_SCHEMA",
    "QUALITY_CONFUSION_VERSION",
    "QUALITY_INITIAL_CYCLES",
    "QUALITY_MARKOV_FORMAT",
    "QUALITY_MARKOV_SCHEMA",
    "QUALITY_MARKOV_VERSION",
    "QualityConfusionModel",
    "QualityMarkovModel",
    "SequencingModelError",
    "parse_quality_confusion",
    "parse_quality_markov",
]
