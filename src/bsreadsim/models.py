"""Pluggable read-processing models.

Every model works on a whole batch at once and draws only from the NumPy
generator it is given, so results depend on the batch, never on scheduling.

* A methylation model realizes one state per fragment site (the molecule's
  pattern). ``BernoulliMethylation`` draws each site independently;
  ``BiLSTMMethylation`` is the interface of a sequence model that conditions
  each site on its neighbours (no trained network ships yet).
* A quality model produces Phred scores per read cycle.
* An error model turns true base calls into observed ones given qualities.

Read arrays are in sequencing orientation: shape (reads, read_length).
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Sequence
from dataclasses import dataclass
import json
from pathlib import Path
from typing import Protocol, TypeVar
import warnings

import numpy as np

from .batch import N, FragmentBatch
from .errors import BSReadSimError

MAX_MODEL_BYTES = 8 * 1024 * 1024
INITIAL_CYCLES = 5
_UINT32_MAX = (1 << 32) - 1
Model = TypeVar("Model")


class MethylationModel(Protocol):
    name: str  # as for --meth-model
    site_kmers: bool  # whether sample reads batch.site_kmers (htsim sends them only if so)

    def sample(self, batch: FragmentBatch, rng: np.random.Generator) -> np.ndarray:
        """Return bool[S]: whether each site of the batch is methylated."""


class QualityModel(Protocol):
    scores: tuple[int, ...]

    def sample(self, mates: np.ndarray, read_length: int,
               rng: np.random.Generator) -> np.ndarray:
        """Return uint8[reads, read_length] Phred scores; mates are 1 or 2."""


class ErrorModel(Protocol):
    def apply(self, bases: np.ndarray, qualities: np.ndarray, mates: np.ndarray,
              rng: np.random.Generator) -> np.ndarray:
        """Return the observed base codes for true ``bases``."""


class BernoulliMethylation:
    """Each site is methylated independently with its profile probability."""

    name = "bernoulli"
    site_kmers = False

    def sample(self, batch: FragmentBatch, rng: np.random.Generator) -> np.ndarray:
        return rng.random(len(batch.site_probabilities)) < batch.site_probabilities


@dataclass(frozen=True)
class SiteFeatures:
    """What a sequence model sees of each site of a batch, sites in fragment order."""

    offsets: np.ndarray  # int64[n + 1], the sites of fragment i are offsets[i]:offsets[i + 1]
    kmers: np.ndarray  # uint16[S], the forward-strand 7-mer centred on the site (A=0, C=1,
    #                    G=2, T=3, first base highest); 65535 if it has an N or passes a contig end
    contexts: np.ndarray  # uint8[S], CG, CHG, or CHH
    probabilities: np.ndarray  # float32[S], the level of the methylation profile
    gaps: np.ndarray  # int64[S], template bases since the previous site (0 for the first)


def site_features(batch: FragmentBatch) -> SiteFeatures:
    gaps = np.diff(batch.site_positions, prepend=0)
    starts = batch.site_offsets[:-1]
    gaps[starts[starts < batch.site_offsets[1:]]] = 0  # the first site of each fragment
    return SiteFeatures(offsets=batch.site_offsets, kmers=batch.site_kmers,
                        contexts=batch.site_contexts, probabilities=batch.site_probabilities,
                        gaps=gaps)


class BiLSTMMethylation:
    """Methylation states that depend on the neighbouring sites of the same molecule.

    A bidirectional LSTM reads the sites of each fragment in order (``SiteFeatures``)
    and gives each site's probability of being methylated, given the whole fragment;
    ``sample`` draws the states from it as ``BernoulliMethylation`` draws them from the
    profile. ``predict`` is the network. None is trained yet, so ``methylation_model``
    falls back to Bernoulli for ``--meth-model bilstm``.
    """

    name = "bilstm"
    site_kmers = True

    def __init__(self, predict: Callable[[SiteFeatures], np.ndarray]):
        self.predict = predict

    def sample(self, batch: FragmentBatch, rng: np.random.Generator) -> np.ndarray:
        probabilities = np.asarray(self.predict(site_features(batch)), dtype=np.float64)
        if probabilities.shape != (len(batch.site_probabilities),):
            raise BSReadSimError("the BiLSTM must give one probability per site")
        return rng.random(len(probabilities)) < probabilities


def methylation_model(name: str) -> MethylationModel:
    """The methylation model of ``--meth-model``."""
    if name == "bilstm":  # no trained network yet: construct BiLSTMMethylation(predict) here
        warnings.warn("BiLSTM methylation model is not available yet; falling back to Bernoulli",
                      RuntimeWarning, stacklevel=2)
    return BernoulliMethylation()


class UniformQuality:
    def __init__(self, phred: int):
        self.phred = phred
        self.scores = (phred,)

    def sample(self, mates: np.ndarray, read_length: int,
               rng: np.random.Generator) -> np.ndarray:
        return np.full((len(mates), read_length), self.phred, dtype=np.uint8)


class MarkovQuality:
    """Five per-cycle initial distributions, then a first-order Markov chain."""

    def __init__(self, scores: Sequence[int], initial_counts: list, transition_counts: list):
        self.scores = tuple(scores)
        self._score_values = np.array(scores, dtype=np.uint8)
        self._initial = _cumulative(initial_counts)  # [mate, cycle, state]
        self._transitions = _cumulative(transition_counts)  # [mate, from, to]

    @classmethod
    def load(cls, path: Path) -> MarkovQuality:
        """The model of a quality-markov JSON file."""
        return _load_file(path, cls.from_json)

    @classmethod
    def from_json(cls, payload: bytes) -> MarkovQuality:
        document = _load(payload, "quality-markov")
        scores = _scores(document)
        initial, transitions = [], []
        for path, mate in _mates(document, ("initial_counts", "transition_counts")):
            initial.append(_matrix(mate["initial_counts"], INITIAL_CYCLES,
                                   len(scores), path + ".initial_counts"))
            transitions.append(_matrix(mate["transition_counts"], len(scores),
                                       len(scores), path + ".transition_counts"))
        return cls(scores, initial, transitions)

    def sample(self, mates: np.ndarray, read_length: int,
               rng: np.random.Generator) -> np.ndarray:
        mate_index = np.asarray(mates, dtype=np.int64) - 1
        states = np.zeros((len(mate_index), read_length), dtype=np.int64)
        for cycle in range(read_length):
            if cycle < INITIAL_CYCLES:
                rows = self._initial[mate_index, cycle]
            else:
                rows = self._transitions[mate_index, states[:, cycle - 1]]
            states[:, cycle] = _draw(rows, rng)
        return self._score_values[states]


class UniformErrors:
    """Each called base is replaced by one of the other three with ``rate``."""

    def __init__(self, rate: float):
        self.rate = rate

    def apply(self, bases: np.ndarray, qualities: np.ndarray, mates: np.ndarray,
              rng: np.random.Generator) -> np.ndarray:
        error = (rng.random(bases.shape) < self.rate) & (bases != N)
        shift = rng.integers(1, 4, size=bases.shape, dtype=np.uint8)
        return np.where(error, (bases + shift) % 4, bases).astype(np.uint8)


class ConfusionErrors:
    """Observed base drawn from a per-mate, per-quality 4x4 transition matrix."""

    def __init__(self, scores: Sequence[int], matrices: list):
        self.scores = tuple(scores)
        self._cumulative = _cumulative(matrices)  # [mate, quality, true, observed]
        self._index = np.full(256, -1, dtype=np.int64)
        self._index[list(scores)] = np.arange(len(scores))

    @classmethod
    def load(cls, path: Path) -> ConfusionErrors:
        """The model of a quality-confusion JSON file."""
        return _load_file(path, cls.from_json)

    @classmethod
    def from_json(cls, payload: bytes) -> ConfusionErrors:
        document = _load(payload, "quality-confusion")
        scores = _scores(document)
        matrices = []
        for path, mate in _mates(document, ("base_transition_counts",)):
            values = mate["base_transition_counts"]
            if not isinstance(values, list) or len(values) != len(scores):
                raise BSReadSimError(path + ".base_transition_counts needs one matrix per score")
            matrices.append([
                _matrix(matrix, 4, 4, f"{path}.base_transition_counts[{q}]")
                for q, matrix in enumerate(values)
            ])
        return cls(scores, matrices)

    def apply(self, bases: np.ndarray, qualities: np.ndarray, mates: np.ndarray,
              rng: np.random.Generator) -> np.ndarray:
        quality_index = self._index[qualities]
        if np.any(quality_index < 0):
            raise BSReadSimError("error model lacks a quality score used by reads")
        mate_index = np.broadcast_to(
            (np.asarray(mates, dtype=np.int64) - 1)[:, None], bases.shape)
        called = bases != N
        rows = self._cumulative[
            mate_index[called], quality_index[called], bases[called]]
        observed = bases.copy()
        observed[called] = _draw(rows, rng)
        return observed


def _cumulative(counts: list) -> np.ndarray:
    values = np.asarray(counts, dtype=np.float64)
    totals = values.sum(axis=-1, keepdims=True)
    cumulative = np.cumsum(values, axis=-1) / totals
    cumulative[..., -1] = 1.0
    return cumulative


def _draw(cumulative_rows: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """One categorical draw per row of cumulative probabilities."""
    u = rng.random(len(cumulative_rows))
    return (u[:, None] >= cumulative_rows).sum(axis=1)


# ----------------------------------------------------------- JSON documents

def _load_file(path: Path, parse: Callable[[bytes], Model]) -> Model:
    try:
        with path.open("rb") as source:
            return parse(source.read(MAX_MODEL_BYTES + 1))
    except (OSError, BSReadSimError) as error:
        raise BSReadSimError(f"cannot load sequencing model {path}: {error}") from error


def _load(payload: bytes, schema: str) -> dict:
    if not 0 < len(payload) <= MAX_MODEL_BYTES:
        raise BSReadSimError("model file must contain 1 byte to 8 MiB")

    def unique(pairs):
        keys = [key for key, _ in pairs]
        if len(keys) != len(set(keys)):
            raise BSReadSimError("duplicate JSON key in model")
        return dict(pairs)

    def reject(value):
        raise BSReadSimError("non-finite number in model: " + value)

    try:
        document = json.loads(payload.decode("utf-8"), object_pairs_hook=unique,
                              parse_constant=reject)
    except (UnicodeDecodeError, json.JSONDecodeError) as error:
        raise BSReadSimError(f"model is not valid UTF-8 JSON: {error}") from error
    if not isinstance(document, dict) or set(document) != {
            "schema", "quality_scores", "mates"}:
        raise BSReadSimError("model must contain exactly schema, quality_scores, mates")
    if document["schema"] != schema:
        raise BSReadSimError("model schema must be " + schema)
    return document


def _scores(document: dict) -> tuple[int, ...]:
    scores = document["quality_scores"]
    if (
        not isinstance(scores, list) or not scores
        or any(type(s) is not int or not 0 <= s <= 93 for s in scores)
        or any(a >= b for a, b in zip(scores, scores[1:], strict=False))
    ):
        raise BSReadSimError("quality_scores must be strictly increasing integers in [0, 93]")
    return tuple(scores)


def _mates(document: dict, keys: tuple[str, ...]) -> Iterator[tuple[str, dict]]:
    mates = document["mates"]
    if not isinstance(mates, list) or len(mates) != 2:
        raise BSReadSimError("mates must list R1 and R2")
    for index, mate in enumerate(mates):
        path = f"$.mates[{index}]"
        if not isinstance(mate, dict) or set(mate) != set(keys):
            raise BSReadSimError(f"{path} must contain exactly {', '.join(keys)}")
        yield path, mate


def _matrix(rows: object, height: int, width: int, path: str) -> list:
    if not isinstance(rows, list) or len(rows) != height:
        raise BSReadSimError(f"{path} must have {height} rows")
    for row in rows:
        if (
            not isinstance(row, list) or len(row) != width
            or any(type(c) is not int or not 0 <= c <= _UINT32_MAX for c in row)
        ):
            raise BSReadSimError(f"{path} rows must hold {width} uint32 counts")
        if sum(row) == 0:
            raise BSReadSimError(f"{path} has a row with zero total")
    return rows
