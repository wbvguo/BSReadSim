"""Methylation and sequencing models: their documents, their sampling, and their loading
for a run."""

import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from bsreadsim.batch import FragmentBatch
from bsreadsim.errors import BSReadSimError
from bsreadsim.models import (
    BernoulliMethylation,
    BiLSTMMethylation,
    ConfusionErrors,
    MarkovQuality,
    UniformErrors,
    UniformQuality,
    methylation_model,
    site_features,
)
from bsreadsim.reads import ReadSimulator
from bsreadsim.settings import Settings


def quality_document(scores=(20, 30)):
    width = len(scores)
    mate = {
        "initial_counts": [[1] * width for _ in range(5)],
        "transition_counts": [[1] * width for _ in range(width)],
    }
    return {"schema": "quality-markov", "quality_scores": list(scores),
            "mates": [mate, mate]}


def error_document(scores=(20, 30)):
    identity = [[int(i == j) for j in range(4)] for i in range(4)]
    mate = {"base_transition_counts": [identity for _ in scores]}
    return {"schema": "quality-confusion", "quality_scores": list(scores),
            "mates": [mate, mate]}


def encode(document) -> bytes:
    return json.dumps(document).encode()


class ModelDocumentTests(unittest.TestCase):
    def test_valid_documents_parse(self):
        quality = MarkovQuality.from_json(encode(quality_document()))
        self.assertEqual(quality.scores, (20, 30))
        errors = ConfusionErrors.from_json(encode(error_document()))
        self.assertEqual(errors.scores, (20, 30))

    def test_malformed_documents_are_rejected(self):
        cases = {
            "duplicate key": b'{"schema": "quality-markov", "schema": "x"}',
            "non-finite": encode(quality_document()).replace(b"[1, 1]", b"[NaN, 1]", 1),
            "not JSON": b"{",
            "empty": b"",
        }
        wrong_schema = quality_document()
        wrong_schema["schema"] = "quality-confusion"
        cases["wrong schema"] = encode(wrong_schema)
        unknown_field = quality_document()
        unknown_field["extra"] = 1
        cases["unknown field"] = encode(unknown_field)
        zero_row = quality_document()
        zero_row["mates"][0]["transition_counts"][1] = [0, 0]
        cases["zero row"] = encode(zero_row)
        unsorted = quality_document(scores=(30, 20))
        cases["unsorted scores"] = encode(unsorted)
        one_mate = quality_document()
        one_mate["mates"] = one_mate["mates"][:1]
        cases["one mate"] = encode(one_mate)
        short_initial = quality_document()
        short_initial["mates"][1]["initial_counts"].pop()
        cases["four initial rows"] = encode(short_initial)
        for name, payload in cases.items():
            with self.subTest(name), self.assertRaises(BSReadSimError):
                MarkovQuality.from_json(payload)

        bad_matrix = error_document()
        bad_matrix["mates"][0]["base_transition_counts"][0].pop()
        with self.assertRaises(BSReadSimError):
            ConfusionErrors.from_json(encode(bad_matrix))


class SamplingTests(unittest.TestCase):
    def test_markov_chain_follows_deterministic_transitions(self):
        document = quality_document()
        for mate in document["mates"]:
            mate["initial_counts"] = [[0, 1]] * 5
            mate["transition_counts"] = [[0, 1], [1, 0]]  # 20->30, 30->20
        quality = MarkovQuality.from_json(encode(document))
        scores = quality.sample(np.array([1, 2]), 8, np.random.default_rng(1))
        expected = [30] * 5 + [20, 30, 20]
        self.assertEqual(scores.tolist(), [expected, expected])

    def test_confusion_matrix_maps_bases_by_quality(self):
        document = error_document()
        rotate = [[int(j == (i + 1) % 4) for j in range(4)] for i in range(4)]
        for mate in document["mates"]:
            mate["base_transition_counts"][0] = rotate  # Q20 rotates A->C->G->T
        errors = ConfusionErrors.from_json(encode(document))
        bases = np.array([[0, 1, 2, 3, 4]], dtype=np.uint8)
        observed = errors.apply(bases, np.array([[20, 20, 30, 20, 20]], dtype=np.uint8),
                                np.array([1]), np.random.default_rng(1))
        self.assertEqual(observed.tolist(), [[1, 2, 2, 0, 4]])

    def test_uniform_errors_substitute_other_bases_and_keep_n(self):
        bases = np.tile(np.array([0, 1, 2, 3, 4], dtype=np.uint8), (2000, 1))
        observed = UniformErrors(0.2).apply(bases, None, None, np.random.default_rng(3))
        changed = observed != bases
        self.assertFalse(changed[:, 4].any())
        self.assertAlmostEqual(changed[:, :4].mean(), 0.2, delta=0.01)
        self.assertTrue((observed[:, :4] < 4).all())


def site_batch():
    """Three fragments with 3, 0, and 2 sites."""
    empty = np.zeros(0, dtype=np.int64)
    return FragmentBatch(
        first_ordinal=0, contig=np.zeros(3, dtype=np.uint32),
        start=np.zeros(3, dtype=np.int64), end=np.full(3, 10, dtype=np.int64),
        haplotype=np.zeros(3, dtype=np.uint8), strand=np.zeros(3, dtype=np.uint8),
        template_offsets=np.array([0, 10, 20, 30]), bases=np.zeros(30, dtype=np.uint8),
        site_offsets=np.array([0, 3, 3, 5]), site_positions=np.array([1, 4, 9, 2, 7]),
        site_probabilities=np.array([0.1, 0.5, 0.9, 0.3, 0.7], dtype=np.float32),
        site_contexts=np.array([1, 2, 3, 1, 1], dtype=np.uint8),
        site_asm=np.zeros(5, dtype=bool), site_kmers=np.array([5, 6, 7, 8, 65535], dtype=np.uint16),
        event_offsets=np.zeros(4, dtype=np.int64), event_starts=empty, event_ends=empty,
        event_kinds=np.zeros(0, dtype=np.uint8), event_deleted=empty)


class MethylationModelTests(unittest.TestCase):
    def test_site_features_restart_at_each_fragment(self):
        features = site_features(site_batch())
        self.assertEqual(features.gaps.tolist(), [0, 3, 5, 0, 5])
        self.assertEqual(features.offsets.tolist(), [0, 3, 3, 5])
        self.assertEqual(features.kmers.tolist(), [5, 6, 7, 8, 65535])

    def test_bilstm_draws_from_its_predictions(self):
        batch = site_batch()
        profile = BiLSTMMethylation(lambda features: features.probabilities)
        self.assertEqual(profile.sample(batch, np.random.default_rng(3)).tolist(),
                         BernoulliMethylation().sample(batch, np.random.default_rng(3)).tolist())
        certain = BiLSTMMethylation(lambda features: (features.contexts == 1).astype(float))
        self.assertEqual(certain.sample(batch, np.random.default_rng(3)).tolist(),
                         [True, False, False, True, True])
        with self.assertRaises(BSReadSimError):
            BiLSTMMethylation(lambda features: np.zeros(2)).sample(batch, np.random.default_rng(3))

    def test_bilstm_falls_back_to_bernoulli(self):
        with self.assertWarns(RuntimeWarning):
            self.assertIsInstance(methylation_model("bilstm"), BernoulliMethylation)
        self.assertIsInstance(methylation_model("bernoulli"), BernoulliMethylation)


class ModelLoadingTests(unittest.TestCase):
    def test_uniform_and_markov_models(self):
        self.assertIsInstance(ReadSimulator.from_settings(Settings()).quality, UniformQuality)
        mate = {"initial_counts": [[1, 1]] * 5, "transition_counts": [[1, 1], [1, 1]]}
        document = {"schema": "quality-markov", "quality_scores": [10, 30],
                    "mates": [mate, mate]}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "quality.json"
            path.write_text(json.dumps(document))
            loaded = ReadSimulator.from_settings(Settings(quality_model=path, phred=None))
            self.assertIsInstance(loaded.quality, MarkovQuality)

    def test_missing_or_invalid_models_fail(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "error.json"
            path.write_text("{}")
            for settings in (Settings(error_model=path, error_rate=None),
                             Settings(quality_model=Path(directory) / "absent.json")):
                with self.assertRaises(BSReadSimError):
                    ReadSimulator.from_settings(settings)


if __name__ == "__main__":
    unittest.main()
