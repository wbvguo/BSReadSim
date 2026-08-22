"""Strict artifacts and deterministic vectors for Python sequencing models."""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import FrozenInstanceError
import copy
import json
from pathlib import Path
import sys
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim.process.sequencing import (  # noqa: E402
    MAX_MODEL_BYTES,
    QUALITY_CONFUSION_FORMAT,
    QUALITY_CONFUSION_SCHEMA,
    QUALITY_CONFUSION_VERSION,
    QUALITY_MARKOV_FORMAT,
    QUALITY_MARKOV_SCHEMA,
    QUALITY_MARKOV_VERSION,
    QualityMarkovModel,
    SequencingModelError,
    parse_quality_confusion,
    parse_quality_markov,
)


def _quality_document() -> dict:
    mate_1 = {
        "initial_counts": [
            [1, 3, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 0],
            [0, 1, 0],
        ],
        "transition_counts": [
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 0],
        ],
    }
    mate_2 = {
        "initial_counts": [
            [0, 0, 1],
            [1, 0, 0],
            [0, 1, 0],
            [0, 0, 1],
            [1, 0, 0],
        ],
        "transition_counts": [
            [0, 0, 1],
            [1, 0, 0],
            [0, 1, 0],
        ],
    }
    return {
        "schema": QUALITY_MARKOV_SCHEMA,
        "quality_scores": [10, 20, 30],
        "mates": [mate_1, mate_2],
    }


def _confusion_document() -> dict:
    cyclic = [
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
    ]
    identity = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
    weighted = [
        [1, 2, 3, 4],
        [4, 3, 2, 1],
        [1, 4, 1, 4],
        [4, 1, 4, 1],
    ]
    reverse_cyclic = [row[::-1] for row in cyclic]
    return {
        "schema": QUALITY_CONFUSION_SCHEMA,
        "quality_scores": [10, 20, 30],
        "mates": [
            {"base_transition_counts": [cyclic, identity, weighted]},
            {
                "base_transition_counts": [
                    reverse_cyclic,
                    identity,
                    weighted,
                ]
            },
        ],
    }


def _payload(document: dict) -> bytes:
    return json.dumps(
        document,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


class QualityMarkovModelTests(unittest.TestCase):
    def test_contract_identifiers_are_frozen(self) -> None:
        self.assertEqual(QUALITY_MARKOV_FORMAT, "json")
        self.assertEqual(QUALITY_MARKOV_VERSION, "quality-markov-v1")
        self.assertEqual(QUALITY_MARKOV_SCHEMA, "quality-markov-v1")

    def test_first_five_cycles_then_transition_have_a_frozen_vector(self) -> None:
        model = parse_quality_markov(_payload(_quality_document()))
        self.assertEqual(model.quality_scores, (10, 20, 30))
        self.assertEqual(
            model.sample(7, 0, 11, 0, 9),
            (20, 20, 30, 10, 20, 30, 10, 20, 30),
        )
        self.assertEqual(
            model.sample(7, 0, 11, 1, 8),
            (30, 10, 20, 30, 10, 30, 20, 10),
        )

    def test_worker_and_completion_order_do_not_change_quality(self) -> None:
        model = parse_quality_markov(_payload(_quality_document()))

        def sample(fragment_ordinal: int):
            return model.sample(19, 1, fragment_ordinal, 0, 17)

        expected = tuple(sample(ordinal) for ordinal in range(64))
        with ThreadPoolExecutor(max_workers=4) as executor:
            observed = tuple(
                reversed(
                    tuple(executor.map(sample, reversed(range(64))))
                )
            )
        self.assertEqual(observed, expected)

    def test_integer_counts_define_the_expected_initial_distribution(self) -> None:
        model = parse_quality_markov(_payload(_quality_document()))
        observed = {10: 0, 20: 0, 30: 0}
        for fragment_ordinal in range(10000):
            quality = model.sample(31, 2, fragment_ordinal, 0, 1)[0]
            observed[quality] += 1
        self.assertEqual(observed[30], 0)
        self.assertAlmostEqual(observed[10] / 10000, 0.25, delta=0.02)
        self.assertAlmostEqual(observed[20] / 10000, 0.75, delta=0.02)

    def test_model_and_sampling_boundaries_fail_closed(self) -> None:
        model = parse_quality_markov(_payload(_quality_document()))
        for arguments in (
            (-1, 0, 0, 0, 10),
            (0, -1, 0, 0, 10),
            (0, 1 << 32, 0, 0, 10),
            (0, 0, -1, 0, 10),
            (0, 0, 0, 2, 10),
            (0, 0, 0, 0, 0),
            (0, 0, 0, 0, 1 << 32),
        ):
            with self.subTest(arguments=arguments):
                with self.assertRaises(SequencingModelError):
                    model.sample(*arguments)

        malformed = _quality_document()
        malformed["quality_scores"] = [10, 10]
        with self.assertRaisesRegex(SequencingModelError, "strictly increasing"):
            parse_quality_markov(_payload(malformed))

        malformed = _quality_document()
        malformed["mates"][0]["initial_counts"] = [[1, 0, 0]] * 4
        with self.assertRaisesRegex(SequencingModelError, "five cycle"):
            parse_quality_markov(_payload(malformed))

        malformed = _quality_document()
        malformed["mates"][0]["transition_counts"][1] = [0, 0, 0]
        with self.assertRaisesRegex(SequencingModelError, "positive total"):
            parse_quality_markov(_payload(malformed))

    def test_models_are_immutable_values(self) -> None:
        model = parse_quality_markov(_payload(_quality_document()))
        with self.assertRaises(FrozenInstanceError):
            model.quality_scores = (40,)  # type: ignore[misc]
        with self.assertRaises(SequencingModelError):
            QualityMarkovModel((), model.mates)


class QualityConfusionModelTests(unittest.TestCase):
    def test_contract_identifiers_are_frozen(self) -> None:
        self.assertEqual(QUALITY_CONFUSION_FORMAT, "json")
        self.assertEqual(QUALITY_CONFUSION_VERSION, "quality-confusion-v1")
        self.assertEqual(
            QUALITY_CONFUSION_SCHEMA,
            "quality-confusion-v1",
        )

    def test_quality_specific_calls_have_a_frozen_vector(self) -> None:
        model = parse_quality_confusion(_payload(_confusion_document()))
        self.assertEqual(
            model.sample(
                7,
                0,
                11,
                0,
                (0, 1, 2, 3, 0, 4),
                (10, 20, 10, 20, 30, 30),
            ),
            (1, 1, 3, 3, 1, 4),
        )
        self.assertEqual(
            model.sample(7, 0, 11, 1, (0, 1, 2, 3), (10, 10, 10, 10)),
            (2, 1, 0, 3),
        )

    def test_missing_quality_shape_and_counts_fail_closed(self) -> None:
        model = parse_quality_confusion(_payload(_confusion_document()))
        with self.assertRaisesRegex(SequencingModelError, "absent"):
            model.sample(0, 0, 0, 0, (0,), (11,))
        with self.assertRaisesRegex(SequencingModelError, "equal length"):
            model.sample(0, 0, 0, 0, (0, 1), (10,))
        with self.assertRaisesRegex(SequencingModelError, "read length"):
            model.sample(0, 0, 0, 0, (), ())
        with self.assertRaisesRegex(SequencingModelError, "protocol encoding"):
            model.sample(0, 0, 0, 0, (5,), (10,))

        malformed = _confusion_document()
        malformed["mates"][0]["base_transition_counts"][0] = [[1, 0, 0, 0]] * 3
        with self.assertRaisesRegex(SequencingModelError, "four source-base"):
            parse_quality_confusion(_payload(malformed))

        malformed = _confusion_document()
        malformed["mates"][0]["base_transition_counts"][0][0][0] = True
        with self.assertRaisesRegex(SequencingModelError, "unsigned 32-bit"):
            parse_quality_confusion(_payload(malformed))

    def test_integer_confusion_counts_define_base_call_frequencies(self) -> None:
        model = parse_quality_confusion(_payload(_confusion_document()))
        observed = [0, 0, 0, 0]
        for fragment_ordinal in range(10000):
            called = model.sample(
                31,
                2,
                fragment_ordinal,
                0,
                (0,),
                (30,),
            )[0]
            observed[called] += 1
        for actual, expected in zip(observed, (0.1, 0.2, 0.3, 0.4)):
            self.assertAlmostEqual(actual / 10000, expected, delta=0.02)

        malformed = _confusion_document()
        malformed["mates"][0]["base_transition_counts"][0][0][0] = 1 << 32
        with self.assertRaisesRegex(SequencingModelError, "unsigned 32-bit"):
            parse_quality_confusion(_payload(malformed))


class StrictJsonTests(unittest.TestCase):
    def test_duplicate_nonfinite_invalid_utf8_and_size_are_rejected(self) -> None:
        invalid_payloads = (
            b'{"schema":"x","schema":"y"}',
            b'{"schema":NaN}',
            b"\xff",
            b"",
            b" " * (MAX_MODEL_BYTES + 1),
        )
        for payload in invalid_payloads:
            with self.subTest(size=len(payload)):
                with self.assertRaises(SequencingModelError):
                    parse_quality_markov(payload)

    def test_unknown_missing_schema_and_nonbytes_are_rejected(self) -> None:
        document = _quality_document()
        document["unknown"] = 1
        with self.assertRaisesRegex(SequencingModelError, "unexpected unknown"):
            parse_quality_markov(_payload(document))

        document = _quality_document()
        del document["mates"]
        with self.assertRaisesRegex(SequencingModelError, "missing mates"):
            parse_quality_markov(_payload(document))

        document = _quality_document()
        document["schema"] = "future"
        with self.assertRaisesRegex(SequencingModelError, "unsupported"):
            parse_quality_markov(_payload(document))

        with self.assertRaisesRegex(SequencingModelError, "must be bytes"):
            parse_quality_markov("{}")  # type: ignore[arg-type]

    def test_input_documents_are_not_mutated(self) -> None:
        quality = _quality_document()
        confusion = _confusion_document()
        expected_quality = copy.deepcopy(quality)
        expected_confusion = copy.deepcopy(confusion)
        parse_quality_markov(_payload(quality))
        parse_quality_confusion(_payload(confusion))
        self.assertEqual(quality, expected_quality)
        self.assertEqual(confusion, expected_confusion)


if __name__ == "__main__":
    unittest.main()
