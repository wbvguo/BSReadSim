"""Quality Markov and quality-specific error models drive base calls."""

import json
import unittest

from tests.semantic import support

READ_LENGTH = 30
LOW, HIGH = 20, 30
# Deterministic models: R1 alternates after five fixed initial cycles,
# R2 always reports Q20. Q20 calls rotate A->C->G->T->A; Q30 calls are exact.
R1_INITIAL = (HIGH, LOW, HIGH, HIGH, LOW)
ROTATE = {"A": "C", "C": "G", "G": "T", "T": "A", "N": "N"}


def one_hot(quality):
    return [int(quality == LOW), int(quality == HIGH)]


QUALITY_MODEL = {
    "schema": "quality-markov",
    "quality_scores": [LOW, HIGH],
    "mates": [
        {
            "initial_counts": [one_hot(q) for q in R1_INITIAL],
            "transition_counts": [one_hot(HIGH), one_hot(LOW)],
        },
        {
            "initial_counts": [one_hot(LOW)] * 5,
            "transition_counts": [one_hot(LOW), one_hot(HIGH)],
        },
    ],
}
IDENTITY = [[int(i == j) for j in range(4)] for i in range(4)]
ROTATION = [[int(j == (i + 1) % 4) for j in range(4)] for i in range(4)]
ERROR_MODEL = {
    "schema": "quality-confusion",
    "quality_scores": [LOW, HIGH],
    "mates": [{"base_transition_counts": [ROTATION, IDENTITY]}] * 2,
}


def expected_qualities(mate):
    if mate == 2:
        return [LOW] * READ_LENGTH
    qualities = list(R1_INITIAL)
    while len(qualities) < READ_LENGTH:
        qualities.append(HIGH if qualities[-1] == LOW else LOW)
    return qualities


def models_run():
    quality = support.write_text("quality.json", json.dumps(QUALITY_MODEL))
    error = support.write_text("error.json", json.dumps(ERROR_MODEL))
    return support.simulate(
        "wgs-models", "wgs", "-n", 4_000, "-l", READ_LENGTH,
        "--insert-mean", 200, "--insert-sd", 20, "--insert-min", 100,
        "--insert-max", 300, "--mutation-rate", 0,
        "--quality-model", quality, "--error-model", error, "--format", "bam",
    )


class SequencingModelTests(unittest.TestCase):
    def test_quality_strings_follow_the_markov_model(self):
        for record in models_run().bam.records:
            cycles = list(record.quality)
            if record.is_reverse:
                cycles.reverse()
            self.assertEqual(cycles, expected_qualities(record.mate_number))

    def test_base_calls_follow_the_quality_specific_matrix(self):
        run = models_run()
        for record in run.bam.records:
            truth = support.GENOME[run.contig_of(record)][
                record.position:record.reference_end]
            called = record.sequence
            states = record.tags["zt"]
            if record.is_reverse:
                truth = support.reverse_complement(truth)
                called = support.reverse_complement(called)
                states = states[::-1]
            for cycle, quality in enumerate(expected_qualities(record.mate_number)):
                error = bool(support.ZT_VALUES[states[cycle]] & 32)
                if quality == LOW and truth[cycle] != "N":
                    self.assertEqual(called[cycle], ROTATE[truth[cycle]])
                    self.assertTrue(error)
                else:
                    self.assertEqual(called[cycle], truth[cycle])
                    self.assertFalse(error)


if __name__ == "__main__":
    unittest.main()
