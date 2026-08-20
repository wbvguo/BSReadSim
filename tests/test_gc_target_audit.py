"""Scientific contracts for the fixed-insert GC-target audit harness."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_ROOT = REPOSITORY_ROOT / "data" / "experiments"
sys.path.insert(0, str(EXPERIMENT_ROOT))

from audit_wgbs_gc_target import (  # noqa: E402
    AuditError,
    _bin_for_gc_counts,
    _candidate_space,
    _distribution_rows,
    _insert_metrics,
    _load_profile,
)


class GcTargetAuditTests(unittest.TestCase):
    def test_single_column_profile_is_strict_and_digest_bound(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            valid_bytes = b"0.1\n0.3\n0.6\n"
            valid = root / "valid.tsv"
            valid.write_bytes(valid_bytes)
            observed = _load_profile(
                valid, hashlib.sha256(valid_bytes).hexdigest()
            )
            np.testing.assert_array_equal(observed, [0.1, 0.3, 0.6])

            legacy_bytes = b"0\t0.1\n1\t1\n"
            legacy = root / "legacy.tsv"
            legacy.write_bytes(legacy_bytes)
            with self.assertRaisesRegex(AuditError, "one strict probability"):
                _load_profile(
                    legacy, hashlib.sha256(legacy_bytes).hexdigest()
                )

            with self.assertRaisesRegex(AuditError, "SHA-256 changed"):
                _load_profile(valid, "0" * 64)

            invalid_sum = root / "invalid-sum.tsv"
            invalid_sum.write_bytes(b"0.1\n0.3\n0.5\n")
            with self.assertRaisesRegex(AuditError, "sum to one"):
                _load_profile(
                    invalid_sum,
                    hashlib.sha256(invalid_sum.read_bytes()).hexdigest(),
                )

    def test_integer_half_up_bins_and_candidate_histogram_match_brute_force(
        self,
    ) -> None:
        counts = np.arange(5, dtype=np.uint32)
        np.testing.assert_array_equal(
            _bin_for_gc_counts(counts, 4, 3),
            [0, 1, 1, 2, 2],
        )

        with tempfile.TemporaryDirectory() as temporary:
            reference = Path(temporary) / "tiny.fa"
            reference.write_bytes(b">chrGC\nAAAACCCCGGGGTTTT\n")
            candidates = _candidate_space(
                reference,
                {
                    "paired_end": True,
                    "read_length_1": 4,
                    "read_length_2": 4,
                    "insert_min": 4,
                    "insert_mean": 4,
                    "insert_max": 4,
                    "insert_stddev": 0,
                    "max_ambiguous_fraction": 0,
                },
                3,
                chunk_starts=3,
            )

        self.assertEqual(candidates.contig_name, "chrGC")
        self.assertEqual(candidates.valid_start_count, 13)
        np.testing.assert_array_equal(candidates.candidate_counts, [2, 4, 7])
        np.testing.assert_array_equal(
            candidates.bin_by_start,
            [0, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 0],
        )

    def test_variable_insert_uses_mean_proxy_and_actual_length_bins(self) -> None:
        counts = np.asarray([0, 2, 4], dtype=np.uint32)
        lengths = np.asarray([4, 5, 8], dtype=np.uint32)
        np.testing.assert_array_equal(
            _bin_for_gc_counts(counts, lengths, 5),
            [0, 2, 2],
        )

        with tempfile.TemporaryDirectory() as temporary:
            reference = Path(temporary) / "tiny.fa"
            reference.write_bytes(b">chrGC\nAAAACCCCGGGGTTTT\n")
            candidates = _candidate_space(
                reference,
                {
                    "paired_end": True,
                    "read_length_1": 3,
                    "read_length_2": 3,
                    "insert_min": 3,
                    "insert_mean": 4,
                    "insert_max": 8,
                    "insert_stddev": 1,
                    "max_ambiguous_fraction": 0,
                },
                3,
                chunk_starts=3,
            )

        self.assertEqual(candidates.calibration_length, 4)
        np.testing.assert_array_equal(candidates.candidate_counts, [2, 4, 7])

        rows, metrics = _insert_metrics(
            3,
            np.asarray([0, 25, 50, 25], dtype=np.uint64),
            np.asarray([0, 20, 60, 20], dtype=np.uint64),
        )
        self.assertEqual([row["insert_length"] for row in rows], [4, 5, 6])
        self.assertAlmostEqual(
            metrics["profile_vs_uniform_insert_overlap"], 0.9
        )
        self.assertAlmostEqual(metrics["profile_insert_mean"], 5.0)
        self.assertAlmostEqual(metrics["uniform_insert_mean"], 5.0)

    def test_target_overlap_and_calibrated_acceptance(self) -> None:
        probabilities = np.asarray([0.2, 0.3, 0.5])
        candidates = np.asarray([20, 40, 70], dtype=np.uint64)
        observed = np.asarray([20, 30, 50], dtype=np.uint64)

        rows, metrics = _distribution_rows(
            probabilities, candidates, observed
        )

        np.testing.assert_allclose(
            [row["expected_fragment_probability"] for row in rows],
            probabilities,
        )
        np.testing.assert_allclose(
            [row["calibrated_acceptance_probability"] for row in rows],
            [1.0, 0.75, 5.0 / 7.0],
        )
        self.assertAlmostEqual(metrics["target_vs_observed_overlap"], 1.0)
        self.assertLess(metrics["candidate_vs_observed_overlap"], 1.0)

        with self.assertRaisesRegex(AuditError, "no eligible fragment start"):
            _distribution_rows(
                np.asarray([0.5, 0.5]),
                np.asarray([10, 0], dtype=np.uint64),
                np.asarray([10, 0], dtype=np.uint64),
            )

        projected_rows, projected_metrics = _distribution_rows(
            np.asarray([0.9, 0.1]),
            np.asarray([10, 0], dtype=np.uint64),
            np.asarray([100, 0], dtype=np.uint64),
            project_unreachable=True,
        )
        self.assertEqual(
            [row["expected_fragment_probability"] for row in projected_rows],
            [1.0, 0.0],
        )
        self.assertAlmostEqual(
            projected_metrics["dropped_target_probability"], 0.1
        )
        self.assertAlmostEqual(
            projected_metrics["target_vs_observed_overlap"], 0.9
        )
        self.assertAlmostEqual(
            projected_metrics["projected_target_vs_observed_overlap"], 1.0
        )

    def test_retained_chr21_result_recomputes_from_all_100_bins(self) -> None:
        profile = np.asarray(
            [
                float(line)
                for line in (
                    EXPERIMENT_ROOT / "wgbs-gc-target-mock.tsv"
                ).read_text(encoding="ascii").splitlines()
            ]
        )
        with (
            EXPERIMENT_ROOT / "chr21-30x-gc-target-v1.tsv"
        ).open("r", encoding="utf-8", newline="") as input_file:
            reader = csv.DictReader(input_file, delimiter="\t")
            retained_rows = list(reader)
        summary = json.loads(
            (
                EXPERIMENT_ROOT / "chr21-30x-gc-target-v1.json"
            ).read_text(encoding="utf-8")
        )

        self.assertEqual(len(profile), 100)
        self.assertAlmostEqual(float(np.sum(profile)), 1.0, places=12)
        self.assertEqual(len(retained_rows), 100)
        self.assertEqual(
            reader.fieldnames,
            [
                "bin",
                "gc_ratio_center",
                "input_target_probability",
                "normalized_target_probability",
                "eligible_start_count",
                "eligible_start_probability",
                "calibrated_acceptance_probability",
                "expected_fragment_probability",
                "observed_fragment_count",
                "observed_fragment_probability",
                "mean_fragments_per_eligible_start",
                "relative_bin_coverage",
                "observed_minus_expected",
            ],
        )
        np.testing.assert_array_equal(
            [int(row["bin"]) for row in retained_rows], np.arange(100)
        )
        np.testing.assert_allclose(
            [float(row["gc_ratio_center"]) for row in retained_rows],
            np.arange(100) / 99.0,
            rtol=0,
            atol=5e-13,
        )
        np.testing.assert_allclose(
            [float(row["input_target_probability"]) for row in retained_rows],
            profile,
            rtol=0,
            atol=0,
        )
        candidates = np.asarray(
            [int(row["eligible_start_count"]) for row in retained_rows],
            dtype=np.uint64,
        )
        observed = np.asarray(
            [int(row["observed_fragment_count"]) for row in retained_rows],
            dtype=np.uint64,
        )
        expected_rows, recomputed = _distribution_rows(
            profile, candidates, observed
        )

        self.assertEqual(int(candidates.sum()), 40075248)
        self.assertEqual(int(observed.sum()), 4670998)
        for retained, expected in zip(retained_rows, expected_rows):
            for field in (
                "normalized_target_probability",
                "eligible_start_probability",
                "calibrated_acceptance_probability",
                "expected_fragment_probability",
                "observed_fragment_probability",
                "observed_minus_expected",
            ):
                self.assertAlmostEqual(
                    float(retained[field]), expected[field], places=11
                )
            if expected["eligible_start_count"] == 0:
                self.assertEqual(
                    retained["mean_fragments_per_eligible_start"], "NA"
                )
                self.assertEqual(retained["relative_bin_coverage"], "NA")
            else:
                self.assertAlmostEqual(
                    float(retained["mean_fragments_per_eligible_start"]),
                    expected["mean_fragments_per_eligible_start"],
                    places=11,
                )
                self.assertAlmostEqual(
                    float(retained["relative_bin_coverage"]),
                    expected["relative_bin_coverage"],
                    places=10,
                )
        self.assertTrue(summary["accounting"]["stream_matches_fastq_run"])
        self.assertEqual(
            summary["accounting"]["stream_sha256"],
            "5ee3ccac0260328cc579c5610277201ad2303838a87b71a8d76c411ef691e12c",
        )
        self.assertEqual(summary["accounting"]["skipped_fragment_count"], 3781335)
        for key, value in recomputed.items():
            self.assertAlmostEqual(summary["metrics"][key], value, places=15)

if __name__ == "__main__":
    unittest.main()
