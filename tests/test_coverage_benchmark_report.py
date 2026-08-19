"""Integrity checks for the retained one-CPU coverage benchmark."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_ROOT = REPOSITORY_ROOT / "data" / "benchmarks"
REPORT = BENCHMARK_ROOT / "chr21-wgbs-uniform-vs-target-1cpu.json"


def _load_harness():
    path = BENCHMARK_ROOT / "compare_wgbs_coverage.py"
    specification = importlib.util.spec_from_file_location(
        "compare_wgbs_coverage", path
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load coverage benchmark harness")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class CoverageBenchmarkReportTests(unittest.TestCase):
    def test_retained_report_recomputes_and_preserves_output_gates(self) -> None:
        harness = _load_harness()
        report = json.loads(REPORT.read_text(encoding="utf-8"))

        self.assertEqual(
            report["schema"],
            "bsreadsim-wgbs-coverage-paired-benchmark-1",
        )
        self.assertEqual(report["environment"]["git_status"], "")
        self.assertEqual(
            report["environment"]["git_commit"],
            "4b215777d276869efcc03cd09baf67db1929f7fd",
        )
        self.assertEqual(report["execution_lane"]["affinity"], [0])
        self.assertEqual(
            report["execution_lane"]["paired_order"],
            ["uniform", "profile", "profile", "uniform", "uniform", "profile"],
        )
        self.assertEqual(report["workload"]["fragments"], 500000)
        self.assertEqual(report["workload"]["reads"], 1000000)
        self.assertEqual(report["workload"]["workers"], 1)
        self.assertEqual(report["workload"]["core_workers"], 1)

        for lane_name in ("core-producer", "FASTQ-E2E"):
            lane = report["lanes"][lane_name]
            measurements = lane["measurements"]
            self.assertEqual(len(measurements), 6)
            recomputed = harness._summarize(measurements)
            for field in (
                "profile_to_uniform_throughput_ratio",
                "profile_throughput_change_percent",
                "profile_wall_overhead_percent",
            ):
                self.assertAlmostEqual(
                    lane["summary"][field], recomputed[field], places=15
                )
            for coverage in ("uniform", "profile"):
                for field in (
                    "median_wall_seconds",
                    "median_fragments_per_second",
                    "minimum_fragments_per_second",
                    "maximum_fragments_per_second",
                ):
                    self.assertAlmostEqual(
                        lane["summary"][coverage][field],
                        recomputed[coverage][field],
                        places=15,
                    )

        e2e = report["lanes"]["FASTQ-E2E"]
        for coverage in ("uniform", "profile"):
            measurements = [
                item
                for item in e2e["measurements"]
                if item["coverage"] == coverage
            ]
            signatures = {
                json.dumps(
                    {
                        "core_counts": item["core_counts"],
                        "python_counts": item["python_counts"],
                        "outputs": item["outputs"],
                    },
                    sort_keys=True,
                )
                for item in measurements
            }
            self.assertEqual(len(signatures), 1)
            self.assertEqual(
                {item["role"] for item in measurements[0]["outputs"]},
                {"read1", "read2"},
            )

        uniform = next(
            item for item in e2e["measurements"] if item["coverage"] == "uniform"
        )
        profile = next(
            item for item in e2e["measurements"] if item["coverage"] == "profile"
        )
        self.assertEqual(uniform["core_counts"]["skipped_fragment_count"], 0)
        self.assertEqual(profile["core_counts"]["skipped_fragment_count"], 405921)
        acceptance = 500000 / (500000 + 405921)
        self.assertAlmostEqual(
            e2e["summary"]["profile_observed_proposal_acceptance"],
            acceptance,
            places=15,
        )


if __name__ == "__main__":
    unittest.main()
