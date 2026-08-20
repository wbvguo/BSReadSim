"""Integrity gates for the current one-million-read evidence."""

from __future__ import annotations

import json
from pathlib import Path
import statistics
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BENCHMARK_ROOT = REPOSITORY_ROOT / "data" / "benchmarks"
SOURCE_COMMIT = "7fd2301d4803fd470607a148402853443c3fd4bd"


def _load(name: str):
    return json.loads((BENCHMARK_ROOT / name).read_text(encoding="utf-8"))


class CurrentBenchmarkReportTests(unittest.TestCase):
    def test_wgbs_component_report_preserves_interleaved_gates(self) -> None:
        report = _load("bsreadsim-0.3.0-wgbs-components-1m-reads-2026-08-19.json")
        self.assertEqual(report["workload"]["reads"], 1000000)
        self.assertEqual(report["environment"]["git_commit"], SOURCE_COMMIT)
        self.assertEqual(report["environment"]["git_status"], "")
        for lane in report["lanes"].values():
            self.assertEqual(len(lane["measurements"]), 6)
            for coverage in ("uniform", "profile"):
                values = [item["wall_seconds"] for item in lane["measurements"] if item["coverage"] == coverage]
                self.assertEqual(len(values), 3)
                self.assertAlmostEqual(lane["summary"][coverage]["median_wall_seconds"], statistics.median(values), places=15)

    def test_output_mode_report_keeps_fastq_identity_and_roles(self) -> None:
        report = _load("bsreadsim-0.3.0-output-modes-1m-reads-2026-08-19.json")
        self.assertTrue(report["gates"]["fastq_hashes_match_across_policies"])
        expected = {
            "production": {"read1", "read2"},
            "debug": {"read1", "read2", "truth"},
            "truth_bam": {"read1", "read2", "truth_bam"},
        }
        fastq_signatures = []
        for policy, roles in expected.items():
            lane = report["policies"][policy]
            self.assertEqual(lane["environment"]["git_commit"], SOURCE_COMMIT)
            outputs = lane["measurements"][0]["outputs"]
            self.assertEqual({item["role"] for item in outputs}, roles)
            fastq_signatures.append([(item["role"], item["record_count"], item["sha256"]) for item in outputs if item["role"] in {"read1", "read2"}])
        self.assertEqual(len({json.dumps(item) for item in fastq_signatures}), 1)

    def test_technology_report_is_balanced_and_reproducible(self) -> None:
        report = _load("bsreadsim-0.3.0-technologies-1m-reads-2026-08-19.json")
        self.assertEqual(report["environment"]["git_commit"], SOURCE_COMMIT)
        self.assertEqual(report["environment"]["git_status"], "")
        self.assertEqual(report["workload"]["reads"], 1000000)
        self.assertEqual(report["workload"]["reference_sha256"], "5d436a0de36d479aa65d751f3fc56435da77a0aec9832116915e3dd87fd79235")
        self.assertEqual(report["workload"]["targets_sha256"], "826f87e42598a6428b52cc1490a8c8902d196c43ac9dfafd5cf68eebf66d8e8d")
        measurements = report["measurements"]
        self.assertEqual(len(measurements), 9)
        self.assertEqual(len({json.dumps(item["core_counts"], sort_keys=True) for item in measurements}), 1)
        for technology in ("WGBS", "RRBS", "TBS"):
            values = [item for item in measurements if item["technology"] == technology]
            self.assertEqual(len(values), 3)
            signatures = {json.dumps({"core": item["core_counts"], "python": item["python_counts"], "outputs": item["outputs"]}, sort_keys=True) for item in values}
            self.assertEqual(len(signatures), 1)
            self.assertAlmostEqual(report["summary"][technology]["median_wall_seconds"], statistics.median(item["wall_seconds"] for item in values), places=15)


if __name__ == "__main__":
    unittest.main()
