"""Focused contracts for the generic end-to-end benchmark driver."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DRIVER = REPOSITORY_ROOT / "data" / "benchmarks" / "benchmark.py"


def _load_driver():
    specification = importlib.util.spec_from_file_location(
        "current_benchmark_driver", DRIVER
    )
    if specification is None or specification.loader is None:
        raise RuntimeError("cannot load benchmark driver")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class BenchmarkDriverTests(unittest.TestCase):
    def test_expected_roles_cover_modes_truth_bam_and_single_end(self) -> None:
        driver = _load_driver()
        paired = {
            "fragments": {"paired_end": True},
            "output": {"truth_bam": False},
        }
        self.assertEqual(driver._expected_output_roles(paired, "production"), {"read1", "read2"})
        self.assertEqual(driver._expected_output_roles(paired, "debug"), {"read1", "read2", "truth"})
        paired["output"]["truth_bam"] = True
        self.assertEqual(driver._expected_output_roles(paired, "production"), {"read1", "read2", "truth_bam"})
        self.assertEqual(driver._expected_output_roles(paired, "debug"), {"read1", "read2", "truth", "truth_bam"})
        single = {"fragments": {"paired_end": False}, "output": {"truth_bam": False}}
        self.assertEqual(driver._expected_output_roles(single, "production"), {"read1"})

    def test_reduced_warmup_changes_only_the_copied_read_pair_count(self) -> None:
        driver = _load_driver()
        template = {
            "fragments": {"paired_end": True, "read_pairs": 500000},
            "output": {"directory": "original", "truth_bam": False},
        }
        observed = driver._document_for_warmup(template, Path("/tmp/warmup"), 50000)
        self.assertEqual(observed["fragments"]["read_pairs"], 50000)
        self.assertEqual(observed["output"]["directory"], "/tmp/warmup")
        self.assertEqual(template["fragments"]["read_pairs"], 500000)
        self.assertEqual(template["output"]["directory"], "original")

        depth_template = {
            "fragments": {"paired_end": True, "depth": 30},
            "output": {"directory": "original", "truth_bam": False},
        }
        with self.assertRaisesRegex(ValueError, "--read-pairs"):
            driver._document_for_warmup(depth_template, Path("/tmp/warmup"), 50000)


if __name__ == "__main__":
    unittest.main()
