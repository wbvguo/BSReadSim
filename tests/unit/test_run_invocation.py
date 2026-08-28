"""Tests for normalized-config to full run-command reconstruction."""

from pathlib import Path
import tempfile
import unittest

from bsreadsim.cli import build_parser, build_run_document
from bsreadsim.run.config import LoadedRunConfig, normalize_run_config
from bsreadsim.run.invocation import FullCommandError, build_full_run_argv
from bsreadsim.run.prepare import materialize_master_seed


class FullRunCommandTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.directory = Path(self.temporary_directory.name).resolve()
        for name, content in (
            ("quality.json", b'{"quality":"model"}\n'),
            ("error.json", b'{"error":"model"}\n'),
            ("coverage.tsv", b"0.5\n"),
        ):
            (self.directory / name).write_bytes(content)

    def effective(self, argv: list[str]) -> LoadedRunConfig:
        arguments = build_parser().parse_args(argv)
        loaded = normalize_run_config(
            build_run_document(arguments, self.directory), self.directory
        )
        return materialize_master_seed(loaded, entropy=lambda bits: 123456789)

    def assert_full_round_trip(self, argv: list[str]) -> list[str]:
        effective = self.effective(argv)
        received = ["bsreadsim", *argv]
        full_argv = build_full_run_argv(effective.normalized, received)
        reparsed = build_parser().parse_args(full_argv[1:])
        rebuilt = normalize_run_config(
            build_run_document(reparsed, self.directory), self.directory
        )
        self.assertEqual(rebuilt.canonical_json, effective.canonical_json)
        self.assertEqual(full_argv.count("--seed"), 1)
        self.assertEqual(
            full_argv[full_argv.index("--seed") + 1], "123456789"
        )
        return full_argv

    def test_every_released_technology_round_trips_through_its_full_command(
        self,
    ) -> None:
        cases = {
            "wgbs": [],
            "wgs": [],
            "rrbs": ["--cut-site", "C|CGG"],
            "tbs": ["--targets", "targets.bed"],
            "wes": ["--targets", "targets.bed"],
            "ts": ["--targets", "targets.bed"],
        }
        for technology, specific in cases.items():
            with self.subTest(technology=technology):
                full_argv = self.assert_full_round_trip(
                    [
                        "run",
                        technology,
                        "--reference",
                        "reference.fa",
                        "--output",
                        "output-{}".format(technology),
                        "--reads",
                        "6",
                        *specific,
                    ]
                )
                self.assertEqual(full_argv[:3], ["bsreadsim", "run", technology])

    def test_full_command_expands_advanced_bisulfite_settings(self) -> None:
        argv = [
            "run",
            "wgbs",
            "--reference",
            "reference.fa",
            "--output",
            "advanced-output",
            "--depth",
            "2.5",
            "--vcf",
            "sample.vcf",
            "--cgmap",
            "sample.cgmap",
            "--asm",
            "sample.asm",
            "--seed-mut",
            "4",
            "--seed-phase",
            "5",
            "--seed-meth",
            "6",
            "--single-end",
            "--read-length",
            "75",
            "--insert-min",
            "80",
            "--insert-mean",
            "100",
            "--insert-max",
            "120",
            "--insert-sd",
            "5.5",
            "--max-ambiguous-fraction",
            "0.02",
            "--indel-fraction",
            "0.2",
            "--indel-extension-probability",
            "0.3",
            "--homozygous-only",
            "--beta-cg",
            "2,3",
            "--beta-chg",
            "4,5",
            "--beta-chh",
            "6,7",
            "--cpg-only",
            "--cgmap-pool",
            "--methylation-model",
            "bilstm",
            "--no-update-variant-boundaries",
            "--conversion-rate",
            "0.97",
            "--undirectional",
            "--quality-model",
            str(self.directory / "quality.json"),
            "--error-model",
            str(self.directory / "error.json"),
            "--threads",
            "12",
            "--prefix",
            "advanced",
            "--format",
            "bam",
            "--gzip-level",
            "4",
            "--fragment-summary",
            "--fragment-realization",
            "--save-truth",
            "--core",
            "custom-core",
        ]

        full_argv = self.assert_full_round_trip(argv)

        for option in (
            "--vcf",
            "--cgmap",
            "--asm",
            "--single-end",
            "--homozygous-only",
            "--cpg-only",
            "--cgmap-pool",
            "--no-update-variant-boundaries",
            "--undirectional",
            "--quality-model",
            "--error-model",
            "--fragment-summary",
            "--fragment-realization",
            "--save-methdb",
            "--save-vcf",
            "--core",
        ):
            self.assertIn(option, full_argv)
        self.assertNotIn("--save-truth", full_argv)
        self.assertEqual(full_argv[full_argv.index("--core") + 1], "custom-core")

    def test_full_command_covers_nonuniform_sampling_modes(self) -> None:
        cases = (
            [
                "run",
                "wgs",
                "-r",
                "reference.fa",
                "-o",
                "wgs-output",
                "-n",
                "2",
                "--gc-profile",
                str(self.directory / "coverage.tsv"),
            ],
            [
                "run",
                "rrbs",
                "-r",
                "reference.fa",
                "-o",
                "rrbs-output",
                "-n",
                "2",
                "--cut-site",
                "C|CGG",
                "--rrbs-candidates",
                "candidates.bed",
                "--sampling",
                "score",
            ],
            [
                "run",
                "tbs",
                "-r",
                "reference.fa",
                "-o",
                "tbs-output",
                "-n",
                "2",
                "--targets",
                "targets.bed",
                "--sampling",
                "score",
            ],
        )
        for argv in cases:
            with self.subTest(technology=argv[1]):
                self.assert_full_round_trip(argv)

    def test_unrepresentable_cli_state_is_rejected(self) -> None:
        effective = self.effective(
            [
                "run",
                "wgbs",
                "-r",
                "reference.fa",
                "-o",
                "output",
                "-n",
                "2",
            ]
        )
        normalized = effective.normalized
        normalized["fragments"]["read_length_2"] = 101

        with self.assertRaisesRegex(FullCommandError, "different R1 and R2"):
            build_full_run_argv(normalized, ["bsreadsim", "run", "wgbs"])


if __name__ == "__main__":
    unittest.main()
