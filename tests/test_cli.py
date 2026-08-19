"""Command-line smoke tests and direct-argument projection contracts."""

import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
sys.path.insert(0, str(SOURCE_ROOT))

from bsreadsim.cli import (  # noqa: E402
    CommandLineError,
    build_parser,
    build_rrbs_catalog_document,
    build_run_document,
)
from bsreadsim import __version__  # noqa: E402
from bsreadsim.config import normalize_run_config  # noqa: E402


def run_module(*arguments: str) -> subprocess.CompletedProcess:
    environment = os.environ.copy()
    existing_pythonpath = environment.get("PYTHONPATH")
    pythonpath_entries = [str(SOURCE_ROOT)]
    if existing_pythonpath:
        pythonpath_entries.append(existing_pythonpath)
    environment["PYTHONPATH"] = os.pathsep.join(pythonpath_entries)

    return subprocess.run(
        [sys.executable, "-m", "bsreadsim", *arguments],
        cwd=REPOSITORY_ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )


def run_module_without_site_packages(*arguments: str) -> subprocess.CompletedProcess:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(SOURCE_ROOT)
    return subprocess.run(
        [sys.executable, "-S", "-m", "bsreadsim", *arguments],
        cwd=REPOSITORY_ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )


class CommandLineTests(unittest.TestCase):
    def test_version(self) -> None:
        result = run_module("--version")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "bsreadsim {}".format(__version__))

    def test_version_does_not_import_runtime_dependencies(self) -> None:
        result = run_module_without_site_packages("--version")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "bsreadsim {}".format(__version__))

    def test_help(self) -> None:
        result = run_module("--help")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("usage: bsreadsim", result.stdout)
        self.assertIn("--version", result.stdout)
        self.assertIn("run", result.stdout)
        self.assertIn("catalog", result.stdout)
        self.assertNotIn("run-config", result.stdout)

    def test_catalog_rrbs_help_exposes_direct_cli_domain(self) -> None:
        result = run_module("catalog", "rrbs", "--help")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--reference", result.stdout)
        self.assertIn("--cut-site", result.stdout)
        self.assertIn("no JSON configuration file", result.stdout)
        self.assertIn("--output", result.stdout)
        self.assertIn("--core", result.stdout)

    def test_rrbs_catalog_is_projected_without_a_json_config(self) -> None:
        arguments = build_parser().parse_args(
            [
                "catalog",
                "rrbs",
                "-r",
                "reference.fa",
                "-o",
                "candidates.bed",
                "--cut-site",
                "C|CGG",
                "--read-length",
                "50",
                "--insert-min",
                "50",
            ]
        )
        self.assertFalse(hasattr(arguments, "config"))

        normalized = normalize_run_config(
            build_rrbs_catalog_document(arguments, REPOSITORY_ROOT),
            REPOSITORY_ROOT,
        ).normalized

        self.assertEqual(normalized["technology"], "RRBS")
        self.assertEqual(normalized["rrbs"], {"cut_sites": ["C|CGG"]})
        self.assertEqual(normalized["coverage"], {"kind": "uniform"})
        self.assertEqual(normalized["mutation"]["rate"], 0)
        self.assertEqual(normalized["fragments"]["read_pairs"], 1)

        variant_arguments = build_parser().parse_args(
            [
                "catalog", "rrbs", "-r", "reference.fa", "-o",
                "candidates.bed", "--cut-site", "C|CGG", "--vcf",
                "variants.vcf",
            ]
        )
        with self.assertRaisesRegex(CommandLineError, "explicit --seed"):
            build_rrbs_catalog_document(variant_arguments, REPOSITORY_ROOT)

    def test_run_help_exposes_the_component_boundary(self) -> None:
        result = run_module("run", "--help")

        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--reference", result.stdout)
        self.assertIn("--read-pairs", result.stdout)
        self.assertIn("--coverage-profile", result.stdout)
        self.assertIn("--bed-methyl", result.stdout)
        self.assertIn("--asm-bed", result.stdout)
        self.assertIn("--core", result.stdout)
        self.assertIn("--mode", result.stdout)

    def test_no_command_is_not_reported_as_a_successful_run(self) -> None:
        result = run_module()

        self.assertEqual(result.returncode, 2)
        self.assertIn("usage: bsreadsim", result.stdout)

    def test_run_defaults_to_production(self) -> None:
        arguments = build_parser().parse_args(
            ["run", "-r", "reference.fa", "-o", "output", "-n", "10"]
        )

        self.assertEqual(arguments.mode, "production")
        self.assertEqual(arguments.technology, "WGBS")
        self.assertFalse(arguments.truth_bam)
        self.assertFalse(hasattr(arguments, "protocol_major"))

    def test_debug_mode_is_explicit(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run",
                "-r",
                "reference.fa",
                "-o",
                "output",
                "-n",
                "10",
                "--mode",
                "debug",
            ]
        )

        self.assertEqual(arguments.mode, "debug")

    def test_truth_bam_flag_projects_into_the_output_contract(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run",
                "-r",
                "reference.fa",
                "-o",
                "output",
                "-n",
                "10",
                "--truth-bam",
            ]
        )
        document = build_run_document(arguments, REPOSITORY_ROOT)
        self.assertTrue(document["output"]["truth_bam"])

    def test_json_run_config_command_is_not_available(self) -> None:
        result = run_module("run-config", "run.json", "--mode", "debug")

        self.assertEqual(result.returncode, 2)
        self.assertIn("invalid choice", result.stderr)

    def test_profile_path_is_hashed_and_projected_without_json(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            profile_bytes = b"0.1\n0.7\n0.2\n"
            (directory / "coverage.tsv").write_bytes(profile_bytes)
            arguments = build_parser().parse_args(
                [
                    "run",
                    "-r",
                    "reference.fa",
                    "-o",
                    "output",
                    "-n",
                    "100",
                    "--coverage-profile",
                    "coverage.tsv",
                    "--insert-size",
                    "350",
                    "--workers",
                    "2",
                ]
            )

            document = build_run_document(arguments, directory)
            normalized = normalize_run_config(document, directory).normalized

        artifact = normalized["coverage"]["artifact"]
        self.assertEqual(normalized["coverage"]["kind"], "profile")
        self.assertEqual(artifact["format"], "tsv")
        self.assertEqual(artifact["version"], "wgbs-gc-target-v1")
        self.assertEqual(
            artifact["sha256"], hashlib.sha256(profile_bytes).hexdigest()
        )
        self.assertEqual(normalized["fragments"]["insert_min"], 350)
        self.assertEqual(normalized["fragments"]["insert_max"], 350)
        self.assertEqual(normalized["fragments"]["insert_stddev"], 0)
        self.assertEqual(normalized["mutation"]["rate"], 0)
        self.assertEqual(normalized["execution"]["workers"], 2)
        self.assertEqual(normalized["sequencing"]["quality"]["phred"], 40)
        self.assertEqual(normalized["sequencing"]["error"]["rate"], 0.005)

    def test_direct_projection_rejects_conflicting_shortcuts(self) -> None:
        common = ["run", "-r", "ref.fa", "-o", "out", "-n", "10"]
        insert_conflict = build_parser().parse_args(
            common + ["--insert-size", "350", "--insert-min", "100"]
        )
        with self.assertRaisesRegex(CommandLineError, "insert range"):
            build_run_document(insert_conflict, REPOSITORY_ROOT)

        profile_on_rrbs = build_parser().parse_args(
            common
            + [
                "--technology",
                "RRBS",
                "--cut-site",
                "C|CGG",
                "--coverage-profile",
                "missing.tsv",
            ]
        )
        with self.assertRaisesRegex(CommandLineError, "WGBS only"):
            build_run_document(profile_on_rrbs, REPOSITORY_ROOT)

    def test_rrbs_arguments_project_the_enzyme_catalog(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run",
                "-r",
                "reference.fa",
                "-o",
                "output",
                "-n",
                "24",
                "--technology",
                "RRBS",
                "--cut-site",
                "C|CGG",
                "--cut-site",
                "A|ATT",
                "--read-length",
                "50",
                "--insert-min",
                "50",
            ]
        )

        normalized = normalize_run_config(
            build_run_document(arguments, REPOSITORY_ROOT),
            REPOSITORY_ROOT,
        ).normalized

        self.assertEqual(normalized["technology"], "RRBS")
        self.assertEqual(normalized["rrbs"]["cut_sites"], ["C|CGG", "A|ATT"])
        self.assertEqual(normalized["coverage"], {"kind": "uniform"})

    def test_rrbs_candidate_scores_are_opt_in(self) -> None:
        common = [
            "run",
            "-r",
            "reference.fa",
            "-o",
            "output",
            "-n",
            "24",
            "--technology",
            "RRBS",
            "--cut-site",
            "C|CGG",
            "--rrbs-candidates",
            "scored.bed",
            "--read-length",
            "50",
            "--insert-min",
            "50",
        ]
        uniform = normalize_run_config(
            build_run_document(build_parser().parse_args(common), REPOSITORY_ROOT),
            REPOSITORY_ROOT,
        ).normalized
        profiled = normalize_run_config(
            build_run_document(
                build_parser().parse_args(common + ["--rrbs-score"]),
                REPOSITORY_ROOT,
            ),
            REPOSITORY_ROOT,
        ).normalized

        self.assertEqual(uniform["coverage"], {"kind": "uniform"})
        self.assertEqual(profiled["coverage"], {"kind": "profile"})
        self.assertEqual(uniform["mutation"]["rate"], 0)
        self.assertTrue(
            uniform["rrbs"]["candidate_bed"].endswith("/scored.bed")
        )

        missing = build_parser().parse_args(
            [
                "run", "-r", "reference.fa", "-o", "output", "-n", "24",
                "--technology", "RRBS", "--cut-site", "C|CGG",
                "--rrbs-score",
            ]
        )
        with self.assertRaisesRegex(CommandLineError, "requires --rrbs-candidates"):
            build_run_document(missing, REPOSITORY_ROOT)

    def test_tbs_arguments_project_targets_and_score_weights(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run",
                "-r",
                "reference.fa",
                "-o",
                "output",
                "-n",
                "24",
                "--technology",
                "TBS",
                "--targets",
                "targets.bed",
                "--target-score",
                "--insert-size",
                "150",
                "--fragment-center-stddev",
                "0",
            ]
        )

        normalized = normalize_run_config(
            build_run_document(arguments, REPOSITORY_ROOT),
            REPOSITORY_ROOT,
        ).normalized

        self.assertEqual(normalized["technology"], "TBS")
        self.assertEqual(normalized["coverage"], {"kind": "target-score"})
        self.assertEqual(normalized["tbs"]["fragment_center_stddev"], 0)
        self.assertEqual(normalized["fragments"]["insert_mean"], 150)

    def test_vcf_automatically_disables_de_novo_mutation(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run",
                "-r",
                "reference.fa",
                "-o",
                "output",
                "-n",
                "10",
                "--vcf",
                "variants.vcf",
            ]
        )

        document = build_run_document(arguments, REPOSITORY_ROOT)

        self.assertEqual(document["mutation"]["rate"], 0)

    def test_bed_methyl_and_asm_bed_project_as_explicit_inputs(self) -> None:
        common = [
            "run", "-r", "reference.fa", "-o", "output", "-n", "10"
        ]
        arguments = build_parser().parse_args(
            common
            + [
                "--vcf", "variants.vcf",
                "--bed-methyl", "levels.bedmethyl.gz",
                "--asm-bed", "alleles.bed.gz",
            ]
        )

        document = build_run_document(arguments, REPOSITORY_ROOT)

        self.assertEqual(
            set(document["inputs"]),
            {"vcf", "bed_methyl", "asm_bed"},
        )
        self.assertEqual(document["mutation"]["rate"], 0)

        without_vcf = build_parser().parse_args(
            common + ["--asm-bed", "alleles.bed"]
        )
        with self.assertRaisesRegex(CommandLineError, "requires --vcf"):
            build_run_document(without_vcf, REPOSITORY_ROOT)

        with self.assertRaises(SystemExit):
            build_parser().parse_args(
                common
                + [
                    "--cgmap", "levels.cgmap",
                    "--bed-methyl", "levels.bedmethyl",
                ]
            )

    def test_retired_protocol_selector_is_not_accepted(self) -> None:
        result = run_module(
            "run",
            "-r",
            "reference.fa",
            "-o",
            "output",
            "-n",
            "10",
            "--protocol-major",
            "1",
        )

        self.assertEqual(result.returncode, 2)
        self.assertIn("unrecognized arguments", result.stderr)


if __name__ == "__main__":
    unittest.main()
