"""Command-line smoke tests and direct-argument projection contracts."""

import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from bsreadsim import __version__
from bsreadsim.cli import (
    CommandLineError,
    build_methdb_document,
    build_parser,
    build_rrbs_document,
    build_run_document,
    build_variant_document,
)
from bsreadsim.run.config import normalize_run_config


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"


def run_module(*arguments: str) -> subprocess.CompletedProcess:
    environment = os.environ.copy()
    entries = [str(SOURCE_ROOT)]
    if environment.get("PYTHONPATH"):
        entries.append(environment["PYTHONPATH"])
    environment["PYTHONPATH"] = os.pathsep.join(entries)
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
    def test_version_and_top_level_help_are_lightweight(self) -> None:
        result = run_module_without_site_packages("--version")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual(result.stdout.strip(), "bsreadsim {}".format(__version__))

        help_result = run_module("--help")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("run", help_result.stdout)
        self.assertIn("build", help_result.stdout)
        self.assertIn("export", help_result.stdout)
        self.assertNotIn("get", help_result.stdout)
        self.assertNotIn("resources", help_result.stdout)
        self.assertNotIn("catalog", help_result.stdout)

    def test_public_help_uses_consistent_prepared_truth_terms(self) -> None:
        wgbs = run_module("run", "wgbs", "--help")
        variants = run_module("build", "variants", "--help")
        methdb = run_module("build", "methdb", "--help")
        for result in (wgbs, variants, methdb):
            self.assertEqual(result.returncode, 0, result.stderr)

        help_text = " ".join(
            "\n".join(result.stdout for result in (wgbs, variants, methdb)).split()
        ).lower()
        self.assertIn("prepared variant set", help_text)
        self.assertIn("prepared methylation profile", help_text)
        self.assertIn("simulation truth artifacts", help_text)
        self.assertNotIn("variant catalog", help_text)
        self.assertNotIn("methylation catalog", help_text)

    def test_export_test_fasta_copies_exact_resource_bytes(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "test.fa"
            copied = run_module_without_site_packages(
                "export", "test-fasta", "-o", str(output)
            )
            self.assertEqual(copied.returncode, 0, copied.stderr)
            self.assertTrue(output.read_bytes().startswith(b">chr10\n"))

    def test_export_requires_a_target_and_its_paths(self) -> None:
        for arguments in (
            ("export",),
            ("export", "test-fasta"),
            ("export", "methdb"),
            ("export", "methdb", "-i", "sample.methdb"),
            ("export", "methdb", "-o", "sample.bed.gz"),
        ):
            result = run_module_without_site_packages(*arguments)
            self.assertEqual(result.returncode, 2)

    def test_export_methdb_parses_input_and_compression_mode(self) -> None:
        compressed = build_parser().parse_args(
            [
                "export", "methdb", "-i", "sample.methdb",
                "-o", "sample.bed.gz",
            ]
        )
        self.assertEqual(compressed.input, Path("sample.methdb"))
        self.assertFalse(compressed.no_compression)

        plain = build_parser().parse_args(
            [
                "export", "methdb", "-i", "sample.methdb",
                "-o", "sample.bed", "--no-compression",
            ]
        )
        self.assertTrue(plain.no_compression)

    def test_build_rrbs_projects_a_fragment_domain(self) -> None:
        arguments = build_parser().parse_args(
            [
                "build", "rrbs", "-r", "reference.fa", "-o", "candidates.bed",
                "--cut-site", "C|CGG", "-l", "50", "--insert-min", "50",
                "--seed-mut", "7", "--seed-phase", "8",
            ]
        )
        normalized = normalize_run_config(
            build_rrbs_document(arguments, REPOSITORY_ROOT), REPOSITORY_ROOT
        ).normalized
        self.assertEqual(normalized["technology"], "RRBS")
        self.assertEqual(normalized["reads"]["count"], 2)
        self.assertEqual(normalized["mutation"]["rate"], 0)
        self.assertEqual(normalized["seeds"]["mutation"], "7")
        self.assertEqual(normalized["seeds"]["phasing"], "8")

    def test_run_stage_seeds_are_derived_but_build_defaults_remain_zero(self) -> None:
        run_arguments = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output",
                "-n", "2",
            ]
        )
        self.assertEqual(
            build_run_document(run_arguments, REPOSITORY_ROOT)["seeds"],
            {"mutation": None, "phasing": None, "methylation": None},
        )

        build_arguments = build_parser().parse_args(
            [
                "build", "methdb", "-r", "reference.fa", "-o", "truth.methdb",
            ]
        )
        self.assertEqual(
            build_methdb_document(build_arguments, REPOSITORY_ROOT)["seeds"],
            {"mutation": "0", "phasing": "0", "methylation": "0"},
        )

    def test_rrbs_cut_sites_have_a_default_and_accept_lists(self) -> None:
        default_run = build_parser().parse_args(
            [
                "run", "rrbs", "-r", "reference.fa", "-o", "output",
                "-n", "2",
            ]
        )
        self.assertEqual(
            build_run_document(default_run, REPOSITORY_ROOT)["rrbs"]["cut_sites"],
            ["C|CGG"],
        )

        default_build = build_parser().parse_args(
            [
                "build", "rrbs", "-r", "reference.fa", "-o", "candidates.bed",
            ]
        )
        self.assertEqual(
            build_rrbs_document(default_build, REPOSITORY_ROOT)["rrbs"]["cut_sites"],
            ["C|CGG"],
        )

        comma_separated = build_parser().parse_args(
            [
                "run", "rrbs", "-r", "reference.fa", "-o", "output",
                "-n", "2", "--cut-site", "C|CGG,G|ANTC",
            ]
        )
        self.assertEqual(
            build_run_document(comma_separated, REPOSITORY_ROOT)["rrbs"]["cut_sites"],
            ["C|CGG", "G|ANTC"],
        )

        mixed_case = build_parser().parse_args(
            [
                "run", "rrbs", "-r", "reference.fa", "-o", "output",
                "-n", "2", "--cut-site", "c|cgg,g|aNtc",
            ]
        )
        self.assertEqual(
            build_run_document(mixed_case, REPOSITORY_ROOT)["rrbs"]["cut_sites"],
            ["C|CGG", "G|ANTC"],
        )

        repeated = build_parser().parse_args(
            [
                "run", "rrbs", "-r", "reference.fa", "-o", "output",
                "-n", "2", "--cut-site", "C|CGG", "--cut-site", "G|ANTC",
            ]
        )
        with self.assertRaisesRegex(
            CommandLineError, "--cut-site may be supplied once"
        ):
            build_run_document(repeated, REPOSITORY_ROOT)

        explicit = build_parser().parse_args(
            [
                "run", "rrbs", "-r", "reference.fa", "-o", "output",
                "-n", "2", "--cut-site", "G|ANTC",
            ]
        )
        self.assertEqual(
            build_run_document(explicit, REPOSITORY_ROOT)["rrbs"]["cut_sites"],
            ["G|ANTC"],
        )

        duplicate = build_parser().parse_args(
            [
                "run", "rrbs", "-r", "reference.fa", "-o", "output",
                "-n", "2", "--cut-site", "C|CGG,c|cgg",
            ]
        )
        with self.assertRaisesRegex(
            CommandLineError, "duplicate --cut-site after case normalization"
        ):
            build_run_document(duplicate, REPOSITORY_ROOT)

        help_result = run_module("run", "rrbs", "--help")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("default: C|CGG", " ".join(help_result.stdout.split()))

    def test_build_variants_accepts_generation_or_vcf_normalization(self) -> None:
        generated = build_parser().parse_args(
            [
                "build", "variants", "-r", "reference.fa", "-o", "truth.vcf.gz",
                "--mutation-rate", "0.002", "--seed-mut", "19",
            ]
        )
        normalized = normalize_run_config(
            build_variant_document(generated, REPOSITORY_ROOT), REPOSITORY_ROOT
        ).normalized
        self.assertEqual(normalized["mutation"]["rate"], 0.002)
        self.assertEqual(normalized["seeds"]["mutation"], "19")
        self.assertEqual(normalized["reads"]["count"], 1)

        from_vcf = build_parser().parse_args(
            [
                "build", "variants", "-r", "reference.fa", "-o", "truth.vcf.gz",
                "--vcf", "input.vcf.gz", "--seed-phase", "23",
            ]
        )
        document = build_variant_document(from_vcf, REPOSITORY_ROOT)
        self.assertEqual(document["mutation"]["rate"], 0)
        self.assertEqual(document["seeds"]["phasing"], "23")
        self.assertEqual(document["inputs"]["vcf"], "input.vcf.gz")

    def test_vcf_and_mutation_rate_are_mutually_exclusive(self) -> None:
        cases = (
            (
                "run", "wgbs", "-r", "reference.fa", "-o", "output", "-n", "1",
                "--vcf", "input.vcf", "--mutation-rate", "0",
            ),
            (
                "build", "variants", "-r", "reference.fa", "-o", "truth.vcf.gz",
                "--vcf", "input.vcf", "--mutation-rate", "0",
            ),
            (
                "build", "methdb", "-r", "reference.fa", "-o", "truth.methdb",
                "--vcf", "input.vcf", "--mutation-rate", "0",
            ),
            (
                "build", "rrbs", "-r", "reference.fa", "-o", "candidates.bed",
                "--cut-site", "C|CGG", "--vcf", "input.vcf",
                "--mutation-rate", "0",
            ),
        )
        for arguments in cases:
            with self.subTest(command=arguments[:2]):
                result = run_module(*arguments)
                self.assertEqual(result.returncode, 2)
                self.assertIn("not allowed with argument", result.stderr)

    def test_build_methdb_projects_split_seeds_and_beta_parameters(self) -> None:
        arguments = build_parser().parse_args(
            [
                "build", "methdb", "-r", "reference.fa", "-o", "truth.methdb",
                "--mutation-rate", "0", "--seed-mut", "3", "--seed-phase", "4",
                "--seed-meth", "5", "--beta-cg", "2,8",
            ]
        )
        document = build_methdb_document(arguments, REPOSITORY_ROOT)
        self.assertEqual(
            document["seeds"],
            {"mutation": "3", "phasing": "4", "methylation": "5"},
        )
        self.assertEqual(document["methylation"]["beta"]["CG"], [2.0, 8.0])

    def test_build_methdb_accepts_all_poolable_text_profiles(self) -> None:
        for option, path, input_name in (
            ("--cgmap", "profile.cgmap", "cgmap"),
            ("--bedmethyl", "profile.bedmethyl", "bed_methyl"),
            ("--methbg", "profile.methbg", "methbg"),
            ("--methbed", "profile.methbed", "methbed"),
        ):
            with self.subTest(option=option):
                arguments = build_parser().parse_args(
                    [
                        "build", "methdb", "-r", "reference.fa",
                        "-o", "truth.methdb", option, path, "--pool-meth",
                    ]
                )
                document = build_methdb_document(arguments, REPOSITORY_ROOT)
                self.assertEqual(document["inputs"][input_name], path)
                self.assertTrue(document["methylation"]["cgmap_pool"])

    def test_build_methdb_accepts_asm_without_vcf(self) -> None:
        arguments = build_parser().parse_args(
            [
                "build", "methdb", "-r", "reference.fa", "-o", "truth.methdb",
                "--asm", "profile.asm", "--seed-phase", "9",
            ]
        )
        document = build_methdb_document(arguments, REPOSITORY_ROOT)
        self.assertEqual(document["inputs"], {"asm": "profile.asm"})
        self.assertEqual(document["mutation"]["rate"], 0)
        self.assertEqual(document["seeds"]["phasing"], "9")

    def test_beta_parameters_require_comma_delimited_pairs(self) -> None:
        result = run_module(
            "build", "methdb", "-r", "reference.fa", "-o", "truth.methdb",
            "--beta-cg", "2", "8",
        )
        self.assertEqual(result.returncode, 2)
        self.assertIn("must be ALPHA,BETA", result.stderr)

        help_result = run_module("build", "methdb", "--help")
        self.assertEqual(help_result.returncode, 0, help_result.stderr)
        self.assertIn("--beta-cg ALPHA,BETA", help_result.stdout)

    def test_run_subcommands_expose_only_relevant_sampling_options(self) -> None:
        wgbs = run_module("run", "wgbs", "--help")
        rrbs = run_module("run", "rrbs", "--help")
        tbs = run_module("run", "tbs", "--help")
        wgs = run_module("run", "wgs", "--help")
        wes = run_module("run", "wes", "--help")
        ts = run_module("run", "ts", "--help")
        for result in (wgbs, rrbs, tbs, wgs, wes, ts):
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("--reads N", result.stdout)
            self.assertNotIn("--fragments", result.stdout)
            self.assertIn("-f {fastq,fastq.gz,bam}", result.stdout)
            self.assertIn("--seed-mut", result.stdout)
            self.assertIn("default: 1000000", result.stdout)
            self.assertIn("omit to derive from --seed", result.stdout)
        self.assertIn("--gc-profile", wgbs.stdout)
        self.assertIn("--gc-profile", wgs.stdout)
        self.assertIn("--sampling {uniform,gc}", wgbs.stdout)
        self.assertIn("--sampling {uniform,gc}", wgs.stdout)
        self.assertIn("--sampling {uniform,score}", rrbs.stdout)
        self.assertIn("--sampling {uniform,score}", tbs.stdout)
        self.assertIn("--sampling {uniform,score}", wes.stdout)
        self.assertIn("--sampling {uniform,score}", ts.stdout)
        for result, input_option in (
            (wgbs, "--gc-profile"),
            (wgs, "--gc-profile"),
            (rrbs, "--cut-site"),
            (tbs, "--targets"),
            (wes, "--targets"),
            (ts, "--targets"),
        ):
            self.assertLess(
                result.stdout.index("--sampling"),
                result.stdout.index(input_option),
            )
        for result in (wgs, wes, ts):
            self.assertNotIn("--conversion-rate", result.stdout)
            self.assertNotIn("--seed-meth", result.stdout)
            self.assertNotIn("--save-methdb", result.stdout)

    def test_standard_run_documents_bypass_bisulfite_chemistry(self) -> None:
        cases = (
            ("wgs", "WGS", ()),
            (
                "wes",
                "WES",
                ("--targets", "exons.bed"),
            ),
            (
                "ts",
                "TS",
                ("--targets", "panel.bed"),
            ),
        )
        for command, technology, extra in cases:
            with self.subTest(command=command):
                arguments = build_parser().parse_args(
                    [
                        "run",
                        command,
                        "-r",
                        "reference.fa",
                        "-o",
                        "output",
                        "-n",
                        "10",
                        "--save-truth",
                        *extra,
                    ]
                )
                document = build_run_document(arguments, REPOSITORY_ROOT)
                normalized = normalize_run_config(
                    document, REPOSITORY_ROOT
                ).normalized

                self.assertEqual(document["technology"], technology)
                self.assertEqual(normalized["technology"], technology)
                self.assertEqual(document["sequencing"]["conversion_rate"], 0.0)
                self.assertFalse(document["output"]["save_methdb"])
                self.assertTrue(document["output"]["save_vcf"])
                self.assertEqual(document["inputs"], {})
                if technology in ("WES", "TS"):
                    self.assertEqual(document["fragments"]["insert_sd"], 25.0)
                    self.assertEqual(document["fragments"]["insert_min"], 100)
                    self.assertEqual(document["fragments"]["insert_max"], 1000)

    def test_methdb_owns_variants_and_defaults_mutation_to_zero(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output",
                "-n", "10", "--methdb", "profile.methdb",
            ]
        )
        document = build_run_document(arguments, REPOSITORY_ROOT)
        self.assertEqual(document["mutation"]["rate"], 0.0)
        self.assertEqual(document["inputs"]["methdb"], "profile.methdb")

        for conflicting in (
            ("--vcf", "variants.vcf"),
            ("--mutation-rate", "0.01"),
        ):
            parsed = build_parser().parse_args(
                [
                    "run", "wgbs", "-r", "reference.fa", "-o", "output",
                    "-n", "10", "--methdb", "profile.methdb", *conflicting,
                ]
            )
            with self.assertRaisesRegex(CommandLineError, "embeds variants"):
                build_run_document(parsed, REPOSITORY_ROOT)

    def test_text_methylation_profiles_allow_pooling(self) -> None:
        for option, path, input_name in (
            ("--cgmap", "profile.cgmap", "cgmap"),
            ("--bedmethyl", "profile.bedmethyl", "bed_methyl"),
            ("--methbg", "profile.methbg", "methbg"),
            ("--methbed", "profile.methbed", "methbed"),
        ):
            with self.subTest(option=option):
                arguments = build_parser().parse_args(
                    [
                        "run", "wgbs", "-r", "reference.fa", "-o", "output",
                        "-n", "10", option, path, "--pool-meth",
                    ]
                )
                document = build_run_document(arguments, REPOSITORY_ROOT)
                self.assertEqual(document["mutation"]["rate"], 0.001)
                self.assertEqual(document["inputs"][input_name], path)
                self.assertTrue(document["methylation"]["cgmap_pool"])

    def test_run_defaults_and_truth_flags_project_cleanly(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output", "-n", "10",
                "--save-truth",
            ]
        )
        document = build_run_document(arguments, REPOSITORY_ROOT)
        self.assertEqual(document["technology"], "WGBS")
        self.assertEqual(document["reads"]["count"], 10)
        self.assertNotIn("count", document["fragments"])
        self.assertEqual(document["output"]["format"], "fastq.gz")
        self.assertTrue(document["output"]["save_methdb"])
        self.assertTrue(document["output"]["save_vcf"])

        empty_vcf = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output", "-n", "10",
                "--mutation-rate", "0", "--save-vcf",
            ]
        )
        self.assertEqual(
            build_run_document(empty_vcf, REPOSITORY_ROOT)["mutation"]["rate"], 0
        )

    def test_bam_format_controls_annotation_output(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output", "-n", "10",
                "-f", "bam", "--fragment-realization",
            ]
        )
        document = build_run_document(arguments, REPOSITORY_ROOT)
        self.assertEqual(document["output"]["format"], "bam")
        self.assertTrue(document["output"]["fragment_summary"])
        self.assertTrue(document["output"]["fragment_realization"])

        invalid = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output", "-n", "10",
                "--fragment-summary",
            ]
        )
        with self.assertRaisesRegex(CommandLineError, "requires --format bam"):
            build_run_document(invalid, REPOSITORY_ROOT)

    def test_gc_profile_is_hashed_and_disables_default_mutations(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            profile = b"0.1\n0.7\n0.2\n"
            (directory / "gc.tsv").write_bytes(profile)
            arguments = build_parser().parse_args(
                [
                    "run", "wgbs", "-r", "reference.fa", "-o", "output", "-n", "100",
                    "--sampling", "gc", "--gc-profile", "gc.tsv", "--threads", "2",
                ]
            )
            normalized = normalize_run_config(
                build_run_document(arguments, directory), directory
            ).normalized
        artifact = normalized["coverage"]["artifact"]
        self.assertEqual(artifact["sha256"], hashlib.sha256(profile).hexdigest())
        self.assertEqual(normalized["fragments"]["insert_sd"], 25)
        self.assertEqual(normalized["mutation"]["rate"], 0)
        self.assertEqual(normalized["execution"]["threads"], 2)

    def test_gc_sampling_requires_an_explicit_profile_pair(self) -> None:
        common = [
            "run", "wgbs", "-r", "reference.fa", "-o", "output", "-n", "10",
        ]
        missing_profile = build_parser().parse_args(
            common + ["--sampling", "gc"]
        )
        with self.assertRaisesRegex(
            CommandLineError, "--sampling gc requires --gc-profile"
        ):
            build_run_document(missing_profile, REPOSITORY_ROOT)

        implicit_profile = build_parser().parse_args(
            common + ["--gc-profile", "gc.tsv"]
        )
        with self.assertRaisesRegex(
            CommandLineError, "--gc-profile requires --sampling gc"
        ):
            build_run_document(implicit_profile, REPOSITORY_ROOT)

    def test_rrbs_and_tbs_sampling_modes_are_unambiguous(self) -> None:
        rrbs_common = [
            "run", "rrbs", "-r", "reference.fa", "-o", "output", "-n", "24",
            "--cut-site", "C|CGG", "--rrbs-candidates", "scored.bed",
            "-l", "50", "--insert-min", "50",
        ]
        uniform = build_run_document(
            build_parser().parse_args(rrbs_common), REPOSITORY_ROOT
        )
        scored = build_run_document(
            build_parser().parse_args(rrbs_common + ["--sampling", "score"]),
            REPOSITORY_ROOT,
        )
        self.assertEqual(uniform["coverage"], {"kind": "uniform"})
        self.assertEqual(scored["coverage"], {"kind": "profile"})

        tbs = build_parser().parse_args(
            [
                "run", "tbs", "-r", "reference.fa", "-o", "output", "-n", "24",
                "--sampling", "score", "--targets", "targets.bed",
                "--insert-mean", "150", "--insert-sd", "0",
            ]
        )
        document = build_run_document(tbs, REPOSITORY_ROOT)
        self.assertEqual(document["coverage"], {"kind": "target-score"})
        self.assertEqual(document["fragments"]["insert_min"], 150)
        self.assertEqual(document["fragments"]["insert_mean"], 150)
        self.assertEqual(document["fragments"]["insert_max"], 150)
        self.assertEqual(document["fragments"]["insert_sd"], 0)

    def test_input_vcf_and_methylation_inputs_project_explicitly(self) -> None:
        arguments = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output", "-n", "10",
                "--vcf", "variants.vcf", "--bedmethyl", "levels.bed.gz",
                "--asm-bed", "alleles.bed.gz", "--seed-phase", "9",
            ]
        )
        document = build_run_document(arguments, REPOSITORY_ROOT)
        self.assertEqual(document["mutation"]["rate"], 0)
        self.assertEqual(set(document["inputs"]), {"vcf", "bed_methyl", "asm_bed"})
        self.assertEqual(document["seeds"]["phasing"], "9")

        asm_only_arguments = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output",
                "-n", "10", "--asm", "alleles.asm", "--seed-phase", "11",
            ]
        )
        asm_only = build_run_document(
            asm_only_arguments, REPOSITORY_ROOT
        )
        self.assertEqual(set(asm_only["inputs"]), {"asm"})
        self.assertEqual(asm_only["mutation"]["rate"], 0)
        self.assertEqual(asm_only["seeds"]["phasing"], "11")

    def test_read_count_defaults_and_requires_complete_pairs(self) -> None:
        default = build_parser().parse_args(
            ["run", "wgbs", "-r", "reference.fa", "-o", "output"]
        )
        self.assertEqual(
            build_run_document(default, REPOSITORY_ROOT)["reads"],
            {"count": 1_000_000},
        )

        paired = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output",
                "--reads", "10",
            ]
        )
        self.assertEqual(
            build_run_document(paired, REPOSITORY_ROOT)["reads"],
            {"count": 10},
        )

        odd_paired = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output",
                "--reads", "9",
            ]
        )
        with self.assertRaisesRegex(CommandLineError, "must be even"):
            build_run_document(odd_paired, REPOSITORY_ROOT)

        single = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output",
                "--reads", "9", "--single-end",
            ]
        )
        self.assertEqual(
            build_run_document(single, REPOSITORY_ROOT)["reads"],
            {"count": 9},
        )

        depth = build_parser().parse_args(
            [
                "run", "wgbs", "-r", "reference.fa", "-o", "output",
                "--depth", "12.5",
            ]
        )
        self.assertEqual(
            build_run_document(depth, REPOSITORY_ROOT)["reads"],
            {"depth": 12.5},
        )

        conflict = run_module(
            "run", "wgbs", "-r", "reference.fa", "-o", "output",
            "--reads", "10", "--depth", "12.5",
        )
        self.assertEqual(conflict.returncode, 2)
        self.assertIn("not allowed with argument", conflict.stderr)

    def test_old_command_and_option_names_are_not_accepted(self) -> None:
        for arguments in (
            ("catalog", "variants", "--help"),
            ("resources", "list"),
            ("run", "wgbs", "-r", "ref.fa", "-o", "out", "-n", "1", "--bam"),
            (
                "run", "wgbs", "-r", "ref.fa", "-o", "out", "-n", "1",
                "--read-pairs", "1",
            ),
            (
                "run", "wgbs", "-r", "ref.fa", "-o", "out", "-n", "1",
                "--insert-size", "150",
            ),
            (
                "run", "wgbs", "-r", "ref.fa", "-o", "out",
                "--fragments", "2",
            ),
            (
                "run", "wgbs", "-r", "ref.fa", "-o", "out", "-n", "2",
                "--bed-methyl", "profile.bedmethyl",
            ),
            (
                "run", "wgbs", "-r", "ref.fa", "-o", "out", "-n", "2",
                "--cgmap", "profile.cgmap", "--pool-methylation-values",
            ),
            (
                "run", "wgbs", "-r", "ref.fa", "-o", "out", "-n", "2",
                "--cgmap", "profile.cgmap", "--cgmap-pool",
            ),
        ):
            result = run_module(*arguments)
            self.assertEqual(result.returncode, 2)


if __name__ == "__main__":
    unittest.main()
