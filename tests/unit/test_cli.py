"""Command-line parsing and the settings it resolves to."""

from contextlib import redirect_stderr
import io
import os
from pathlib import Path
import subprocess
import sys
import unittest

from bsreadsim import __version__
from bsreadsim.cli import build_parser, from_arguments
from bsreadsim.errors import BSReadSimError
from bsreadsim.htsim import HtsimCore
from bsreadsim.output import full_command
from bsreadsim.settings import Settings, derive_stage_seed

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]


def run_module(*arguments: str, site: bool = True) -> subprocess.CompletedProcess:
    environment = dict(os.environ, PYTHONPATH=str(REPOSITORY_ROOT / "src"))
    flags = [] if site else ["-S"]
    return subprocess.run([sys.executable, *flags, "-m", "bsreadsim", *arguments],
                          cwd=REPOSITORY_ROOT, env=environment, check=False,
                          capture_output=True, text=True)


def settings(*arguments: str) -> Settings:
    return from_arguments(build_parser().parse_args(list(arguments)), entropy=lambda bits: 7)


def run_settings(assay: str, *arguments: str) -> Settings:
    return settings("run", assay, "-r", "ref.fa", "-o", "out", *arguments)


class CommandLineTests(unittest.TestCase):
    def test_version_and_help_do_not_need_site_packages(self):
        for option in ("-v", "--version"):
            result = run_module(option, site=False)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertEqual(result.stdout.strip(), "bsreadsim {}".format(__version__))
        result = run_module("--help", site=False)
        self.assertEqual(result.returncode, 0, result.stderr)
        for command in ("run", "validate", "build", "export"):
            self.assertIn(command, result.stdout)

    def test_help_uses_prepared_truth_terms(self):
        text = " ".join(" ".join(
            run_module(*command, "--help").stdout
            for command in (("run", "wgbs"), ("build", "variants"), ("build", "methdb"))
        ).split()).lower()
        for phrase in ("prepared variant set", "prepared methylation profile",
                       "simulation truth artifacts"):
            self.assertIn(phrase, text)

    def test_options_are_neither_abbreviated_nor_shortened(self):
        parser = build_parser()
        for arguments in (["-s", "1"], ["--see", "1"], ["--no-update-variant-boundaries"]):
            with self.assertRaises(SystemExit), redirect_stderr(io.StringIO()):
                parser.parse_args(["run", "wgbs", "-r", "r", "-o", "o", *arguments])

    def test_export_requires_its_paths(self):
        for arguments in (("export",), ("export", "methdb"),
                          ("export", "methdb", "-i", "x.methdb")):
            self.assertEqual(run_module(*arguments, site=False).returncode, 2)


class RunSettingsTests(unittest.TestCase):
    def test_defaults(self):
        resolved = run_settings("wgbs")
        self.assertEqual((resolved.reads, resolved.paired_end, resolved.read_length),
                         (1_000_000, True, (100,)))
        self.assertEqual((resolved.insert_min, resolved.insert_mean, resolved.insert_max,
                          resolved.insert_sd), (100, 400, 1000, 25.0))
        self.assertEqual((resolved.mutation_rate, resolved.conversion_rate, resolved.phred,
                          resolved.error_rate), (0.001, 0.998, 40, 0.005))
        self.assertEqual(resolved.reference, Path("ref.fa").resolve())
        self.assertEqual(resolved.output, Path("out").resolve())

    def test_omitted_seeds_are_generated_and_derived(self):
        resolved = run_settings("wgbs")
        self.assertEqual(resolved.seed, 7)
        self.assertEqual(resolved.seed_mut, derive_stage_seed(7, "mutation"))
        self.assertEqual(resolved.seed_phase, derive_stage_seed(7, "phasing"))
        self.assertEqual(resolved.seed_meth, derive_stage_seed(7, "methylation"))
        explicit = run_settings("wgbs", "--seed", "5", "--seed-mut", "9")
        self.assertEqual((explicit.seed, explicit.seed_mut), (5, 9))
        self.assertEqual(explicit.seed_phase, derive_stage_seed(5, "phasing"))
        with self.assertRaises(SystemExit), redirect_stderr(io.StringIO()):
            run_settings("wgbs", "--seed", str(1 << 64))

    def test_fixed_inserts_need_no_bounds(self):
        fixed = run_settings("wgbs", "--insert-sd", "0", "--insert-mean", "250")
        self.assertEqual((fixed.insert_min, fixed.insert_max), (250, 250))
        rrbs = run_settings("rrbs", "--insert-sd", "0", "--insert-mean", "250")
        self.assertEqual((rrbs.insert_min, rrbs.insert_max), (100, 1000))

    def test_mates_may_have_different_lengths(self):
        mates = run_settings("wgbs", "-l", "150,100")
        self.assertEqual((mates.read_length, mates.read_lengths), ((150, 100), (150, 100)))
        self.assertIn("150,100", HtsimCore(Path("htsim")).argv(mates))
        same = run_settings("wgbs", "-l", "120,120")
        self.assertEqual((same.read_length, same.read_lengths), ((120,), (120, 120)))
        for value in ("150,", "1,2,3", "x", "-5"):
            with self.assertRaises(SystemExit, msg=value), redirect_stderr(io.StringIO()):
                run_settings("wgbs", "-l", value)

    def test_fixed_inputs_disable_default_mutations(self):
        for option in ("--vcf", "--methdb", "--gc-profile"):
            self.assertEqual(run_settings("wgbs", option, "x").mutation_rate, 0.0)
        self.assertEqual(run_settings("rrbs", "--rrbs-candidates", "x").mutation_rate, 0.0)
        self.assertEqual(run_settings("wgbs", "--mutation-rate", "0.01").mutation_rate, 0.01)

    def test_reads_must_form_whole_fragments(self):
        self.assertEqual(run_settings("wgbs", "-n", "4").fragments, 2)
        self.assertEqual(run_settings("wgbs", "-n", "3", "--single-end").fragments, 3)
        with self.assertRaises(BSReadSimError):
            run_settings("wgbs", "-n", "3")
        depth = run_settings("wgbs", "-d", "2.5")
        self.assertEqual((depth.reads, depth.depth, depth.fragments), (None, 2.5, None))

    def test_standard_assays_have_no_bisulfite_chemistry(self):
        resolved = run_settings("wgs", "--save-truth")
        self.assertEqual((resolved.conversion_rate, resolved.directional), (0.0, True))
        self.assertEqual((resolved.save_vcf, resolved.save_methdb), (True, False))
        with self.assertRaises(SystemExit), redirect_stderr(io.StringIO()):
            run_settings("wgs", "--cgmap", "x")

    def test_python_owned_options_are_checked(self):
        for arguments in (("--fragment-summary",), ("--prefix", "a/b"), ("--threads", "0"),
                          ("--conversion-rate", "2"), ("--phred", "94")):
            with self.assertRaises(BSReadSimError, msg=arguments):
                run_settings("wgbs", *arguments)
        realization = run_settings("wgbs", "--format", "bam", "--fragment-realization")
        self.assertTrue(realization.fragment_summary)
        model = run_settings("wgbs", "--quality-model", "q.json")
        self.assertIsNone(model.phred)

    def test_cut_sites(self):
        self.assertEqual(run_settings("rrbs").cut_sites, ("C|CGG",))
        self.assertEqual(run_settings("rrbs", "--cut-site", "c|cgg,T|CGA").cut_sites,
                         ("C|CGG", "T|CGA"))
        with self.assertRaises(BSReadSimError):
            run_settings("rrbs", "--cut-site", "C|CGG,c|cgg")
        self.assertEqual(run_settings("wgbs").cut_sites, ())


class BuildSettingsTests(unittest.TestCase):
    def test_build_seeds_default_to_zero(self):
        resolved = settings("build", "methdb", "-r", "ref.fa", "-o", "x.methdb")
        self.assertEqual((resolved.seed, resolved.seed_mut, resolved.seed_phase,
                          resolved.seed_meth), (0, 0, 0, 0))
        self.assertEqual(resolved.mutation_rate, 0.001)

    def test_build_rrbs_and_variants(self):
        rrbs = settings("build", "rrbs", "-r", "ref.fa", "-o", "c.bed", "-l", "36")
        self.assertEqual((rrbs.technology, rrbs.read_length, rrbs.mutation_rate),
                         ("RRBS", (36,), 0.0))
        variants = settings("build", "variants", "-r", "ref.fa", "-o", "v.vcf.gz",
                            "--vcf", "in.vcf")
        self.assertEqual((variants.mutation_rate, variants.vcf),
                         (0.0, Path("in.vcf").resolve()))


class CoreArgumentTests(unittest.TestCase):
    def test_core_options_mirror_the_settings(self):
        resolved = run_settings("rrbs", "--seed", "3", "--cut-site", "C|CGG,T|CGA",
                                "--single-end", "-n", "10")
        argv = HtsimCore(Path("htsim")).argv(resolved, threads=4)
        options = dict(zip(argv[1::2], argv[2::2], strict=True))
        self.assertEqual(argv[0], "htsim")
        self.assertEqual(options["--technology"], "RRBS")
        self.assertEqual(options["--cut-sites"], "C|CGG,T|CGA")
        self.assertEqual(options["--paired-end"], "false")
        self.assertEqual(options["--fragments"], "10")
        self.assertEqual(options["--threads"], "4")
        self.assertNotIn("--config-sha256", options)
        self.assertNotIn("--vcf", options)
        self.assertNotIn("--depth", options)
        self.assertNotIn("--conversion-rate", options)

    def test_full_command_makes_seeds_explicit(self):
        argv = ["bsreadsim", "run", "wgs", "-r", "ref.fa", "-o", "out", "--seed-mut", "4"]
        resolved = from_arguments(build_parser().parse_args(argv[1:]), entropy=lambda bits: 7)
        self.assertEqual(full_command(resolved, argv), " ".join(argv) + " --seed 7 --seed-phase {}"
                         .format(derive_stage_seed(7, "phasing")))


if __name__ == "__main__":
    unittest.main()
