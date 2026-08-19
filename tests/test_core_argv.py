"""Tests for the normalized-config to C++ core argv boundary."""

from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim.config import ConfigValidationError, normalize_run_config  # noqa: E402
from bsreadsim.core_argv import (  # noqa: E402
    CoreArgvError,
    build_core_argv,
)
from bsreadsim.preparation import FileDigest, PreparedRun, prepare_run  # noqa: E402


RUN_ID = "00000000-0000-4000-8000-0000000000ab"


def base_config(technology="WGBS"):
    config = {
        "schema_version": "1.0",
        "seed": "0",
        "reference": "reference.fa",
        "inputs": {},
        "technology": technology,
        "mutation": {},
        "fragments": {
            "paired_end": False,
            "read_length_1": 4,
            "insert_min": 4,
            "insert_mean": 6,
            "insert_max": 10,
            "read_pairs": 2,
        },
        "methylation": {
            "beta": {
                "CG": [0.5, 0.5],
                "CHG": [0.05, 0.10],
                "CHH": [0.01, 0.20],
            },
        },
        "sequencing": {
            "conversion_rate": 0.998,
            "directional": True,
            "quality": {"kind": "uniform", "phred": 40},
            "error": {"kind": "uniform", "rate": 0.0},
        },
        "execution": {},
        "output": {
            "directory": "output",
            "prefix": "sample",
            "compression": "gzip",
        },
    }
    if technology == "RRBS":
        config["rrbs"] = {"cut_sites": ["C|CGG", "CCTN|AGG"]}
    elif technology == "TBS":
        config["tbs"] = {
            "bed": "targets.bed",
            "fragment_center_stddev": 12.5,
        }
    return config


def option_value(argv, option):
    index = argv.index(option)
    return argv[index + 1]


class CoreArgvTests(unittest.TestCase):
    def setUp(self):
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.directory = Path(self.temporary_directory.name)
        self.reference_bytes = b">chr1\nACGT\n"
        (self.directory / "reference.fa").write_bytes(self.reference_bytes)

    def prepared(self, document=None):
        config = base_config() if document is None else document
        loaded = normalize_run_config(config, self.directory)
        return prepare_run(loaded)

    def artifact(self, filename, artifact_format, version):
        payload = ("artifact:" + filename).encode("utf-8")
        (self.directory / filename).write_bytes(payload)
        return {
            "path": filename,
            "format": artifact_format,
            "version": version,
            "sha256": hashlib.sha256(payload).hexdigest(),
        }

    def test_wgbs_single_end_projection_is_exact_and_deterministic(self):
        prepared = self.prepared()
        reference_digest = hashlib.sha256(self.reference_bytes).hexdigest()

        observed = build_core_argv(prepared, RUN_ID, Path("/opt/bin/htsim-core"))
        expected = (
            "/opt/bin/htsim-core",
            "--truth-columns",
            "none",
            "--run-id",
            RUN_ID,
            "--config-sha256",
            prepared.config.sha256,
            "--seed",
            "0",
            "--protocol-batch-fragments",
            "64",
            "--reference",
            str((self.directory / "reference.fa").resolve()),
            "--reference-sha256",
            reference_digest,
            "--technology",
            "WGBS",
            "--paired-end",
            "false",
            "--read-length-1",
            "4",
            "--insert-min",
            "4",
            "--insert-mean",
            "6",
            "--insert-max",
            "10",
            "--insert-stddev",
            "25",
            "--read-pairs",
            "2",
            "--max-ambiguous-fraction",
            "0.05",
            "--chunk-size",
            "10000",
            "--core-workers",
            "1",
            "--coverage",
            "uniform",
            "--mutation-rate",
            "0.001",
            "--indel-fraction",
            "0.15",
            "--indel-extension-probability",
            "0.15",
            "--homozygous-only",
            "false",
            "--collect-non-cpg",
            "true",
            "--cgmap-pool",
            "false",
            "--update-variant-boundaries",
            "true",
            "--beta-cg",
            "0.5,0.5",
            "--beta-chg",
            "0.05,0.1",
            "--beta-chh",
            "0.01,0.2",
        )

        self.assertEqual(observed, expected)
        self.assertEqual(build_core_argv(prepared, RUN_ID, expected[0]), expected)
        self.assertNotIn("--read-length-2", observed)
        self.assertNotIn("--rrbs-cut-site", observed)
        self.assertNotIn("--tbs-bed", observed)
        self.assertNotIn("--coverage-profile", observed)

    def test_truth_policy_and_batch_size_are_runtime_options(self):
        prepared = self.prepared()

        argv = build_core_argv(
            prepared,
            RUN_ID,
            "htsim-core",
            truth_columns="none",
        )

        self.assertNotIn("--protocol-major", argv)
        self.assertEqual(option_value(argv, "--truth-columns"), "none")
        self.assertEqual(option_value(argv, "--protocol-batch-fragments"), "64")
        self.assertNotIn("protocol", prepared.config.normalized)

        bounded = build_core_argv(
            prepared,
            RUN_ID,
            "htsim-core",
            protocol_batch_fragments=7,
        )
        self.assertEqual(option_value(bounded, "--protocol-batch-fragments"), "7")

        for truth_columns in ("", "FULL", None):
            with self.subTest(truth_columns=truth_columns):
                with self.assertRaisesRegex(CoreArgvError, "truth_columns"):
                    build_core_argv(
                        prepared,
                        RUN_ID,
                        "htsim-core",
                        truth_columns=truth_columns,
                    )
        for batch_fragments in (0, 65, True, "7"):
            with self.subTest(protocol_batch_fragments=batch_fragments):
                with self.assertRaisesRegex(
                    CoreArgvError, "protocol_batch_fragments"
                ):
                    build_core_argv(
                        prepared,
                        RUN_ID,
                        "htsim-core",
                        protocol_batch_fragments=batch_fragments,
                    )

    def test_core_worker_count_crosses_only_the_core_boundary(self):
        document = base_config()
        document["execution"]["core_workers"] = 4
        document["execution"]["workers"] = 12
        prepared = self.prepared(document)

        argv = build_core_argv(prepared, RUN_ID, "htsim-core")

        self.assertEqual(option_value(argv, "--core-workers"), "4")
        self.assertNotIn("--workers", argv)

    def test_tbs_profile_projects_all_core_inputs_and_hides_python_models(self):
        document = base_config("TBS")
        document["fragments"] = {
            "paired_end": True,
            "read_length_1": 4,
            "read_length_2": 4,
            "insert_min": 4,
            "insert_mean": 6,
            "insert_max": 10,
            "insert_stddev": 2.5,
            "depth": 3.25,
        }
        for filename in (
            "targets.bed",
            "sample; literal.vcf",
            "sample.cgmap",
            "sample.asm",
        ):
            (self.directory / filename).write_bytes(filename.encode("utf-8"))
        document["inputs"] = {
            "vcf": "sample; literal.vcf",
            "cgmap": "sample.cgmap",
            "asm": "sample.asm",
        }
        coverage_artifact = self.artifact("coverage.json", "json", "coverage-v1")
        quality_artifact = self.artifact("quality.npz", "npz", "quality-v8")
        error_artifact = self.artifact("error.json", "json", "error-v7")
        document["coverage"] = {
            "kind": "profile",
            "artifact": coverage_artifact,
        }
        document["sequencing"]["quality"] = {
            "kind": "markov",
            "artifact": quality_artifact,
        }
        document["sequencing"]["error"] = {
            "kind": "quality-confusion",
            "artifact": error_artifact,
        }
        prepared = self.prepared(document)

        argv = build_core_argv(prepared, RUN_ID, "htsim-core")

        self.assertEqual(option_value(argv, "--technology"), "TBS")
        self.assertEqual(option_value(argv, "--paired-end"), "true")
        self.assertEqual(option_value(argv, "--read-length-2"), "4")
        self.assertEqual(option_value(argv, "--depth"), "3.25")
        self.assertNotIn("--read-pairs", argv)
        for name in ("vcf", "cgmap", "asm"):
            role = prepared.file_for_role("input." + name)
            self.assertEqual(option_value(argv, "--" + name), str(role.path))
            self.assertEqual(option_value(argv, "--" + name + "-sha256"), role.sha256)
        tbs = prepared.file_for_role("input.tbs-bed")
        self.assertEqual(option_value(argv, "--tbs-bed"), str(tbs.path))
        self.assertEqual(option_value(argv, "--tbs-bed-sha256"), tbs.sha256)
        self.assertEqual(option_value(argv, "--tbs-center-stddev"), "12.5")
        profile = prepared.file_for_role("model.coverage")
        self.assertEqual(option_value(argv, "--coverage-profile"), str(profile.path))
        self.assertEqual(option_value(argv, "--coverage-profile-format"), "json")
        self.assertEqual(
            option_value(argv, "--coverage-profile-version"), "coverage-v1"
        )

        self.assertEqual(option_value(argv, "--coverage-profile-sha256"), profile.sha256)

        python_only_options = {
            "--conversion-rate",
            "--directional",
            "--quality",
            "--error",
            "--workers",
            "--max-in-flight-fragments",
            "--output",
            "--compression",
        }
        self.assertTrue(python_only_options.isdisjoint(argv))
        for artifact in (quality_artifact, error_artifact):
            self.assertNotIn(str((self.directory / artifact["path"]).resolve()), argv)
            self.assertNotIn(artifact["version"], argv)
        self.assertIn(
            str((self.directory / "sample; literal.vcf").resolve()), argv
        )
        self.assertTrue(all(isinstance(argument, str) for argument in argv))

    def test_bed_methyl_inputs_project_distinct_core_options_and_digests(self):
        for filename in ("sample.vcf", "levels.bedmethyl", "levels.asm.bed"):
            (self.directory / filename).write_bytes(filename.encode("utf-8"))
        document = base_config()
        document["mutation"]["rate"] = 0
        document["inputs"] = {
            "vcf": "sample.vcf",
            "bed_methyl": "levels.bedmethyl",
            "asm_bed": "levels.asm.bed",
        }
        prepared = self.prepared(document)

        argv = build_core_argv(prepared, RUN_ID, "htsim-core")

        for name in ("vcf", "bed_methyl", "asm_bed"):
            role = prepared.file_for_role("input." + name)
            option = "--" + name.replace("_", "-")
            self.assertEqual(option_value(argv, option), str(role.path))
            self.assertEqual(option_value(argv, option + "-sha256"), role.sha256)
        self.assertNotIn("--cgmap", argv)
        self.assertNotIn("--asm", argv)

    def test_rrbs_cut_sites_repeat_in_normalized_order(self):
        prepared = self.prepared(base_config("RRBS"))

        argv = build_core_argv(prepared, RUN_ID, "htsim-core")

        sites = [
            argv[index + 1]
            for index, argument in enumerate(argv[:-1])
            if argument == "--rrbs-cut-site"
        ]
        self.assertEqual(sites, ["C|CGG", "CCTN|AGG"])
        self.assertNotIn("--tbs-bed", argv)

    def test_rrbs_candidate_profile_projects_path_without_hash_option(self):
        (self.directory / "rrbs-candidates.bed").write_text(
            "chr1\t0\t4\tchr1:0-4\t0.5\t.\t3\t4\t2\t2\n",
            encoding="utf-8",
        )
        document = base_config("RRBS")
        document["mutation"]["rate"] = 0
        document["rrbs"]["candidate_bed"] = "rrbs-candidates.bed"
        document["coverage"] = {"kind": "profile"}
        prepared = self.prepared(document)

        argv = build_core_argv(prepared, RUN_ID, "htsim-core")

        self.assertEqual(
            option_value(argv, "--rrbs-candidate-bed"),
            str((self.directory / "rrbs-candidates.bed").resolve()),
        )
        self.assertEqual(option_value(argv, "--coverage"), "profile")
        self.assertNotIn("--rrbs-candidate-bed-sha256", argv)
        self.assertNotIn("--coverage-profile", argv)
        self.assertTrue(
            {"input.rrbs-candidate-bed", "model.coverage"}.isdisjoint(
                {item.role for item in prepared.files}
            )
        )

    def test_tbs_target_score_projects_without_a_second_artifact(self):
        (self.directory / "targets.bed").write_text(
            "chr1\t1\t2\ttarget\t3\t+\n", encoding="utf-8"
        )
        document = base_config("TBS")
        document["coverage"] = {"kind": "target-score"}
        prepared = self.prepared(document)

        argv = build_core_argv(prepared, RUN_ID, "htsim-core")

        self.assertEqual(option_value(argv, "--technology"), "TBS")
        self.assertEqual(option_value(argv, "--coverage"), "target-score")
        self.assertNotIn("--coverage-profile", argv)
        self.assertNotIn("model.coverage", {
            item.role for item in prepared.files
        })

    def test_json_number_and_boolean_spellings_are_stable(self):
        document = base_config()
        document["fragments"].pop("read_pairs")
        document["fragments"]["depth"] = 2.5
        document["fragments"]["insert_stddev"] = -0.0
        document["mutation"]["rate"] = 1e-7
        document["mutation"]["indel_fraction"] = 1.0
        document["mutation"]["homozygous_only"] = True
        document["methylation"]["collect_non_cpg"] = False
        prepared = self.prepared(document)

        argv = build_core_argv(prepared, RUN_ID, "htsim-core")

        self.assertEqual(option_value(argv, "--depth"), json.dumps(2.5))
        self.assertEqual(option_value(argv, "--insert-stddev"), json.dumps(-0.0))
        self.assertEqual(option_value(argv, "--mutation-rate"), json.dumps(1e-7))
        self.assertEqual(option_value(argv, "--indel-fraction"), json.dumps(1.0))
        self.assertEqual(option_value(argv, "--homozygous-only"), "true")
        self.assertEqual(option_value(argv, "--collect-non-cpg"), "false")

    def test_invalid_identity_and_executable_fail_closed(self):
        prepared = self.prepared()
        for run_id in (
            "not-a-uuid",
            RUN_ID.upper(),
            "{" + RUN_ID + "}",
            None,
        ):
            with self.subTest(run_id=run_id):
                with self.assertRaisesRegex(CoreArgvError, "canonical lowercase UUID"):
                    build_core_argv(prepared, run_id, "htsim-core")

        for executable in ("", "bad\x00path", b"htsim-core", 3):
            with self.subTest(executable=executable):
                with self.assertRaises(CoreArgvError):
                    build_core_argv(prepared, RUN_ID, executable)

        with self.assertRaisesRegex(CoreArgvError, "PreparedRun"):
            build_core_argv(object(), RUN_ID, "htsim-core")

    def test_missing_duplicate_and_unexpected_roles_fail_closed(self):
        prepared = self.prepared()
        reference = prepared.file_for_role("reference")
        cases = (
            (
                replace(prepared, files=()),
                "missing prepared file roles",
            ),
            (
                replace(prepared, files=prepared.files + (reference,)),
                "duplicate prepared file role",
            ),
            (
                replace(
                    prepared,
                    files=prepared.files
                    + (
                        FileDigest(
                            role="unexpected",
                            path=self.directory / "unexpected",
                            size_bytes=0,
                            sha256="0" * 64,
                        ),
                    ),
                ),
                "unexpected prepared file roles",
            ),
            (
                PreparedRun(config=prepared.config, files=list(prepared.files)),
                "files must be a tuple",
            ),
        )
        for invalid, message in cases:
            with self.subTest(message=message):
                with self.assertRaisesRegex(CoreArgvError, message):
                    build_core_argv(invalid, RUN_ID, "htsim-core")

    def test_role_paths_digests_and_config_identity_cannot_be_forged(self):
        prepared = self.prepared()
        reference = prepared.file_for_role("reference")
        forged_files = (
            replace(reference, path=self.directory / "other.fa"),
            replace(reference, sha256="A" * 64),
            replace(reference, size_bytes=-1),
            replace(reference, declared_sha256="0" * 64),
        )
        for forged in forged_files:
            invalid = replace(prepared, files=(forged,))
            with self.subTest(forged=forged):
                with self.assertRaises(CoreArgvError):
                    build_core_argv(invalid, RUN_ID, "htsim-core")

        with self.assertRaisesRegex(ConfigValidationError, "SHA-256"):
            replace(prepared.config, sha256="0" * 64)

        unmaterialized_document = base_config()
        unmaterialized_document.pop("seed")
        unmaterialized = normalize_run_config(
            unmaterialized_document, self.directory
        )
        with self.assertRaises(CoreArgvError):
            build_core_argv(
                PreparedRun(config=unmaterialized, files=()),
                RUN_ID,
                "htsim-core",
            )


if __name__ == "__main__":
    unittest.main()
