"""Tests for CLI document normalization and stable run identity."""

import copy
import hashlib
from pathlib import Path
import tempfile
import unittest


from bsreadsim.run.config import (
    ConfigValidationError,
    UINT64_MAX,
    normalize_run_config,
)


MODEL_SHA = "1" * 64


def base_config(technology: str = "WGBS") -> dict:
    config = {
        "reference": "inputs/reference.fa",
        "inputs": {},
        "technology": technology,
        "mutation": {},
        "seeds": {"mutation": "0", "phasing": "0", "methylation": "0"},
        "fragments": {
            "paired_end": True,
            "read_length_1": 100,
            "read_length_2": 100,
            "insert_min": 100,
            "insert_max": 1000,
            "depth": 20,
        },
        "methylation": {
            "beta": {
                "CG": [0.5, 0.5],
                "CHG": [0.01, 0.05],
                "CHH": [0.01, 0.05],
            },
        },
        "sequencing": {
            "conversion_rate": 0.998,
            "directional": True,
            "quality": {"kind": "uniform", "phred": 40},
            "error": {"kind": "uniform", "rate": 0.005},
        },
        "execution": {},
        "output": {"directory": "results", "prefix": "sim"},
    }
    if technology == "RRBS":
        config["rrbs"] = {"cut_sites": ["C|CGG"]}
    elif technology in ("TBS", "WES", "TS"):
        config["tbs"] = {"bed": "inputs/targets.bed"}
    if technology in ("WGS", "WES", "TS"):
        config["sequencing"]["conversion_rate"] = 0.0
    return config


class NormalizedConfigTests(unittest.TestCase):
    def test_json_run_file_apis_are_not_available(self) -> None:
        from bsreadsim.run import config as config_module
        from bsreadsim.run import execute as execute_module

        self.assertFalse(hasattr(config_module, "load_run_config"))
        self.assertFalse(hasattr(execute_module, "run_config"))

    def test_wgbs_is_the_normalized_default(self) -> None:
        document = base_config()
        document.pop("technology")

        loaded = normalize_run_config(document, self.base_directory)

        self.assertEqual(loaded.normalized["technology"], "WGBS")
        self.assertNotIn("site_model", loaded.normalized["methylation"])

    def test_unimplemented_site_model_config_is_rejected(self) -> None:
        document = base_config()
        document["methylation"]["site_model"] = {"kind": "bernoulli"}

        with self.assertRaises(ConfigValidationError):
            normalize_run_config(document, self.base_directory)

    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.base_directory = Path(self.temporary_directory.name)

    def test_valid_wgbs_materializes_all_defaults(self) -> None:
        loaded = normalize_run_config(base_config(), self.base_directory)

        self.assertEqual(loaded.normalized["coverage"], {"kind": "uniform"})
        self.assertEqual(
            loaded.normalized["mutation"],
            {
                "rate": 0.001,
                "indel_fraction": 0.15,
                "indel_extension_probability": 0.15,
                "homozygous_only": False,
            },
        )
        self.assertEqual(loaded.normalized["fragments"]["insert_mean"], 400)
        self.assertEqual(loaded.normalized["fragments"]["insert_sd"], 25)
        self.assertEqual(loaded.normalized["execution"]["workers"], 1)
        self.assertEqual(loaded.normalized["execution"]["core_workers"], 1)
        self.assertEqual(loaded.normalized["execution"]["chunk_size"], 10000)
        self.assertEqual(
            loaded.normalized["execution"]["max_in_flight_fragments"], 4096
        )
        self.assertEqual(loaded.normalized["output"]["format"], "fastq.gz")
        self.assertEqual(loaded.normalized["output"]["gzip_level"], 6)
        self.assertFalse(loaded.normalized["output"]["save_methdb"])
        self.assertFalse(loaded.normalized["output"]["save_vcf"])
        self.assertNotIn("details", loaded.normalized["output"])
        self.assertNotIn("shuffle", loaded.normalized["output"])
        self.assertNotIn("overwrite", loaded.normalized["output"])

    def test_valid_rrbs_and_tbs_sections(self) -> None:
        rrbs = normalize_run_config(base_config("RRBS"), self.base_directory)
        tbs_config = base_config("TBS")
        tbs_config["coverage"] = {"kind": "target-score"}
        tbs = normalize_run_config(tbs_config, self.base_directory)

        self.assertEqual(rrbs.normalized["rrbs"]["cut_sites"], ["C|CGG"])
        self.assertEqual(tbs.normalized["tbs"]["fragment_center_stddev"], 50)
        self.assertEqual(tbs.normalized["coverage"], {"kind": "target-score"})
        self.assertEqual(
            tbs.normalized["tbs"]["bed"],
            str((self.base_directory / "inputs/targets.bed").resolve()),
        )

    def test_rrbs_candidate_bed_profile_has_no_model_artifact(self) -> None:
        config = base_config("RRBS")
        config["seed"] = "7"
        config["mutation"]["rate"] = 0
        config["rrbs"]["candidate_bed"] = "profiles/rrbs-candidates.bed"
        config["coverage"] = {"kind": "profile"}

        loaded = normalize_run_config(config, self.base_directory)

        self.assertEqual(loaded.normalized["coverage"], {"kind": "profile"})
        self.assertEqual(
            loaded.normalized["rrbs"]["candidate_bed"],
            str((self.base_directory / "profiles/rrbs-candidates.bed").resolve()),
        )

        missing = base_config("RRBS")
        missing["coverage"] = {"kind": "profile"}
        with self.assertRaisesRegex(ConfigValidationError, "candidate BED"):
            normalize_run_config(missing, self.base_directory)

        wgbs = base_config("WGBS")
        wgbs["coverage"] = {"kind": "profile"}
        with self.assertRaisesRegex(ConfigValidationError, "requires RRBS"):
            normalize_run_config(wgbs, self.base_directory)

    def test_rrbs_variant_candidate_bed_preserves_split_seed_defaults(self) -> None:
        config = base_config("RRBS")
        config["rrbs"]["candidate_bed"] = "profiles/rrbs-candidates.bed"
        config["inputs"]["vcf"] = "inputs/sample.vcf"
        config["mutation"]["rate"] = 0

        normalized = normalize_run_config(config, self.base_directory).normalized
        self.assertEqual(
            normalized["seeds"],
            {"mutation": "0", "phasing": "0", "methylation": "0"},
        )

    def test_paths_resolve_against_cli_invocation_directory(self) -> None:
        config = base_config()
        config["inputs"] = {
            "vcf": "inputs/sample.vcf.gz",
            "cgmap": "../profiles/sample.CGmap.gz",
            "asm": "~/literal-path/sample.asm.gz",
        }
        config["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "models/coverage.json",
                "sha256": MODEL_SHA,
            },
        }
        invocation_directory = self.base_directory / "invocation"
        invocation_directory.mkdir()
        loaded = normalize_run_config(config, invocation_directory)

        self.assertEqual(
            loaded.normalized["reference"],
            str((invocation_directory / "inputs/reference.fa").resolve()),
        )
        self.assertEqual(
            loaded.normalized["inputs"]["vcf"],
            str((invocation_directory / "inputs/sample.vcf.gz").resolve()),
        )
        self.assertEqual(
            loaded.normalized["inputs"]["cgmap"],
            str((invocation_directory / "../profiles/sample.CGmap.gz").resolve()),
        )
        self.assertEqual(
            loaded.normalized["inputs"]["asm"],
            str((invocation_directory / "~/literal-path/sample.asm.gz").resolve()),
        )
        self.assertEqual(
            loaded.normalized["coverage"]["artifact"]["path"],
            str((invocation_directory / "models/coverage.json").resolve()),
        )
        self.assertEqual(
            loaded.normalized["output"]["directory"],
            str((invocation_directory / "results").resolve()),
        )

    def test_bed_methyl_inputs_resolve_and_remain_mutually_exclusive(self) -> None:
        config = base_config()
        config["inputs"] = {
            "vcf": "inputs/sample.vcf",
            "bed_methyl": "profiles/sample.bedmethyl.gz",
            "asm_bed": "profiles/sample.asm.bed.gz",
        }

        loaded = normalize_run_config(config, self.base_directory)

        self.assertEqual(
            loaded.normalized["inputs"]["bed_methyl"],
            str((self.base_directory / "profiles/sample.bedmethyl.gz").resolve()),
        )
        self.assertEqual(
            loaded.normalized["inputs"]["asm_bed"],
            str((self.base_directory / "profiles/sample.asm.bed.gz").resolve()),
        )

        both_levels = base_config()
        both_levels["inputs"] = {
            "cgmap": "levels.cgmap",
            "bed_methyl": "levels.bedmethyl",
        }
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(both_levels, self.base_directory)

        both_asm = base_config()
        both_asm["inputs"] = {
            "vcf": "sample.vcf",
            "asm": "levels.asm",
            "asm_bed": "levels.asm.bed",
        }
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(both_asm, self.base_directory)

        missing_vcf = base_config()
        missing_vcf["inputs"] = {"asm_bed": "levels.asm.bed"}
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(missing_vcf, self.base_directory)

    def test_seed_zero_and_maximum_remain_decimal_strings(self) -> None:
        for seed in (0, UINT64_MAX):
            with self.subTest(seed=seed):
                config = base_config()
                config["seed"] = str(seed)
                loaded = normalize_run_config(config, self.base_directory)
                self.assertEqual(loaded.master_seed, seed)
                self.assertEqual(loaded.normalized["seed"], str(seed))

    def test_seed_above_u64_is_rejected(self) -> None:
        config = base_config()
        config["seed"] = str(UINT64_MAX + 1)
        with self.assertRaisesRegex(ConfigValidationError, "unsigned 64-bit"):
            normalize_run_config(config, self.base_directory)

    def test_output_session_annotation_is_not_a_user_config_field(self) -> None:
        config = base_config()
        config["output"]["details"] = "full"
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(config, self.base_directory)

    def test_gzip_level_is_normalized_and_schema_validated(self) -> None:
        default_loaded = normalize_run_config(base_config(), self.base_directory)
        config = base_config()
        config["output"]["gzip_level"] = 1
        loaded = normalize_run_config(config, self.base_directory)
        self.assertEqual(loaded.normalized["output"]["gzip_level"], 1)
        self.assertNotEqual(loaded.sha256, default_loaded.sha256)

        explicit_default = base_config()
        explicit_default["output"]["gzip_level"] = 6
        self.assertEqual(
            normalize_run_config(explicit_default, self.base_directory).sha256,
            default_loaded.sha256,
        )
        for value in (-1, 10, True):
            with self.subTest(value=value):
                invalid = base_config()
                invalid["output"]["gzip_level"] = value
                with self.assertRaises(ConfigValidationError):
                    normalize_run_config(invalid, self.base_directory)

    def test_format_is_an_explicit_output_policy(self) -> None:
        enabled = base_config()
        enabled["output"]["format"] = "bam"
        loaded = normalize_run_config(enabled, self.base_directory)
        self.assertEqual(loaded.normalized["output"]["format"], "bam")
        self.assertNotIn("details", loaded.normalized["output"])

        invalid = base_config()
        invalid["output"]["format"] = "yes"
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(invalid, self.base_directory)

    def test_user_annotation_field_is_rejected(self) -> None:
        config = base_config()
        config["output"]["details"] = "full"
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(config, self.base_directory)

    def test_core_worker_bounds_are_schema_validated(self) -> None:
        for value in (0, 65):
            with self.subTest(value=value):
                config = base_config()
                config["execution"]["core_workers"] = value
                with self.assertRaises(ConfigValidationError):
                    normalize_run_config(config, self.base_directory)

    def test_depth_and_fragment_count_are_mutually_exclusive(self) -> None:
        config = base_config()
        config["fragments"]["count"] = 100
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(config, self.base_directory)

    def test_insert_relationship_and_read_lengths_are_checked(self) -> None:
        invalid_relationship = base_config()
        invalid_relationship["fragments"]["insert_min"] = 500
        invalid_relationship["fragments"]["insert_max"] = 600
        invalid_read_length = base_config()
        invalid_read_length["fragments"]["insert_min"] = 10
        invalid_read_length["fragments"]["insert_mean"] = 100
        invalid_read_length["fragments"]["insert_max"] = 200

        with self.assertRaisesRegex(ConfigValidationError, "insert_min"):
            normalize_run_config(invalid_relationship, self.base_directory)
        with self.assertRaisesRegex(ConfigValidationError, "insert_min"):
            normalize_run_config(invalid_read_length, self.base_directory)

        asymmetric = base_config()
        asymmetric["fragments"]["read_length_2"] = 99
        normalized = normalize_run_config(asymmetric, self.base_directory)
        self.assertEqual(normalized.normalized["fragments"]["read_length_2"], 99)

        fixed_tbs = base_config("TBS")
        fixed_tbs["fragments"].update(
            {
                "read_length_1": 150,
                "read_length_2": 150,
                "insert_mean": 300,
                "insert_sd": 0,
            }
        )
        fixed = normalize_run_config(fixed_tbs, self.base_directory).normalized
        self.assertEqual(fixed["fragments"]["insert_mean"], 300)
        self.assertEqual(fixed["fragments"]["insert_min"], 100)

    def test_cpp_numeric_projection_bounds_are_checked(self) -> None:
        oversized_insert = base_config()
        oversized_insert["fragments"]["insert_max"] = 1 << 32
        oversized_number = base_config()
        oversized_number["fragments"]["insert_sd"] = 10**400

        for config in (oversized_insert, oversized_number):
            with self.subTest(config=config):
                with self.assertRaises(ConfigValidationError):
                    normalize_run_config(config, self.base_directory)

    def test_technology_sections_are_required_and_exclusive(self) -> None:
        rrbs_without_section = base_config()
        rrbs_without_section["technology"] = "RRBS"
        wgbs_with_tbs = base_config()
        wgbs_with_tbs["tbs"] = {"bed": "targets.bed"}

        for config in (rrbs_without_section, wgbs_with_tbs):
            with self.subTest(config=config):
                with self.assertRaises(ConfigValidationError):
                    normalize_run_config(config, self.base_directory)

        oversized_rrbs = base_config("RRBS")
        oversized_rrbs["rrbs"]["cut_sites"] = ["A" * 1025 + "|A"]
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(oversized_rrbs, self.base_directory)

    def test_standard_technologies_use_whole_genome_or_target_sections(self) -> None:
        for technology in ("WGS", "WES", "TS"):
            with self.subTest(technology=technology):
                normalized = normalize_run_config(
                    base_config(technology), self.base_directory
                ).normalized
                self.assertEqual(normalized["technology"], technology)
                self.assertEqual("tbs" in normalized, technology in ("WES", "TS"))

    def test_standard_technologies_reject_methylation_inputs_and_truth(self) -> None:
        methylation_input = base_config("WGS")
        methylation_input["inputs"]["cgmap"] = "sample.cgmap"
        methdb_truth = base_config("WGS")
        methdb_truth["output"]["save_methdb"] = True
        realization = base_config("WGS")
        realization["output"].update(
            {"format": "bam", "fragment_summary": True, "fragment_realization": True}
        )
        conversion = base_config("WGS")
        conversion["sequencing"]["conversion_rate"] = 0.5
        orientation = base_config("WGS")
        orientation["sequencing"]["directional"] = False

        for document in (
            methylation_input,
            methdb_truth,
            realization,
            conversion,
            orientation,
        ):
            with self.subTest(document=document):
                with self.assertRaisesRegex(
                    ConfigValidationError, "standard sequencing"
                ):
                    normalize_run_config(document, self.base_directory)

    def test_methdb_is_the_only_variant_authority(self) -> None:
        methdb = base_config()
        methdb["inputs"]["methdb"] = "inputs/profile.methdb"
        methdb["mutation"]["rate"] = 0
        normalized = normalize_run_config(
            methdb, self.base_directory
        ).normalized
        self.assertIn("methdb", normalized["inputs"])

        external_vcf = copy.deepcopy(methdb)
        external_vcf["inputs"]["vcf"] = "inputs/variants.vcf"
        de_novo = copy.deepcopy(methdb)
        de_novo["mutation"]["rate"] = 0.01
        for document in (external_vcf, de_novo):
            with self.subTest(document=document):
                with self.assertRaisesRegex(
                    ConfigValidationError, "embedded variants|VCF or overlays"
                ):
                    normalize_run_config(document, self.base_directory)

    def test_invalid_model_sha_and_conflicting_declarations_are_rejected(self) -> None:
        invalid_sha = base_config()
        invalid_sha["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "models/shared.json",
                "sha256": "ABC",
            },
        }
        with self.assertRaises(ConfigValidationError):
            normalize_run_config(invalid_sha, self.base_directory)

        conflicting = base_config()
        conflicting["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "models/shared.json",
                "sha256": "1" * 64,
            },
        }
        conflicting["sequencing"]["quality"] = {
            "kind": "markov",
            "artifact": {
                "path": "models/shared.json",
                "sha256": "2" * 64,
            },
        }
        with self.assertRaisesRegex(ConfigValidationError, "conflicting sha256"):
            normalize_run_config(conflicting, self.base_directory)

    def test_artifact_format_and_version_fields_are_not_accepted(self) -> None:
        config = base_config()
        config["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "models/coverage.tsv",
                "format": "tsv",
                "version": "obsolete",
                "sha256": MODEL_SHA,
            },
        }

        with self.assertRaises(ConfigValidationError):
            normalize_run_config(config, self.base_directory)

    def test_unknown_key_and_invalid_schema_type_are_rejected(self) -> None:
        unknown = base_config()
        unknown["unexpected"] = True
        wrong_type = base_config()
        wrong_type["execution"]["workers"] = "four"

        for config in (unknown, wrong_type):
            with self.subTest(config=config):
                with self.assertRaises(ConfigValidationError):
                    normalize_run_config(config, self.base_directory)

    def test_canonical_json_and_digest_are_stable(self) -> None:
        first = base_config()
        first["seed"] = "42"
        second = {key: copy.deepcopy(first[key]) for key in reversed(first)}

        loaded_first = normalize_run_config(first, self.base_directory)
        loaded_second = normalize_run_config(second, self.base_directory)

        self.assertEqual(loaded_first.canonical_json, loaded_second.canonical_json)
        self.assertEqual(loaded_first.sha256, loaded_second.sha256)
        self.assertEqual(
            loaded_first.sha256,
            hashlib.sha256(loaded_first.canonical_json.encode("utf-8")).hexdigest(),
        )
        self.assertNotIn(" ", loaded_first.canonical_json)
        self.assertNotIn("\n", loaded_first.canonical_json)

    def test_input_mapping_is_not_mutated(self) -> None:
        config = base_config()
        before = copy.deepcopy(config)
        normalize_run_config(config, self.base_directory)
        self.assertEqual(config, before)

if __name__ == "__main__":
    unittest.main()
