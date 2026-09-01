"""Tests for run seed materialization and immutable input hashing."""

import hashlib
from pathlib import Path
import tempfile
import unittest


from bsreadsim.run.config import normalize_run_config
from bsreadsim.run.prepare import (
    FileDigest,
    PreparationError,
    derive_stage_seed,
    materialize_master_seed,
    materialize_run_seeds,
    prepare_run,
    snapshot_prepared_file,
)


def base_config() -> dict:
    return {
        "reference": "reference.fa",
        "inputs": {},
        "technology": "WGBS",
        "mutation": {},
        "seeds": {"mutation": "0", "phasing": "0", "methylation": "0"},
        "reads": {"count": 2},
        "fragments": {
            "paired_end": False,
            "read_length_1": 4,
            "insert_min": 4,
            "insert_mean": 6,
            "insert_max": 10,
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
        "output": {"directory": "output", "prefix": "sample"},
    }


class PreparationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.base_directory = Path(self.temporary_directory.name)
        self.reference_bytes = b">chr1\nACGT\n"
        (self.base_directory / "reference.fa").write_bytes(self.reference_bytes)

    def normalized(self, document=None):
        return normalize_run_config(
            base_config() if document is None else document,
            self.base_directory,
        )

    def test_omitted_seed_is_materialized_before_identity_is_frozen(self) -> None:
        loaded = self.normalized()
        self.assertIsNone(loaded.master_seed)

        effective = materialize_master_seed(loaded, entropy=lambda bits: 42)

        self.assertEqual(effective.master_seed, 42)
        self.assertEqual(effective.normalized["seed"], "42")
        self.assertNotEqual(effective.sha256, loaded.sha256)
        self.assertIn('"seed":"42"', effective.canonical_json)

    def test_explicit_seed_never_consumes_entropy(self) -> None:
        document = base_config()
        document["seed"] = "0"
        loaded = self.normalized(document)

        def forbidden_entropy(bits):
            raise AssertionError("entropy was consumed for an explicit seed")

        self.assertIs(
            materialize_master_seed(loaded, entropy=forbidden_entropy), loaded
        )

    def test_omitted_stage_seeds_are_derived_from_the_master_seed(self) -> None:
        document = base_config()
        document["seeds"] = {
            "mutation": None,
            "phasing": None,
            "methylation": None,
        }

        effective = materialize_run_seeds(
            self.normalized(document), entropy=lambda bits: 42
        )

        self.assertEqual(effective.master_seed, 42)
        self.assertEqual(
            effective.normalized["seeds"],
            {
                "mutation": "7049704659189676143",
                "phasing": "15840450471693409671",
                "methylation": "7747208672942775485",
            },
        )

    def test_explicit_stage_seeds_override_derivation(self) -> None:
        document = base_config()
        document["seed"] = "42"
        document["seeds"] = {
            "mutation": "7",
            "phasing": None,
            "methylation": "11",
        }

        def forbidden_entropy(bits):
            raise AssertionError("entropy was consumed for an explicit seed")

        effective = materialize_run_seeds(
            self.normalized(document), entropy=forbidden_entropy
        )

        self.assertEqual(
            effective.normalized["seeds"],
            {
                "mutation": "7",
                "phasing": "15840450471693409671",
                "methylation": "11",
            },
        )

    def test_stage_seed_derivation_has_frozen_domain_vectors(self) -> None:
        self.assertEqual(
            derive_stage_seed(0, "mutation"), 9839267760111386065
        )
        self.assertEqual(
            derive_stage_seed((1 << 64) - 1, "methylation"),
            370557882438090985,
        )
        self.assertNotEqual(
            derive_stage_seed(7, "mutation"),
            derive_stage_seed(7, "phasing"),
        )
        with self.assertRaisesRegex(PreparationError, "domain"):
            derive_stage_seed(7, "reads")

    def test_invalid_entropy_result_is_rejected(self) -> None:
        loaded = self.normalized()
        for value in (-1, 1 << 64, True, "12"):
            with self.subTest(value=value):
                with self.assertRaisesRegex(PreparationError, "unsigned 64-bit"):
                    materialize_master_seed(
                        loaded, entropy=lambda bits, value=value: value
                    )

    def test_prepare_hashes_inputs_and_verifies_model_artifact(self) -> None:
        vcf_bytes = b"##fileformat=VCFv4.2\n"
        model_bytes = b'{"kind":"coverage"}\n'
        (self.base_directory / "sample.vcf").write_bytes(vcf_bytes)
        (self.base_directory / "coverage.json").write_bytes(model_bytes)
        document = base_config()
        document["seed"] = "7"
        document["inputs"] = {"vcf": "sample.vcf"}
        document["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "coverage.json",
                "sha256": hashlib.sha256(model_bytes).hexdigest(),
            },
        }

        prepared = prepare_run(self.normalized(document), hash_chunk_size=3)

        self.assertEqual(
            tuple(file.role for file in prepared.files),
            ("reference", "input.vcf", "model.coverage"),
        )
        reference = prepared.file_for_role("reference")
        self.assertEqual(reference.size_bytes, len(self.reference_bytes))
        self.assertEqual(
            reference.sha256, hashlib.sha256(self.reference_bytes).hexdigest()
        )
        model = prepared.file_for_role("model.coverage")
        self.assertEqual(model.declared_sha256, model.sha256)

    def test_declared_model_digest_mismatch_fails_before_launch(self) -> None:
        (self.base_directory / "coverage.json").write_bytes(b"actual")
        document = base_config()
        document["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "coverage.json",
                "sha256": "0" * 64,
            },
        }

        with self.assertRaisesRegex(PreparationError, "SHA-256 mismatch"):
            prepare_run(self.normalized(document), entropy=lambda bits: 1)

    def test_missing_and_non_regular_inputs_fail_closed(self) -> None:
        missing = base_config()
        missing["reference"] = "missing.fa"
        with self.assertRaisesRegex(PreparationError, "cannot hash input"):
            prepare_run(self.normalized(missing), entropy=lambda bits: 1)

        directory = base_config()
        directory["reference"] = "."
        with self.assertRaisesRegex(PreparationError, "not a regular file"):
            prepare_run(self.normalized(directory), entropy=lambda bits: 1)

    def test_hash_chunk_size_is_validated(self) -> None:
        loaded = self.normalized()
        for value in (0, -1, True, 1.5):
            with self.subTest(value=value):
                with self.assertRaisesRegex(PreparationError, "positive integer"):
                    prepare_run(loaded, entropy=lambda bits: 1, hash_chunk_size=value)

    def test_prepared_model_can_be_revalidated_into_an_exact_snapshot(self) -> None:
        model_bytes = b'{"schema":"model"}\n'
        model_path = self.base_directory / "quality.json"
        model_path.write_bytes(model_bytes)
        document = base_config()
        document["sequencing"]["quality"] = {
            "kind": "markov",
            "artifact": {
                "path": "quality.json",
                "sha256": hashlib.sha256(model_bytes).hexdigest(),
            },
        }
        prepared = prepare_run(self.normalized(document), entropy=lambda bits: 7)
        identity = prepared.file_for_role("model.quality")

        self.assertEqual(
            snapshot_prepared_file(identity, maximum_size=1024, chunk_size=3),
            model_bytes,
        )

        model_path.write_bytes(b'{"schema":"changed!"}\n')
        with self.assertRaisesRegex(PreparationError, "changed before use"):
            snapshot_prepared_file(identity, maximum_size=1024)

    def test_snapshot_rejects_forged_metadata_and_size_limits(self) -> None:
        identity = FileDigest(
            role="model.quality",
            path=self.base_directory / "reference.fa",
            size_bytes=len(self.reference_bytes),
            sha256=hashlib.sha256(self.reference_bytes).hexdigest(),
            declared_sha256="0" * 64,
        )
        with self.assertRaisesRegex(PreparationError, "declaration"):
            snapshot_prepared_file(identity, maximum_size=1024)

        identity = FileDigest(
            role="model.quality",
            path=self.base_directory / "reference.fa",
            size_bytes=len(self.reference_bytes),
            sha256=hashlib.sha256(self.reference_bytes).hexdigest(),
            declared_sha256=hashlib.sha256(self.reference_bytes).hexdigest(),
        )
        with self.assertRaisesRegex(PreparationError, "size limit"):
            snapshot_prepared_file(identity, maximum_size=1)
        for maximum_size in (0, True):
            with self.subTest(maximum_size=maximum_size):
                with self.assertRaisesRegex(PreparationError, "positive integer"):
                    snapshot_prepared_file(identity, maximum_size=maximum_size)

        invalid_path = FileDigest(
            role="model.quality",
            path=Path("relative.json"),
            size_bytes=0,
            sha256=hashlib.sha256(b"").hexdigest(),
            declared_sha256=hashlib.sha256(b"").hexdigest(),
        )
        with self.assertRaisesRegex(PreparationError, "path must be absolute"):
            snapshot_prepared_file(invalid_path, maximum_size=1024)


if __name__ == "__main__":
    unittest.main()
