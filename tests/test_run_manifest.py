"""Tests for the final run manifest and its publication boundary."""

import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim import __version__  # noqa: E402
from bsreadsim.run.config import normalize_run_config  # noqa: E402
from bsreadsim.run.manifest import (  # noqa: E402
    ManifestError,
    build_complete_manifest,
    validate_header_projection,
    verify_complete_manifest,
)
from bsreadsim.output import OutputConfig, OutputSession  # noqa: E402
from process_support import UniformProcessConfig, process_fragment  # noqa: E402
from bsreadsim.run.prepare import prepare_run  # noqa: E402
from bsreadsim.native.protocol import (  # noqa: E402
    AmbiguityPolicy,
    BaseEncoding,
    Contig,
    Header,
    Technology,
    Trailer,
)
from bsreadsim.rng import RNG_CONTRACT  # noqa: E402
from bsreadsim.process.batch import READ_NAME_CONTRACT  # noqa: E402
from tests.test_process_stages import make_fragment  # noqa: E402
from tests.test_run_prepare import base_config  # noqa: E402


RUN_ID = "12345678-1234-4234-8234-123456789abc"


class ManifestTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.directory = Path(self.temporary_directory.name).resolve()
        self.reference = self.directory / "reference.fa"
        self.reference_sequence = b"ACGT" * 50
        self.reference.write_bytes(b">chr1\n" + self.reference_sequence + b"\n")

        document = base_config()
        document["reference"] = "reference.fa"
        document["seed"] = "7"
        document["output"] = {
            "directory": "output",
            "prefix": "sample",
            "compression": "none",
        }
        self.prepared = prepare_run(
            normalize_run_config(document, self.directory)
        )
        self.header = Header(
            run_id=RUN_ID,
            core_version=__version__,
            config_schema_version="1.1",
            rng_contract=RNG_CONTRACT,
            master_seed=7,
            normalized_config_sha256=bytes.fromhex(self.prepared.config.sha256),
            technology=Technology.WGBS,
            has_details=False,
            mates_per_fragment=1,
            base_encoding=BaseEncoding.ACGTN_U8,
            ambiguity_policy=AmbiguityPolicy.PRESERVE_N,
            read_length_r1=4,
            read_length_r2=0,
            contigs=(
                Contig(
                    name="chr1",
                    length=len(self.reference_sequence),
                    reference_sha256=hashlib.sha256(
                        self.reference_sequence
                    ).digest(),
                ),
            ),
        )
        self.trailer = Trailer(
            fragment_count=1,
            fragment_batch_count=1,
            mate_count=1,
            template_base_count=8,
            methylation_site_count=2,
            skipped_fragment_count=0,
            per_contig_fragment_counts=(1,),
            stream_sha256=hashlib.sha256(b"stream").digest(),
        )
        self.output_config = OutputConfig(
            directory=self.directory / "output",
            prefix="sample",
            paired_end=False,
            compression="none",
        )

    def staged_outputs(self):
        transaction = OutputSession(self.output_config)
        transaction.write_fragment(
            process_fragment(
                make_fragment(paired_end=False),
                "chr1",
                UniformProcessConfig(
                    master_seed=7,
                    conversion_rate=1.0,
                    error_rate=0.0,
                ),
            )
        )
        return transaction, transaction.finalize()

    def test_run_manifest_cross_checks_and_commits_the_run(self) -> None:
        transaction, outputs = self.staged_outputs()
        manifest = build_complete_manifest(
            self.prepared, self.header, self.trailer, outputs
        )

        verify_complete_manifest(manifest.document)
        self.assertEqual(manifest.document["status"], "complete")
        self.assertEqual(manifest.document["randomness"]["master_seed"], "7")
        self.assertEqual(
            manifest.document["versions"]["read_name"], READ_NAME_CONTRACT
        )
        self.assertEqual(manifest.document["counts"]["core"]["fragment_count"], 1)
        self.assertFalse((self.directory / "output/sample.R1.fastq").exists())

        transaction.commit(manifest.canonical_json)

        manifest_path = self.directory / "output/sample.manifest.json"
        observed = json.loads(manifest_path.read_text(encoding="utf-8"))
        verify_complete_manifest(observed)
        self.assertTrue((self.directory / "output/sample.R1.fastq").is_file())

    def test_fastq_only_manifest_accepts_exactly_the_read_roles(self) -> None:
        document = base_config()
        document["reference"] = "reference.fa"
        document["seed"] = "7"
        document["output"] = {
            "directory": "fastq-output",
            "prefix": "sample",
            "compression": "none",
        }
        prepared = prepare_run(
            normalize_run_config(document, self.directory)
        )
        header = copy.deepcopy(self.header)
        object.__setattr__(
            header,
            "normalized_config_sha256",
            bytes.fromhex(prepared.config.sha256),
        )
        object.__setattr__(header, "has_details", False)
        output_config = OutputConfig(
            directory=self.directory / "fastq-output",
            prefix="sample",
            paired_end=False,
            compression="none",
        )
        transaction = OutputSession(output_config)
        transaction.write_fragment(
            process_fragment(
                make_fragment(paired_end=False),
                "chr1",
                UniformProcessConfig(
                    master_seed=7,
                    conversion_rate=1.0,
                    error_rate=0.0,
                ),
                include_details=False,
            )
        )
        outputs = transaction.finalize()
        manifest = build_complete_manifest(
            prepared,
            header,
            self.trailer,
            outputs,
        )

        self.assertEqual(
            tuple(item["role"] for item in manifest.document["outputs"]),
            ("read1",),
        )
        transaction.commit(manifest.canonical_json)

    def test_run_manifest_digest_detects_any_stable_field_change(self) -> None:
        transaction, outputs = self.staged_outputs()
        self.addCleanup(transaction.abort)
        manifest = build_complete_manifest(
            self.prepared, self.header, self.trailer, outputs
        )
        tampered = manifest.as_dict()
        tampered["counts"]["core"]["skipped_fragment_count"] = 99

        with self.assertRaisesRegex(ManifestError, "mismatch"):
            verify_complete_manifest(tampered)

    def test_header_projection_mismatches_fail_before_publication(self) -> None:
        transaction, outputs = self.staged_outputs()
        self.addCleanup(transaction.abort)
        invalid_header = copy.deepcopy(self.header)
        object.__setattr__(invalid_header, "master_seed", 8)

        with self.assertRaisesRegex(ManifestError, "master seed"):
            build_complete_manifest(
                self.prepared, invalid_header, self.trailer, outputs
            )
        self.assertFalse((self.directory / "output/sample.manifest.json").exists())

    def test_header_projection_checks_execution_shape_and_annotation_mode(self) -> None:
        header = Header(
            run_id=RUN_ID,
            core_version=__version__,
            config_schema_version="1.1",
            rng_contract=RNG_CONTRACT,
            master_seed=7,
            normalized_config_sha256=bytes.fromhex(self.prepared.config.sha256),
            technology=Technology.WGBS,
            has_details=False,
            mates_per_fragment=1,
            base_encoding=BaseEncoding.ACGTN_U8,
            ambiguity_policy=AmbiguityPolicy.PRESERVE_N,
            read_length_r1=4,
            read_length_r2=0,
            contigs=(
                Contig(
                    name="chr1",
                    length=len(self.reference_sequence),
                    reference_sha256=hashlib.sha256(
                        self.reference_sequence
                    ).digest(),
                ),
            ),
        )

        validate_header_projection(self.prepared, header)

        failures = (
            (replace(header, technology=Technology.RRBS), "technology"),
            (replace(header, mates_per_fragment=2, read_length_r2=4), "mate"),
            (replace(header, read_length_r1=5), "R1"),
            (replace(header, has_details=True), "detail policy"),
            (
                replace(header, ambiguity_policy=AmbiguityPolicy.RESOLVE_ONCE),
                "ambiguity",
            ),
        )
        for invalid, pattern in failures:
            with self.subTest(pattern=pattern):
                with self.assertRaisesRegex(ManifestError, pattern):
                    validate_header_projection(self.prepared, invalid)

    def test_bam_requires_full_projection(self) -> None:
        document = base_config()
        document["reference"] = "reference.fa"
        document["seed"] = "7"
        document["output"] = {
            "directory": "bam-output",
            "prefix": "sample",
            "compression": "none",
            "bam": True,
        }
        prepared = prepare_run(normalize_run_config(document, self.directory))
        header = replace(
            self.header,
            normalized_config_sha256=bytes.fromhex(prepared.config.sha256),
            has_details=True,
        )

        validate_header_projection(prepared, header)
        with self.assertRaisesRegex(ManifestError, "detail policy"):
            validate_header_projection(
                prepared,
                replace(header, has_details=False),
            )

    def test_count_and_output_role_mismatches_fail_closed(self) -> None:
        transaction, outputs = self.staged_outputs()
        self.addCleanup(transaction.abort)
        invalid_trailer = copy.deepcopy(self.trailer)
        object.__setattr__(invalid_trailer, "mate_count", 2)
        with self.assertRaisesRegex(ManifestError, "mate count"):
            build_complete_manifest(
                self.prepared, self.header, invalid_trailer, outputs
            )

        invalid_outputs = copy.deepcopy(outputs)
        object.__setattr__(invalid_outputs, "fragment_count", 2)
        with self.assertRaisesRegex(ManifestError, "fragment counts"):
            build_complete_manifest(
                self.prepared, self.header, self.trailer, invalid_outputs
            )


if __name__ == "__main__":
    unittest.main()
