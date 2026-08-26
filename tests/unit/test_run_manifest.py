"""Tests for the final run manifest and its publication boundary."""

import copy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import shlex
import tempfile
import unittest


from bsreadsim import __version__
from bsreadsim.cli import build_parser, build_run_document
from bsreadsim.run.config import normalize_run_config
from bsreadsim.run.manifest import (
    MANIFEST_VERSION,
    ManifestError,
    build_complete_manifest,
    validate_header_projection,
    verify_complete_manifest,
)
from bsreadsim.output import OutputConfig, OutputSession
from tests.helpers.process_support import UniformProcessConfig, process_fragment
from bsreadsim.run.prepare import prepare_run
from bsreadsim.htsim.protocol import (
    AmbiguityPolicy,
    BaseEncoding,
    Contig,
    Header,
    Technology,
    Trailer,
)
from bsreadsim.rng import RNG_CONTRACT
from bsreadsim.process.batch import READ_NAME_CONTRACT
from tests.unit.test_process_stages import make_fragment
from tests.unit.test_run_prepare import base_config


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
        document["fragments"]["insert_sd"] = 25.0
        document["output"] = {
            "directory": "output",
            "prefix": "sample",
            "format": "fastq",
        }
        self.prepared = prepare_run(
            normalize_run_config(document, self.directory)
        )
        self.header = Header(
            run_id=RUN_ID,
            core_version=__version__,
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
            format="fastq",
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
        self.assertEqual(manifest.document["version"], MANIFEST_VERSION)
        details = manifest.document["details"]
        self.assertEqual(details["randomness"]["master_seed"], "7")
        self.assertEqual(
            details["randomness"]["methylation_seed"], "0"
        )
        self.assertEqual(details["randomness"]["mutation_seed"], "0")
        self.assertEqual(details["randomness"]["phasing_seed"], "0")
        effective = details["configuration"]
        self.assertNotIn("seed", effective)
        self.assertNotIn("seeds", effective)
        self.assertEqual(effective["coverage"]["type"], "uniform")
        self.assertNotIn("kind", effective["coverage"])
        self.assertEqual(effective["methylation"]["beta"]["CG"], "0.5,0.5")
        self.assertEqual(
            details["contracts"]["read_name"], READ_NAME_CONTRACT
        )
        self.assertEqual(manifest.document["summary"]["fragment_count"], 1)
        self.assertEqual(manifest.document["summary"]["read_count"], 1)
        self.assertEqual(
            manifest.document["command"], {"interface": "python-api"}
        )
        self.assertFalse((self.directory / "output/sample.R1.fastq").exists())

        transaction.commit(manifest.canonical_json)

        manifest_path = self.directory / "output/sample.manifest.json"
        observed = json.loads(manifest_path.read_text(encoding="utf-8"))
        verify_complete_manifest(observed)
        manifest_text = manifest_path.read_text(encoding="utf-8")
        self.assertTrue(manifest_text.startswith("{\n  \"command\": {\n"))
        self.assertTrue(manifest_text.endswith("\n"))
        self.assertTrue((self.directory / "output/sample.R1.fastq").is_file())

    def test_run_manifest_records_user_and_full_cli_commands(self) -> None:
        transaction, outputs = self.staged_outputs()
        self.addCleanup(transaction.abort)
        spaced_reference = self.directory / "reference with spaces.fa"
        spaced_reference.write_bytes(self.reference.read_bytes())
        argv = (
            "bsreadsim",
            "run",
            "wgbs",
            "-r",
            "reference with spaces.fa",
            "-o",
            "output",
            "-n",
            "1",
            "--single-end",
            "--read-length",
            "4",
            "--prefix",
            "sample",
            "--format",
            "fastq",
        )
        arguments = build_parser().parse_args(argv[1:])
        prepared = prepare_run(
            normalize_run_config(
                build_run_document(arguments, self.directory), self.directory
            ),
            entropy=lambda bits: 7,
        )
        header = replace(
            self.header,
            normalized_config_sha256=bytes.fromhex(prepared.config.sha256),
        )
        manifest = build_complete_manifest(
            prepared,
            header,
            self.trailer,
            outputs,
            invocation_argv=argv,
        )

        command = manifest.document["command"]
        self.assertEqual(command["interface"], "cli")
        self.assertEqual(shlex.split(command["user_command"]), list(argv))
        self.assertEqual(
            command["user_command"],
            "bsreadsim run wgbs -r 'reference with spaces.fa' -o output -n 1 "
            "--single-end --read-length 4 --prefix sample --format fastq",
        )
        full_argv = shlex.split(command["full_command"])
        self.assertEqual(command["full_command"], shlex.join(full_argv))
        self.assertEqual(full_argv[:3], ["bsreadsim", "run", "wgbs"])
        self.assertEqual(full_argv[full_argv.index("--seed") + 1], "7")
        for option in (
            "--reference",
            "--output",
            "--fragments",
            "--seed-mut",
            "--seed-phase",
            "--seed-meth",
            "--read-length",
            "--insert-min",
            "--insert-mean",
            "--insert-max",
            "--insert-sd",
            "--max-ambiguous-fraction",
            "--mutation-rate",
            "--indel-fraction",
            "--indel-extension-probability",
            "--beta-cg",
            "--beta-chg",
            "--beta-chh",
            "--methylation-model",
            "--conversion-rate",
            "--phred",
            "--error-rate",
            "--workers",
            "--core-workers",
            "--chunk-size",
            "--max-in-flight-fragments",
            "--prefix",
            "--format",
            "--gzip-level",
        ):
            self.assertIn(option, full_argv)

        full_arguments = build_parser().parse_args(full_argv[1:])
        round_trip = normalize_run_config(
            build_run_document(full_arguments, self.directory), self.directory
        )
        self.assertEqual(round_trip.canonical_json, prepared.config.canonical_json)
        verify_complete_manifest(manifest.document)

    def test_full_cli_command_canonicalizes_the_effective_seed(self) -> None:
        transaction, outputs = self.staged_outputs()
        self.addCleanup(transaction.abort)
        argv = (
            "bsreadsim",
            "run",
            "wgbs",
            "-r",
            "reference.fa",
            "-n",
            "1",
            "-s7",
        )
        manifest = build_complete_manifest(
            self.prepared,
            self.header,
            self.trailer,
            outputs,
            invocation_argv=argv,
        )

        command = manifest.document["command"]
        self.assertEqual(shlex.split(command["user_command"]), list(argv))
        full_argv = shlex.split(command["full_command"])
        self.assertNotIn("-s7", full_argv)
        self.assertEqual(full_argv.count("--seed"), 1)
        self.assertEqual(full_argv[full_argv.index("--seed") + 1], "7")
        verify_complete_manifest(manifest.document)

        tampered = manifest.as_dict()
        tampered_full_argv = shlex.split(tampered["command"]["full_command"])
        seed_index = tampered_full_argv.index("--seed") + 1
        tampered_full_argv[seed_index] = "8"
        tampered["command"]["full_command"] = shlex.join(tampered_full_argv)
        with self.assertRaisesRegex(ManifestError, "command"):
            verify_complete_manifest(tampered)

    def test_standard_manifest_disables_methylation_model_and_sites(self) -> None:
        document = base_config()
        document["reference"] = "reference.fa"
        document["seed"] = "7"
        document["technology"] = "WGS"
        document["sequencing"]["conversion_rate"] = 0.0
        document["output"] = {
            "directory": "output",
            "prefix": "sample",
            "format": "fastq",
        }
        prepared = prepare_run(normalize_run_config(document, self.directory))
        header = replace(
            self.header,
            normalized_config_sha256=bytes.fromhex(prepared.config.sha256),
            technology=Technology.WGS,
        )
        trailer = replace(self.trailer, methylation_site_count=0)
        transaction, outputs = self.staged_outputs()
        self.addCleanup(transaction.abort)

        manifest = build_complete_manifest(prepared, header, trailer, outputs)

        self.assertEqual(
            manifest.document["details"]["models"]["methylation_state"],
            {"effective": "disabled", "requested": "disabled"},
        )
        with self.assertRaisesRegex(ManifestError, "methylation sites"):
            build_complete_manifest(
                prepared,
                header,
                replace(trailer, methylation_site_count=1),
                outputs,
            )

    def test_fastq_only_manifest_accepts_exactly_the_read_roles(self) -> None:
        document = base_config()
        document["reference"] = "reference.fa"
        document["seed"] = "7"
        document["output"] = {
            "directory": "fastq-output",
            "prefix": "sample",
            "format": "fastq",
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
            format="fastq",
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
        tampered["summary"]["skipped_fragment_count"] = 99

        with self.assertRaisesRegex(ManifestError, "mismatch"):
            verify_complete_manifest(tampered)

    def test_manifest_version_is_strict(self) -> None:
        transaction, outputs = self.staged_outputs()
        self.addCleanup(transaction.abort)
        manifest = build_complete_manifest(
            self.prepared, self.header, self.trailer, outputs
        ).as_dict()

        manifest["version"] = True
        with self.assertRaisesRegex(ManifestError, "version"):
            verify_complete_manifest(manifest)

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
                replace(header, ambiguity_policy=1),
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
            "format": "bam",
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
