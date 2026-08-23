"""Tests for transactional, deterministic simulation output."""

import gzip
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest


from bsreadsim.output import (
    OutputConfig,
    OutputError,
    OutputSession,
    BamConfig,
)
from bsreadsim.output.fastq import (
    format_fragment_records_trusted,
)
from tests.helpers.process_support import UniformProcessConfig, process_fragment
from bsreadsim.run.manifest import MANIFEST_SCHEMA_VERSION
from tests.unit.test_process_stages import make_fragment


def manifest_json(summary):
    document = {
        "counts": {
            "python": {
                "fragment_count": summary.fragment_count,
                "mate_count": summary.mate_count,
                "records_by_role": {
                    item.role: item.record_count for item in summary.files
                },
            }
        },
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
        "outputs": [
            {
                "path": str(item.path),
                "record_count": item.record_count,
                "role": item.role,
                "sha256": item.sha256,
                "size_bytes": item.size_bytes,
            }
            for item in summary.files
        ],
        "reproducibility": {"scope": "output component test"},
        "status": "complete",
    }
    identity = json.dumps(
        document, allow_nan=False, separators=(",", ":"), sort_keys=True
    )
    document["reproducibility"]["sha256"] = hashlib.sha256(
        identity.encode("utf-8")
    ).hexdigest()
    return json.dumps(
        document, allow_nan=False, separators=(",", ":"), sort_keys=True
    )


def processed(
    ordinal=0,
    paired_end=False,
    compact_base_states=False,
    include_details=True,
):
    return process_fragment(
        make_fragment(paired_end=paired_end, ordinal=ordinal),
        "chr1",
        UniformProcessConfig(
            master_seed=7,
            conversion_rate=1.0,
            error_rate=0.0,
            quality_phred=30,
        ),
        compact_base_states=compact_base_states,
        include_details=include_details,
    )


class OutputSessionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.directory = Path(self.temporary_directory.name).resolve()

    def config(self, **overrides):
        values = {
            "directory": self.directory,
            "prefix": "sample",
            "paired_end": False,
            "compression": "none",
        }
        values.update(overrides)
        return OutputConfig(**values)

    def bam_config(self, program: str) -> BamConfig:
        return BamConfig(
            writer_argv=(sys.executable, "-c", program),
            sam_header=(
                b"@HD\tVN:1.6\tSO:unsorted\n"
                b"@SQ\tSN:chr1\tLN:1000\n"
                b"@RG\tID:run\tSM:sample\n"
                b"@PG\tID:bsreadsim\tPN:bsreadsim\tVN:test\n"
            ),
            references=(("chr1", 1000),),
            read_group_id="run",
        )

    def test_single_end_files_are_invisible_until_commit(self) -> None:
        transaction = OutputSession(self.config())
        transaction.write_fragment(processed())
        self.assertFalse((self.directory / "sample.R1.fastq").exists())
        self.assertTrue(transaction.staging_directory.exists())

        summary = transaction.finalize()
        self.assertFalse((self.directory / "sample.R1.fastq").exists())
        committed_manifest = manifest_json(summary)
        transaction.commit(committed_manifest)

        fastq = self.directory / "sample.R1.fastq"
        self.assertEqual(
            fastq.read_text(encoding="utf-8"),
            "@chr1:101-108:0/1\nACTGT\n+\n?????\n",
        )
        self.assertEqual(summary.fragment_count, 1)
        self.assertEqual(summary.mate_count, 1)
        self.assertEqual(
            (self.directory / "sample.manifest.json").read_text(encoding="utf-8"),
            committed_manifest + "\n",
        )
        for item in summary.files:
            self.assertEqual(item.size_bytes, item.path.stat().st_size)
            self.assertEqual(
                item.sha256, hashlib.sha256(item.path.read_bytes()).hexdigest()
            )
            self.assertEqual(item.record_count, 1)

    def test_fastq_output_accepts_no_annotations(self) -> None:
        config = self.config(paired_end=True)
        with OutputSession(config) as transaction:
            transaction.write_fragment(
                processed(
                    paired_end=True,
                    )
            )
            summary = transaction.finalize()
            transaction.commit(manifest_json(summary))

        self.assertEqual(
            tuple(item.role for item in summary.files),
            ("read1", "read2"),
        )
        self.assertTrue((self.directory / "sample.R1.fastq").is_file())
        self.assertTrue((self.directory / "sample.R2.fastq").is_file())

    def test_fastq_only_formatted_batch_is_supported(self) -> None:
        config = self.config(paired_end=True)
        fragments = tuple(
            processed(
                ordinal=ordinal,
                paired_end=True,
            )
            for ordinal in range(3)
        )
        formatted = tuple(
            format_fragment_records_trusted(
                fragment,
                paired_end=True,
            )
            for fragment in fragments
        )
        with OutputSession(config) as transaction:
            transaction.write_formatted_batch(
                0,
                tuple((len(read1), len(read2)) for read1, read2 in formatted),
                b"".join(item[0] for item in formatted),
                b"".join(item[1] for item in formatted),
            )
            summary = transaction.finalize()
            transaction.commit(manifest_json(summary))

        self.assertEqual({item.role for item in summary.files}, {"read1", "read2"})
        self.assertEqual(summary.mate_count, 6)
        self.assertTrue(all(item.record_count == 3 for item in summary.files))

    def test_paired_mates_are_published_together(self) -> None:
        with OutputSession(self.config(paired_end=True)) as transaction:
            transaction.write_fragment(processed(paired_end=True))
            summary = transaction.finalize()
            transaction.commit(manifest_json(summary))

        self.assertEqual(summary.mate_count, 2)
        self.assertTrue((self.directory / "sample.R1.fastq").is_file())
        read2 = (self.directory / "sample.R2.fastq").read_text(encoding="utf-8")
        self.assertTrue(read2.startswith("@chr1:101-108:0/2\nCGACA\n"))

    def test_bam_is_streamed_and_counted_per_mate(self) -> None:
        program = (
            "import sys; value=sys.stdin.buffer.read(); "
            "sys.stdout.buffer.write(b'fake-bam\\n' + value)"
        )
        config = self.config(
            paired_end=True,
            bam=self.bam_config(program),
        )
        with OutputSession(config) as transaction:
            transaction.write_fragment(processed(paired_end=True))
            summary = transaction.finalize()
            transaction.commit(manifest_json(summary))

        bam = self.directory / "sample.bam"
        self.assertTrue(bam.read_bytes().startswith(b"fake-bam\n@HD\tVN:1.6"))
        self.assertEqual(summary.file_for_role("bam").record_count, 2)
        self.assertEqual(
            tuple(item.role for item in summary.files),
            ("bam",),
        )

    def test_bam_writer_failure_aborts_every_output(self) -> None:
        program = "import sys; sys.stdin.buffer.read(); raise SystemExit(9)"
        transaction = OutputSession(
            self.config(
                    bam=self.bam_config(program),
            )
        )
        transaction.write_fragment(processed())
        with self.assertRaisesRegex(OutputError, "status 9"):
            transaction.finalize()

        self.assertFalse(transaction.staging_directory.exists())
        self.assertFalse((self.directory / "sample.R1.fastq").exists())
        self.assertFalse((self.directory / "sample.bam").exists())

    def test_fragment_order_and_mate_mode_fail_closed(self) -> None:
        transaction = OutputSession(self.config())
        with self.assertRaisesRegex(OutputError, "start at zero"):
            transaction.write_fragment(processed(ordinal=1))
        self.assertFalse(transaction.staging_directory.exists())
        self.assertFalse((self.directory / "sample.R1.fastq").exists())

        transaction = OutputSession(self.config())
        with self.assertRaisesRegex(OutputError, "cardinality"):
            transaction.write_fragment(processed(paired_end=True))
        self.assertFalse((self.directory / "sample.R1.fastq").exists())

    def test_context_exception_aborts_without_final_outputs(self) -> None:
        staging = None
        with self.assertRaisesRegex(RuntimeError, "downstream failure"):
            with OutputSession(self.config()) as transaction:
                staging = transaction.staging_directory
                transaction.write_fragment(processed())
                raise RuntimeError("downstream failure")

        self.assertIsNotNone(staging)
        self.assertFalse(staging.exists())
        self.assertFalse((self.directory / "sample.R1.fastq").exists())

    def test_existing_destination_is_never_overwritten(self) -> None:
        existing = self.directory / "sample.R1.fastq"
        existing.write_text("keep\n", encoding="utf-8")
        with self.assertRaisesRegex(OutputError, "already exists"):
            OutputSession(self.config())
        self.assertEqual(existing.read_text(encoding="utf-8"), "keep\n")

    def test_gzip_bytes_are_deterministic_and_have_zero_mtime(self) -> None:
        outputs = []
        for index in range(2):
            directory = self.directory / str(index)
            directory.mkdir()
            config = OutputConfig(
                directory=directory,
                prefix="same",
                paired_end=False,
                compression="gzip",
            )
            with OutputSession(config) as transaction:
                transaction.write_fragment(processed())
                summary = transaction.finalize()
                transaction.commit(manifest_json(summary))
            outputs.append((directory / "same.R1.fastq.gz").read_bytes())

        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(outputs[0][4:8], b"\x00\x00\x00\x00")
        with gzip.open(self.directory / "0/same.R1.fastq.gz", "rt") as stream:
            self.assertEqual(
                stream.read(),
                "@chr1:101-108:0/1\nACTGT\n+\n?????\n",
            )

    def test_gzip_level_changes_only_the_compressed_representation(self) -> None:
        outputs = []
        for level in (1, 6):
            directory = self.directory / str(level)
            directory.mkdir()
            config = OutputConfig(
                directory=directory,
                prefix="same",
                paired_end=False,
                compression="gzip",
                gzip_level=level,
                )
            with OutputSession(config) as transaction:
                for ordinal in range(100):
                    transaction.write_fragment(
                        processed(ordinal=ordinal, include_details=False)
                    )
                summary = transaction.finalize()
                transaction.commit(manifest_json(summary))
            compressed = (directory / "same.R1.fastq.gz").read_bytes()
            outputs.append((compressed, gzip.decompress(compressed)))

        self.assertEqual(outputs[0][1], outputs[1][1])
        self.assertNotEqual(outputs[0][0], outputs[1][0])

    def test_formatted_batch_is_byte_identical_to_fragment_writes(self) -> None:
        fragments = tuple(
            processed(ordinal=ordinal, paired_end=True) for ordinal in range(3)
        )
        directories = (self.directory / "fragment", self.directory / "batch")
        for directory in directories:
            directory.mkdir()

        with OutputSession(
            OutputConfig(
                directory=directories[0],
                prefix="same",
                paired_end=True,
                compression="gzip",
            )
        ) as transaction:
            for fragment in fragments:
                transaction.write_fragment(fragment)
            fragment_summary = transaction.finalize()
            transaction.commit(manifest_json(fragment_summary))

        formatted = tuple(
            format_fragment_records_trusted(fragment, paired_end=True)
            for fragment in fragments
        )
        lengths = tuple(
            (len(read1), len(read2))
            for read1, read2 in formatted
        )
        with OutputSession(
            OutputConfig(
                directory=directories[1],
                prefix="same",
                paired_end=True,
                compression="gzip",
            )
        ) as transaction:
            transaction.write_formatted_batch(
                0,
                lengths,
                memoryview(b"".join(item[0] for item in formatted)),
                memoryview(b"".join(item[1] for item in formatted)),
            )
            batch_summary = transaction.finalize()
            transaction.commit(manifest_json(batch_summary))

        for name in (
            "same.R1.fastq.gz",
            "same.R2.fastq.gz",
        ):
            self.assertEqual(
                (directories[0] / name).read_bytes(),
                (directories[1] / name).read_bytes(),
            )
        self.assertEqual(fragment_summary.fragment_count, 3)
        self.assertEqual(batch_summary.mate_count, 6)

    def test_invalid_formatted_batch_aborts_the_transaction(self) -> None:
        transaction = OutputSession(self.config())
        staging = transaction.staging_directory

        with self.assertRaisesRegex(OutputError, "region sizes"):
            transaction.write_formatted_batch(
                0,
                ((4, 0),),
                b"abc",
                None,
            )

        self.assertFalse(staging.exists())
        self.assertFalse((self.directory / "sample.R1.fastq").exists())

    def test_invalid_policy_is_rejected_before_staging(self) -> None:
        with self.assertRaises(OutputError):
            self.config(prefix="../escape")
        with self.assertRaisesRegex(OutputError, "absolute"):
            OutputConfig(
                directory=Path("relative"),
                prefix="sample",
                paired_end=False,
            )
        for value in (-1, 10, True, 1.5):
            with self.subTest(gzip_level=value):
                with self.assertRaisesRegex(OutputError, "gzip_level"):
                    self.config(gzip_level=value)

    def test_run_manifest_is_required_and_must_be_canonical_complete_json(self) -> None:
        for value in ("{}", '{"status": "complete"}', "not-json"):
            with self.subTest(value=value):
                config = self.config(prefix="case{}".format(len(value)))
                with OutputSession(config) as transaction:
                    transaction.write_fragment(processed())
                    transaction.finalize()
                    with self.assertRaises(OutputError):
                        transaction.commit(value)
                self.assertFalse(
                    (self.directory / "{}.manifest.json".format(config.prefix)).exists()
                )

    def test_run_manifest_must_identify_the_exact_staged_bytes(self) -> None:
        transaction = OutputSession(self.config())
        transaction.write_fragment(processed())
        summary = transaction.finalize()
        document = json.loads(manifest_json(summary))
        document["outputs"][0]["sha256"] = "0" * 64
        document["reproducibility"].pop("sha256")
        identity = json.dumps(document, separators=(",", ":"), sort_keys=True)
        document["reproducibility"]["sha256"] = hashlib.sha256(
            identity.encode("utf-8")
        ).hexdigest()
        tampered = json.dumps(document, separators=(",", ":"), sort_keys=True)

        with self.assertRaisesRegex(OutputError, "staged outputs"):
            transaction.commit(tampered)
        self.assertFalse((self.directory / "sample.R1.fastq").exists())
        self.assertFalse((self.directory / "sample.manifest.json").exists())


if __name__ == "__main__":
    unittest.main()
