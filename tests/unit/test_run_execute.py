"""Tests for bounded C++/Python orchestration boundaries."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import tempfile
import unittest


from bsreadsim.run.config import (
    ConfigValidationError,
    normalize_run_config,
)
from bsreadsim.run.execute import (
    PipelineError,
    _build_execution_plan,
    _protocol_batch_fragment_limit,
    run_prepared,
)
from bsreadsim.run.prepare import prepare_run
from tests.unit.test_run_prepare import base_config


RUN_ID = "12345678-1234-4234-8234-123456789abc"


def baseline_config() -> dict:
    document = base_config()
    document["seed"] = "7"
    document["mutation"]["rate"] = 0
    document["fragments"].update(
        {
            "insert_min": 5,
            "insert_mean": 5,
            "insert_max": 5,
            "insert_sd": 0,
        }
    )
    return document


def quality_model_bytes(scores=(10, 30)) -> bytes:
    width = len(scores)
    initial = []
    for cycle in range(5):
        row = [0] * width
        row[cycle % width] = 1
        initial.append(row)
    transitions = []
    for state in range(width):
        row = [0] * width
        row[(state + 1) % width] = 1
        transitions.append(row)
    mate = {
        "initial_counts": initial,
        "transition_counts": transitions,
    }
    return json.dumps(
        {
            "schema": "quality-markov",
            "quality_scores": list(scores),
            "mates": [mate, mate],
        },
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def error_model_bytes(scores=(10, 30, 40)) -> bytes:
    identity = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
    mate = {"base_transition_counts": [identity for _ in scores]}
    return json.dumps(
        {
            "schema": "quality-confusion",
            "quality_scores": list(scores),
            "mates": [mate, mate],
        },
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


class PipelineTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.directory = Path(self.temporary_directory.name).resolve()
        (self.directory / "reference.fa").write_bytes(b">chr1\nACGTCGTAA\n")

    def prepared(self, document=None):
        return prepare_run(
            normalize_run_config(
                baseline_config() if document is None else document,
                self.directory,
            )
        )


    def test_full_annotation_bounds_fragment_local_object_lifetimes(self) -> None:
        self.assertEqual(
            _protocol_batch_fragment_limit(
                emit_details=True,
                max_in_flight=2048,
            ),
            128,
        )
        self.assertEqual(
            _protocol_batch_fragment_limit(
                emit_details=False,
                max_in_flight=2048,
            ),
            1024,
        )
        self.assertEqual(
            _protocol_batch_fragment_limit(
                emit_details=True,
                max_in_flight=2048,
                materializes_detail_objects=False,
            ),
            1024,
        )
        self.assertEqual(
            _protocol_batch_fragment_limit(
                emit_details=True,
                max_in_flight=3,
            ),
            3,
        )
        self.assertEqual(
            _protocol_batch_fragment_limit(
                emit_details=True,
                max_in_flight=3,
                materializes_detail_objects=False,
            ),
            3,
        )

    def test_global_thread_budget_is_split_without_stage_multiplication(self) -> None:
        fastq = _build_execution_plan(16, "fastq.gz")
        self.assertEqual(fastq.core_workers, 8)
        self.assertEqual(fastq.process_workers, 8)
        self.assertEqual(fastq.bam_compression_threads, 0)
        self.assertTrue(fastq.use_process_pool)
        self.assertEqual(fastq.protocol_batch_fragments, 1024)
        self.assertEqual(fastq.max_in_flight_fragments, 16384)

        bam = _build_execution_plan(16, "bam")
        self.assertEqual(
            bam.core_workers
            + bam.process_workers
            + bam.bam_compression_threads
            + 1,
            16,
        )
        self.assertEqual(bam.bam_compression_threads, 2)

        serial = _build_execution_plan(1, "fastq.gz")
        self.assertFalse(serial.use_process_pool)
        self.assertEqual(serial.core_workers, 1)
        self.assertEqual(serial.process_workers, 1)

    def test_advanced_sequencing_models_load_before_process_launch(self) -> None:
        quality_bytes = quality_model_bytes()
        error_bytes = error_model_bytes()
        quality_path = self.directory / "quality.json"
        error_path = self.directory / "error.json"
        quality_path.write_bytes(quality_bytes)
        error_path.write_bytes(error_bytes)
        document = baseline_config()
        document["sequencing"]["quality"] = {
            "kind": "markov",
            "artifact": {
                "path": "quality.json",
                "sha256": hashlib.sha256(quality_bytes).hexdigest(),
            },
        }
        document["sequencing"]["error"] = {
            "kind": "quality-confusion",
            "artifact": {
                "path": "error.json",
                "sha256": hashlib.sha256(error_bytes).hexdigest(),
            },
        }
        prepared = self.prepared(document)

        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                prepared,
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        quality_path.write_bytes(quality_bytes.replace(b"10", b"11", 1))
        with self.assertRaisesRegex(PipelineError, "changed before use"):
            run_prepared(
                prepared,
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )
        self.assertFalse((self.directory / "output").exists())

    def test_sequencing_model_contract_and_quality_domain_fail_closed(self) -> None:
        quality_bytes = quality_model_bytes()
        error_bytes = error_model_bytes(scores=(10,))
        (self.directory / "quality.json").write_bytes(quality_bytes)
        (self.directory / "error.json").write_bytes(error_bytes)
        document = baseline_config()
        document["sequencing"]["quality"] = {
            "kind": "markov",
            "artifact": {
                "path": "quality.json",
                "sha256": hashlib.sha256(quality_bytes).hexdigest(),
            },
        }
        document["sequencing"]["error"] = {
            "kind": "quality-confusion",
            "artifact": {
                "path": "error.json",
                "sha256": hashlib.sha256(error_bytes).hexdigest(),
            },
        }
        with self.assertRaisesRegex(PipelineError, "missing quality score.*30"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        uniform = baseline_config()
        uniform["sequencing"]["error"] = document["sequencing"]["error"]
        with self.assertRaisesRegex(PipelineError, "missing quality score.*40"):
            run_prepared(
                self.prepared(uniform),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )
        self.assertFalse((self.directory / "output").exists())

    def test_de_novo_mutation_passes_all_technology_gates(self) -> None:
        document = baseline_config()
        document["mutation"]["rate"] = 0.1
        prepared = self.prepared(document)

        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                prepared,
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        profile_bytes = b"0.5\n0.5\n"
        (self.directory / "coverage.tsv").write_bytes(profile_bytes)
        document["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "coverage.tsv",
                "sha256": hashlib.sha256(profile_bytes).hexdigest(),
            },
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document["technology"] = "RRBS"
        document["rrbs"] = {"cut_sites": ["C|CGG"]}
        document["coverage"] = {"kind": "uniform"}
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        (self.directory / "variants.vcf").write_text("fixture\n", encoding="utf-8")
        document = baseline_config()
        document["mutation"]["rate"] = 0.1
        document["inputs"] = {"vcf": "variants.vcf"}
        with self.assertRaisesRegex(PipelineError, "mutually exclusive"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document["inputs"] = {}
        document["methylation"]["update_variant_boundaries"] = False
        with self.assertRaisesRegex(PipelineError, "variant generation requires"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        self.assertFalse((self.directory / "output").exists())

    def test_depth_passes_all_technology_capability_gates(self) -> None:
        document = baseline_config()
        document["reads"] = {"depth": 2}

        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document["technology"] = "RRBS"
        document["rrbs"] = {"cut_sites": ["C|CGG"]}
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        (self.directory / "targets.bed").write_text(
            "chr1\t2\t8\tprobe\t1\t+\n", encoding="utf-8"
        )
        document["technology"] = "TBS"
        document.pop("rrbs")
        document["tbs"] = {
            "bed": "targets.bed",
            "center_sd": 0,
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        self.assertFalse((self.directory / "output").exists())

    def test_variable_insert_capability_boundary(self) -> None:
        document = baseline_config()
        document["fragments"].update(
            {
                "insert_mean": 9,
                "insert_max": 16,
                "insert_sd": 4,
            }
        )
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        profile_bytes = b"0\n0.5\n0.5\n"
        (self.directory / "coverage.tsv").write_bytes(profile_bytes)
        profiled = baseline_config()
        profiled["fragments"].update(document["fragments"])
        profiled["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "coverage.tsv",
                "sha256": hashlib.sha256(profile_bytes).hexdigest(),
            },
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(profiled),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        (self.directory / "variants.vcf").write_text("fixture\n", encoding="utf-8")
        with_vcf = baseline_config()
        with_vcf["fragments"].update(document["fragments"])
        with_vcf["inputs"] = {"vcf": "variants.vcf"}
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(with_vcf),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        (self.directory / "levels.asm").write_text("fixture\n", encoding="utf-8")
        with_asm = baseline_config()
        with_asm["fragments"].update(document["fragments"])
        with_asm["inputs"] = {
            "vcf": "variants.vcf",
            "asm": "levels.asm",
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(with_asm),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        with_mutation = baseline_config()
        with_mutation["fragments"].update(document["fragments"])
        with_mutation["mutation"]["rate"] = 0.1
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(with_mutation),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        (self.directory / "targets.bed").write_text(
            "chr1\t3\t4\ttarget\t0\t+\n", encoding="utf-8"
        )
        tbs = baseline_config()
        tbs["technology"] = "TBS"
        tbs["fragments"].update(document["fragments"])
        tbs["tbs"] = {
            "bed": "targets.bed",
            "center_sd": 0,
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(tbs),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        self.assertFalse((self.directory / "output").exists())

    def test_tbs_invalid_center_distribution_gate_precedes_side_effects(self) -> None:
        (self.directory / "targets.bed").write_text(
            "chr1\t3\t4\ttarget\t0\t+\n", encoding="utf-8"
        )
        document = baseline_config()
        document["technology"] = "TBS"
        document["tbs"] = {
            "bed": "targets.bed",
            "center_sd": -1,
        }
        with self.assertRaises(ConfigValidationError):
            self.prepared(document)

        self.assertFalse((self.directory / "output").exists())

    def test_tbs_nonzero_center_distribution_passes_capability_gate(self) -> None:
        (self.directory / "targets.bed").write_text(
            "chr1\t3\t4\ttarget\t0\t+\n", encoding="utf-8"
        )
        document = baseline_config()
        document["technology"] = "TBS"
        document["fragments"].update(
            {
                "insert_min": 2,
                "insert_mean": 5,
                "insert_max": 8,
                "insert_sd": 0,
            }
        )
        document["tbs"] = {
            "bed": "targets.bed",
            "center_sd": 1,
        }
        prepared = self.prepared(document)

        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                prepared,
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        self.assertFalse((self.directory / "output").exists())

    def test_wgbs_gc_profile_contract_precedes_process_side_effects(self) -> None:
        profile_bytes = b"0\n0.5\n0.5\n"
        (self.directory / "coverage.tsv").write_bytes(profile_bytes)
        document = baseline_config()
        document["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "coverage.tsv",
                "sha256": hashlib.sha256(profile_bytes).hexdigest(),
            },
        }
        prepared = self.prepared(document)
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                prepared,
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        wrong_technology = baseline_config()
        wrong_technology["technology"] = "RRBS"
        wrong_technology["rrbs"] = {"cut_sites": ["C|CGG"]}
        wrong_technology["coverage"] = document["coverage"]
        with self.assertRaisesRegex(PipelineError, "WGBS only"):
            run_prepared(
                self.prepared(wrong_technology),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        self.assertFalse((self.directory / "output").exists())

    def test_tbs_target_score_passes_only_its_technology_gate(self) -> None:
        (self.directory / "targets.bed").write_text(
            "chr1\t3\t4\ttarget\t3\t+\n", encoding="utf-8"
        )
        document = baseline_config()
        document["technology"] = "TBS"
        document["tbs"] = {
            "bed": "targets.bed",
            "center_sd": 0,
        }
        document["coverage"] = {"kind": "target-score"}
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        wrong_technology = baseline_config()
        wrong_technology["coverage"] = {"kind": "target-score"}
        with self.assertRaisesRegex(PipelineError, "TBS only"):
            run_prepared(
                self.prepared(wrong_technology),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

    def test_vcf_passes_wgbs_rrbs_and_tbs_capability_gates(self) -> None:
        (self.directory / "variants.vcf").write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "chr1\t4\t.\tT\tC\t.\tPASS\t.\tGT\t1|0\n",
            encoding="utf-8",
        )
        document = baseline_config()
        document["inputs"] = {"vcf": "variants.vcf"}
        prepared = self.prepared(document)

        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                prepared,
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        profile_bytes = b"0.5\n0.5\n"
        (self.directory / "coverage.tsv").write_bytes(profile_bytes)
        document["coverage"] = {
            "kind": "profile",
            "artifact": {
                "path": "coverage.tsv",
                "sha256": hashlib.sha256(profile_bytes).hexdigest(),
            },
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document["technology"] = "RRBS"
        document["rrbs"] = {"cut_sites": ["C|CGG"]}
        document["coverage"] = {"kind": "uniform"}
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        (self.directory / "targets.bed").write_text(
            "chr1\t3\t4\tprobe\t1\t+\n", encoding="utf-8"
        )
        document = baseline_config()
        document["inputs"] = {"vcf": "variants.vcf"}
        document["technology"] = "TBS"
        document["tbs"] = {
            "bed": "targets.bed",
            "center_sd": 0,
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document = baseline_config()
        document["inputs"] = {"vcf": "variants.vcf"}
        document["methylation"]["update_variant_boundaries"] = False
        with self.assertRaisesRegex(PipelineError, "update_variant_boundaries=true"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )
        self.assertFalse((self.directory / "output").exists())

    def test_methylation_profiles_and_asm_formats_pass_capability_gate(self) -> None:
        (self.directory / "levels.cgmap").write_text(
            "chr1\tC\t2\tCG\tCG\t1\t4\t4\n",
            encoding="utf-8",
        )
        (self.directory / "variants.vcf").write_text("fixture\n", encoding="utf-8")
        (self.directory / "levels.asm").write_text("fixture\n", encoding="utf-8")
        (self.directory / "levels.bedmethyl").write_text(
            "fixture\n", encoding="utf-8"
        )
        (self.directory / "levels.asm.bed").write_text(
            "fixture\n", encoding="utf-8"
        )
        document = baseline_config()
        document["inputs"] = {"cgmap": "levels.cgmap"}
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document["methylation"]["cgmap_pool"] = True
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        missing_profile = baseline_config()
        missing_profile["methylation"]["cgmap_pool"] = True
        with self.assertRaisesRegex(PipelineError, "text methylation profile"):
            run_prepared(
                self.prepared(missing_profile),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document["inputs"] = {
            "vcf": "variants.vcf",
            "cgmap": "levels.cgmap",
            "asm": "levels.asm",
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document["inputs"] = {
            "cgmap": "levels.cgmap",
            "asm": "levels.asm",
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        document["inputs"] = {
            "vcf": "variants.vcf",
            "bed_methyl": "levels.bedmethyl",
            "asm_bed": "levels.asm.bed",
        }
        with self.assertRaisesRegex(PipelineError, "cannot resolve"):
            run_prepared(
                self.prepared(document),
                core_executable=self.directory / "missing-core",
                run_id=RUN_ID,
            )

        self.assertFalse((self.directory / "output").exists())

if __name__ == "__main__":
    unittest.main()
