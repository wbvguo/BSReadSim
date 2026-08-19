"""Validate the public RRBS/TBS CLI examples through the real core boundary."""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
from types import SimpleNamespace

from bsreadsim.cli import main as cli_main
from bsreadsim.manifest import verify_complete_manifest
from bsreadsim.runtime import resolve_core_executable


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXAMPLE_ROOT = REPOSITORY_ROOT / "data" / "experiments"


def _copy_examples(directory: Path) -> None:
    for name in (
        "mock-reference.fa",
        "mock-targets.bed",
        "mock-variants.vcf",
    ):
        shutil.copy2(str(EXAMPLE_ROOT / name), str(directory / name))


def _truth_rows(output_directory: Path, prefix: str) -> list[dict]:
    truth = output_directory / "{}.truth.jsonl".format(prefix)
    return [
        json.loads(line)
        for line in truth.read_text(encoding="utf-8").splitlines()
    ]


def _validate_manifest(result, technology: str) -> dict:
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    verify_complete_manifest(manifest)
    normalized = manifest["config"]["normalized"]
    if normalized["technology"] != technology:
        raise SystemExit("{} manifest lost its technology".format(technology))
    if manifest["counts"]["core"]["fragment_count"] != 24:
        raise SystemExit("{} emitted the wrong fragment count".format(technology))
    if manifest["counts"]["core"]["skipped_fragment_count"] != 0:
        raise SystemExit("{} unexpectedly skipped a fragment".format(technology))
    return manifest


def _validate_variant_run(result, technology: str, rows: list[dict]) -> None:
    manifest = _validate_manifest(result, technology)
    roles = {item["role"] for item in manifest["inputs"]}
    if "input.vcf" not in roles:
        raise SystemExit("{} manifest lost its verified VCF snapshot".format(technology))
    event_rows = [row for row in rows if row["variant_events"]]
    if not event_rows:
        raise SystemExit("{} mock never emitted typed variant truth".format(technology))
    for row in event_rows:
        if row["haplotype"] != 0:
            raise SystemExit("{} emitted a haplotype-0 event on haplotype 1".format(technology))
        if any(
            event["kind"] != "INSERTION" or event["phased_haplotype"] != 0
            for event in row["variant_events"]
        ):
            raise SystemExit("{} variant truth lost event kind or phase".format(technology))


def _run_cli(arguments: list[str]) -> Path:
    stdout = io.StringIO()
    stderr = io.StringIO()
    with contextlib.redirect_stdout(stdout), contextlib.redirect_stderr(stderr):
        status = cli_main(arguments)
    if status != 0:
        raise SystemExit(
            "direct CLI failed: {}".format(stderr.getvalue().strip())
        )
    lines = [line for line in stdout.getvalue().splitlines() if line]
    if len(lines) != 1:
        raise SystemExit("direct CLI did not report exactly one output path")
    return Path(lines[0]).resolve(strict=True)


def main(argv: list[str]) -> int:
    if len(argv) > 2:
        raise SystemExit("usage: check_technology_examples.py [CORE_EXECUTABLE]")
    core = (
        resolve_core_executable()
        if len(argv) == 1
        else Path(argv[1]).resolve(strict=True)
    )

    with tempfile.TemporaryDirectory(
        prefix="bsreadsim-technology-examples-"
    ) as temporary:
        directory = Path(temporary).resolve()
        _copy_examples(directory)

        rrbs_domain_arguments = [
            "--seed", "20260813",
            "--cut-site", "C|CGG",
            "--read-length", "4",
            "--insert-min", "8",
            "--insert-mean", "16",
            "--insert-max", "24",
            "--insert-stddev", "0",
            "--max-ambiguous-fraction", "0",
            "--mutation-rate", "0",
        ]
        direct_run_arguments = [
            "-r", str(directory / "mock-reference.fa"),
            "-n", "24",
            "--technology", "RRBS",
            *rrbs_domain_arguments,
            "--beta-cg", "2", "5",
            "--beta-chg", "1", "9",
            "--beta-chh", "1", "19",
            "--conversion-rate", "1",
            "--phred", "35",
            "--error-rate", "0",
            "--workers", "2",
            "--chunk-size", "5",
            "--max-in-flight-fragments", "4",
            "--compression", "none",
            "--mode", "debug",
            "--core", str(core),
        ]
        rrbs_manifest = _run_cli(
            [
                "run",
                *direct_run_arguments,
                "-o", str(directory / "runs" / "rrbs-mock"),
                "-p", "rrbs_mock",
            ]
        )
        rrbs = SimpleNamespace(manifest_path=rrbs_manifest)
        _validate_manifest(rrbs, "RRBS")
        rrbs_rows = _truth_rows(directory / "runs" / "rrbs-mock", "rrbs_mock")
        cut_boundaries = {3, 11, 19, 27, 35, 43, 51}
        for row in rrbs_rows:
            mates = row["mates"]
            fragment_start = mates[0]["reference_start"]
            fragment_end = mates[1]["reference_end"]
            if (
                fragment_start not in cut_boundaries
                or fragment_end not in cut_boundaries
                or fragment_end - fragment_start not in {8, 16, 24}
            ):
                raise SystemExit("RRBS emitted a non-MspI-bounded fragment")
            if row["fragment_conversion_mode"] != "C_TO_T":
                raise SystemExit("RRBS directional conversion mode changed")

        candidate_bed = directory / "rrbs-candidates.bed"
        exported = _run_cli(
            [
                "catalog", "rrbs",
                "-r", str(directory / "mock-reference.fa"),
                "-o", str(candidate_bed),
                *rrbs_domain_arguments,
                "--core", str(core),
            ]
        )
        if exported != candidate_bed or not candidate_bed.is_file():
            raise SystemExit("RRBS candidate command did not publish its BED")
        candidate_lines = candidate_bed.read_text(encoding="utf-8").splitlines()
        expected_header = (
            "#chrom\tstart\tend\tcandidate_id\tscore\tstrand"
            "\thaplotype_mask\ttemplate_length\tgc_count"
            "\trestriction_site_count"
        )
        if not candidate_lines or candidate_lines[0] != expected_header:
            raise SystemExit("RRBS candidate BED header changed")
        if any(
            token in candidate_lines[0].lower()
            for token in ("sha", "hash", "fingerprint")
        ):
            raise SystemExit("RRBS candidate BED unexpectedly contains a hash field")
        candidate_rows = [line.split("\t") for line in candidate_lines[1:]]
        if not candidate_rows or any(len(row) != 10 for row in candidate_rows):
            raise SystemExit("RRBS candidate BED did not emit ten-field rows")
        if any(
            "@h" in row[3] or len(row[3]) > len(row[0]) + 32
            for row in candidate_rows
        ):
            raise SystemExit("RRBS candidate ID is not short and haplotype-neutral")

        # Uniform mode consumes the same exchange BED, accepts absent scores,
        # and ignores both score values and Python row order.
        uniform_bed = directory / "rrbs-candidates-uniform.bed"
        uniform_rows = []
        for row in reversed(candidate_rows):
            changed = list(row)
            changed[4] = "."
            uniform_rows.append("\t".join(changed))
        uniform_bed.write_text(
            expected_header + "\n" + "\n".join(uniform_rows) + "\n",
            encoding="utf-8",
        )
        uniform_manifest = _run_cli(
            [
                "run",
                *direct_run_arguments,
                "-o", str(directory / "runs" / "rrbs-bed-uniform"),
                "-p", "rrbs_bed_uniform",
                "--rrbs-candidates", str(uniform_bed),
            ]
        )
        uniform_result = SimpleNamespace(manifest_path=uniform_manifest)
        _validate_manifest(uniform_result, "RRBS")
        uniform_truth = _truth_rows(
            directory / "runs" / "rrbs-bed-uniform", "rrbs_bed_uniform"
        )
        if uniform_truth != rrbs_rows:
            raise SystemExit("RRBS uniform BED round-trip changed fixed-seed output")

        # Profile mode gives all probability mass to one known candidate, so
        # every emitted fragment must use that reference envelope.
        selected = candidate_rows[0]
        selected_envelope = (int(selected[1]), int(selected[2]))
        profile_rows = []
        for row in reversed(candidate_rows):
            changed = list(row)
            changed[4] = "1" if row[3] == selected[3] else "0"
            profile_rows.append("\t".join(changed))
        profile_bed = directory / "rrbs-candidates-profile.bed"
        profile_bed.write_text(
            expected_header + "\n" + "\n".join(profile_rows) + "\n",
            encoding="utf-8",
        )
        profile_manifest_path = _run_cli(
            [
                "run",
                *direct_run_arguments,
                "-o", str(directory / "runs" / "rrbs-profile"),
                "-p", "rrbs_profile",
                "--rrbs-candidates", str(profile_bed),
                "--rrbs-score",
            ]
        )
        profile_result = SimpleNamespace(manifest_path=profile_manifest_path)
        profile_manifest = _validate_manifest(profile_result, "RRBS")
        if profile_manifest["config"]["normalized"]["coverage"] != {
            "kind": "profile"
        }:
            raise SystemExit("RRBS candidate profile mode was not recorded")
        if any(
            item["role"] == "input.rrbs-candidate-bed"
            for item in profile_manifest["inputs"]
        ):
            raise SystemExit("RRBS profile unexpectedly hashed its candidate BED")
        if not profile_manifest["config"]["normalized"]["rrbs"][
            "candidate_bed"
        ].endswith("rrbs-candidates-profile.bed"):
            raise SystemExit("RRBS profile manifest lost its candidate BED path")
        profile_truth = _truth_rows(
            directory / "runs" / "rrbs-profile", "rrbs_profile"
        )
        for row in profile_truth:
            envelope = (
                row["mates"][0]["reference_start"],
                row["mates"][1]["reference_end"],
            )
            if envelope != selected_envelope:
                raise SystemExit("RRBS profile selected a zero-weight candidate")

        tbs_run_arguments = [
            "-r", str(directory / "mock-reference.fa"),
            "-n", "24",
            "--technology", "TBS",
            "--targets", str(directory / "mock-targets.bed"),
            "--target-score",
            "--fragment-center-stddev", "0",
            "--seed", "20260813",
            "--read-length", "4",
            "--insert-size", "12",
            "--max-ambiguous-fraction", "0",
            "--mutation-rate", "0",
            "--beta-cg", "2", "5",
            "--beta-chg", "1", "9",
            "--beta-chh", "1", "19",
            "--conversion-rate", "1",
            "--phred", "35",
            "--error-rate", "0",
            "--workers", "2",
            "--chunk-size", "5",
            "--max-in-flight-fragments", "4",
            "--compression", "none",
            "--mode", "debug",
            "--core", str(core),
        ]
        tbs_manifest = _run_cli(
            [
                "run",
                *tbs_run_arguments,
                "-o", str(directory / "runs" / "tbs-mock"),
                "-p", "tbs_mock",
            ]
        )
        tbs = SimpleNamespace(manifest_path=tbs_manifest)
        manifest = _validate_manifest(tbs, "TBS")
        roles = {item["role"] for item in manifest["inputs"]}
        if "input.tbs-bed" not in roles:
            raise SystemExit("TBS manifest lost its verified BED snapshot")
        if manifest["config"]["normalized"]["coverage"] != {
            "kind": "target-score"
        }:
            raise SystemExit("TBS manifest lost target-score coverage")

        tbs_rows = _truth_rows(directory / "runs" / "tbs-mock", "tbs_mock")
        expected = {
            (10, 22): "C_TO_T",
            (26, 38): "G_TO_A",
        }
        observed = set()
        for row in tbs_rows:
            mates = row["mates"]
            envelope = (
                mates[0]["reference_start"],
                mates[1]["reference_end"],
            )
            if envelope not in expected:
                raise SystemExit("TBS selected the zero-weight or invalid target")
            if row["fragment_conversion_mode"] != expected[envelope]:
                raise SystemExit("TBS BED strand did not control conversion mode")
            observed.add(envelope)
        if observed != set(expected):
            raise SystemExit("TBS weighted mock did not exercise both positive targets")

        rrbs_vcf_manifest = _run_cli(
            [
                "run",
                *direct_run_arguments,
                "--vcf", str(directory / "mock-variants.vcf"),
                "-o", str(directory / "runs" / "rrbs-vcf-mock"),
                "-p", "rrbs_vcf_mock",
            ]
        )
        rrbs_vcf = SimpleNamespace(manifest_path=rrbs_vcf_manifest)
        rrbs_vcf_rows = _truth_rows(
            directory / "runs" / "rrbs-vcf-mock", "rrbs_vcf_mock"
        )
        _validate_variant_run(rrbs_vcf, "RRBS", rrbs_vcf_rows)

        tbs_vcf_manifest = _run_cli(
            [
                "run",
                *tbs_run_arguments,
                "--vcf", str(directory / "mock-variants.vcf"),
                "-o", str(directory / "runs" / "tbs-vcf-mock"),
                "-p", "tbs_vcf_mock",
            ]
        )
        tbs_vcf = SimpleNamespace(manifest_path=tbs_vcf_manifest)
        tbs_vcf_rows = _truth_rows(
            directory / "runs" / "tbs-vcf-mock", "tbs_vcf_mock"
        )
        _validate_variant_run(tbs_vcf, "TBS", tbs_vcf_rows)

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
