"""Exercise an installed wheel without a source-tree core override."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import List, Optional

from bsreadsim import __version__
from bsreadsim.manifest import verify_complete_manifest
from bsreadsim.runtime import packaged_core_candidate, resolve_core_executable


def _baseline_arguments(output_directory: str) -> List[str]:
    return [
        "-r",
        "tiny.fa",
        "-o",
        output_directory,
        "-n",
        "4",
        "--seed",
        "17",
        "--mutation-rate",
        "0",
        "--read-length",
        "3",
        "--insert-size",
        "5",
        "--max-ambiguous-fraction",
        "0",
        "--beta-cg",
        "2",
        "5",
        "--beta-chg",
        "3",
        "4",
        "--beta-chh",
        "5",
        "2",
        "--conversion-rate",
        "0.998",
        "--phred",
        "37",
        "--error-rate",
        "0.01",
        "--coverage-profile",
        "coverage.tsv",
        "--workers",
        "2",
        "--chunk-size",
        "2",
        "--max-in-flight-fragments",
        "2",
        "--prefix",
        "sample",
        "--compression",
        "gzip",
    ]


def _run_cli(
    directory: Path,
    output_directory: str,
    *,
    mode: Optional[str] = None,
):
    arguments = [
        sys.executable,
        "-m",
        "bsreadsim",
        "run",
        *_baseline_arguments(output_directory),
    ]
    if mode is not None:
        arguments.extend(("--mode", mode))
    completed = subprocess.run(
        arguments,
        cwd=str(directory),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0 or completed.stderr:
        raise SystemExit(
            "installed CLI failed: status={} stderr={!r}".format(
                completed.returncode, completed.stderr
            )
        )
    manifest_path = Path(completed.stdout.strip())
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verify_complete_manifest(manifest)
    if manifest["versions"]["protocol"] != "2.0":
        raise SystemExit("installed manifest recorded the wrong protocol")
    return manifest_path, manifest


def main() -> int:
    packaged = packaged_core_candidate().resolve(strict=True)
    if resolve_core_executable() != packaged:
        raise SystemExit("the installed package did not select its bundled core")
    version = subprocess.run(
        [str(packaged), "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    expected_version = "htsim-core {}".format(__version__)
    if version.returncode != 0 or version.stdout.strip() != expected_version:
        raise SystemExit("the bundled core version check failed: {!r}".format(version))

    with tempfile.TemporaryDirectory(prefix="bsreadsim-installed-") as temporary:
        directory = Path(temporary).resolve()
        (directory / "tiny.fa").write_bytes(b">chr1\nACGTCGTAA\n")
        profile_bytes = b"0.5\n0.5\n"
        (directory / "coverage.tsv").write_bytes(profile_bytes)
        manifest_path, manifest = _run_cli(directory, "output")
        expected_manifest_path = directory / "output" / "sample.manifest.json"
        if manifest_path != expected_manifest_path:
            raise SystemExit("installed CLI reported the wrong manifest")
        normalized = manifest["config"]["normalized"]
        if normalized["technology"] != "WGBS":
            raise SystemExit("installed WGBS default was not materialized")
        artifact = normalized["coverage"]["artifact"]
        if artifact["sha256"] != hashlib.sha256(profile_bytes).hexdigest():
            raise SystemExit("installed direct CLI recorded the wrong profile digest")
        if artifact["version"] != "wgbs-gc-target-v1":
            raise SystemExit("installed direct CLI recorded the wrong profile version")
        if manifest["counts"]["core"]["fragment_count"] != 4:
            raise SystemExit("installed pipeline emitted the wrong fragment count")

        if {item["role"] for item in manifest["outputs"]} != {"read1", "read2"}:
            raise SystemExit("installed default mode did not remain FASTQ-only")
        if (directory / "output" / "sample.truth.jsonl.gz").exists():
            raise SystemExit("installed default mode emitted a Truth artifact")

        _, debug_manifest = _run_cli(directory, "debug-output", mode="debug")
        if {item["role"] for item in debug_manifest["outputs"]} != {
            "read1", "read2", "truth"
        }:
            raise SystemExit("installed debug mode omitted Truth")
        if debug_manifest["counts"]["core"] != manifest["counts"]["core"]:
            raise SystemExit("installed debug mode changed core scientific counts")
        for field in ("fragment_count", "mate_count"):
            if (
                debug_manifest["counts"]["python"][field]
                != manifest["counts"]["python"][field]
            ):
                raise SystemExit(
                    "installed debug mode changed Python {}".format(field)
                )

        fastq_output = directory / "output"
        debug_output = directory / "debug-output"
        for name in ("sample.R1.fastq.gz", "sample.R2.fastq.gz"):
            if (fastq_output / name).read_bytes() != (debug_output / name).read_bytes():
                raise SystemExit(
                    "installed debug mode changed {} bytes".format(name)
                )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
