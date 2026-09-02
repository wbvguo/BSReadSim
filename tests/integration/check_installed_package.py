"""Exercise an installed wheel without a source-tree core override."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

from bsreadsim import __version__
from bsreadsim.htsim.protocol import PROTOCOL_VERSION
from bsreadsim.run.manifest import verify_complete_manifest
from bsreadsim.htsim.launch import packaged_core_candidate, resolve_core_executable


def _baseline_arguments(output_directory: str) -> list[str]:
    return [
        "-r",
        "tiny.fa",
        "-o",
        output_directory,
        "-n",
        "8",
        "--seed",
        "17",
        "--mutation-rate",
        "0",
        "--read-length",
        "3",
        "--insert-mean",
        "5",
        "--insert-sd",
        "0",
        "--max-ambiguous-fraction",
        "0",
        "--beta-cg",
        "2,5",
        "--beta-chg",
        "3,4",
        "--beta-chh",
        "5,2",
        "--conversion-rate",
        "0.998",
        "--phred",
        "37",
        "--error-rate",
        "0.01",
        "--sampling",
        "gc",
        "--gc-profile",
        "coverage.tsv",
        "--threads",
        "2",
        "--prefix",
        "sample",
        "--format",
        "fastq.gz",
    ]


def _run_cli(
    directory: Path,
    output_directory: str,
):
    arguments = [
        sys.executable,
        "-m",
        "bsreadsim",
        "run",
        "wgbs",
        *_baseline_arguments(output_directory),
    ]
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
    manifest_text = manifest_path.read_text(encoding="utf-8")
    manifest = json.loads(manifest_text)
    verify_complete_manifest(manifest)
    expected_argv = [
        "bsreadsim",
        "run",
        "wgbs",
        *_baseline_arguments(output_directory),
    ]
    if shlex.split(manifest["command"]["user_command"]) != expected_argv:
        raise SystemExit("installed CLI manifest did not preserve argv")
    if not manifest_text.startswith("{\n  \"command\": {"):
        raise SystemExit("installed CLI manifest is not pretty-printed JSON")
    if manifest["details"]["protocol_version"] != PROTOCOL_VERSION:
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
        effective = manifest["details"]["configuration"]
        if effective["technology"] != "WGBS":
            raise SystemExit("installed WGBS default was not materialized")
        artifact = effective["coverage"]["artifact"]
        if artifact["sha256"] != hashlib.sha256(profile_bytes).hexdigest():
            raise SystemExit("installed direct CLI recorded the wrong profile digest")
        if set(artifact) != {"path", "sha256"}:
            raise SystemExit("installed direct CLI recorded obsolete artifact metadata")
        if manifest["summary"]["fragment_count"] != 4:
            raise SystemExit("installed pipeline emitted the wrong fragment count")

        if {item["role"] for item in manifest["outputs"]} != {"read1", "read2"}:
            raise SystemExit("installed run did not emit exactly paired FASTQ")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
