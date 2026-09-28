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
from bsreadsim import htsim
from bsreadsim.htsim import CORE_FILENAME, HtsimCore



def require_complete(document: dict) -> None:
    if document.get("status") != "complete":
        raise SystemExit("manifest is not complete")

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
    require_complete(manifest)
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
    if manifest["details"]["software_versions"] != {"core": __version__, "python": __version__}:
        raise SystemExit("installed manifest recorded the wrong versions")
    return manifest_path, manifest


def main() -> int:
    packaged = (Path(htsim.__file__).parent / "bin" / CORE_FILENAME).resolve(strict=True)
    if HtsimCore.find().path != packaged:
        raise SystemExit("the installed package did not select its bundled core")
    version = subprocess.run(
        [str(packaged), "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    expected_version = "htsim {}".format(__version__)
    if version.returncode != 0 or version.stdout.strip() != expected_version:
        raise SystemExit("the bundled core version check failed: {!r}".format(version))
    command = subprocess.run(  # the htsim command runs the bundled core
        [str(Path(sys.executable).parent / "htsim"), "--version"],
        check=False,
        capture_output=True,
        text=True,
    )
    if command.returncode != 0 or command.stdout.strip() != expected_version:
        raise SystemExit("the installed htsim command failed: {!r}".format(command))

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
        profile = [item for item in manifest["inputs"] if item["role"] == "gc_profile"]
        if [item["sha256"] for item in profile] != [hashlib.sha256(profile_bytes).hexdigest()]:
            raise SystemExit("installed direct CLI recorded the wrong profile digest")
        if manifest["summary"]["fragment_count"] != 4:
            raise SystemExit("installed pipeline emitted the wrong fragment count")

        if {item["role"] for item in manifest["outputs"]} != {"read1", "read2"}:
            raise SystemExit("installed run did not emit exactly paired FASTQ")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
