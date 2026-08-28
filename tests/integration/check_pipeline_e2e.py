"""Run the sole protocol through the real inline and shared pipeline paths."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import io
import json
import hashlib
from pathlib import Path
import shlex
import sys
import tempfile
from types import SimpleNamespace

from bsreadsim.cli import main as cli_main
from bsreadsim.htsim.protocol import PROTOCOL_VERSION
from bsreadsim.run.manifest import verify_complete_manifest


def quality_model() -> bytes:
    mate = {
        "initial_counts": [[1, 3], [3, 1], [1, 3], [3, 1], [1, 3]],
        "transition_counts": [[1, 3], [3, 1]],
    }
    return json.dumps(
        {
            "schema": "quality-markov",
            "quality_scores": [10, 30],
            "mates": [mate, mate],
        },
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def error_model() -> bytes:
    rotate = [
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
    ]
    mate = {"base_transition_counts": [rotate, rotate]}
    return json.dumps(
        {
            "schema": "quality-confusion",
            "quality_scores": [10, 30],
            "mates": [mate, mate],
        },
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")


def data_files(directory: Path) -> dict:
    return {
        path.name: path.read_bytes()
        for path in sorted(directory.iterdir())
        if not path.name.endswith(".manifest.json")
    }


def _run_cli(
    core: Path,
    root: Path,
    name: str,
    *,
    workers: int,
    advanced: bool = False,
    output_format: str = "fastq",
):
    arguments = [
        "run",
        "wgbs",
        "-r", str(root / "tiny.fa"),
        "-o", str(root / name),
        "-n", "14",
        "--seed", "1311768467463790320",
        "--mutation-rate", "0",
        "--indel-fraction", "0.15",
        "--indel-extension-probability", "0.3",
        "--read-length", "3",
        "--insert-mean", "5",
        "--insert-sd", "0",
        "--max-ambiguous-fraction", "0",
        "--beta-cg", "2,5",
        "--beta-chg", "3,4",
        "--beta-chh", "5,2",
        "--conversion-rate", "0.73",
        "--undirectional",
        "--threads", str(workers),
        "--prefix", "sample",
        "--format", output_format,
        "--core", str(core),
    ]
    if advanced:
        arguments.extend(
            [
                "--quality-model", str(root / "quality.json"),
                "--error-model", str(root / "error.json"),
            ]
        )
    else:
        arguments.extend(["--phred", "31", "--error-rate", "0.19"])
    stdout = io.StringIO()
    stderr = io.StringIO()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        status = cli_main(arguments)
    if status != 0:
        raise SystemExit("direct CLI failed: {}".format(stderr.getvalue().strip()))
    manifest_path = Path(stdout.getvalue().strip()).resolve(strict=True)
    manifest_text = manifest_path.read_text(encoding="utf-8")
    document = json.loads(manifest_text)
    verify_complete_manifest(document)
    if not manifest_text.startswith("{\n  \"command\": {"):
        raise SystemExit("CLI manifest is not pretty-printed JSON")
    if shlex.split(document["command"]["user_command"]) != [
        "bsreadsim",
        *arguments,
    ]:
        raise SystemExit("CLI manifest did not preserve argv")
    if document["details"]["protocol_version"] != PROTOCOL_VERSION:
        raise SystemExit("manifest recorded the wrong observed protocol")
    summary = document["summary"]
    if summary["fragment_count"] != 7 or summary["read_count"] != 14:
        raise SystemExit("pipeline accounting changed")
    return SimpleNamespace(manifest_path=manifest_path)


def run(
    core: Path,
    root: Path,
    name: str,
    *,
    workers: int,
    output_format: str = "fastq",
):
    return _run_cli(
        core,
        root,
        name,
        workers=workers,
        output_format=output_format,
    )


def run_advanced(core: Path, root: Path, name: str, *, workers: int):
    return _run_cli(
        core,
        root,
        name,
        workers=workers,
        advanced=True,
    )


def run_direct_profile_cli(core: Path, root: Path) -> None:
    profile_bytes = b"0\n0.5\n0.5\n"
    profile_path = root / "coverage.tsv"
    profile_path.write_bytes(profile_bytes)
    output = root / "direct-profile"
    stdout = io.StringIO()
    with redirect_stdout(stdout):
        status = cli_main(
            [
                "run",
                "wgbs",
                "-r", str(root / "tiny.fa"),
                "-o", str(output),
                "-n", "14",
                "--seed", "1311768467463790320",
                "--mutation-rate", "0",
                "--read-length", "3",
                "--insert-min", "3",
                "--insert-mean", "5",
                "--insert-max", "8",
                "--insert-sd", "1",
                "--max-ambiguous-fraction", "0",
                "--gc-profile", str(profile_path),
                "--error-rate", "0",
                "--threads", "1",
                "--core", str(core),
            ]
        )
    if status != 0:
        raise SystemExit("direct profile CLI failed")
    manifest_path = Path(stdout.getvalue().strip())
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    verify_complete_manifest(document)
    effective = document["details"]["configuration"]
    artifact = effective["coverage"]["artifact"]
    if effective["coverage"]["type"] != "profile":
        raise SystemExit("direct CLI lost profile coverage")
    if artifact["sha256"] != hashlib.sha256(profile_bytes).hexdigest():
        raise SystemExit("direct CLI recorded the wrong profile digest")
    fragments = effective["fragments"]
    if (
        fragments["insert_min"] != 3
        or fragments["insert_mean"] != 5
        or fragments["insert_max"] != 8
        or fragments["insert_sd"] != 1.0
    ):
        raise SystemExit("direct profile CLI changed variable insert parameters")
    if document["summary"]["fragment_count"] != 7:
        raise SystemExit("direct profile CLI emitted the wrong fragment count")
    if {item["role"] for item in document["outputs"]} != {"read1", "read2"}:
        raise SystemExit("direct profile CLI did not remain FASTQ-only")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_pipeline_e2e.py CORE")
    core = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="bsreadsim-pipeline-") as temporary:
        root = Path(temporary)
        (root / "tiny.fa").write_bytes(
            b">chrTiny\nACGTNCGCGTACGTCGCGTACGTNCGCGTACGT\n"
        )
        (root / "quality.json").write_bytes(quality_model())
        (root / "error.json").write_bytes(error_model())

        run_direct_profile_cli(core, root)

        run(core, root, "fastq-inline", workers=1)
        run(core, root, "fastq-pool", workers=2)
        if data_files(root / "fastq-inline") != data_files(root / "fastq-pool"):
            raise SystemExit("worker count changed ordered FASTQ bytes")

        run(
            core,
            root,
            "fastq-gzip-inline",
            workers=1,
            output_format="fastq.gz",
        )
        run(
            core,
            root,
            "fastq-gzip-pool",
            workers=2,
            output_format="fastq.gz",
        )
        if data_files(root / "fastq-gzip-inline") != data_files(
            root / "fastq-gzip-pool"
        ):
            raise SystemExit("thread count changed ordered FASTQ.gz bytes")

        # Non-uniform sequencing policies use the general typed path, whose
        # worker count must still preserve exact bytes.
        inline_advanced = run_advanced(
            core, root, "advanced-inline", workers=1
        )
        pooled_advanced = run_advanced(
            core, root, "advanced-pool", workers=2
        )
        for result in (inline_advanced, pooled_advanced):
            document = json.loads(
                result.manifest_path.read_text(encoding="utf-8")
            )
            verify_complete_manifest(document)
            if document["details"]["protocol_version"] != PROTOCOL_VERSION:
                raise SystemExit("advanced manifest protocol is wrong")
        if data_files(root / "advanced-inline") != data_files(root / "advanced-pool"):
            raise SystemExit("general typed path changed advanced FASTQ bytes")


if __name__ == "__main__":
    main()
