"""Run the sole protocol through the real inline and shared pipeline paths."""

from __future__ import annotations

from contextlib import redirect_stderr, redirect_stdout
import io
import json
import hashlib
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace

from bsreadsim.cli import main as cli_main
from bsreadsim.manifest import verify_complete_manifest


def quality_model() -> bytes:
    mate = {
        "initial_counts": [[1, 3], [3, 1], [1, 3], [3, 1], [1, 3]],
        "transition_counts": [[1, 3], [3, 1]],
    }
    return json.dumps(
        {
            "schema": "bsreadsim-quality-markov-v1",
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
            "schema": "bsreadsim-quality-confusion-v1",
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
    truth: str,
    advanced: bool = False,
):
    arguments = [
        "run",
        "-r", str(root / "tiny.fa"),
        "-o", str(root / name),
        "-n", "7",
        "--seed", "1311768467463790320",
        "--mutation-rate", "0",
        "--indel-fraction", "0.15",
        "--indel-extension-probability", "0.3",
        "--read-length", "3",
        "--insert-size", "5",
        "--max-ambiguous-fraction", "0",
        "--beta-cg", "2", "5",
        "--beta-chg", "3", "4",
        "--beta-chh", "5", "2",
        "--conversion-rate", "0.73",
        "--undirectional",
        "--workers", str(workers),
        "--core-workers", "2",
        "--chunk-size", "3",
        "--max-in-flight-fragments", "2",
        "--prefix", "sample",
        "--compression", "none",
        "--mode", "debug" if truth == "full" else "production",
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
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    verify_complete_manifest(document)
    if document["versions"]["protocol"] != "2.0":
        raise SystemExit("manifest recorded the wrong observed protocol")
    counts = document["counts"]["python"]
    if counts["fragment_count"] != 7 or counts["mate_count"] != 14:
        raise SystemExit("v2 pipeline accounting changed")
    return SimpleNamespace(manifest_path=manifest_path)


def run(core: Path, root: Path, name: str, *, workers: int, truth: str):
    return _run_cli(core, root, name, workers=workers, truth=truth)


def run_advanced(core: Path, root: Path, name: str, *, workers: int):
    return _run_cli(
        core,
        root,
        name,
        workers=workers,
        truth="none",
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
                "-r", str(root / "tiny.fa"),
                "-o", str(output),
                "-n", "7",
                "--seed", "1311768467463790320",
                "--mutation-rate", "0",
                "--read-length", "3",
                "--insert-min", "3",
                "--insert-mean", "5",
                "--insert-max", "8",
                "--insert-stddev", "1",
                "--max-ambiguous-fraction", "0",
                "--coverage-profile", str(profile_path),
                "--error-rate", "0",
                "--workers", "1",
                "--core", str(core),
            ]
        )
    if status != 0:
        raise SystemExit("direct profile CLI failed")
    manifest_path = Path(stdout.getvalue().strip())
    document = json.loads(manifest_path.read_text(encoding="utf-8"))
    verify_complete_manifest(document)
    normalized = document["config"]["normalized"]
    artifact = normalized["coverage"]["artifact"]
    if normalized["coverage"]["kind"] != "profile":
        raise SystemExit("direct CLI lost profile coverage")
    if artifact["sha256"] != hashlib.sha256(profile_bytes).hexdigest():
        raise SystemExit("direct CLI recorded the wrong profile digest")
    fragments = normalized["fragments"]
    if (
        fragments["insert_min"] != 3
        or fragments["insert_mean"] != 5
        or fragments["insert_max"] != 8
        or fragments["insert_stddev"] != 1.0
    ):
        raise SystemExit("direct profile CLI changed variable insert parameters")
    if document["counts"]["core"]["fragment_count"] != 7:
        raise SystemExit("direct profile CLI emitted the wrong fragment count")
    if {item["role"] for item in document["outputs"]} != {"read1", "read2"}:
        raise SystemExit("direct profile CLI did not remain FASTQ-only")


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_pipeline_e2e.py CORE")
    core = Path(sys.argv[1]).resolve()
    with tempfile.TemporaryDirectory(prefix="bsreadsim-pipeline-v2-") as temporary:
        root = Path(temporary)
        (root / "tiny.fa").write_bytes(
            b">chrTiny\nACGTNCGCGTACGTCGCGTACGTNCGCGTACGT\n"
        )
        (root / "quality.json").write_bytes(quality_model())
        (root / "error.json").write_bytes(error_model())

        run_direct_profile_cli(core, root)

        run(core, root, "truth-inline", workers=1, truth="full")
        run(core, root, "truth-pool", workers=2, truth="full")
        expected_truth = data_files(root / "truth-inline")
        if expected_truth != data_files(root / "truth-pool"):
            raise SystemExit("worker count changed ordered FASTQ/Truth bytes")

        run(core, root, "fastq-inline", workers=1, truth="none")
        run(core, root, "fastq-pool", workers=2, truth="none")
        expected_fastq = data_files(root / "fastq-inline")
        if expected_fastq != data_files(root / "fastq-pool"):
            raise SystemExit("worker count changed ordered FASTQ bytes")
        for role in ("sample.R1.fastq", "sample.R2.fastq"):
            if expected_fastq[role] != expected_truth[role]:
                raise SystemExit("Truth policy changed FASTQ bytes")

        # Non-uniform sequencing policies use the compact typed fallback when
        # truth is disabled; worker count must still preserve exact bytes.
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
            if document["versions"]["protocol"] != "2.0":
                raise SystemExit("advanced manifest protocol is wrong")
        if data_files(root / "advanced-inline") != data_files(root / "advanced-pool"):
            raise SystemExit("compact typed fallback changed advanced FASTQ bytes")


if __name__ == "__main__":
    main()
