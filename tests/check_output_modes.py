"""Exercise the public production/debug CLI contract with the real core."""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

from bsreadsim.manifest import verify_complete_manifest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def _environment() -> dict:
    environment = os.environ.copy()
    python_paths = [str(REPOSITORY_ROOT / "src")]
    if environment.get("PYTHONPATH"):
        python_paths.append(environment["PYTHONPATH"])
    environment["PYTHONPATH"] = os.pathsep.join(python_paths)
    return environment


def _run_direct(root: Path, core: Path, output: str, *arguments: str):
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            "run",
            "-r",
            "tiny.fa",
            "-o",
            output,
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
            "--workers",
            "1",
            "--core-workers",
            "1",
            "--chunk-size",
            "4",
            "--max-in-flight-fragments",
            "4",
            "--prefix",
            "sample",
            "--compression",
            "gzip",
            "--core",
            str(core),
            *arguments,
        ],
        cwd=str(root),
        env=_environment(),
        check=False,
        capture_output=True,
        text=True,
    )


def _manifest(completed: subprocess.CompletedProcess) -> dict:
    if completed.returncode != 0 or completed.stderr:
        raise SystemExit(
            "CLI run failed: status={} stderr={!r}".format(
                completed.returncode,
                completed.stderr,
            )
        )
    path = Path(completed.stdout.strip())
    document = json.loads(path.read_text(encoding="utf-8"))
    verify_complete_manifest(document)
    return document


def _fastq_bytes(directory: Path, mate: int) -> bytes:
    with gzip.open(directory / "sample.R{}.fastq.gz".format(mate), "rb") as stream:
        return stream.read()


def _require_standard_comment_lines(value: bytes) -> None:
    lines = value.splitlines()
    if not lines or len(lines) % 4 != 0:
        raise SystemExit("FASTQ output does not contain complete four-line records")
    if any(lines[index + 2] != b"+" for index in range(0, len(lines), 4)):
        raise SystemExit("production/debug FASTQ comment line is not exactly '+'")


def _require_read_name_contract(value: bytes, mate: int) -> None:
    lines = value.splitlines()
    pattern = re.compile(
        rb"@chr1:([1-9][0-9]*)-([1-9][0-9]*):([0-9]+)/" + str(mate).encode()
    )
    for ordinal, offset in enumerate(range(0, len(lines), 4)):
        match = pattern.fullmatch(lines[offset])
        if match is None:
            raise SystemExit("FASTQ identifier violates read-name v1")
        left, right, observed_ordinal = map(int, match.groups())
        if left > right or observed_ordinal != ordinal:
            raise SystemExit("FASTQ identifier fields violate read-name v1")


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_output_modes.py CORE_EXECUTABLE")
    core = Path(sys.argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-output-modes-") as temporary:
        root = Path(temporary).resolve()
        (root / "tiny.fa").write_bytes(b">chr1\nACGTCGTAA\n")

        production = _manifest(_run_direct(root, core, "production"))
        debug = _manifest(
            _run_direct(root, core, "debug", "--mode", "debug")
        )

        for name, document, truth, roles in (
            ("production", production, "none", {"read1", "read2"}),
            ("debug", debug, "full", {"read1", "read2", "truth"}),
        ):
            if document["versions"]["protocol"] != "2.0":
                raise SystemExit("{} recorded the wrong protocol".format(name))
            if document["versions"].get("read_name") != "bsreadsim-read-name-v1":
                raise SystemExit("{} omitted read-name v1".format(name))
            if document["config"]["normalized"]["output"]["truth"] != truth:
                raise SystemExit("{} normalized the wrong Truth policy".format(name))
            if {item["role"] for item in document["outputs"]} != roles:
                raise SystemExit("{} emitted the wrong artifact roles".format(name))
            if document["counts"]["core"]["fragment_count"] != 4:
                raise SystemExit("{} changed core fragment accounting".format(name))
            if document["counts"]["python"]["fragment_count"] != 4:
                raise SystemExit("{} changed Python fragment accounting".format(name))

        production_directory = root / "production"
        debug_directory = root / "debug"
        if (production_directory / "sample.truth.jsonl.gz").exists():
            raise SystemExit("production emitted a Truth artifact")
        if not (debug_directory / "sample.truth.jsonl.gz").is_file():
            raise SystemExit("debug did not emit a Truth artifact")

        for mate in (1, 2):
            production_fastq = _fastq_bytes(production_directory, mate)
            debug_fastq = _fastq_bytes(debug_directory, mate)
            if production_fastq != debug_fastq:
                raise SystemExit("production/debug R{} FASTQ bytes differ".format(mate))
            _require_standard_comment_lines(production_fastq)
            _require_read_name_contract(production_fastq, mate)

        rejected = _run_direct(root, core, "conflict", "--truth", "full")
        if rejected.returncode == 0 or "unrecognized arguments" not in rejected.stderr:
            raise SystemExit("CLI silently accepted a user Truth setting")
        if (root / "conflict").exists():
            raise SystemExit("invalid Truth argument created output side effects")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
