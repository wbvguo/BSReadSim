"""Exercise fixed MethDB export, reload, binding, and model fallback."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from bsreadsim.run.manifest import verify_complete_manifest


def _run(root: Path, core: Path, *arguments: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "bsreadsim", "run", *arguments, "--core", str(core)],
        cwd=str(root),
        env=os.environ.copy(),
        check=False,
        capture_output=True,
        text=True,
    )


def _manifest(completed: subprocess.CompletedProcess) -> dict:
    if completed.returncode != 0:
        raise SystemExit(
            "MethDB run failed: status={} stderr={!r}".format(
                completed.returncode, completed.stderr
            )
        )
    path = Path(completed.stdout.strip())
    document = json.loads(path.read_text(encoding="utf-8"))
    verify_complete_manifest(document)
    return document


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_methdb.py CORE_EXECUTABLE")
    core = Path(sys.argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-methdb-") as temporary:
        root = Path(temporary).resolve()
        (root / "reference.fa").write_bytes(
            b">chr1\n" + b"ACGT" * 300 + b"\n"
        )
        common = (
            "-r", "reference.fa",
            "-n", "8",
            "--methdb-seed", "23",
            "--mutation-rate", "0",
            "--read-length", "50",
            "--insert-size", "150",
            "--max-ambiguous-fraction", "0",
            "--fragment-realization",
        )
        exported = _run(
            root,
            core,
            *common,
            "-o", "exported",
            "--seed", "17",
            "--save-methdb", "catalog.methdb",
        )
        if exported.stderr:
            raise SystemExit(f"MethDB export emitted unexpected stderr: {exported.stderr!r}")
        first = _manifest(exported)
        snapshot = root / "catalog.methdb"
        if not snapshot.is_file() or snapshot.stat().st_size == 0:
            raise SystemExit("MethDB export did not publish a snapshot")
        snapshot_sha256 = hashlib.sha256(snapshot.read_bytes()).hexdigest()
        if {item["role"] for item in first["outputs"]} != {"bam"}:
            raise SystemExit("rich BAM MethDB run emitted forbidden sidecars")
        if first["randomness"]["catalog_seed"] != "23":
            raise SystemExit("MethDB manifest lost its catalog seed")
        first_methdb = [
            item for item in first["inputs"] if item["role"] == "input.methdb"
        ]
        if len(first_methdb) != 1 or first_methdb[0]["sha256"] != snapshot_sha256:
            raise SystemExit("MethDB export run did not reload its exact snapshot")

        loaded = _run(
            root,
            core,
            *common,
            "-o", "loaded",
            "--seed", "999",
            "--methdb", "catalog.methdb",
            "--methylation-model", "bilstm",
        )
        if "falling back to BernoulliStateModel" not in loaded.stderr:
            raise SystemExit("BiLSTM request omitted its explicit fallback warning")
        second = _manifest(loaded)
        state_model = second["models"]["methylation_state"]
        if state_model != {
            "contract": "bernoulli-site-v1",
            "effective": "bernoulli",
            "requested": "bilstm",
        }:
            raise SystemExit("BiLSTM fallback is not auditable in the manifest")
        if second["randomness"]["master_seed"] != "999":
            raise SystemExit("loaded MethDB incorrectly froze the simulation seed")
        second_methdb = [
            item for item in second["inputs"] if item["role"] == "input.methdb"
        ]
        if len(second_methdb) != 1 or second_methdb[0]["sha256"] != snapshot_sha256:
            raise SystemExit("loaded run did not authenticate the fixed MethDB")

        rejected = _run(
            root,
            core,
            "-r", "reference.fa",
            "-o", "rejected",
            "-n", "1",
            "--seed", "999",
            "--methdb-seed", "24",
            "--mutation-rate", "0",
            "--read-length", "50",
            "--insert-size", "150",
            "--bam",
            "--methdb", "catalog.methdb",
        )
        if rejected.returncode == 0 or "catalog binding is incompatible" not in rejected.stderr:
            raise SystemExit("MethDB accepted a mismatched catalog seed")
        rejected_root = root / "rejected"
        if rejected_root.exists() and any(rejected_root.iterdir()):
            raise SystemExit("rejected MethDB run published output files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
