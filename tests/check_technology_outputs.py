"""Exercise RRBS and TBS through the real annotated BAM pipeline."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

from bsreadsim.run.manifest import verify_complete_manifest
from check_bam import _parse_bam


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
EXPERIMENT_ROOT = REPOSITORY_ROOT / "data" / "experiments"


def _environment() -> dict:
    return os.environ.copy()


def _qname_envelope(query_name: str):
    locus, ordinal_text = query_name.rsplit(":", 1)
    contig, interval = locus.rsplit(":", 1)
    left_text, right_text = interval.split("-", 1)
    return contig, int(left_text), int(right_text), int(ordinal_text, 16)


def _run(core: Path, root: Path, name: str, technology_arguments):
    output = root / name
    command = [
        sys.executable,
        "-m",
        "bsreadsim",
        "run",
        "-r",
        str(root / "mock-reference.fa"),
        "-o",
        str(output),
        "-n",
        "12",
        "--seed",
        "20260813",
        "--read-length",
        "4",
        "--max-ambiguous-fraction",
        "0",
        "--mutation-rate",
        "0",
        "--conversion-rate",
        "1",
        "--phred",
        "35",
        "--error-rate",
        "0",
        "--workers",
        "2",
        "--chunk-size",
        "5",
        "--max-in-flight-fragments",
        "4",
        "--prefix",
        "sample",
        "--bam",
        "--fragment-summary",
        "--core",
        str(core),
        *technology_arguments,
    ]
    completed = subprocess.run(
        command,
        cwd=str(root),
        env=_environment(),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0 or completed.stderr:
        raise SystemExit(
            "{} pipeline failed: status={} stderr={!r}".format(
                name, completed.returncode, completed.stderr
            )
        )
    manifest = json.loads(
        Path(completed.stdout.strip()).read_text(encoding="utf-8")
    )
    verify_complete_manifest(manifest)
    if {item["role"] for item in manifest["outputs"]} != {"bam"}:
        raise SystemExit("{} emitted the wrong artifact roles".format(name))
    _, _, records, _ = _parse_bam(output / "sample.bam")
    if len(records) != 24:
        raise SystemExit("{} emitted the wrong BAM record count".format(name))
    return records


def main(argv) -> int:
    if len(argv) != 2:
        raise SystemExit("usage: check_technology_outputs.py CORE_EXECUTABLE")
    core = Path(argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-technologies-") as value:
        root = Path(value).resolve()
        shutil.copy2(EXPERIMENT_ROOT / "mock-reference.fa", root)
        shutil.copy2(EXPERIMENT_ROOT / "mock-targets.bed", root)

        rrbs = _run(
            core,
            root,
            "rrbs",
            ("--technology", "RRBS", "--cut-site", "C|CGG", "--insert-size", "8"),
        )
        for record in rrbs[::2]:
            contig, left, right, _ = _qname_envelope(record["query_name"])
            if contig != "chrMock" or right - left + 1 != 8:
                raise SystemExit("RRBS emitted a non-MspI fragment envelope")
            if (record["aux"]["zf"][1][0] >> 4) & 0x7:
                raise SystemExit("RRBS changed its directional conversion mode")

        tbs = _run(
            core,
            root,
            "tbs",
            (
                "--technology",
                "TBS",
                "--targets",
                str(root / "mock-targets.bed"),
                "--target-score",
                "--fragment-center-stddev",
                "0",
                "--insert-size",
                "12",
            ),
        )
        expected_modes = {(11, 22): 0, (27, 38): 1}
        observed = set()
        for record in tbs[::2]:
            contig, left, right, _ = _qname_envelope(record["query_name"])
            envelope = (left, right)
            if contig != "chrMock" or envelope not in expected_modes:
                raise SystemExit("TBS selected an invalid or zero-weight target")
            mode = (record["aux"]["zf"][1][0] >> 4) & 0x7
            if mode != expected_modes[envelope]:
                raise SystemExit("TBS strand did not control conversion mode")
            observed.add(envelope)
        if observed != set(expected_modes):
            raise SystemExit("TBS weighted fixture did not exercise both targets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
