"""Verify production/debug scientific parity through the real pipeline."""

from __future__ import annotations

import io
import json
from pathlib import Path
import sys
import tempfile
from contextlib import redirect_stderr, redirect_stdout

from bsreadsim.cli import main as cli_main
from bsreadsim.manifest import verify_complete_manifest


def _run(core: Path, root: Path, name: str, *, mode: str, workers: int):
    stdout = io.StringIO()
    stderr = io.StringIO()
    arguments = [
        "run",
        "-r", str(root / "reference.fa"),
        "-o", str(root / name),
        "-n", "257",
        "--seed", "81985529216486895",
        "--vcf", str(root / "variants.vcf"),
        "--mutation-rate", "0",
        "--indel-fraction", "0.15",
        "--indel-extension-probability", "0.3",
        "--read-length", "12",
        "--insert-size", "24",
        "--max-ambiguous-fraction", "1",
        "--beta-cg", "2", "5",
        "--beta-chg", "3", "4",
        "--beta-chh", "5", "2",
        "--conversion-rate", "0.73",
        "--undirectional",
        "--phred", "31",
        "--error-rate", "0.19",
        "--workers", str(workers),
        "--core-workers", "2",
        "--chunk-size", "17",
        "--max-in-flight-fragments", "64",
        "--prefix", "sample",
        "--compression", "none",
        "--mode", mode,
        "--core", str(core),
    ]
    with redirect_stdout(stdout), redirect_stderr(stderr):
        status = cli_main(arguments)
    if status != 0:
        raise SystemExit("CLI run failed: {}".format(stderr.getvalue().strip()))
    manifest_path = Path(stdout.getvalue().strip()).resolve(strict=True)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    verify_complete_manifest(manifest)
    if manifest["counts"]["core"]["fragment_count"] != 257:
        raise SystemExit("pipeline changed the requested fragment count")
    return manifest_path, manifest


def _fastq_bytes(root: Path, name: str) -> tuple[bytes, bytes]:
    output = root / name
    return (
        (output / "sample.R1.fastq").read_bytes(),
        (output / "sample.R2.fastq").read_bytes(),
    )


def _validate_truth(root: Path, name: str) -> None:
    rows = [
        json.loads(line)
        for line in (root / name / "sample.truth.jsonl")
        .read_text(encoding="utf-8")
        .splitlines()
    ]
    if len(rows) != 257:
        raise SystemExit("debug mode emitted the wrong Truth row count")
    event_kinds = set()
    for expected_ordinal, row in enumerate(rows):
        if row["fragment_ordinal"] != expected_ordinal:
            raise SystemExit("Truth fragment ordinals are not consecutive")
        for event in row["variant_events"]:
            event_kinds.add(event["kind"])
        for mate in row["mates"]:
            annotations = mate["annotations"]
            if [item["read_offset"] for item in annotations] != list(
                range(len(annotations))
            ):
                raise SystemExit("Truth base annotations lost read offsets")
    if event_kinds != {"SNV", "INSERTION", "DELETION"}:
        raise SystemExit("debug Truth did not preserve all variant event kinds")


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit("usage: check_protocol_scientific_parity.py CORE_EXECUTABLE")
    core = Path(argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-scientific-parity-") as tmp:
        root = Path(tmp)
        sequence = bytearray(b"ACGT" * 40)
        sequence[90] = ord("N")
        reference = root / "reference.fa"
        reference.write_bytes(b">chrParity\n" + bytes(sequence) + b"\n")
        vcf = root / "variants.vcf"
        vcf.write_bytes(
            b"##fileformat=VCFv4.3\n"
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            b"chrParity\t25\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\n"
            b"chrParity\t41\t.\tA\tAGG\t.\tPASS\t.\tGT\t1|1\n"
            b"chrParity\t61\t.\tAC\tA\t.\tPASS\t.\tGT\t1|1\n"
        )

        runs = (
            ("production-w1", "production", 1),
            ("production-w2", "production", 2),
            ("debug-w1", "debug", 1),
            ("debug-w2", "debug", 2),
        )
        manifests = {}
        for name, mode, workers in runs:
            _, manifests[name] = _run(core, root, name, mode=mode, workers=workers)

        expected = _fastq_bytes(root, "production-w1")
        for name, _, _ in runs[1:]:
            if _fastq_bytes(root, name) != expected:
                raise SystemExit("mode or worker count changed FASTQ bytes")
        if (root / "production-w1" / "sample.truth.jsonl").exists():
            raise SystemExit("production mode emitted a Truth artifact")
        if {item["role"] for item in manifests["production-w1"]["outputs"]} != {
            "read1", "read2"
        }:
            raise SystemExit("production manifest recorded unexpected output roles")
        if {item["role"] for item in manifests["debug-w1"]["outputs"]} != {
            "read1", "read2", "truth"
        }:
            raise SystemExit("debug manifest omitted its Truth role")
        _validate_truth(root, "debug-w1")
        _validate_truth(root, "debug-w2")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
