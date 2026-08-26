"""Exercise deterministic VCF.gz export and core round-trip loading."""

from __future__ import annotations

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import uuid

from bsreadsim.cli import build_parser, build_run_document
from bsreadsim.htsim.launch import build_core_argv
from bsreadsim.run.config import normalize_run_config
from bsreadsim.run.prepare import prepare_run


def _export(root: Path, core: Path, output: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            "build",
            "variants",
            "--reference",
            "reference.fa",
            "--output",
            output,
            "--mutation-rate",
            "0.2",
            "--seed-mut",
            "71",
            "--core",
            str(core),
        ],
        cwd=root,
        env=os.environ.copy(),
        check=False,
        capture_output=True,
        text=True,
    )


def _run_with_saved_vcf(root: Path, core: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            "run",
            "wgbs",
            "--reference",
            "reference.fa",
            "--output",
            "saved-run",
            "--fragments",
            "1",
            "--single-end",
            "--read-length",
            "25",
            "--insert-mean",
            "25",
            "--insert-sd",
            "0",
            "--max-ambiguous-fraction",
            "0",
            "--mutation-rate",
            "0.2",
            "--seed-mut",
            "71",
            "--seed",
            "5",
            "--save-vcf",
            "--core",
            str(core),
        ],
        cwd=root,
        env=os.environ.copy(),
        check=False,
        capture_output=True,
        text=True,
    )


def _build_empty_vcf(root: Path, core: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            "build",
            "variants",
            "--reference",
            "reference.fa",
            "--output",
            "empty.vcf.gz",
            "--mutation-rate",
            "0",
            "--core",
            str(core),
        ],
        cwd=root,
        env=os.environ.copy(),
        check=False,
        capture_output=True,
        text=True,
    )


def _build_from_vcf(root: Path, core: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            "build",
            "variants",
            "--reference",
            "reference.fa",
            "--output",
            "normalized.vcf.gz",
            "--vcf",
            "input.vcf",
            "--seed-phase",
            "13",
            "--core",
            str(core),
        ],
        cwd=root,
        env=os.environ.copy(),
        check=False,
        capture_output=True,
        text=True,
    )


def _run_with_input_vcf_truth(
    root: Path, core: Path
) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            "run",
            "wgbs",
            "--reference",
            "reference.fa",
            "--output",
            "input-vcf-run",
            "--fragments",
            "1",
            "--single-end",
            "--read-length",
            "25",
            "--insert-mean",
            "25",
            "--insert-sd",
            "0",
            "--max-ambiguous-fraction",
            "0",
            "--vcf",
            "input.vcf",
            "--seed-phase",
            "13",
            "--seed",
            "5",
            "--save-vcf",
            "--core",
            str(core),
        ],
        cwd=root,
        env=os.environ.copy(),
        check=False,
        capture_output=True,
        text=True,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_variant_catalog.py CORE_EXECUTABLE")
    core = Path(sys.argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-variant-catalog-") as temporary:
        root = Path(temporary).resolve()
        (root / "reference.fa").write_text(
            ">chr1\n" + "ACGT" * 300 + "\n",
            encoding="ascii",
        )

        first = _export(root, core, "first.vcf.gz")
        second = _export(root, core, "second.vcf.gz")
        if first.returncode != 0 or second.returncode != 0:
            raise SystemExit(
                "variant export failed: first={!r} second={!r}".format(
                    first.stderr, second.stderr
                )
            )
        first_path = root / "first.vcf.gz"
        second_path = root / "second.vcf.gz"
        first_bytes = first_path.read_bytes()
        if first_bytes != second_path.read_bytes():
            raise SystemExit("fixed-input variant VCF.gz export is not byte-stable")
        bgzf_eof = bytes.fromhex(
            "1f8b08040000000000ff0600424302001b0003000000000000000000"
        )
        if first_bytes[12:14] != b"BC" or not first_bytes.endswith(bgzf_eof):
            raise SystemExit("variant export is gzip but not canonical BGZF")
        with gzip.open(first_path, "rt", encoding="ascii", newline="") as input_file:
            text = input_file.read()
        if not text.startswith("##fileformat=VCFv4.3\n"):
            raise SystemExit("variant export omitted its VCF header")
        rows = [line for line in text.splitlines() if line and line[0] != "#"]
        if not rows or any(
            row.split("\t")[9] not in {"1|0", "0|1", "1|1"}
            for row in rows
        ):
            raise SystemExit("variant export omitted records or phased genotypes")

        empty = _build_empty_vcf(root, core)
        if empty.returncode != 0:
            raise SystemExit(
                "zero-rate variant export failed: {!r}".format(empty.stderr)
            )
        with gzip.open(
            root / "empty.vcf.gz", "rt", encoding="ascii", newline=""
        ) as input_file:
            empty_text = input_file.read()
        if not empty_text.startswith("##fileformat=VCFv4.3\n") or any(
            line and line[0] != "#" for line in empty_text.splitlines()
        ):
            raise SystemExit("zero-rate export is not a header-only VCF.gz")

        (root / "input.vcf").write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0/1\n",
            encoding="ascii",
        )
        normalized = _build_from_vcf(root, core)
        if normalized.returncode != 0:
            raise SystemExit(
                "input VCF normalization failed: {!r}".format(normalized.stderr)
            )
        normalized_path = root / "normalized.vcf.gz"
        with gzip.open(
            normalized_path, "rt", encoding="ascii", newline=""
        ) as input_file:
            normalized_text = input_file.read()
        normalized_rows = [
            line
            for line in normalized_text.splitlines()
            if line and line[0] != "#"
        ]
        if len(normalized_rows) != 1 or normalized_rows[0].split("\t")[9] not in {
            "0|1",
            "1|0",
        }:
            raise SystemExit("build variants did not phase the input VCF")

        input_truth_run = _run_with_input_vcf_truth(root, core)
        if input_truth_run.returncode != 0:
            raise SystemExit(
                "run --vcf --save-vcf failed: {!r}".format(
                    input_truth_run.stderr
                )
            )
        input_truth_path = (
            root / "input-vcf-run" / "truth" / "sim.variants.vcf.gz"
        )
        if input_truth_path.read_bytes() != normalized_path.read_bytes():
            raise SystemExit(
                "run --vcf --save-vcf diverged from build variants --vcf"
            )

        saved_run = _run_with_saved_vcf(root, core)
        if saved_run.returncode != 0:
            raise SystemExit(
                "run --save-vcf failed: {!r}".format(saved_run.stderr)
            )
        saved_vcf = root / "saved-run" / "truth" / "sim.variants.vcf.gz"
        if saved_vcf.read_bytes() != first_bytes:
            raise SystemExit("run --save-vcf diverged from build variants")

        arguments = build_parser().parse_args(
            [
                "run",
                "wgbs",
                "--reference",
                "reference.fa",
                "--output",
                "roundtrip-run",
                "--fragments",
                "1",
                "--single-end",
                "--read-length",
                "25",
                "--insert-mean",
                "25",
                "--insert-sd",
                "0",
                "--max-ambiguous-fraction",
                "0",
                "--vcf",
                "first.vcf.gz",
                "--seed-phase",
                "71",
            ]
        )
        loaded = normalize_run_config(build_run_document(arguments, root), root)
        prepared = prepare_run(loaded.with_master_seed(5))
        argv = build_core_argv(
            prepared,
            str(uuid.uuid4()),
            core,
            emit_details=False,
        )
        roundtrip = subprocess.run(
            argv,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            check=False,
        )
        if roundtrip.returncode != 0:
            raise SystemExit(
                "exported VCF.gz was not accepted as input: {!r}".format(
                    roundtrip.stderr.decode("utf-8", errors="replace")
                )
            )

        collision = _export(root, core, "first.vcf.gz")
        if collision.returncode == 0 or "already exists" not in collision.stderr:
            raise SystemExit("variant export overwrote an existing destination")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
