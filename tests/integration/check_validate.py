"""Exercise public input validation, including sparse reference coverage."""

from __future__ import annotations

import gzip
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def _run(core: Path, *arguments: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            "validate",
            *arguments,
            "--core",
            str(core),
        ],
        check=False,
        capture_output=True,
        env=os.environ.copy(),
        text=True,
    )


def _write(path: Path, text: str) -> None:
    path.write_text(text, encoding="ascii", newline="")


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_validate.py CORE_EXECUTABLE")
    core = Path(sys.argv[1]).resolve(strict=True)

    with tempfile.TemporaryDirectory(prefix="bsreadsim-validate-") as temporary:
        root = Path(temporary)
        reference = root / "reference.fa.gz"
        with gzip.open(reference, "wt", encoding="ascii", newline="") as output:
            for index in range(1, 26):
                output.write(">chr{:02d}\nACGTCGTA\n".format(index))

        vcf = root / "subset.vcf"
        _write(
            vcf,
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            "chr03\t1\tkept-snv\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
            "chr03\t2\treference\tC\tG\t.\tPASS\t.\tGT\t0/0\n"
            "chr20\t1\tlong-insertion\tA\tAAAAAA\t.\tPASS\t.\tGT\t0/1\n"
            "chr20\t4\tkept-insertion\tT\tTA\t.\tPASS\t.\tGT\t1/1\n"
            "chr20\t5\tcomplex\tC\tAG\t.\tPASS\t.\tGT\t0/1\n"
            "chr20\t7\tmnp\tTA\tCG\t.\tPASS\t.\tGT\t1/1\n",
        )
        cgmap = root / "subset.CGmap"
        _write(
            cgmap,
            "chr03\tC\t2\tCG\tCG\t0.5\t1\t2\n"
            "chr20\tC\t2\tCG\tCG\t0.75\t3\t4\n",
        )

        valid = _run(
            core,
            "-r",
            str(reference),
            "--vcf",
            str(vcf),
            "--cgmap",
            str(cgmap),
            "--json",
        )
        if valid.returncode != 0:
            raise SystemExit(
                "valid sparse inputs were rejected:\n{}\n{}".format(
                    valid.stdout, valid.stderr
                )
            )
        summary = json.loads(valid.stdout)
        if (
            summary["reference"]["contigs"] != 25
            or summary["vcf"]["contigs"] != 2
            or summary["methylation"]["contigs"] != 2
            or summary["vcf"]["rows"] != 6
            or summary["vcf"]["retained"] != 2
            or summary["vcf"]["reference_genotypes"] != 1
            or summary["vcf"]["skipped"]
            != {
                "total": 3,
                "mnp": 1,
                "complex_replacement": 1,
                "long_indel": 1,
            }
        ):
            raise SystemExit("sparse-input validation summary is incorrect")

        strict = _run(
            core,
            "-r",
            str(reference),
            "--vcf",
            str(vcf),
            "--json",
            "--strict",
        )
        if strict.returncode != 1 or json.loads(strict.stdout) != {
            **summary,
            "methylation": None,
        }:
            raise SystemExit("strict validation did not fail with usable JSON")

        reversed_vcf = root / "reversed.vcf"
        _write(
            reversed_vcf,
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            "chr20\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n"
            "chr03\t1\t.\tA\tT\t.\tPASS\t.\tGT\t0/1\n",
        )
        unsorted_vcf = _run(
            core, "-r", str(reference), "--vcf", str(reversed_vcf)
        )
        if unsorted_vcf.returncode == 0 or "sorted in reference order" not in unsorted_vcf.stderr:
            raise SystemExit("out-of-order subset VCF was not rejected clearly")

        duplicate_cgmap = root / "duplicate.CGmap"
        _write(
            duplicate_cgmap,
            "chr03\tC\t2\tCG\tCG\t0.5\t1\t2\n"
            "chr03\tC\t2\tCG\tCG\t0.5\t1\t2\n",
        )
        unsorted_cgmap = _run(
            core, "-r", str(reference), "--cgmap", str(duplicate_cgmap)
        )
        if (
            unsorted_cgmap.returncode == 0
            or "FASTA order with unique positions" not in unsorted_cgmap.stderr
        ):
            raise SystemExit("duplicate CGmap position was not rejected clearly")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
