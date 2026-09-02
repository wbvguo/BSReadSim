"""Validate every documented source-distribution example with the real core."""

from __future__ import annotations

import gzip
import os
from pathlib import Path
import subprocess
import sys
import tempfile


def _run(core: Path, *arguments: str) -> None:
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            *arguments,
            "--core",
            str(core),
        ],
        check=False,
        capture_output=True,
        env=os.environ.copy(),
        text=True,
    )
    if completed.returncode != 0:
        raise SystemExit(
            "example-data command failed: {}\nstdout:\n{}\nstderr:\n{}".format(
                " ".join(arguments), completed.stdout, completed.stderr
            )
        )


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_example_data.py CORE_EXECUTABLE")
    core = Path(sys.argv[1]).resolve(strict=True)
    repository = Path(__file__).resolve().parents[2]
    examples = repository / "data" / "examples"

    with gzip.open(
        examples / "test.vcf.gz", "rt", encoding="ascii", newline=""
    ) as input_file:
        input_vcf_rows = [
            line.rstrip("\n").split("\t")
            for line in input_file
            if line and line[0] != "#"
        ]
    input_records = {
        (row[0], int(row[1]), row[3], row[4]) for row in input_vcf_rows
    }
    input_het_snvs = {
        (row[0], int(row[1]), row[3], row[4])
        for row in input_vcf_rows
        if len(row[3]) == 1
        and len(row[4]) == 1
        and row[3] in "ACGT"
        and row[4] in "ACGT"
        and row[9].split(":")[row[8].split(":").index("GT")]
        .replace("|", "/")
        in {"0/1", "1/0"}
    }
    directed_substitution_counts: dict[tuple[str, str], int] = {}
    for _chrom, _pos, ref, alt in input_het_snvs:
        key = (ref, alt)
        directed_substitution_counts[key] = (
            directed_substitution_counts.get(key, 0) + 1
        )
    if (
        len(input_het_snvs) != 227
        or len(directed_substitution_counts) != 12
        or max(directed_substitution_counts.values())
        - min(directed_substitution_counts.values())
        > 1
    ):
        raise SystemExit(
            "example VCF must contain 227 heterozygous SNVs balanced across "
            "all 12 directed substitutions"
        )
    expected_complex_records = {
        ("chr10", 5160, "C", "TCGA"),
        ("chr10", 6083, "CCAC", "A"),
    }
    if (
        len(input_vcf_rows) != 1557
        or not expected_complex_records.issubset(input_records)
    ):
        raise SystemExit("example VCF lost its two skipped complex records")
    expected_chr11_indels = {
        (44957, "A", "AT"),
        (102318, "T", "TAT"),
        (176408, "A", "ATAT"),
        (216504, "A", "AATAT"),
        (239187, "TT", "T"),
        (276207, "TTA", "T"),
        (338258, "TTTT", "T"),
        (410865, "TTTTT", "T"),
    }
    chr11_indels = {
        (int(row[1]), row[3], row[4])
        for row in input_vcf_rows
        if row[0] == "chr11" and len(row[3]) != len(row[4])
    }
    if chr11_indels != expected_chr11_indels:
        raise SystemExit("example VCF lost its chr11 1-4 bp indel matrix")

    reference_sequences: dict[str, list[str]] = {}
    reference_name = ""
    with (examples / "test.fa").open("rt", encoding="ascii") as input_file:
        for line in input_file:
            line = line.strip()
            if line.startswith(">"):
                reference_name = line[1:].split()[0]
                reference_sequences[reference_name] = []
            elif line:
                reference_sequences[reference_name].append(line.upper())
    expected_cgmap_sites: set[tuple[str, str, str]] = set()
    for chromosome, chunks in reference_sequences.items():
        sequence = "".join(chunks)
        for offset in range(len(sequence) - 1):
            if sequence[offset : offset + 2] == "CG":
                expected_cgmap_sites.add((chromosome, "C", str(offset + 1)))
                expected_cgmap_sites.add((chromosome, "G", str(offset + 2)))

    with (examples / "test.CGmap").open(
        "rt", encoding="ascii", newline=""
    ) as input_file:
        plain_cgmap_rows = [
            line.rstrip("\n").split("\t") for line in input_file
        ]
    with gzip.open(
        examples / "test.CGmap.gz", "rt", encoding="ascii", newline=""
    ) as input_file:
        cgmap_rows = [line.rstrip("\n").split("\t") for line in input_file]
    if cgmap_rows != plain_cgmap_rows:
        raise SystemExit("plain and compressed example CGmaps must match")
    observed_cgmap_sites = {
        (row[0], row[1], row[2]) for row in cgmap_rows if len(row) >= 3
    }
    if (
        len(cgmap_rows) != len(expected_cgmap_sites)
        or observed_cgmap_sites != expected_cgmap_sites
        or any(
            len(row) != 8 or row[3:5] != ["CG", "CG"]
            for row in cgmap_rows
        )
    ):
        raise SystemExit(
            "example CGmap must contain all {:,} CG-only rows".format(
                len(expected_cgmap_sites)
            )
        )

    with gzip.open(
        examples / "test.asm.gz", "rt", encoding="ascii", newline=""
    ) as input_file:
        asm_rows = [line.rstrip("\n").split("\t") for line in input_file]
    if len(asm_rows) != 1001 or any(len(row) != 13 for row in asm_rows):
        raise SystemExit("example ASS must contain a header and 1000 data rows")
    asm_data_rows = asm_rows[1:]
    asm_targets = {(row[0], int(row[5])) for row in asm_data_rows}
    asm_linked_snvs = {
        (row[0], int(row[1]), row[2], row[3], row[4])
        for row in asm_data_rows
    }
    substitution_counts: dict[tuple[str, str], int] = {}
    for row in asm_data_rows:
        substitution_class = tuple(sorted((row[3], row[4])))
        substitution_counts[substitution_class] = (
            substitution_counts.get(substitution_class, 0) + 1
        )
    if (
        len(asm_targets) != 1000
        or asm_linked_snvs
        != {
            (chrom, pos, ref, ref, alt)
            for chrom, pos, ref, alt in input_het_snvs
        }
        or len(substitution_counts) != 6
        or max(substitution_counts.values()) - min(substitution_counts.values()) > 1
        or max(abs(int(row[1]) - int(row[5])) for row in asm_data_rows) > 99
    ):
        raise SystemExit(
            "example ASS must have 1000 unique targets, use all 227 linked "
            "heterozygous SNVs, balance the six substitution classes, and "
            "keep every target within 99 bp of its SNP"
        )

    with tempfile.TemporaryDirectory(prefix="bsreadsim-examples-") as temporary:
        output = Path(temporary)
        _run(
            core,
            "build",
            "variants",
            "--reference",
            str(examples / "test.fa"),
            "--output",
            str(output / "variants.vcf.gz"),
            "--vcf",
            str(examples / "test.vcf.gz"),
            "--seed-phase",
            "7",
        )
        _run(
            core,
            "build",
            "methdb",
            "--reference",
            str(examples / "test.fa"),
            "--output",
            str(output / "profile.methdb"),
            "--cgmap",
            str(examples / "test.CGmap.gz"),
            "--asm",
            str(examples / "test.asm.gz"),
            "--vcf",
            str(examples / "test.vcf.gz"),
            "--seed-phase",
            "7",
        )

        for name in ("variants.vcf.gz", "profile.methdb"):
            path = output / name
            if not path.is_file() or path.stat().st_size == 0:
                raise SystemExit(f"example-data command did not create {name}")

        with gzip.open(
            output / "variants.vcf.gz", "rt", encoding="ascii", newline=""
        ) as input_file:
            saved_rows = [
                line.rstrip("\n").split("\t")
                for line in input_file
                if line and line[0] != "#"
            ]
        if len(saved_rows) != 1555:
            raise SystemExit("saved truth did not skip exactly two complex records")
        saved_chr11_indels = {
            (int(row[1]), row[3], row[4])
            for row in saved_rows
            if row[0] == "chr11" and row[2].startswith("example-")
        }
        if saved_chr11_indels != expected_chr11_indels:
            raise SystemExit("saved truth lost a chr11 1-4 bp indel")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
