#!/usr/bin/env python3
"""Generate the deterministic WGBS/RRBS/TBS throughput corpus."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


CONTIG = "chrRrbsSynthetic"
REFERENCE_LENGTH = 10_000_000
FRAGMENT_LENGTH = 300
TARGET_START = 149
REFERENCE_SHA256 = "5d436a0de36d479aa65d751f3fc56435da77a0aec9832116915e3dd87fd79235"
TARGETS_SHA256 = "826f87e42598a6428b52cc1490a8c8902d196c43ac9dfafd5cf68eebf66d8e8d"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-directory", type=Path, required=True)
    arguments = parser.parse_args()

    output = arguments.output_directory.expanduser().resolve(strict=False)
    output.mkdir(parents=True, exist_ok=True)
    reference = output / "technology-synthetic-10m.fa"
    targets = output / "tbs-targets.bed"
    manifest = output / "corpus-manifest.json"
    for path in (reference, targets, manifest):
        if path.exists():
            raise SystemExit("refusing to overwrite {}".format(path))

    block = b"CCGG" + b"A" * (FRAGMENT_LENGTH - 4)
    repeats, remainder = divmod(REFERENCE_LENGTH, len(block))
    sequence = block * repeats + block[:remainder]
    reference.write_bytes(
        b">" + CONTIG.encode("ascii") + b"\n" + sequence + b"\n"
    )

    with targets.open("w", encoding="ascii", newline="\n") as destination:
        for start in range(TARGET_START, REFERENCE_LENGTH - 1, FRAGMENT_LENGTH):
            destination.write(
                "{}\t{}\t{}\ttarget-{}\t1\t.\n".format(
                    CONTIG,
                    start,
                    start + 1,
                    start,
                )
            )

    reference_sha256 = _sha256(reference)
    targets_sha256 = _sha256(targets)
    if reference_sha256 != REFERENCE_SHA256:
        raise RuntimeError("generated reference digest changed")
    if targets_sha256 != TARGETS_SHA256:
        raise RuntimeError("generated target digest changed")

    document = {
        "schema": "bsreadsim-technology-corpus-1",
        "contig": CONTIG,
        "reference": reference.name,
        "reference_length": REFERENCE_LENGTH,
        "reference_sha256": reference_sha256,
        "rrbs_cut_site": "C|CGG",
        "fragment_length": FRAGMENT_LENGTH,
        "targets": targets.name,
        "target_count": sum(
            1 for _ in range(TARGET_START, REFERENCE_LENGTH - 1, FRAGMENT_LENGTH)
        ),
        "targets_sha256": targets_sha256,
    }
    manifest.write_text(
        json.dumps(document, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
