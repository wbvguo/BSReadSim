# Output policy

`bsreadsim run` has one output decision:

- Without `--bam`, publish R1, optional R2, and the reproducibility
  manifest.
- With `--bam`, publish one annotated BAM and the manifest. FASTQ
  sidecars are omitted because reads and qualities are recoverable with
  `samtools fastq`.

There is no production/debug mode, no `output.details` field, and no
per-fragment Details JSONL artifact. Full Details protocol columns are an internal
transport requirement only when annotated BAM is requested.

`zt` and `zr` are fixed on every BAM record. `zf` is emitted on every
record only when `--fragment-summary` is selected. `zx` is emitted on both
mates when `--fragment-realization` is selected; that option implies `zf` and
BAM output. The normative field layout is [bam-v3.md](bam-v3.md).

## FASTQ contract

FASTQ identifiers use
`@<contig>:<start>-<end>:<ordinal-hex>/<pair-number>`. Coordinates are
one-based and inclusive, the ordinal is variable-width lowercase hexadecimal,
and the third line is exactly `+`. See
[read-name-v2.md](read-name-v2.md).

## Transaction contract

All data artifacts are staged privately. Their record counts, sizes, and
SHA-256 digests are reconciled with the protocol trailer before the complete
manifest is published last. Existing destinations are never overwritten.

## Required gates

1. Default output has exactly R1, optional R2, and manifest roles.
2. BAM output has exactly one `bam` data role plus the manifest.
3. The protocol Details-column policy matches `--bam`.
4. BAM records remain SAM/BAM 1.6 compatible and usable by samtools.
5. Core and Python fragment, mate, site, base, and per-contig counts reconcile.
6. Worker count does not change fixed-input, fixed-seed output semantics.
