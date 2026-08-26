# CGmap methylation profile

Use `--cgmap PATH` for a position-specific CGmap profile. CGmap and bedMethyl
inputs are mutually exclusive.

## Accepted text format

Input may be plain text or gzip-compressed. Empty lines and lines beginning
with `#` are ignored. Every data row has exactly eight tab-separated fields:

```text
CHR  NUC  POS  CONTEXT  DINUC  METH  MC  NC
```

- `CHR` is an exact FASTA contig name. Rows follow FASTA contig order and
  strictly increasing positions within each contig.
- `NUC` is `C` or `G`; `POS` is a positive one-based coordinate.
- `CONTEXT` is `CG`, `CHG`, or `CHH`.
- `DINUC` is `CA`, `CC`, `CG`, or `CT` in cytosine orientation.
- `METH` is a finite decimal in `[0,1]` or lowercase `na`.
- `MC` and `NC` are decimal counts from `0` through `4,294,967,295` satisfying
  `MC <= NC`.

Coordinates, strand, context, and dinucleotide are checked against the full
reference. Unknown contigs, duplicate positions, unresolved boundary/N
contexts, and reference mismatches are rejected. CHG and CHH rows are still
validated when `--cpg-only` is selected, but are not used.

## Methylation behavior

A defined `METH` value replaces the generated Beta probability at the same
reference position and context. Missing positions and `na` use the Beta
fallback.

With a diploid VCF, a profile value applies only where that haplotype remains
reference-equivalent. Variant-created, inserted, or context-changed sites use
the fallback. Add `--cgmap-pool` only when values should form empirical
context distributions rather than remain position-specific; see [CGmap
context pooling](../advanced/cgmap-pooling.md).
