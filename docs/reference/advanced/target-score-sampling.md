# TBS target-score sampling

TBS normally gives every eligible target equal sampling mass. Select
`--sampling score` when column 5 of the BED6 input should instead define each
target's relative expected output.

```bash
bsreadsim run tbs \
  --reference reference.fa \
  --output runs/tbs-weighted \
  --reads 100000 \
  --targets targets.bed \
  --sampling score \
  --insert-mean 300 \
  --insert-sd 0 \
  --mutation-rate 0 \
  --seed 42
```

In this mode every score must be an exact integer from `0` through
`4,294,967,295`. Decimal forms such as `1.0` and `1e3` are accepted when they
represent an integer exactly. A score of zero keeps the row valid but prevents
it from being selected. The eligible score total on each contig must also fit
that range.

Scores are relative and need not sum to one. Multiplying all positive scores
by the same constant leaves the distribution unchanged. At least one eligible
target must have a positive score.

The score represents aggregate output weight, which may include capture,
library, amplification, bisulfite, and mapping effects. It should not be
interpreted automatically as a pure molecular capture probability. Target
coordinates and eligibility still follow the [TBS target
format](../formats/tbs-catalog.md).
