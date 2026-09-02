# Target-score sampling for targeted assays

Targeted bisulfite sequencing (TBS), whole-exome sequencing (WES), and targeted
sequencing (TS) sample eligible targets uniformly by default. Select
`--sampling score` to use BED column 5 as each target's relative expected
output. This example uses TBS:

```bash
bsreadsim run tbs \
  --reference reference.fa \
  --output runs/tbs-weighted \
  --reads 100000 \
  --sampling score \
  --targets targets.bed \
  --insert-mean 300 \
  --center-sd 50 \
  --mutation-rate 0 \
  --seed 42
```

## Generation

The BED6 file defines eligible targets. Target selection happens before the
actual fragment length and center displacement are drawn, so those generation
controls do not change a BED row's relative score mass.

## Sampling

In this mode every score must represent a uint32 value. Decimal forms such as
`1.0` and `1e3` are accepted when they represent an integer exactly. A score of
zero keeps the row valid but prevents it from being selected. The eligible
score total on each contig must also fit uint32.

Scores are relative and need not sum to one. Multiplying all positive scores
by the same constant leaves the distribution unchanged. At least one eligible
target must have a positive score.

## Interpretation

The score represents aggregate output weight, which may include capture,
library, amplification, bisulfite, and mapping effects. It should not be
interpreted automatically as a pure molecular capture probability. Target
coordinates and eligibility still follow the
[capture target BED contract](../formats.md#capture-target-bed).
