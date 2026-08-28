# RRBS candidate BED

RRBS can sample restriction fragments directly without an exchange file. Use
the candidate BED workflow only when an external model must attach relative
scores to the exact fragments discovered by BSReadSim.

## Workflow

Export candidates:

```bash
bsreadsim build rrbs \
  --reference reference.fa \
  --output candidates.bed \
  --cut-site 'C|CGG' \
  --mutation-rate 0
```

An external program may change only the `score` column and must preserve every
row and all other fields. Run with the resulting file:

```bash
bsreadsim run rrbs \
  --reference reference.fa \
  --output runs/rrbs \
  --reads 100000 \
  --cut-site 'C|CGG' \
  --rrbs-candidates candidates.scored.bed \
  --sampling score \
  --mutation-rate 0 \
  --seed 42
```

Reference, cut-site, fragment, ambiguity, VCF, mutation, and focused-seed
options must match the export command. When variants or mutations affect the
candidates, keep `--seed-mut` and `--seed-phase` identical. The run's `--seed`
controls later fragment sampling and is not accepted by `build rrbs`.

## Ten-column format

The file begins with a comment header and tab-separated rows:

```text
#chrom start end candidate_id score strand haplotype_mask template_length gc_count restriction_site_count
```

| Column | Meaning |
| --- | --- |
| `chrom` | exact reference contig name |
| `start`, `end` | zero-based, half-open reference envelope |
| `candidate_id` | deterministic fragment identity |
| `score` | `.` or a finite non-negative relative weight |
| `strand` | `.` for RRBS candidates |
| `haplotype_mask` | `1`, `2`, or `3` |
| `template_length` | physical haplotype fragment length |
| `gc_count` | number of C and G bases in the template |
| `restriction_site_count` | endpoint and internal recognized cut sites |

The exporter initializes scores to `1`. Scores are relative and need not sum
to one. `.` is accepted only when `--sampling score` is not selected. Input order
does not matter, but missing, extra, duplicate, or modified candidates are
rejected when BSReadSim regenerates and verifies the candidate set.
