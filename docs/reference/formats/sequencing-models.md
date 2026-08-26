# Sequencing quality and error models

Use `--quality-model PATH` and `--error-model PATH` to supply empirical JSON
models. Quality is assigned after optional bisulfite conversion; the error model then
samples the final base call conditional on quality.

## Common JSON structure

Both files are strict UTF-8 JSON no larger than 8 MiB. Duplicate keys,
non-finite numbers, unknown fields, malformed dimensions, and zero-total rows
are rejected.

| Model | File format | Content `schema` |
| --- | --- | --- |
| quality | JSON | `quality-markov` |
| error | JSON | `quality-confusion` |

Each document has exactly `schema`, `quality_scores`, and `mates`.
`quality_scores` is a strictly increasing, nonempty list of integer Phred
values from 0 through 93. `mates` always contains exactly two entries in R1,
R2 order, including for single-end runs.

Models store non-negative integer counts rather than normalized probabilities.
Every row must have positive total mass; sampling is proportional to its
counts.

## Quality Markov model

Each mate object contains:

- `initial_counts`: exactly five rows, one per cycle 0–4, with one count per
  quality state; and
- `transition_counts`: one row per leading quality state, with one count per
  following quality state.

For reads shorter than five bases, only the required initial rows are used. At
cycle 5 and later, the row for the preceding sampled quality state is used.

## Quality-specific error model

Each mate object contains only `base_transition_counts`: one 4×4 matrix for
every declared quality score. Rows are source A/C/G/T and columns are final
A/C/G/T. Diagonal entries represent unchanged bases. N remains N because it
has no A/C/G/T confusion row.

The error model must contain every quality value the selected quality policy
can emit. Both model file identities are recorded in the run manifest.
