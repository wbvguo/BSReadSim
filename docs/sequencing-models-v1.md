# Python sequencing quality and error models v1

Status: released normative pure-Python model and pipeline component.

The paper defines quality and sequencing-error assignment after bisulfite
conversion. The advanced quality path draws the first five cycles from
empirical distributions, then follows separate R1/R2 quality transition
matrices. Once quality is known, a quality-specific base transition matrix
samples the final base call. This component freezes that boundary without
pickle, mutable RNG state, or filesystem I/O.

## Artifact envelope

Both artifacts are strict UTF-8 JSON no larger than 8 MiB. Duplicate object
keys, non-finite numbers, unknown fields, malformed dimensions, and zero-total
rows fail closed. The config declaration and content identifier are:

| Model | `format` | `version` | content `schema` |
| --- | --- | --- | --- |
| quality | `json` | `quality-markov-v1` | `quality-markov-v1` |
| error | `json` | `quality-confusion-v1` | `quality-confusion-v1` |

Each document has exactly `schema`, `quality_scores`, and `mates`. The quality
scores are a strictly increasing, nonempty subset of integer Phred values
`[0, 93]`. `mates` has exactly two entries in R1, R2 order even for a
single-end run, so one artifact has an invariant schema.

Counts, rather than normalized floating-point probabilities, are stored in the
artifact. Every count is `uint32`; a row total is positive and checked
`uint64`. Sampling is proportional to these exact counts. This preserves raw
empirical histograms and avoids platform-dependent normalization.

## Quality Markov document

Each mate object has exactly:

- `initial_counts`: exactly five rows, one per cycle 0--4, each with one count
  per quality state;
- `transition_counts`: one row per leading quality state, with one count per
  following quality state.

For a read shorter than five bases, only the required initial rows are sampled.
At cycle 5 and later, the row indexed by the preceding sampled quality state is
used. The RNG key uses stage `quality` and the protocol `contig_index`; the
entity is the global `fragment_ordinal`; the local index is
`(mate_index << 32) | read_offset`.

## Quality-specific error document

Each mate object contains only `base_transition_counts`. It has one 4x4 matrix
per declared quality score. Rows are source A/C/G/T and columns are final
A/C/G/T, both in protocol base order `0,1,2,3`. The diagonal therefore means no
sequencing error. A quality absent from the artifact fails closed. Protocol N
is retained as N because it has no empirical A/C/G/T confusion row.

The RNG key uses stage `sequencing-error` and the same `contig_index`, entity,
and packed local index as quality. The matrix is applied to post-conversion
bases, never to the reference or pre-conversion template.

## Width and ownership boundary

Quality values, base values, and model indices remain small typed values;
individual empirical counts and read offsets are `uint32`. Row totals, packed
mate/base local indices, and fragment ordinals are `uint64`. No u96 coordinate,
counter, identifier, genotype word, `rseq`, or `posidx` field is introduced.

The parser accepts already snapshotted bytes and performs no I/O. Before core
process launch, the pipeline reopens each prepared model, verifies its exact
size and SHA-256 while reading a stable regular file, and retains only the
decoded immutable model in memory. It rejects a quality-confusion model unless
its quality domain contains every state the selected quality policy can emit.
The runtime preserves `conversion -> quality -> sequencing-error`, and the
manifest records both prepared model identities.
