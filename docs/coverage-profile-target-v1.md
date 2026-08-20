# WGBS target GC distribution v1

Status: normative component contract `wgbs-gc-target-v1`.

This component makes the profile describe the requested output distribution,
not a raw rejection probability. If line `i` contains `p_i`, the expected
fraction of emitted fragments in GC bin `i` is `p_i` after calibration against
the supplied reference.

## Artifact boundary

The normalized config records format `tsv`, version `wgbs-gc-target-v1`, and
the SHA-256 of the exact raw bytes. The C++ core verifies the digest before
decoding strict plain text or gzip.

Every physical line contains exactly one finite probability in `[0,1]`. A file
with `B` numbers defines `B` bins, requires `B >= 2`, and the probabilities
must sum to one within absolute tolerance `1e-9`:

```text
0.10
0.25
0.45
0.20
```

Empty lines, comments, explicit bin indices, surrounding whitespace, and extra
fields fail closed. For GC count `g`, that candidate's fragment length `F`, and `B` bins,
the bin is calculated exactly with integer round-half-up:

```text
bin = round_half_up(g * (B - 1) / F)
```

## Fixed-insert exact calibration

Eligibility is identical to uniform fixed-insert WGBS: the insert must fit and
each emitted mate must satisfy the configured ambiguous-base threshold. Let
`N_ci` be the number of eligible starts in contig `c` and bin `i`,
`N_i = sum_c N_ci`, and `p_i` the target output probability.

The core computes

```text
r_i = p_i / N_i
a_i = r_i / max_j(r_j)
W_c = sum_i N_ci * a_i
```

`a_i` is the calibrated rejection-sampling acceptance probability. `W_c` is
the contig allocation weight. Hamilton largest-remainder allocation over the
real-valued `W_c` values fixes each contig's integer fragment count. Combining
that allocation with within-contig rejection yields global output mass
proportional to `N_i*a_i = p_i/max_j(r_j)`, hence the requested `p_i` apart
from finite sampling error and integer contig allocation.

A positive `p_i` with `N_i = 0` is impossible and fails before the protocol
preamble. A zero target bin receives zero acceptance. For depth-based runs,
only contigs with positive calibrated allocation mass contribute to the
effective reference length.

## Variable-insert approximate calibration

Variable inserts keep the addressed clamped-normal proposal from
[insert-length-v1.md](insert-length-v1.md). The core does not replace it with a
fixed length. Instead, it scans eligible starts once at `insert_mean` and uses
those counts as a deterministic plug-in estimate of the mixed-length proposal:

```text
N_ci* = eligible starts in contig c and bin i at insert_mean
N_i*  = sum_c N_ci*
a_i   = (p_i / N_i*) / max_j(p_j / N_j*)
W_c*  = sum_i N_ci* * a_i
```

Hamilton allocation uses `W_c*`. During generation, every proposal first draws
its actual insert length, then its start, and maps the complete actual fragment
to a GC bin. The proposal is retained with probability `a_i` for that bin.

If a positive input bin has no eligible `insert_mean` start, approximate mode
projects the target onto the proxy's reachable support: those bins receive zero
acceptance and the remaining target mass is implicitly renormalized. If no
positive target mass remains, generation fails before the protocol preamble.
The audit reports the dropped input mass. Fixed-insert exact mode continues to
reject any positive unreachable bin rather than projecting it.

This is exact when each contig's variable-insert GC opportunity distribution is
the same as its `insert_mean` proxy. Otherwise the output is an approximation;
the error normally grows as the insert distribution becomes wider or GC varies
rapidly at the insert-length scale. The implementation therefore makes no exact
target claim for variable inserts. Fixed-insert behavior and its exact
calibration are unchanged.

GC rejection can also reweight insert lengths. Conditional on validity, the
emitted insert distribution is proportional to the configured insert proposal
times the expected GC acceptance at that length. Consequently
`minimum/mean/maximum/standard_deviation` define the proposal distribution, not
an independent guarantee on the final histogram when target-GC coverage is
enabled. Accuracy and insert drift must be checked on the intended reference
and fitted profile.

## Sampling and reproducibility

For each per-contig output ordinal `entity` and retry `a`, the core:

- derives the `fragment` key from `(master_seed, contig_index)`;
- draws a valid-start rank at local index `1 + 2*a`;
- draws an independent acceptance variate at `2 + 2*a`;
- accepts exactly when the variate is below `a_i` for that start's bin.

For variable inserts, each per-contig candidate ordinal uses local index `0`
for insert length, `1` for start, and `2` for GC acceptance. An invalid mate or
GC rejection advances the candidate ordinal; the next chunk resumes at the
returned ordinal. Thus fixed and variable paths are both chunk-independent.
Each rejected proposal increments protocol `skipped_fragment_count`; 100,000
attempts for one requested fragment without an acceptance fails closed.

## Released capability boundary

The contract supports fixed or variable insert WGBS on the unmodified reference
sequence. Fixed insert output has the exact expectation described above;
variable insert output uses the explicitly approximate mean-insert plug-in.
VCF input, de novo mutation, and ASM remain rejected because their haplotype
proposal domain is not represented by the reference-only calibration.
