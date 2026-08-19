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
fields fail closed. For GC count `g`, fixed fragment length `F`, and `B` bins,
the bin is calculated exactly with integer round-half-up:

```text
bin = round_half_up(g * (B - 1) / F)
```

## Global calibration

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

## Sampling and reproducibility

For each per-contig output ordinal `entity` and retry `a`, the core:

- derives the `fragment` key from `(master_seed, contig_index)`;
- draws a valid-start rank at local index `1 + 2*a`;
- draws an independent acceptance variate at `2 + 2*a`;
- accepts exactly when the variate is below `a_i` for that start's bin.

Retries do not advance the output ordinal. The next chunk resumes at the next
output ordinal, so chunk size and worker scheduling cannot change the result.
Each rejected proposal increments protocol `skipped_fragment_count`; 100,000
attempts without an acceptance fail closed.

## Released capability boundary

The contract supports WGBS with one fixed insert length and the unmodified
reference sequence. VCF input, de novo mutation, ASM, and variable inserts are
rejected. This boundary prevents an approximate reference histogram from
being presented as an exact target guarantee. Those combinations require a
future target calibration over their actual proposal domains.
