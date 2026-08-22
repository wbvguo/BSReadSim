# WGBS target GC distribution v2

Status: normative component contract `wgbs-gc-target-v2`.

The supplied profile is the requested output distribution. BSReadSim does not
reinterpret it as a proposal or infer how the probabilities were fitted.

## Input

The verified strict TSV contains one finite probability in `[0,1]` per line,
at least two lines, summing to one within `1e-9`. Blank lines, comments, extra
fields, and whitespace fail closed. For GC count `g`, actual physical fragment
length `F`, and `B` bins:

```text
bin = round_half_up(g * (B - 1) / F)
```

## Exact fixed-insert calibration

Eligibility is evaluated after haplotype construction. A physical opportunity
is `(haplotype, reference_start)`, and GC is counted on its actual projected
haplotype fragment. Categories are `(haplotype, gc_bin)`. For target bin mass
`p_i`, each haplotype category receives target `p_i / 2`.

With opportunity count `N_hi`:

```text
r_hi = (p_i / 2) / N_hi
a_hi = r_hi / max(r)
W_c  = sum_hi N_chi * a_hi
```

Hamilton allocation over `W_c`, followed by calibrated rejection sampling,
gives the requested global GC expectation while balancing haplotypes. Positive
target mass with no physical support fails before output. VCF, fixed de novo
mutation, and ASM are supported because calibration owns the post-haplotype
proposal domain.

## Variable insert boundary

Reference-only variable inserts retain the v1 mean-insert proxy approximation
and report unreachable projected mass. Variable insert together with VCF or de
novo variants remains fail-closed: exact category opportunity depends jointly
on insert length, haplotype projection, validity, and GC, and has not yet been
implemented. This restriction concerns the sampler, not the scientific meaning
of the user profile.

## Reproducibility

Sampling uses addressed Philox domains. Retry ordinals, candidate ranks, and
acceptance draws are independent of chunking and worker count. GC rejection is
included in `skipped_fragment_count`.
