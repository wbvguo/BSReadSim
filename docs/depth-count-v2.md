# Target-region depth-to-count contract v2

Status: normative C++ fragment-count boundary for WGBS, RRBS, and TBS.

The C++ generator converts requested depth into biological fragments after
complete eligibility planning. Python never estimates this count.

## Denominator

- WGBS: union length of reference contigs with positive sampling mass.
- TBS: union length of normalized BED target intervals; overlaps count once.
- RRBS: union length of eligible restriction-fragment reference intervals;
  overlaps count once.

The denominator is target reference space, not emitted flank bases, candidate
multiplicity, haplotype-expanded length, or the sum of mate footprints. Variant
indels therefore do not alter it.

For read length `R` and `mate_count` one or two:

```text
fragments = ceil((binary64(target_reference_bases) * binary64(depth))
                 / binary64(R * mate_count))
```

The operation order is normative under round-to-nearest binary64. A zero
denominator, non-finite intermediate, or result above `UINT32_MAX` fails before
the protocol preamble. A positive depth over a nonempty target emits at least
one fragment.

This definition makes requested depth the nominal mean number of sequenced
bases divided by selected target-region bases. Reads may extend into flanks or
overlap; those emitted bases remain in the numerator and are not silently
removed.
