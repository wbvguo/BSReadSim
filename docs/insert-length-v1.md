# Addressed insert-length contract v1

Status: normative C++ fragment-generation component.

Insert length is owned by C++. Version 1 preserves the historical simulator's
clamped-normal semantics while replacing mutable random state with the shared
counter-addressed RNG. For parameters `minimum`, `mean`, `maximum`, and
`standard_deviation`, the result is:

```text
Z = standard_normal(fragment_key(contig), candidate_ordinal, local_index=0)
insert = clamp(mean + trunc_toward_zero(standard_deviation * Z),
               minimum, maximum)
```

Mathematically equivalent comparisons against the lower and upper normalized
bounds occur before floating-point-to-integer conversion. This keeps extreme
finite standard deviations defined without an overflowing cast. Zero standard
deviation returns `mean`; equal minimum/mean/maximum returns that fixed length.

The constraints are `0 < minimum <= mean <= maximum` with a finite,
non-negative standard deviation. Parameters and the returned insert length are
checked `uint32`. The per-contig candidate ordinal and Philox key are `uint64`.
No u96/u128 coordinate, length, or identifier is introduced.

The `fragment` RNG stage reserves `local_index=0` for insert length and
`local_index=1` for start selection. Global fragment ordinal is not used here:
changing allocation or worker/chunk boundaries therefore cannot silently move
the insert-length draw to another address.
