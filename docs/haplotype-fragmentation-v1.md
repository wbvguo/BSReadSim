# Haplotype-first fragmentation v1

This contract freezes the shared boundary used when RRBS or TBS receives a VCF
or deterministic de novo mutation catalog. The order is biological and
architectural: C++ constructs both haplotypes and the diploid MethDB first,
then the selected technology discovers or places physical fragments on those
haplotypes. Python receives only typed fragment records and never reconstructs
this coordinate layer.

## Width and mask boundary

Contig length, reference and haplotype boundary, physical template length,
candidate count, and per-contig requested count are checked `uint32_t` values.
Global fragment ordinals, skipped counts, and aggregate counts are `uint64_t`.
No wider integer is introduced.

Availability uses two bits: bit 0 (numeric 1) denotes zero-based haplotype 0,
bit 1 (numeric 2) denotes haplotype 1, and numeric 3 denotes both. A mask is
validated before conversion to the separate zero-based wire field. It is never
shifted into a coordinate or overloaded with another flag.

## Coordinate layout

For one selected haplotype, the compact layout retains variant boundary
exceptions and an ambiguity bit vector. It does not retain a second whole
reference-position array. Reference runs map linearly. An insertion produces
two distinct physical boundaries at one reference anchor: before all ALT bases
and after all ALT bases. A deletion collapses one physical boundary whose left
fragment ends at the deletion start and whose right fragment starts at the
deletion end.

A physical haplotype slice maps back to one reference truth envelope and uses
the existing typed SNV, insertion, and deletion events. Physical insert length
is therefore independent of reference-envelope width. Insertions are
indivisible wire events: a candidate boundary strictly inside ALT bases
is unavailable. The generator omits a fixed candidate or retries a displaced
candidate; it never truncates the event or silently substitutes reference
bases.

BED reference boundaries map immediately before an insertion at the same
anchor. A BED center strictly inside a deletion has no coordinate on that
haplotype and removes only that haplotype choice.

## Provider composition

- RRBS materializes one haplotype sequence at a time for restriction-motif
  discovery. Candidate length and mate ambiguity are evaluated in physical
  haplotype coordinates. Candidates on a variant-bearing contig carry one mask
  bit and retain their complete typed truth projection.
- TBS keeps BED6 as verified reference-anchored input and computes each
  target's two-bit haplotype eligibility. Sampling is target-first: one BED row
  contributes one uniform choice or one aggregate target output weight, then a
  conditional half/half draw selects between two eligible haplotypes. Fixed
  physical insert length, displacement, mate ambiguity, and capture-strand
  rules are applied on that selected haplotype. The target weight is not
  multiplied by haplotype popcount.
- A contig with no typed events uses the original reference catalog and the
  `both` mask, preserving event-free protocol bytes.

VCF and de novo catalogs are mutually exclusive, and either requires
`update_variant_boundaries=true`. Candidate planning, payload preflight, and
input validation complete before the protocol preamble. Chunk size does not
enter any RNG address or candidate identity. Payload preflight is physical-
window local: for each haplotype it sweeps the maximum event-record bytes that
can coexist inside `insert_max` (RRBS) or the fixed insert (TBS). Variants that
cannot share one fragment are not summed into an artificial whole-contig
limit; the sweep uses `uint32_t` physical boundaries and a `uint64_t` byte
total.
