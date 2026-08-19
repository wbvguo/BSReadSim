# TBS target catalog v1

This document freezes the initial targeted-bisulfite fragment boundary between
the C++ generator and Python postprocessor.  C++ exclusively owns BED parsing,
target validation, fragment placement, sampling, and the protocol
`capture_strand`.  Python never reparses BED coordinates.

## Verified BED6 input

The input is a plain or gzip-compressed BED6 text snapshot whose raw file bytes
must match the prepared SHA-256 digest.  Coordinates are zero-based,
half-open, unsigned 32-bit values and must satisfy
`0 <= start < end <= contig_length`.  Contig names must exist in the verified
reference catalog.

Rows contain exactly six tab-separated fields:

1. reference contig name;
2. interval start;
3. interval end;
4. non-empty target name without control bytes;
5. finite non-negative numeric score;
6. strand `+`, `-`, or `.`.

Empty lines, `#` comments, and UCSC `track`/`browser` metadata lines are
ignored.  LF and CRLF are accepted; bare CR and lines longer than 1 MiB are
rejected.  Targets are canonically sorted by reference order, coordinate,
strand, name, and score.  Two rows with the same contig, interval, and strand
are duplicates and fail closed. The score does not affect uniform sampling;
its optional target-level `output_weight` interpretation is defined by
[tbs-target-score-v1.md](tbs-target-score-v1.md).

Strand maps directly to the wire enum: `+` is `forward`, `-` is `reverse`, and
`.` is `unknown`.  Python uses that typed value to select the bisulfite
conversion orientation.

## Center and candidate rule

The TBS baseline accepts one fixed insert length `F`. A target `[start,end)`
has the integer center

```
center = start + (end - start) / 2
```

using truncating unsigned division. For each selected haplotype, that reference
center is projected to the physical boundary immediately before an insertion
at the same anchor. With `fragment_center_stddev=0`, the fragment starts at
`haplotype_center - floor(F/2)` and spans exactly `F` haplotype bases. Therefore
the target center is at template offset `floor(F/2)` for both odd and even `F`.

With positive center standard deviation `sigma`, each selected target draws
`Z ~ Normal(0,1)` using the versioned `box-muller-normal-v1` component. The
continuous displacement `sigma*Z` is truncated toward zero to a signed 64-bit
integer, then added to the target center. No
coordinate wider than `uint32_t` is emitted or retained.

A proposed fragment is rejected when it crosses its selected haplotype boundary
or when either emitted mate exceeds the configured ambiguous `N` threshold.
Bases in an unsequenced insert interior do not affect this test. At `sigma=0`,
a BED target is eligible when its fixed fragment is valid on at least one
haplotype. At `sigma>0`, it is eligible when at least one haplotype retains the
target center and has at least one sequenceable fixed-length start;
displacement attempts can move an edge target inward. Distinct overlapping
targets remain distinct choices even when they resolve to identical fragment
coordinates.

Uniform coverage is target-first. Every eligible BED row contributes exactly
one target choice, independent of whether one or two haplotypes can materialize
it. Target selection uses the frozen `fragment` RNG stage, the
contig-index-derived key, the per-contig output candidate ordinal as entity,
and local index 1. After target selection, a one-bit eligibility mask fixes the
only available haplotype; mask `both` selects haplotype 0 or 1 with probability
one half using the `haplotype` RNG stage, the global fragment ordinal, and
local index 0. Haplotype availability never changes the target's allocation
mass.

At positive `sigma`, normal attempt zero uses `fragment` local index 2, attempt
one uses 3, and so on. Retries keep the same selected target and haplotype, and
the first valid proposal wins. Every rejected proposal increments protocol
`skipped_fragment_count`; 100,000 failed proposals for one output candidate
fail closed. Sampling is independent of chunk size. Event-free targets retain
`HaplotypeMask::both` until the shared core selects a haplotype. With VCF or de
novo events, the diploid catalog computes the two-bit eligibility mask, selects
conditionally within it, and emits a physical candidate carrying exactly mask
bit 1 or 2. A center deleted from one haplotype therefore changes only the
conditional haplotype choice, not the BED target's sampling mass. Complete
indels are projected back into typed reference truth; a physical boundary
inside inserted ALT bases is retried or omitted rather than truncating its
event.

Target-score coverage replaces only the uniform target-selection step. It uses
the same RNG address and retains the same displacement and retry rules; see
[tbs-target-score-v1.md](tbs-target-score-v1.md).

All contig-local coordinates, candidate counts, and target ordinals use
`uint32_t`.  Global target counts, sample ordinals, and fragment ordinals use
`uint64_t`.  No wider integer representation is part of this contract.
