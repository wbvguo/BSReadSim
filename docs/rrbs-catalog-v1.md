# RRBS candidate catalog v1

This document freezes the C++ ownership boundary for the first RRBS generator.
It follows the paper's method: identify restriction sites on haplotypes, form
all possible fragments in the configured length range, and sample candidates
using uniform or profile-derived probabilities. The Python component never
parses cut motifs or reconstructs this catalog.

## In-memory ownership and optional BED exchange

The runtime does not invoke the historical standalone `rrcut` program. For
each contig, C++ finds cut boundaries and constructs a `CandidateCatalog` (or a
haplotype-aware `DiploidCandidateCatalog`). Uniform simulation can sample that
catalog directly without BED I/O.

For an external probability model, the same catalog can be serialized and
scored through [RRBS candidate BED exchange v1](rrbs-candidate-bed-v1.md).
Python never reconstructs candidates. On import, C++ regenerates the catalog,
matches every fixed row, and uses only the external score. Thus BED is an
exchange table, not a second candidate implementation or an authoritative
source of fragment geometry.

## Cut declarations

A cut declaration is uppercase `LEFT|RIGHT`; `RIGHT` is non-empty and the
combined motif is at most 1024 bases. The `|` is a boundary, not a base. For
example, `C|CGG` recognizes `CCGG` and cuts at `motif_start + 1`.

`A`, `C`, `G`, and `T` match exactly. Motif `N` matches any concrete haplotype
base but does not match an unknown `N`. Matches may overlap. The v1 declaration
is searched exactly on each forward haplotype;
non-palindromic reverse-complement behavior must be requested with another
explicit declaration until a later enzyme schema freezes asymmetric cuts.

Cut coordinates are 0-based physical haplotype boundaries in
`[0, haplotype_length]`. Matches from
different declarations at the same boundary are coalesced while retaining a
recognition count.

## Candidate construction

For each contig, every ordered pair of distinct cut boundaries `(left, right)`
is a candidate when:

- `insert_min <= right - left <= insert_max`;
- R1 and, for paired-end data, R2 each satisfy the configured maximum ambiguous
  base fraction; and
- the read length fits every candidate, as required by the config schema.

The restriction-site count includes both endpoint sites and every internal
site, including multiple recognition declarations coalesced at one boundary.
Synthetic contig-edge cuts are not added: a candidate must be bounded by two
observed restriction sites. Event-free candidates are ordered by `(left,
right)`. Once a contig has typed events, haplotype-0 candidates precede
haplotype-1 candidates and each subcatalog retains `(left, right)` order.
`insert_mean` and `insert_stddev` do not alter an enzyme-bounded RRBS catalog;
they remain in the shared config schema for WGBS/TBS sampling and future model
metadata. RRBS v1 uses only the inclusive minimum/maximum length range.

When a contig has no typed events, every candidate carries the
haplotype-availability mask `both` (bits `1|2`, numeric value 3), preserving the
event-free stream. With VCF or de novo events, motif discovery, physical length,
and mate ambiguity are evaluated independently on the two constructed
haplotypes. A physical candidate then carries exactly bit 1 or bit 2. The mask
is availability metadata, never a zero-based protocol haplotype; conversion to
protocol 0/1 occurs only after validation.

An indel can make the physical template length differ from its reference details
envelope. Complete insertions and deletions are retained as typed events. A cut
strictly inside an inserted ALT string is omitted because the protocol identifies
the insertion as one indivisible event; it is never silently truncated.

## Uniform sampling

Uniform fairness is defined over physical haplotype copies. A mask-1 or mask-2
row contributes allocation mass one and a mask-3 row contributes mass two.
This prevents event-free mask-3 contigs from receiving half the per-copy mass
of contigs represented by separate haplotype rows.

Uniform sampling draws a physical-copy-equivalent catalog rank with replacement using RNG stage
`fragment`, the contig-index-derived key, the per-contig candidate ordinal as
`entity_ordinal`, and local index 1. Candidate ordinals advance by emitted
count across chunks. Therefore chunk size does not affect candidate ranks, and
the uniform catalog reports zero skipped fragments.

Profile scores and the exact optional BED workflow are specified in
[rrbs-candidate-bed-v1.md](rrbs-candidate-bed-v1.md).

All contig coordinates, candidate counts, and per-contig requested counts are
checked `uint32_t`. Global ordinals and aggregate counts are `uint64_t`; no
wider integer is part of this contract.
