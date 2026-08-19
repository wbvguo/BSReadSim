# Variable-span haplotype WGBS contract v1

Status: released normative C++ candidate component.

This component extends [variable-wgbs-v1.md](variable-wgbs-v1.md) to a typed
VCF or de novo mutation catalog. Candidate coordinates and spans stay in
original-reference
coordinates and are checked `uint32` values. A typed VCF catalog may retain
the CGmap and ASM overlays defined by
[variant-methylation-catalog-v1.md](variant-methylation-catalog-v1.md).
Uniform coverage consumes accepted candidates directly. Target-GC coverage
does not accept VCF or de novo mutation input in its released contract.

Allocation uses the fixed-span two-bit index from
[variant-start-index-v1.md](variant-start-index-v1.md) at `insert_max`. This is
a conservative pre-header gate: each allocated contig can emit at least one
maximum-span fragment on at least one haplotype. Permanent eligibility still
uses exactly two adjacent bits per reference start. Bit 0 means protocol
haplotype 0/design `HaplotypeMask` 1, bit 1 means protocol haplotype 1/design
mask 2, both means mask 3, and zero is represented as no candidate rather than
as a `HaplotypeMask` value.

Each attempt uses one per-contig `uint64 candidate_ordinal`. The addressed
insert draw uses `fragment/local_index=0`; the uniform reference-start draw
uses `fragment/local_index=1`. The chosen half-open reference interval is then
projected independently through both complete typed haplotypes. A haplotype is
eligible only when neither boundary cuts an active deletion, the projected
template contains a complete read, and every emitted mate satisfies the N
threshold. Insertions and deletions may therefore change mate contents and the
two-bit eligibility mask without changing the sampled reference span.

The component exposes the same operation as a single addressed proposal:
`candidate_at(candidate_ordinal)` returns either one typed candidate or the
ineligible zero-bit state. Batch sampling is defined only as repeated calls to
that operation and resumes by the returned candidate ordinal.

An accepted candidate carries its nonzero `HaplotypeMask`; the existing
global-fragment-ordinal haplotype draw chooses only within that mask. Rejected
attempts increment the candidate ordinal and `uint64` skipped count. Chunk
boundaries cannot change output, and the production attempt cap is 100,000 per
requested fragment. Exhaustion terminates the core stream so Python can roll
back its private output transaction.

ASM REF/ALT ownership follows the linked VCF event's same mask; no fragment
span, packed `posidx` flag, or sampled haplotype number is allowed to redefine
the linkage.

Reference coordinates, spans, possible-start counts, allocation weights, and
per-contig output counts are `uint32`. Candidate/global ordinals, accumulated
skips, and size sums are `uint64`. No u96 coordinate, identifier, or packed
genotype representation is introduced.
