# Variant-aware WGBS start index v1

Status: normative pre-header eligibility and rank-sampling boundary for fixed
original-reference-span WGBS fragments.

Every possible reference start owns exactly two adjacent bits. Bit 0 means
protocol haplotype 0 is eligible (the design-deck `HaplotypeMask` value 1), bit
1 means protocol haplotype 1 is eligible (mask value 2), and both bits give
mask value 3. A zero pair is not a candidate. Thirty-two starts fit in one
`uint64_t`; the structure does not allocate a coordinate or integer per base.

Eligibility is computed on the final selected haplotype before any protocol
bytes are written. A reference interval is eligible for a haplotype only when:

- neither half-open boundary cuts through an active deletion;
- its projected template is at least `read_length` bases;
- the first `read_length` projected bases pass the N threshold; and
- in paired-end mode, the last `read_length` projected bases also pass it.

The implementation scans each haplotype in `O(L + V)` time. A sparse event
boundary cursor maps increasing reference boundaries into haplotype offsets.
One bit per haplotype base plus a rank prefix counts N in each emitted mate;
inserted/SNV ALT bases are A/C/G/T and therefore contribute zero N. A deletion
may pull a formerly interior N into a mate, and an insertion may displace an N
and make a start eligible; both effects are included. The temporary N-rank
index is discarded after each haplotype. Permanent availability uses only two
bits per possible start plus a sparse rank prefix.

Rank sampling is uniform over reference starts having at least one eligible
haplotype, not over `(start, haplotype)` pairs. Once a start is selected, its
mask constrains the existing global-fragment-ordinal haplotype draw. Thus a
region is not double-weighted merely because both haplotypes are viable. The
fragment RNG address remains `(stage=fragment, contig_index, candidate ordinal,
local=1)`. With no events, validity, rank mapping, and sampled starts are
exactly identical to the reference-only `ValidStartIndex`.

Reference coordinates and start/rank counts remain `uint32`. Haplotype offset
sums and RNG ordinals use `uint64`; no u96 representation is used.
