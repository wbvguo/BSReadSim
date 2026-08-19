# Variable-insert WGBS sampling contract v1

Status: released normative C++ fragment-candidate component.

This reference-only proposal component supplies uniform-coverage WGBS.
Target-GC coverage is intentionally restricted to fixed inserts, so it does
not consume this proposal stream. Typed VCF and de novo catalogs use the sibling
[variable-haplotype-wgbs-v1.md](variable-haplotype-wgbs-v1.md) boundary.
Fixed-insert runs retain the existing valid-rank sampler and byte stream.

For each contig, allocation weight is the number of ambiguity-valid starts at
`insert_max`. This conservative gate guarantees that every allocated contig
can fit every configured insert length and preserves a checked `uint32`
per-contig weight.

Each sampling attempt uses one per-contig `uint64 candidate_ordinal`:

1. `fragment/local_index=0` selects the insert length according to
   [insert-length-v1.md](insert-length-v1.md).
2. If that insert fits, `fragment/local_index=1` selects a start uniformly from
   every fitting start position.
3. The attempt is accepted only when each emitted mate has at most
   `floor(max_ambiguous_fraction * read_length)` ambiguous bases. Bases in the
   unsequenced insert interior do not affect eligibility.

A rejected attempt increments both the candidate ordinal and the `uint64`
skipped count. The next chunk starts at the returned next candidate ordinal,
so chunk size cannot change any length, start, or rejection. Each requested
fragment has a fixed 100,000-attempt cap in production. Exhaustion and ordinal
overflow fail closed.

An integration that enables this component MUST validate the complete
maximum-span allocation domain before writing the protocol preamble. A runtime
attempt-cap exhaustion terminates the core stream and the Python transaction
without publishing outputs. Coordinates,
lengths, possible-start counts, allocation weights, and per-contig output
counts are `uint32`; candidate ordinals and skipped counts are `uint64`. No
u96/u128 coordinate or identifier is used.
