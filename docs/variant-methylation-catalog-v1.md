# Variant-aware diploid methylation catalog v1

Status: normative Beta-baseline MethDB contract after VCF haplotype projection.

The C++ core discovers methylatable sites on each complete haplotype, not on a
fragment slice. A streaming haplotype cursor applies typed SNV/insertion/deletion
events and feeds a five-base sliding window. It therefore classifies C/G bases
correctly when a causal variant lies just outside a fragment, when an insertion
creates a site, or when a deletion joins formerly separated flanks.

The implementation does not materialize whole-contig `reference_positions` or
`base_variant_indices` arrays. It walks each haplotype in `O(L + V)` time, where `L`
is contig length and `V` is event count. Classified site streams are merged by
base origin:

- a mapped origin is a `uint32` reference position;
- an inserted origin is `(uint32 event_ordinal, uint2 insertion_offset)`;
- identical origin/context pairs on both haplotypes are stored once;
- differing or one-sided sites are stored only in the corresponding haplotype
  overlay.

For compact lookup, reference origin id is the zero-extended position. Inserted
origin id sets bit 63 and stores `(event_ordinal << 2) | offset` in the low 34
bits. This lookup key is distinct from the RNG entity in
[methylation-site-entity-v1.md](methylation-site-entity-v1.md). Both fit in one
`uint64_t`; neither uses u96.

## Probability and allele ownership

When a mapped base and context equal the original reference, its Beta draw uses
tag-0 reference identity and is bit-identical to the reference-only catalog.
When both haplotypes share the same variant-created context, they share one
tag-1 mapped or tag-4 insertion probability. A haplotype-specific changed site
uses tag 2/3 or 5/6. Shape selection always follows the final haplotype context.

Protocol allele is decided by the diploid comparison:

- identical sites on both haplotypes are `shared`;
- an unchanged reference site absent or different on the other haplotype is
  `reference_haplotype` while retaining its tag-0 probability;
- a changed or inserted haplotype-specific site is `alternate_haplotype`.

Projection into a fragment matches sites by typed base origin, assigns
consecutive `site_index` values in template order, and uses `reference_pos=-1`
for inserted sites. The projection carries its contig and haplotype identity;
the catalog validates both along with insertion ALT order and event phase. The
same catalog probability is reused by every overlapping fragment.

## CGmap overlay

When a verified CGmap or bedMethyl profile is present, the core first builds
the deterministic Beta catalog above and then applies
[CGmap profile v1](cgmap-profile-v1.md) or
[bedMethyl profile v1](bed-methyl-profile-v1.md). A defined profile value can
replace a shared or `reference_haplotype` site only
when its mapped origin and final context equal the reference record. It changes
the protocol source to `CGMAP`. An `na` row, an `alternate_haplotype` context,
or an inserted origin retains its addressed Beta value and `BETA` source.

## ASM overlay

After any CGmap-compatible overlay, the core applies
[HTSIM ASM profile v1](asm-profile-v1.md) or
[ASM BED profile v1](asm-bed-profile-v1.md). Each row must link to one exact
heterozygous typed VCF SNV. The existing design-deck `HaplotypeMask` determines
which zero-based protocol haplotype carries REF and ALT; no packed `geno_int`
or fragment-wide `posidx` flag is consulted.

An ASM target is eligible only while it is a shared mapped origin with the same
reference-equivalent context on both complete haplotypes. The overlay removes
that shared site and creates one site in each haplotype catalog. The
reference-allele haplotype receives `REF_METH`, allele
`reference_haplotype`, and source `ASM`; the alternate-allele haplotype receives
`ALT_METH`, allele `alternate_haplotype`, and source `ASM`. Thus `ASM >
CGmap/bedMethyl > Beta`, while phasing remains owned exclusively by the VCF
catalog. Missing,
homozygous, non-SNV, allele-mismatched, or context-divergent links fail closed.
