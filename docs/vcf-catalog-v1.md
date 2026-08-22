# VCF event catalog v1

Status: normative contract for predefined genetic variants.

The C++ core is the sole owner of VCF parsing, deterministic phasing, reference
validation, and haplotype construction. Python receives typed variant events
and, when annotated BAM is requested, their per-base provenance through the
Full Details protocol projection.

## Accepted VCF subset

Input is strict plain-text or gzip VCF 4.2/4.3, verified by SHA-256 through the
same opened regular-file snapshot used for parsing. Version 1 requires:

- exactly one named sample and ten tab-separated columns;
- rows strictly sorted by FASTA contig order and increasing one-based POS;
- one ALT allele, diploid `GT`, and alleles only `0` or `1`;
- uppercase A/C/G/T REF and ALT;
- SNVs or pure insertions/deletions of at most four bases after normalization.

Missing genotypes, multiallelic records, symbolic alleles, MNPs, complex
replacements, duplicate positions, and overlapping normalized events fail
closed. A `0/0` record is valid but creates no event. FILTER and INFO are
retained only as input bytes; v1 does not silently filter rows by those fields.

Normalization removes the maximal common prefix and then maximal remaining
common suffix. Coordinates are 0-based half-open in the typed catalog.
Insertions have `reference_start == reference_end` and empty REF; deletions have
empty ALT. Exact REF bases are checked against the independently verified
reference snapshot before the event catalog can be used.

## Haplotype assignment

The catalog retains the frozen two-bit `HaplotypeMask`: 1 is haplotype 1, 2 is
haplotype 2, and 3 is both. These values are availability/applicability bits,
not the protocol's zero-based sampled haplotype.

Phased `0|1` and `1|0` map directly. `1|1` maps to both. For an unphased
heterozygote, the core derives the `haplotype` key from `(catalog_seed,
contig_index)` and draws one Bernoulli value at `(entity=per-contig event
ordinal, local=1)`. True assigns ALT to haplotype
1; false assigns ALT to haplotype 2. Reference-only records do not consume an
event ordinal. The result is independent of chunking and worker scheduling.

All contig coordinates, per-contig event counts, and event ordinals checked at
the file boundary use `uint32_t`; total events and RNG entity ordinals use
`uint64_t`. No u96 value is used. The new catalog never recreates the
underspecified signed `geno_int` bit packing.

## Released generator checkpoint

The released end-to-end generator enables all three typed event kinds for
fixed or variable-insert WGBS with uniform coverage. It validates
the complete VCF, every REF
allele, and a sliding upper bound for every possible protocol payload before
emitting the protocol preamble. Its pre-header
[variant-start-index-v1](variant-start-index-v1.md) eligibility pass rejects a
haplotype/start combination when a boundary cuts a deletion, the projected
template is too short, or either emitted mate exceeds the N threshold. Uniform
fixed-insert sampling ranks only eligible reference starts and reports zero
skips. Variable spans use the addressed proposal stream. Fixed-insert target-GC coverage calibrates physical opportunities after haplotype construction. Variable-insert target-GC with variants remains gated.

Python carries typed `variants` through the internal Full Details transport when BAM is selected. BAM alignment, CIGAR, and summary flags retain deletion effects without a JSON sidecar. RRBS motif discovery and TBS targeting consume
the two constructed haplotypes before fragmentation, as frozen by
[haplotype-fragmentation-v1](haplotype-fragmentation-v1.md). A physical cut or
fragment boundary inside inserted ALT sequence is omitted because the wire
event is indivisible; this is a local candidate
restriction, not a technology-wide VCF restriction.
