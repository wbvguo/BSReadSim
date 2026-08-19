# HTSIM ASM profile v1

Status: normative allele-linked methylation input boundary.

For the equivalent BED6+6 representation, use
[HTSIM ASM BED profile v1](asm-bed-profile-v1.md). The two ASM serializations
are mutually exclusive and share the same VCF/haplotype validation boundary.

This is a BSReadSim simulation profile, not the raw output of the CGmapTools
`asm` command. CGmapTools reports allele-specific sites or regions relative to
a heterozygous variant; those reports must be converted into the site-linked
profile below before simulation. The format uses fourteen columns and replaces
an implicit fragment-wide mutation flag with an exact typed VCF link.

## File contract

The input is plain text or gzip, verified by SHA-256 from the same opened file
snapshot used for parsing. Empty lines and lines whose first byte is `#` are
ignored. Every data row has exactly fourteen tab-separated fields:

| Column | Name | v1 requirement |
| ---: | --- | --- |
| 1 | `CHR` | exact FASTA contig name |
| 2 | `NUC` | `C` or `G` at the methylation target |
| 3 | `POS` | one-based `uint32` methylation-target position |
| 4 | `CONTEXT` | `CG`, `CHG`, or `CHH` |
| 5 | `DINUC` | `CA`, `CC`, `CG`, or `CT` in CGmap orientation |
| 6 | `METH` | finite probability in `[0,1]`, or lowercase `na` |
| 7 | `SNP_POS` | one-based `uint32` position of the linked VCF SNV |
| 8 | `REF` | one uppercase A/C/G/T reference allele |
| 9 | `ALT` | one different uppercase A/C/G/T alternate allele |
| 10 | `REF_METH` | required finite probability in `[0,1]` |
| 11 | `ALT_METH` | required finite probability in `[0,1]` |
| 12 | `FOLD_CHANGE` | finite decimal or lowercase `na`; provenance only |
| 13 | `P_VALUE` | finite probability in `[0,1]` or lowercase `na`; provenance only |
| 14 | `COMMENT` | non-empty caller/provenance text |

Rows follow FASTA contig order and strictly increasing `POS` within a contig;
one target position may occur only once. The core validates `NUC`, `CONTEXT`,
`DINUC`, linked `REF`, and both coordinates against the independently verified
reference before writing protocol output. `METH`, `FOLD_CHANGE`, `P_VALUE`, and
`COMMENT` are validated but do not replace the allele-specific probabilities.

The normalized spool record is 20 bytes: two zero-based `uint32` positions,
two binary32 probabilities, the existing typed context byte, one dinucleotide
base, and the linked REF/ALT bases. Per-contig row counts and coordinates are
`uint32`; total rows and spool offsets are `uint64`. No u96 representation is
used.

## VCF and haplotype ownership

Every row must resolve to exactly one typed, heterozygous SNV in the verified
VCF by `(CHR, SNP_POS, REF, ALT)`. Its frozen `HaplotypeMask` retains the design
meaning: value 1 carries ALT on haplotype 1, value 2 on haplotype 2, and value 3
on both. ASM v1 rejects mask 3 because no reference-allele haplotype remains.
Protocol haplotypes continue to use zero-based values; the conversion from the
mask happens only at this typed boundary.

The methylation target must exist at the same reference origin and context on
both complete haplotypes. The shared MethDB site is then split into two typed
sites: the haplotype carrying the linked REF allele receives `REF_METH` and
`reference_haplotype`; the one carrying ALT receives `ALT_METH` and
`alternate_haplotype`. Both receive source `ASM`. A variant-created, deleted,
inserted, or context-divergent target fails closed instead of being inferred
from `posidx` bit 0 or from whether a fragment covers any mutation.

When multiple level sources are configured, precedence is `ASM > CGmap >
Beta`. The linked VCF event determines allele ownership; source precedence
never changes haplotype phasing.

Fixed and variable WGBS insert spans use the same catalog-level linkage. The
fragment sampler first returns a nonzero two-bit eligibility mask, the core
selects protocol haplotype 0 or 1, and only then does projection select the
corresponding ASM site. Insert length and fragment overlap never alter REF/ALT
ownership.
