# HTSIM ASM BED profile v1

Status: normative BED serialization of the allele-linked ASM boundary.

Use this format through `bsreadsim run --asm-bed PATH` or the paired native
options `--asm-bed PATH --asm-bed-sha256 HEX64`. It is mutually exclusive with
the fourteen-column HTSIM ASM profile and still requires a verified VCF.

This is a project-defined BED6+6 schema, not an arbitrary BED3 file and not a
raw CGmapTools ASM report. The BED prefix makes target intervals compatible
with ordinary BED tooling; the six typed extension fields preserve the exact
SNV link and two allele-specific probabilities required by simulation.

## File contract

Input is plain text or gzip and is verified by SHA-256. Empty lines, `#`
comments, and standard `track ` or `browser ` lines are ignored. Each row has
exactly twelve tab-separated fields:

| Column | Name | v1 requirement |
| ---: | --- | --- |
| 1 | `chrom` | exact FASTA contig name for target and linked SNV |
| 2 | `chromStart` | zero-based target position |
| 3 | `chromEnd` | exactly `chromStart + 1` |
| 4 | `name` | non-empty provenance label |
| 5 | `score` | integer in `[0,1000]`; provenance only |
| 6 | `strand` | `+`, `-`, or `.` for the target C/G |
| 7 | `linkedStart` | zero-based linked VCF SNV position |
| 8 | `linkedEnd` | exactly `linkedStart + 1` |
| 9 | `REF` | one uppercase A/C/G/T reference allele |
| 10 | `ALT` | one different uppercase A/C/G/T alternate allele |
| 11 | `REF_METH` | required finite probability in `[0,1]` |
| 12 | `ALT_METH` | required finite probability in `[0,1]` |

Both intervals must be in range on the same contig. Rows follow FASTA order and
strictly increasing `chromStart`; target positions are unique. `+` requires a
reference C, `-` requires G, and `.` is resolved from FASTA. The complete
CG/CHG/CHH context, dinucleotide, and linked REF base are independently checked
against the verified reference before protocol output.

## VCF and MethDB semantics

The linked interval and REF/ALT alleles must resolve to exactly one typed,
heterozygous VCF SNV. Haplotype ownership, shared-target requirements, failure
on variant-created/deleted/context-divergent targets, and probability assignment
are exactly those in [HTSIM ASM profile v1](asm-profile-v1.md). The BED form
omits the TSV profile's total methylation, fold-change, p-value, and comment
fields because they are provenance-only and never determine the two simulated
allele probabilities.

Source precedence remains `ASM > CGmap/bedMethyl > Beta`. Equivalent ASM TSV
and ASM BED inputs must produce identical typed protocol streams.
