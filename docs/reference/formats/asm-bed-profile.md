# ASM BED profile

Use `--asm-bed PATH` with a VCF for the BED6+6 representation of
allele-specific methylation. It is mutually exclusive with the fourteen-column
[ASM profile](asm-profile.md).

Input may be plain text or gzip-compressed. Empty lines, `#` comments, and
`track` or `browser` lines are ignored. Every row has twelve tab-separated
fields:

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `chrom` | exact FASTA contig name |
| 2 | `chromStart` | zero-based target position |
| 3 | `chromEnd` | exactly `chromStart + 1` |
| 4 | `name` | non-empty provenance label |
| 5 | `score` | integer in `[0,1000]` |
| 6 | `strand` | `+`, `-`, or `.` |
| 7 | `linkedStart` | zero-based linked VCF SNV position |
| 8 | `linkedEnd` | exactly `linkedStart + 1` |
| 9 | `REF` | one uppercase A/C/G/T base |
| 10 | `ALT` | one different uppercase A/C/G/T base |
| 11 | `REF_METH` | probability in `[0,1]` |
| 12 | `ALT_METH` | probability in `[0,1]` |

Both intervals must lie on the same contig. Rows follow FASTA order and
strictly increasing target positions. Target orientation, context, linked REF,
and coordinates are checked against the reference.

The linked variant must be exactly one heterozygous VCF SNV. The REF
haplotype receives `REF_METH` and the ALT haplotype receives `ALT_METH`.
Source precedence is `ASM > CGmap/bedMethyl > Beta`.
