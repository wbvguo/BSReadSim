# ASM profile

Use `--asm PATH` with a VCF to define allele-specific methylation. This is a
BSReadSim input format, not raw CGmapTools `asm` output. It is mutually
exclusive with [ASM BED](asm-bed-profile.md).

## Fourteen-column format

Input may be plain text or gzip-compressed. Empty lines and `#` comments are
ignored. Every row has exactly fourteen tab-separated fields:

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `CHR` | exact FASTA contig name |
| 2 | `NUC` | `C` or `G` at the methylation target |
| 3 | `POS` | one-based target position |
| 4 | `CONTEXT` | `CG`, `CHG`, or `CHH` |
| 5 | `DINUC` | `CA`, `CC`, `CG`, or `CT` |
| 6 | `METH` | probability in `[0,1]` or lowercase `na` |
| 7 | `SNP_POS` | one-based linked VCF SNV position |
| 8 | `REF` | one uppercase A/C/G/T base |
| 9 | `ALT` | one different uppercase A/C/G/T base |
| 10 | `REF_METH` | probability in `[0,1]` |
| 11 | `ALT_METH` | probability in `[0,1]` |
| 12 | `FOLD_CHANGE` | finite decimal or lowercase `na` |
| 13 | `P_VALUE` | probability in `[0,1]` or lowercase `na` |
| 14 | `COMMENT` | non-empty provenance text |

Rows follow FASTA order and strictly increasing target positions. Target and
linked-variant coordinates, bases, and contexts are checked against the
reference.

## VCF relationship

Every row must resolve to exactly one heterozygous VCF SNV by
`(CHR, SNP_POS, REF, ALT)`. The haplotype carrying REF receives `REF_METH`; the
haplotype carrying ALT receives `ALT_METH`. The target must retain the same
reference origin and context on both haplotypes. Variant-created, deleted,
inserted, or context-divergent targets are rejected.

When multiple methylation sources are present, precedence is
`ASM > CGmap/bedMethyl > Beta`.
