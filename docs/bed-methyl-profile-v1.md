# bedMethyl methylation profile v1

Status: normative alternate serialization of the position-specific MethDB
overlay.

Use this format through `bsreadsim run --bed-methyl PATH` or the paired native
options `--bed-methyl PATH --bed-methyl-sha256 HEX64`. It is mutually exclusive
with CGmap input. File suffixes never select a parser.

The accepted schema follows the UCSC/ENCODE bedMethyl BED9+2 contract and its
BED9+9 extension. The first three fields therefore use ordinary BED coordinates:
zero-based `chromStart` and half-open `chromEnd`. See the
[UCSC bedMethyl specification](https://www.genome.ucsc.edu/goldenPath/help/bedMethyl.html)
and the
[modkit column description](https://nanoporetech.github.io/modkit/intro_pileup.html#bedmethyl-column-descriptions).

## File contract

Input is a plain or gzip-compressed regular file whose raw bytes match the
declared SHA-256. Empty lines, `#` comments, and standard `track ` or `browser `
header lines are ignored. Every data row has either 11 fields (BED9+2) or 18
fields (BED9+9), and one file cannot mix the two widths:

| Column | Name | v1 requirement |
| ---: | --- | --- |
| 1 | `chrom` | exact FASTA contig name |
| 2 | `chromStart` | zero-based `uint32` target position |
| 3 | `chromEnd` | exactly `chromStart + 1` |
| 4 | `name` | non-empty modification/motif label; provenance only |
| 5 | `score` | `uint32`; provenance only |
| 6 | `strand` | `+`, `-`, or `.` |
| 7 | `thickStart` | exactly `chromStart` |
| 8 | `thickEnd` | exactly `chromEnd` |
| 9 | `itemRgb` | `0` or an `R,G,B` triple with components in `[0,255]` |
| 10 | `validCoverage` | `uint32`; provenance only |
| 11 | `percentModified` | finite decimal in `[0,100]` |
| 12-18 | extended counts | optional BED9+9 `uint32` counts |

For BED9+9, modified, canonical, and other-modification counts must sum exactly
to `validCoverage`. The remaining deletion, failed, different-base, and no-call
counts are syntax-checked but do not alter the simulation probability.

Rows follow FASTA contig order and strictly increasing `chromStart` within each
contig. Duplicate or returning targets fail. A row must select one reference C
or G whose complete CG/CHG/CHH context is resolvable. `+` requires C, `-`
requires G, and `.` asks the core to infer the orientation from FASTA. The core
derives context and dinucleotide from the verified reference rather than
trusting the display-oriented `name` field.

## MethDB semantics

`percentModified / 100` is the site's methylation probability. Missing genomic
positions retain their deterministic Beta value. Unlike CGmap, bedMethyl v1 has
no `na` value; omit a position when it should use the Beta fallback.

After reference validation, bedMethyl rows normalize to the same typed records
as [CGmap profile v1](cgmap-profile-v1.md). Position-specific overlay,
variant-aware reference equivalence, source precedence, bounded per-contig
spooling, and `cgmap_pool=true` therefore behave identically. For protocol
compatibility the normalized provenance remains `CGMAP` (or `POOLED_CGMAP`);
the manifest records the actual input format as `bedMethyl`.
