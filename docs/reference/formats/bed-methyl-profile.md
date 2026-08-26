# bedMethyl methylation profile

Use `--bed-methyl PATH` for a UCSC/ENCODE bedMethyl profile. It is mutually
exclusive with CGmap input. Input may be plain text or gzip-compressed.

BSReadSim accepts BED9+2 and BED9+9. See the [UCSC bedMethyl
specification](https://www.genome.ucsc.edu/goldenPath/help/bedMethyl.html) and
[modkit column description](https://nanoporetech.github.io/modkit/intro_pileup.html#bedmethyl-column-descriptions).

## Accepted fields

Empty lines, `#` comments, and `track` or `browser` lines are ignored. Every
data row has either 11 fields or 18 fields; one file cannot mix the two widths.

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `chrom` | exact FASTA contig name |
| 2 | `chromStart` | zero-based target position |
| 3 | `chromEnd` | exactly `chromStart + 1` |
| 4 | `name` | non-empty label |
| 5 | `score` | non-negative integer |
| 6 | `strand` | `+`, `-`, or `.` |
| 7 | `thickStart` | exactly `chromStart` |
| 8 | `thickEnd` | exactly `chromEnd` |
| 9 | `itemRgb` | `0` or an `R,G,B` triple |
| 10 | `validCoverage` | non-negative integer |
| 11 | `percentModified` | finite decimal in `[0,100]` |
| 12–18 | extended counts | optional BED9+9 non-negative counts |

For BED9+9, modified, canonical, and other-modification counts must sum to
`validCoverage`. Rows follow FASTA order and strictly increasing positions.
The selected C or G, strand, coordinate, and complete context are checked
against the reference.

`percentModified / 100` is the methylation probability. Omit a position when
it should use the generated Beta fallback. Variant-aware matching, source
precedence, and `--cgmap-pool` behavior match [CGmap](cgmap-profile.md).
