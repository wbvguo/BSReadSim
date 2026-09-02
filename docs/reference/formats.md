# File formats

This page defines the exact contracts for structured files used by BSReadSim,
including inputs and reusable VCF, MethDB, and reduced representation
bisulfite sequencing (RRBS) candidate artifacts. For input files, the
command-line option selects the parser; the input filename suffix does not.
For example, a `.bed` file may be MethBED, bedMethyl, ASM BED, an RRBS
candidate set, or a capture-target file. Output commands still enforce the
suffixes documented in the [CLI reference](cli.md).

For rules about combining inputs, see
[CLI option combination rules](cli.md#option-combination-rules).
For how those inputs affect a simulation, see
[Customize](../simulation/customize.md). Read outputs are documented under
[Outputs](../outputs/index.md).

<div class="reference-jump" markdown>

**Jump to:** [genome and variants](#genome-and-variant-inputs) ·
[methylation profiles](#methylation-profile-inputs) ·
[ASM](#allele-specific-methylation-inputs) ·
[fragment sampling](#fragment-sampling-inputs) ·
[sequencing models](#sequencing-model-inputs)

</div>

## Common rules

- Contig names are case-sensitive and must match the first
  whitespace-delimited token in the FASTA header.
- VCF and CGmap positions are one-based. BED-derived formats, including
  MethBG, use zero-based, half-open intervals.
- Coordinate-bearing inputs must use the same reference build. Matching
  contig names alone do not make coordinates from different builds valid.
- Each VCF, methylation profile, or ASM file may cover a subset of FASTA
  contigs. Omitted contigs are valid. Depending on the input type, they
  contribute no predefined variants, use generated methylation fallback, or
  receive no ASM overrides. Contig blocks that are present must retain their
  relative FASTA order.
- Row-oriented formats use tabs between fields unless a section explicitly
  says that the file contains one value per line. Parsers reject malformed or
  ambiguous rows rather than guessing from the filename suffix.
- VCF, methylation and ASM text profiles, target-GC profiles, capture-target
  BED, and RRBS candidate BED use a 1 MiB decoded-line limit. LF and CRLF line
  endings are accepted; a bare carriage return is rejected.

Compression is part of the selected format contract:

| Format | Accepted storage |
| --- | --- |
| Reference FASTA | Plain text or gzip-compressed text |
| VCF | Plain text or gzip-compressed text |
| MethBG, MethBED, bedMethyl, and CGmap | Plain text or gzip-compressed text |
| ASM and ASM BED | Plain text or gzip-compressed text |
| Target-GC profile and capture-target BED | Plain text or gzip-compressed text |
| RRBS candidate BED | Plain text only |
| Quality and error models | Uncompressed JSON only |
| MethDB | Opaque, internally compressed binary; do not wrap it in gzip |

Plain versus gzip text is detected from the file content, not from a `.gz`
suffix.

`bsreadsim validate` directly checks reference FASTA, VCF, the four text
methylation profiles, ASM, and ASM BED without creating simulation output. It
does not accept MethDB, target-GC profiles, RRBS candidate BED, capture-target
BED, or sequencing-model JSON; those formats are checked when the applicable
`run` or `export` command consumes them. Add `--json` for a machine-readable
coverage and row-count summary, or `--strict` to reject VCF records that would
otherwise be skipped.

## Genome and variants { #genome-and-variant-inputs }

### Reference FASTA { #fasta }

Use `--reference PATH` to supply the reference sequence against which every
coordinate-bearing input is validated. See NCBI's
[nucleotide FASTA format](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/)
for the external format. Plain and gzip-compressed FASTA are accepted, with
the BSReadSim-specific validation rules below.

```text { .no-copy }
>chr1 optional description
ACGTCGATCGATCGATCG
>chr2
TTNACCGGTA
```

The contig names above are `chr1` and `chr2`. Names must be non-empty, valid
UTF-8, unique, and no longer than 1 MiB. Text after the first space or tab is
an ignored description. Except for a tab separator, a header cannot contain
ASCII control characters or DEL. Sequence may contain only A, C, G, T, and N,
in either case. Every contig must contain at least one base. Wrapped sequence
lines are accepted; empty lines and sequence before the first header are
rejected. A reference used for simulation or artifact building must have each
contig length representable as uint32.

### VCF { #vcf }

Use `--vcf PATH` for predefined diploid variants. The external format is
defined by the GA4GH-maintained [VCF 4.2](https://samtools.github.io/hts-specs/VCFv4.2.pdf)
and [VCF 4.3](https://samtools.github.io/hts-specs/VCFv4.3.pdf)
specifications. BSReadSim accepts a restricted textual subset of those
versions, either plain or gzip-compressed; BCF is not accepted.

```text { .no-copy }
##fileformat=VCFv4.3
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	sample
chr1	101	rs1	C	T	.	PASS	.	GT	0|1
chr1	250	.	A	AG	.	PASS	.	GT	1|0
```

The file must contain exactly one `##fileformat=VCFv4.2` or
`##fileformat=VCFv4.3` declaration before one standard ten-column header.
After that header, every physical line must be a data row; blank lines and
additional header or comment lines are rejected. A header-only VCF is valid
and represents an empty variant set.

BSReadSim accepts:

- exactly one named sample and ten tab-separated columns;
- rows in FASTA contig order and nondecreasing one-based `POS`;
- one ALT allele and a diploid `GT` containing only `0` and `1`; `FORMAT` and
  the sample must have matching field counts and contain exactly one `GT`;
- uppercase A, C, G, or T in the normalized event, although a shared indel
  anchor may be `N`; and
- an SNV or a pure insertion or deletion of at most four bases after
  normalization.

MNPs, complex replacements, and pure indels longer than four bases are
syntax-checked and participate in input-order validation, but they are skipped
and do not appear in saved variant truth. They must use valid alleles and
genotypes and have an in-range REF interval, but they are skipped before exact
REF matching and retained-event conflict checks.

Among retained non-reference events, duplicate or overlapping normalized
events are rejected. Distinct, non-overlapping normalized events may share a
source `POS`. Their exact REF bases and coordinates are checked against the
FASTA. Missing genotypes, multiallelic records, and symbolic alleles are always
rejected.

A `0/0` row is syntax-checked and participates in input-order validation, but
creates no variant. It is discarded before event normalization and REF-allele
matching against the FASTA. FILTER and INFO values do not silently remove
rows. Empty or `.` IDs receive deterministic generated IDs; repeated IDs are
disambiguated in saved variant truth.

Phased `0|1` and `1|0` assignments are preserved, while `1|1` applies ALT
to both haplotypes. `--seed-phase` assigns unphased heterozygous variants.
Saved variant truth is a normalized, phased, BGZF-compressed VCF. Reuse it
with `--vcf`; see [Outputs](../outputs/index.md#saved-truth-artifacts) for the
saved artifact.

## Methylation profiles { #methylation-profile-inputs }

Choose the parser that matches the actual input:

| Profile source | Select with | Purpose |
| --- | --- | --- |
| Site-level methylation values | `--cgmap`, `--bedmethyl`, `--methbed`, or `--methbg` | Use values at listed sites, or use them as context pools with `--pool-meth`; otherwise unlisted sites use fallback generation |
| Prepared snapshot | `--methdb` | Reuse a complete prepared methylation profile and its embedded prepared variant set |

See
[methylation profiles](../simulation/customize.md#predefined-methylation) for
fallback and pooling behavior.

All four text-profile parsers accept plain or gzip-compressed input and ignore
empty lines and lines beginning with `#`. MethBG, MethBED, and bedMethyl also
ignore UCSC directives beginning with `track ` or `browser `. A CGmap header
is ignored only when it begins with `#`; an unprefixed header is rejected.
Imported probabilities are stored as uint16 UNORM values, using
`round(probability * 65535)`. Each file must contain at least one data row.

### MethBG { #methbg }

Use `--methbg PATH` for BSReadSim's minimal, methylation-specific bedGraph
table. Every data row has exactly four tab-separated fields.

```text { .no-copy }
#chrom	chromStart	chromEnd	probability
chr1	100	101	0.8234
chr1	205	206	0.14
```

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `chrom` | Exact FASTA contig name |
| 2 | `chromStart` | Zero-based target position |
| 3 | `chromEnd` | Exactly `chromStart + 1` |
| 4 | `probability` | Methylation level in `[0,1]`, used as the site's methylation probability |

Rows follow FASTA contig order and strictly increasing positions within each
contig. Every interval must resolve to an eligible cytosine context; base,
strand, and CG, CHG, or CHH context are derived from the reference.

### MethBED { #methbed }

Use `--methbed PATH` for BSReadSim's BED-compatible methylation format.

```text { .no-copy }
#chrom	chromStart	chromEnd	name	score	strand	meth_count	total_count	base	context
chr1	100	101	.	823	+	82	100	C	CG
```

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `chrom` | Exact FASTA contig name |
| 2 | `chromStart` | Zero-based target position |
| 3 | `chromEnd` | Exactly `chromStart + 1` |
| 4 | `name` | Non-empty site name, or `.` |
| 5 | `score` | Unsigned decimal integer in `[0,1000]`; probability is `score / 1000` |
| 6 | `strand` | `+`, `-`, or `.` |
| 7 | `meth_count` | Optional uint32 methylated count, or `.` |
| 8 | `total_count` | Optional uint32 total count, or `.` |
| 9 | `base` | Optional `C`, `G`, or `.` |
| 10 | `context` | Optional `CG`, `CHG`, `CHH`, or `.` |

Every row must use the same supported width: 6, 8, 9, or 10 fields. Fields 7
and 8 must therefore appear together; when defined, they satisfy
`meth_count <= total_count`. Counts are metadata and do not override `score`.
Producers encode a probability as `floor(probability * 1000 + 0.5)`.

Rows follow FASTA order and strictly increasing positions. Coordinates,
strand, base, and context are checked against the reference whenever supplied.

### bedMethyl { #bedmethyl }

Use `--bedmethyl PATH` for UCSC/ENCODE bedMethyl. BSReadSim accepts BED9+2 and
BED9+9, but one file cannot mix the two widths. See the [UCSC bedMethyl
specification](https://genome.ucsc.edu/goldenPath/help/bedMethyl.html).

```text { .no-copy }
chr1	100	101	m6C	82	+	100	101	0	100	82
```

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `chrom` | Exact FASTA contig name |
| 2 | `chromStart` | Zero-based target position |
| 3 | `chromEnd` | Exactly `chromStart + 1` |
| 4 | `name` | Non-empty label |
| 5 | `score` | uint32 |
| 6 | `strand` | `+`, `-`, or `.` |
| 7 | `thickStart` | Exactly `chromStart` |
| 8 | `thickEnd` | Exactly `chromEnd` |
| 9 | `itemRgb` | `0` or an `R,G,B` triple with components in `[0,255]` |
| 10 | `validCoverage` | uint32 |
| 11 | `percentModified` | Finite decimal in `[0,100]`; probability is `percentModified / 100` |
| 12–18 | BED9+9 counts | Seven uint32 values: `N_mod`, `N_canonical`, `N_other_mod`, `N_delete`, `N_fail`, `N_diff`, and `N_nocall` |

For BED9+9, modified, canonical, and other-modification counts must sum to
`validCoverage`. Only `percentModified` determines the imported methylation
probability; `score`, coverage, and count fields are validated metadata. Rows
follow FASTA order and strictly increasing positions; coordinates, strand,
and cytosine context are checked against the reference.

### CGmap { #cgmap }

The external CGmap format is documented by
[CGmapTools](https://cgmaptools.github.io/cgmaptools_documentation/file-formats.html#cgmap-format).
Use `--cgmap PATH` for its eight-column representation. Every data row has
exactly eight tab-separated fields.

```text { .no-copy }
#CHR	NUC	POS	CONTEXT	DINUC	METH	MC	NC
chr1	C	101	CG	CG	0.82	82	100
chr1	G	102	CG	CG	na	0	0
```

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `CHR` | Exact FASTA contig name |
| 2 | `NUC` | `C` or `G` |
| 3 | `POS` | Positive one-based coordinate |
| 4 | `CONTEXT` | `CG`, `CHG`, or `CHH` |
| 5 | `DINUC` | `CA`, `CC`, `CG`, or `CT` in cytosine orientation |
| 6 | `METH` | Methylation level in `[0,1]`, or lowercase `na` for fallback generation |
| 7 | `MC` | uint32 methylated count |
| 8 | `NC` | uint32 total count satisfying `MC <= NC` |

Rows follow FASTA contig order and strictly increasing positions. Coordinates,
strand, context, and dinucleotide are checked against the complete reference.
`METH` is the authoritative imported probability; `MC` and `NC` are validated
metadata, and their ratio is not compared with `METH`.

### MethDB { #methdb }

MethDB v2 is BSReadSim's immutable snapshot of a prepared methylation profile
and its embedded prepared variant set. It preserves generated or imported
probabilities, ASM overrides, contexts, sources, alleles, and normalized
variants.

**Save:** Create a snapshot with `build methdb`, `--save-methdb`, or
`--save-truth`.

**Reuse:** Load the original snapshot with `--methdb` and the same reference.
Reloading uses the stored profile and variants rather than regenerating them.
Explicit VCF, baseline-profile, ASM, and pooling inputs, as well as nonzero de
novo mutation, are rejected; profile-generation settings do not modify the
snapshot. Build-time mutation, phasing, and methylation seeds are not reload
requirements. See
[CLI option combination rules](cli.md#option-combination-rules).

**Inspect:** Use `bsreadsim export methdb` for a human-readable extended BED
projection. Its `methdb-bed-v2` layout is not the `--methbed` input format. The
export is for inspection and cannot replace the original snapshot or be loaded
with `--methbed`.

MethDB is an opaque binary artifact and must not be edited or concatenated.
The file begins with the `methdb` magic and version byte `2`; version 2 is
the only accepted representation.

??? info "MethDB v2 binary details"

    Each contig has an independently compressed baseline section. Diploid
    contigs add prepared events, reference overlays, and insertion overlays.
    ASM probabilities occupy a separate section.

    | Logical part | Physical section | Purpose |
    | --- | --- | --- |
    | Reference baseline | `BASELINE` | Sorted reference sites and default probabilities |
    | Variant layer | `EVENTS` | Normalized events and haplotype masks |
    | Variant layer | `REF_OVERLAY` | Reference sites changed or removed by events |
    | Variant layer | `INS_OVERLAY` | Methylation sites created on inserted bases |
    | ASM layer | `ASM` | Linked allele-specific probability overrides |

    Sections store uncompressed length and SHA-256 identity and use independent
    zlib frames. Probabilities use uint16 UNORM as
    `round(probability * 65535)`, with maximum quantization error
    `1/131070`. At runtime only the active contig is materialized.

## Allele-specific methylation { #allele-specific-methylation-inputs }

Choose either ASM or ASM BED. Each retained row links one
methylation site to one heterozygous SNV and provides allele-specific
methylation probabilities. A VCF can optionally be provided to preserve its
complete variant set and existing phasing. See
[Customize](../simulation/customize.md#allele-specific-methylation) for how
the values are applied.

### ASM { #asm }

The ASM format is defined by
[CGmapTools](https://cgmaptools.github.io/cgmaptools_documentation/methylation-analysis.html#asm).
Use `--asm PATH` for site-level `cgmaptools asm -m ass` output. Region-level
`-m asr` output is not accepted. Empty lines, lines beginning with `#`, and the
exact standard header shown below are ignored. Every other line must have
exactly thirteen tab-separated fields.

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `Chr` | Exact FASTA contig name |
| 2 | `SNP_Pos` | One-based linked-SNV position |
| 3 | `Ref` | One uppercase A/C/G/T base |
| 4 | `Allele1` | First heterozygous allele in CGmapTools genotype order |
| 5 | `Allele2` | Second, different heterozygous allele |
| 6 | `C_Pos` | One-based methylation-target position |
| 7 | `Allele1_linked_C` | `methylated-unmethylated` uint32 counts for `Allele1` |
| 8 | `Allele2_linked_C` | `methylated-unmethylated` uint32 counts for `Allele2` |
| 9 | `Allele1_linked_C_met` | Methylation probability in `[0,1]` for `Allele1` |
| 10 | `Allele2_linked_C_met` | Methylation probability in `[0,1]` for `Allele2` |
| 11 | `pvalue` | Probability in `[0,1]` |
| 12 | `fdr` | Adjusted probability in `[0,1]` |
| 13 | `ASM` | Uppercase `TRUE` or `FALSE` |

```text { .no-copy }
Chr	SNP_Pos	Ref	Allele1	Allele2	C_Pos	Allele1_linked_C	Allele2_linked_C	Allele1_linked_C_met	Allele2_linked_C_met	pvalue	fdr	ASM
chr1	150	A	T	A	101	8-2	2-8	0.80	0.20	0.001	0.005	TRUE
```

Each support field must contain at least one read in total. Only `TRUE` rows
become overrides; a file with no `TRUE` rows is rejected. `FALSE` rows still
undergo field, range, coordinate-boundary, and contig-order checks, but they
are discarded before reference/context and haplotype-link validation.
Exactly one of `Allele1` and `Allele2` must equal `Ref`; `Allele1` is not
assumed to be REF. BSReadSim maps both values to REF and ALT. Evidence counts,
p-value, and FDR do not alter the supplied probabilities.

### ASM BED { #asm-bed }

Use `--asm-bed PATH` for BSReadSim's editable BED representation. Empty lines,
lines beginning with `#`, and UCSC directives beginning with `track ` or
`browser ` are ignored. Use BED6+6 throughout the file, or BED6+10 when
retaining analysis evidence; the two widths cannot be mixed. At least one
data row is required.

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `chrom` | Exact FASTA contig name |
| 2 | `chromStart` | Zero-based target position |
| 3 | `chromEnd` | Exactly `chromStart + 1` |
| 4 | `name` | Non-empty provenance label |
| 5 | `score` | Integer in `[0,1000]`; display only |
| 6 | `strand` | `+`, `-`, or `.` |
| 7 | `linkedStart` | Zero-based linked-SNV position |
| 8 | `linkedEnd` | Exactly `linkedStart + 1` |
| 9 | `REF` | One uppercase A/C/G/T base |
| 10 | `ALT` | One different uppercase A/C/G/T base |
| 11 | `REF_METH` | REF-allele methylation probability in `[0,1]` |
| 12 | `ALT_METH` | ALT-allele methylation probability in `[0,1]` |
| 13 | `REF_SUPPORT` | BED6+10 only: `methylated-unmethylated` uint32 counts with positive total |
| 14 | `ALT_SUPPORT` | BED6+10 only: `methylated-unmethylated` uint32 counts with positive total |
| 15 | `P_VALUE` | BED6+10 only: finite probability in `[0,1]` |
| 16 | `Q_VALUE` | BED6+10 only: finite probability in `[0,1]` |

```text { .no-copy }
chr1	100	101	asm-1	600	+	149	150	A	T	0.2	0.8	2-8	8-2	0.001	0.005
```

BED6+10 includes all four evidence fields as a complete group; none accepts
`.` or `na`. They do not override `REF_METH` or `ALT_METH`. For consistent
browser display, producers may set `score` to
`round(abs(REF_METH - ALT_METH) * 1000)`; BSReadSim does not use it in the
simulation.

### ASM validation

The validation below applies to `TRUE` CGmapTools rows and every ASM BED data
row. Without a VCF, retained links define the minimal heterozygous SNV set;
repeated identical links are collapsed and conflicting alleles are rejected.
With a VCF, every link must resolve to exactly one matching heterozygous SNV.
The linked SNV and target must be on the same contig, and the REF allele,
coordinates, strand, and target context must agree with the reference. Target
positions must be unique and remain the same CG, CHG, or CHH site on both
haplotypes. Duplicate targets, inserted or deleted targets, and
context-divergent targets are rejected. Contig blocks follow FASTA order.

## Fragment generation and sampling { #fragment-sampling-inputs }

### Target-GC profile { #target-gc-profile }

Use `--sampling gc --gc-profile PATH` with whole-genome bisulfite sequencing
(WGBS) or whole-genome sequencing (WGS) to load BSReadSim's target fragment-GC
profile. The plain or gzip-compressed file contains one finite probability in
`[0,1]` per physical line, at least two lines, and must sum to one within
`1e-9`. Blank lines, comments, extra fields, and surrounding whitespace are
rejected.

```text { .no-copy }
0.05
0.15
0.60
0.15
0.05
```

The first line is the lowest GC bin and the last is the highest. For GC count
`g`, fragment length `F`, and `B` bins, BSReadSim assigns:

```text
bin = floor(g * (B - 1) / F + 0.5)
```

For fixed-insert runs, positive probability in a bin with no eligible fragment
is rejected. Reference-only variable-insert runs instead drop unreachable
positive probability and renormalize the reachable bins; they fail if no
positive-probability bin is reachable. Supported variant and insert-length
combinations are listed under
[CLI option combination rules](cli.md#option-combination-rules).

### RRBS candidate BED { #rrbs-candidate-bed }

Use `build rrbs` to export the exact candidate domain for reduced
representation bisulfite sequencing (RRBS), let an external model change only
its `score` column, then load it with `--sampling score --rrbs-candidates
PATH`. See the
[RRBS tutorial](../simulation/tutorials.md#rrbs) for the commands.

This BSReadSim-defined exchange file is uncompressed plain text with exactly
ten tab-separated fields per data row. Empty lines and lines beginning with
`#` are ignored; at least one candidate row is required.

```text { .no-copy }
#chrom start end candidate_id score strand haplotype_mask template_length gc_count restriction_site_count
chrR	1	5	chrR:1-5~0	1	.	1	4	2	2
chrR	1	5	chrR:1-5~1	1	.	2	4	2	2
```

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `chrom` | Exact FASTA contig name |
| 2–3 | `start`, `end` | uint32 zero-based, half-open reference envelope satisfying `start <= end <= contig length` |
| 4 | `candidate_id` | Non-empty deterministic fragment identity without whitespace or control bytes |
| 5 | `score` | `.` or a finite non-negative relative weight |
| 6 | `strand` | `.` |
| 7 | `haplotype_mask` | `1`, `2`, or `3` |
| 8 | `template_length` | Positive uint32 physical haplotype fragment length |
| 9 | `gc_count` | uint32 number of C and G bases, no greater than `template_length` |
| 10 | `restriction_site_count` | uint32 count of recognized endpoint and internal cut sites |

Scores are initialized to `1`, are relative, and need not sum to one. `.` is
accepted only outside score-weighted sampling. Input order may change, but
missing, extra, duplicate, or modified candidate identities are rejected. The
build and run must use matching domain-defining settings.

### Capture target BED { #capture-target-bed }

Use `--targets PATH` with targeted bisulfite sequencing (TBS), whole-exome
sequencing (WES), or targeted sequencing (TS). Input uses the six-column
layout of [UCSC BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1), with
zero-based, half-open coordinates and the BSReadSim-specific score rules
below.

| Column | Name | Requirement |
| ---: | --- | --- |
| 1 | `chrom` | Exact FASTA contig name |
| 2 | `start` | Target start satisfying `0 <= start < end` |
| 3 | `end` | Target end no greater than the contig length |
| 4 | `name` | Non-empty target name without control bytes |
| 5 | `score` | Finite non-negative relative weight |
| 6 | `strand` | `+`, `-`, or `.` |

```text { .no-copy }
chr1	1000	1100	target-1	1	.
chr1	2500	2550	target-2	3	+
```

Empty lines, `#` comments, and UCSC `track` or `browser` lines are
ignored, but at least one target row is required. Rows may appear in any order;
BSReadSim groups them by FASTA contig and sorts them by coordinate. Duplicate
contig, interval, and strand combinations are rejected. Score is metadata
under uniform sampling. Under `--sampling score`, it must be uint32; see
[target-score sampling](advanced/target-score-sampling.md).

## Sequencing models { #sequencing-model-inputs }

Use `--quality-model PATH` and `--error-model PATH` for empirically estimated
models defined by BSReadSim's JSON schemas below. Both are uncompressed,
strict UTF-8 JSON: the raw file must contain between 1 byte and 8 MiB
inclusive. Duplicate keys, non-finite numbers, unknown fields, malformed
dimensions, and zero-total rows are rejected.

| Model | Content `schema` | Mate fields |
| --- | --- | --- |
| Quality | `quality-markov` | `initial_counts`, `transition_counts` |
| Error | `quality-confusion` | `base_transition_counts` |

Each document contains exactly `schema`, `quality_scores`, and `mates`.
`quality_scores` is a non-empty, strictly increasing list of integer Phred
values from 0 through 93. `mates` has exactly two entries in R1, R2 order,
including for single-end runs. Every count is uint32;
every row must have positive total mass, and sampling is proportional to the
counts.

### Quality Markov model

Each mate has `initial_counts` with exactly five rows and one count per
quality state, plus `transition_counts` with one row for each possible
preceding quality state and one count per following state. The first five
cycles use their corresponding initial row; later cycles use the row for the
previously sampled quality.

??? example "Minimal one-state quality model"

    ```json
    {
      "schema": "quality-markov",
      "quality_scores": [30],
      "mates": [
        {
          "initial_counts": [
            [1],
            [1],
            [1],
            [1],
            [1]
          ],
          "transition_counts": [[1]]
        },
        {
          "initial_counts": [
            [1],
            [1],
            [1],
            [1],
            [1]
          ],
          "transition_counts": [[1]]
        }
      ]
    }
    ```

### Quality-specific error model

Each mate contains `base_transition_counts`: one 4×4 matrix per declared
quality score. Rows are source A/C/G/T and columns are final A/C/G/T; diagonal
entries leave the base unchanged. N remains N. The error model must cover
every quality value emitted by the selected quality policy.

??? example "Minimal one-state error model"

    ```json
    {
      "schema": "quality-confusion",
      "quality_scores": [30],
      "mates": [
        {
          "base_transition_counts": [
            [
              [1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]
            ]
          ]
        },
        {
          "base_transition_counts": [
            [
              [1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]
            ]
          ]
        }
      ]
    }
    ```
