# CLI parameters and defaults

This page is the complete reference for BSReadSim's public command-line
interface. It lists every option, its accepted value, its effective default,
and the commands on which it is available. BSReadSim uses explicit subcommands
for simulations and reusable artifact construction. The interface is
path-based: provide input paths and BSReadSim computes, verifies, and records
file identities internally. SHA-256 digests are not command-line inputs.

In the tables below, **required** means that the option must be supplied,
**none** means that no input or model is selected, and **off** means that a
flag is disabled unless it is present. Conditional defaults are called out
next to the option and explained immediately below the table.

## Synopsis

```text
bsreadsim run wgbs [OPTIONS]
bsreadsim run rrbs [OPTIONS]
bsreadsim run tbs  [OPTIONS]
bsreadsim run wgs  [OPTIONS]
bsreadsim run wes  [OPTIONS]
bsreadsim run ts   [OPTIONS]

bsreadsim build variants [OPTIONS]
bsreadsim build methdb   [OPTIONS]
bsreadsim build rrbs     [OPTIONS]

bsreadsim export test-fasta -o PATH
bsreadsim export methdb -i INPUT.methdb -o OUTPUT.bed.gz
bsreadsim export methdb -i INPUT.methdb -o OUTPUT.bed --no-compression
```

Use `bsreadsim --help`, `bsreadsim COMMAND --help`, or
`bsreadsim COMMAND SUBCOMMAND --help` to inspect the interface installed in
the current environment.

## Choose a command

| Command | Use it for | Primary result |
| --- | --- | --- |
| `run wgbs` | whole-genome fragment sampling | FASTQ or annotated BAM plus a manifest |
| `run rrbs` | restriction-fragment sampling | FASTQ or annotated BAM plus a manifest |
| `run tbs` | BED-targeted fragment sampling | FASTQ or annotated BAM plus a manifest |
| `run wgs` | ordinary whole-genome sequencing | FASTQ or annotated BAM plus a manifest |
| `run wes` | whole-exome capture from a target BED | FASTQ or annotated BAM plus a manifest |
| `run ts` | ordinary targeted sequencing from a target BED | FASTQ or annotated BAM plus a manifest |
| `build variants` | build a reusable prepared variant set | one BGZF-compressed `.vcf.gz` |
| `build methdb` | build a reusable prepared methylation profile | one MethDB file |
| `build rrbs` | export the exact RRBS candidate domain | one candidate BED |
| `export test-fasta` | copy the bundled synthetic reference | one FASTA file |
| `export methdb` | decode a MethDB methylation profile | one BED or BED.gz file |

`run` commands always require a technology subcommand. `build` commands do
not generate reads or a run manifest.

## Value conventions

- A **probability** is a finite number in `[0, 1]`.
- A **uint64 seed** is a decimal integer from `0` through
  `18446744073709551615`.
- A **flag** takes no value; its presence enables the behavior.
- Relative paths are resolved from the working directory. Manifests record
  resolved paths and input SHA-256 identities.
- Run output directories are created as needed. Existing directories are
  allowed, but an existing destination artifact is never overwritten.
- Standalone `build` outputs require an existing parent directory and a new,
  non-existing destination path.

### Global options

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| `-h` | `--help` | flag | show help and exit |
| — | `--version` | flag | print the installed version and exit |

## `run`: shared options

The general options in this section are shared by all six `run` commands.
Methylation and bisulfite-chemistry options are exposed only by `run wgbs`,
`run rrbs`, and `run tbs`; technology-specific options are documented below.

### Required inputs and dataset size

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| `-r` | `--reference` | FASTA path | required |
| `-o` | `--output` | directory path | required |
| `-n` | `--reads` | SE: integer `1..4294967295`; PE: even integer `2..8589934590` | required unless `--depth` is used |
| `-d` | `--depth` | finite number `> 0` | required unless `--reads` is used |

Choose exactly one of `--reads` and `--depth`.

`--reads` counts total read records across all output mates. In paired-end
mode the value must be even: `--reads 1000` emits 500 R1/R2 pairs. In
single-end mode, the same option emits 1000 R1 records.

`--depth D` follows the same read-count contract. It first resolves the raw
number of reads, rounds up to a complete fragment bundle, and then derives the
internal fragment count:

```text
raw_reads = effective_reference_bases * D / read_length
resolved_reads = emitted_mates * ceil(raw_reads / emitted_mates)
resolved_fragments = resolved_reads / emitted_mates
```

The effective region is technology-specific: eligible contigs for WGBS/WGS,
the union of eligible restriction-fragment envelopes for RRBS, and the union
of target intervals for TBS/WES/TS. The requested depth is retained in the
effective configuration; resolved read and fragment counts are retained in the
manifest summary.

### Random seeds

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| `-s` | `--seed` | uint64 | generated from OS entropy and recorded |
| — | `--seed-mut` | uint64 | `0` |
| — | `--seed-phase` | uint64 | `0` |
| — | `--seed-meth` | uint64 | `0`; bisulfite modes only |

The seeds address independent random domains:

| Seed | Controls |
| --- | --- |
| `--seed-mut` | de novo mutation generation |
| `--seed-phase` | deterministic assignment of unphased VCF heterozygotes |
| `--seed-meth` | generated methylation probabilities and empirical-pool draws in bisulfite modes |
| `--seed` | fragment geometry and selection, qualities, and sequencing errors; also methylation states and conversion in bisulfite modes |

Fix all relevant seeds for a byte-reproducible run. Keep the focused
biological seeds fixed while changing only `--seed` to draw a different read
sample from the same prepared genome and, for bisulfite modes, methylome.

### Genetic variation

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--vcf` | one-sample diploid VCF path | none |
| — | `--mutation-rate` | probability | conditional; see below |
| — | `--indel-fraction` | probability | `0.15` |
| — | `--indel-extension-probability` | probability | `0.15` |
| — | `--homozygous-only` | flag | off |

`--mutation-rate` is the total de novo mutation-event rate. Among generated
events, `--indel-fraction` selects the fraction that are insertions or
deletions; the remainder are SNVs. `--indel-extension-probability` controls
extension of indels up to the supported four-base limit.

The run default is `0.001`, except that it becomes `0` when `--gc-profile`,
`--rrbs-candidates`, or `--methdb` is supplied. `--vcf` selects predefined
variants and automatically records an internal mutation rate of `0`; it is
mutually exclusive with `--mutation-rate`. A MethDB already embeds the prepared
variant set, so it excludes `--vcf` and accepts only an explicitly supplied
`--mutation-rate 0`. For auditable experiments without a VCF, specify
`--mutation-rate` explicitly when the inferred default is not desired.

### Methylation inputs

These options are available only for WGBS, RRBS, and TBS. WGS, WES, and TS do
not construct or scan a methylome.

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--cgmap` | CGmap path | none |
| — | `--bed-methyl` | BED9+2 or BED9+9 path | none |
| — | `--methdb` | MethDB path | none |
| — | `--asm` | BSReadSim ASM path | none; requires `--vcf` |
| — | `--asm-bed` | ASM BED6+6 path | none; requires `--vcf` |

`--cgmap`, `--bed-methyl`, and `--methdb` are mutually exclusive. `--asm` and
`--asm-bed` are mutually exclusive. A loaded MethDB cannot be combined with a
CGmap, bedMethyl, ASM overlay, or `--cgmap-pool`.

### Methylation generation

These options are available only for WGBS, RRBS, and TBS.

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--beta-cg` | positive floats as `ALPHA,BETA` | `0.5,0.5` |
| — | `--beta-chg` | positive floats as `ALPHA,BETA` | `0.01,0.05` |
| — | `--beta-chh` | positive floats as `ALPHA,BETA` | `0.01,0.05` |
| — | `--cpg-only` | flag | off; retain CG, CHG, and CHH |
| — | `--cgmap-pool` | flag | off |
| — | `--methylation-model` | `bernoulli` or `bilstm` | `bernoulli` |
| — | `--no-update-variant-boundaries` | flag | off; update affected contexts |

Without a methylation input, probabilities are generated from the three Beta
distributions. Missing CGmap or bedMethyl positions also use these fallbacks.
`--cgmap-pool` requires `--cgmap` or `--bed-methyl` and samples values by
contig and context instead of preserving their input coordinates.

Variant-aware runs require boundary updates, so
`--no-update-variant-boundaries` is rejected with `--vcf` or a positive
mutation rate. `bilstm` is reserved but currently warns and executes the
released Bernoulli state model.

### Fragment geometry

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--single-end` | flag | paired-end |
| `-l` | `--read-length` | integer `1..10000` | `100` |
| — | `--insert-min` | positive integer | `100` |
| — | `--insert-mean` | positive integer | `400` |
| — | `--insert-max` | positive integer | `1000` |
| — | `--insert-sd` | non-negative float | `25`; `0` for TBS/WES/TS |
| — | `--max-ambiguous-fraction` | probability | `0.05` |

Set `--insert-mean N --insert-sd 0` to use one fixed insert length for WGBS,
WGS, TBS, WES, or TS. If min/max are omitted in that mode, both resolve to the
mean. All
normalized configurations satisfy:

```text
insert_min <= insert_mean <= insert_max
```

For fixed whole-genome or targeted sampling, each read must fit `insert_mean`.
A positive WGBS/WGS SD uses a clamped-normal distribution and each read must
fit `insert_min`. For RRBS,
restriction sites determine realized fragment lengths and `--insert-min` /
`--insert-max` define the retained size window; mean and SD do not resample
restriction fragments. TBS, WES, and TS use and default to `--insert-sd 0`.

`--max-ambiguous-fraction` is checked independently for each emitted mate.

### Chemistry, quality, and error

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--conversion-rate` | probability | `0.998`; bisulfite modes only |
| — | `--undirectional` | flag | off; directional OT/OB library; bisulfite modes only |
| `-q` | `--phred` | integer `0..93` | `40` |
| — | `--quality-model` | quality Markov JSON path | none |
| `-e` | `--error-rate` | probability | `0.005` |
| — | `--error-model` | quality-confusion JSON path | none |

`--conversion-rate` is the probability that an unmethylated convertible base
is converted. It and `--undirectional` apply only to WGBS, RRBS, and TBS.
Directional simulation samples independent original-top (OT/Watson) and
original-bottom (OB/Crick) fragments with equal probability. `--undirectional`
adds the complementary CTOT and CTOB orientations, sampling all four with
equal probability. Target-strand constraints retain their Watson/Crick source
while directionality controls whether the original or complementary molecule
is sequenced.
WGS, WES, and TS bypass both methylation realization and bisulfite conversion.
`--quality-model` replaces uniform `--phred`; the two cannot be combined.
`--error-model` replaces uniform `--error-rate`; those two also cannot be
combined. Model files are hashed automatically and their identities are
recorded in the manifest.

### Read output and simulation truth

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| `-p` | `--prefix` | `[A-Za-z0-9._-]+`, at most 128 characters | `sim` |
| `-f` | `--format` | `fastq`, `fastq.gz`, or `bam` | `fastq.gz` |
| — | `--gzip-level` | integer `0..9` | `6` |
| — | `--fragment-summary` | flag | off; BAM only |
| — | `--fragment-realization` | flag | off; BAM and bisulfite modes only, implies summary |
| — | `--save-methdb` | flag | off; bisulfite modes only |
| — | `--save-vcf` | flag | off |
| — | `--save-truth` | flag | off; enables all applicable truth flags |

`--format` selects exactly one read representation. BAM replaces FASTQ rather
than accompanying it. `--gzip-level` controls compressed FASTQ or BAM output;
it has no byte-level effect on uncompressed FASTQ. Compressed FASTQ is a
standards-compliant concatenation of deterministic gzip members, which lets
processing workers compress ordered batches in parallel.

The reusable biological state has two consistent names throughout the user
documentation: the **prepared variant set** and, in bisulfite modes, the
**prepared methylation profile**. The save flags publish these as simulation
truth artifacts with predictable names:

```text
OUTPUT/
├── PREFIX.R1.fastq.gz       # or .fastq, or PREFIX.bam
├── PREFIX.R2.fastq.gz       # paired FASTQ only
├── truth/
│   ├── PREFIX.methdb
│   └── PREFIX.variants.vcf.gz
└── PREFIX.manifest.json
```

With no variants, `--save-vcf` writes the prepared variant set as a valid
header-only VCF.gz. With an input VCF, it writes the normalized and
deterministically phased set used by the run. With an input MethDB,
`--save-methdb` verifies and copies that exact methylation profile into the new
truth directory. In WGS, WES, and TS, `--save-truth` publishes only the saved
variant set because those modes have no methylation profile.

### Execution controls

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| `-t` | `--threads` | integer `1..256` | `1`; one CPU budget for generation, processing, and compression |
| — | `--core` | executable path | bundled `htsim-core` |

`--threads` is a total CPU budget, not a per-stage multiplier. BSReadSim splits
it between ordered fragment construction and read processing; BAM runs also
reserve capacity for the HTSlib writer and BGZF compression. Protocol batches
and in-flight buffers are derived automatically. The same inputs and seeds
produce identical output bytes at every thread count. `--core` is primarily a
development and testing override.

## Technology-specific `run` options

### `run wgbs`

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--gc-profile` | target-GC TSV path | uniform sampling |

Supplying `--gc-profile` requests an output fragment-GC distribution. A fixed
insert supports reference-only, VCF, and de novo-variant runs. A variable
insert with a GC profile currently supports only a reference-only genome and
uses mean-insert calibration.

```bash
bsreadsim run wgbs \
  -r reference.fa \
  -o runs/wgbs \
  -n 1000000 \
  -l 150 \
  --insert-mean 300 \
  --insert-sd 0 \
  --mutation-rate 0.001 \
  --seed-mut 11 \
  --seed-meth 12 \
  -s 42 \
  -f bam \
  --save-truth
```

### `run wgs`

WGS uses the same whole-genome opportunity and optional GC-profile sampler as
WGBS, but emits ordinary reads: no methylome is constructed and no bisulfite
conversion is applied.

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--gc-profile` | target-GC TSV path | uniform sampling |

```bash
bsreadsim run wgs \
  -r reference.fa \
  -o runs/wgs \
  -n 1000000 \
  -l 150 \
  --insert-mean 300 \
  --insert-sd 0 \
  --mutation-rate 0 \
  -s 42 \
  -f bam \
  --save-truth
```

### `run rrbs`

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--cut-site` | motif; repeatable | required |
| — | `--rrbs-candidates` | candidate BED path | discover candidates internally |
| — | `--sampling` | `uniform` or `score` | `uniform` |

A motif places one `|` at the cut boundary, for example `C|CGG`. Repeat
`--cut-site` for multi-enzyme libraries. `--sampling score` requires
`--rrbs-candidates`; BSReadSim regenerates the candidate domain and permits an
external tool to change only the score column.

When changing the RRBS size window, keep `--insert-mean` inside that window
because the shared geometry contract still requires
`min <= mean <= max`. Candidate eligibility itself depends on min/max, not SD.

```bash
bsreadsim run rrbs \
  -r reference.fa \
  -o runs/rrbs \
  -n 1000000 \
  --cut-site 'C|CGG' \
  -l 100 \
  --insert-min 100 \
  --insert-mean 250 \
  --insert-max 500 \
  --insert-sd 25 \
  --mutation-rate 0 \
  -s 42
```

### `run tbs`

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--targets` | BED6 path | required |
| — | `--sampling` | `uniform` or `score` | `uniform` |
| — | `--fragment-center-stddev` | non-negative float | `50` |

TBS uses and defaults to `--insert-sd 0`, with `--insert-mean` as its fixed length.
`--fragment-center-stddev 0` centers every fragment exactly; larger values add
deterministic normal displacement. With `--sampling score`, BED column 5
supplies relative target mass.

```bash
bsreadsim run tbs \
  -r reference.fa \
  -o runs/tbs \
  -n 1000000 \
  --targets targets.bed \
  --sampling score \
  -l 150 \
  --insert-mean 300 \
  --insert-sd 0 \
  --fragment-center-stddev 50 \
  --mutation-rate 0 \
  -s 42
```

### `run wes`

WES uses the BED capture sampler without methylation or bisulfite chemistry.
Provide exon or other exome intervals as strand-aware BED6 targets. Target
scores and placement controls have the same meaning as in TBS.

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--targets` | BED6 path | required |
| — | `--sampling` | `uniform` or `score` | `uniform` |
| — | `--fragment-center-stddev` | non-negative float | `50` |

```bash
bsreadsim run wes \
  -r reference.fa \
  -o runs/wes \
  -n 1000000 \
  --targets exome.bed \
  --sampling uniform \
  -l 150 \
  --insert-mean 300 \
  --insert-sd 0 \
  --mutation-rate 0 \
  -s 42
```

### `run ts`

TS is the general ordinary targeted-sequencing mode. It shares the WES/TBS
BED capture sampler and differs from WES by its recorded technology identity
and intended target set.

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| — | `--targets` | BED6 path | required |
| — | `--sampling` | `uniform` or `score` | `uniform` |
| — | `--fragment-center-stddev` | non-negative float | `50` |

```bash
bsreadsim run ts \
  -r reference.fa \
  -o runs/ts \
  -n 1000000 \
  --targets panel.bed \
  --sampling score \
  -l 150 \
  --insert-mean 300 \
  --insert-sd 0 \
  --fragment-center-stddev 50 \
  --mutation-rate 0 \
  -s 42
```

## Compatibility rules worth checking first

- `--vcf` and `--mutation-rate` are mutually exclusive; `--vcf` automatically
  selects an internal mutation rate of `0`.
- `--asm` and `--asm-bed` require `--vcf`.
- `--methdb` embeds the prepared variants and therefore excludes `--vcf`,
  positive de novo mutation rates, methylation overlays, and pooling.
- `--cgmap-pool` requires `--cgmap` or `--bed-methyl`.
- `--quality-model` excludes `--phred`; `--error-model` excludes
  `--error-rate`.
- `--fragment-summary` and `--fragment-realization` require `--format bam`.
- WGBS variable-insert GC-profile sampling does not support variants.
- RRBS score sampling requires a matching candidate BED.
- TBS, WES, and TS default to and require `--insert-sd 0`; `--insert-mean` is their fixed length.
- WGS, WES, and TS reject methylation inputs and bisulfite-only truth or
  realization controls.

Invalid combinations fail before the final read artifacts and manifest are
published.

## `build variants`

Build a deterministic de novo variant set, a header-only no-variant set, or a
normalized and phased copy of an input VCF.

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| `-r` | `--reference` | FASTA path | required |
| `-o` | `--output` | new `.vcf.gz` path | required |
| — | `--vcf` | one-sample diploid VCF path | none |
| — | `--mutation-rate` | probability | `0.001`; mutually exclusive with `--vcf` |
| — | `--indel-fraction` | probability | `0.15` |
| — | `--indel-extension-probability` | probability | `0.15` |
| — | `--homozygous-only` | flag | off |
| — | `--seed-mut` | uint64 | `0` |
| — | `--seed-phase` | uint64 | `0` |
| — | `--core` | executable path | bundled `htsim-core` |

Generated variant set:

```bash
bsreadsim build variants \
  -r reference.fa \
  -o truth.variants.vcf.gz \
  --mutation-rate 0.001 \
  --indel-fraction 0.15 \
  --seed-mut 11
```

Normalize and phase an input VCF:

```bash
bsreadsim build variants \
  -r reference.fa \
  -o sample.phased.vcf.gz \
  --vcf sample.vcf.gz \
  --seed-phase 12
```

Supplying `--vcf` automatically selects an internal mutation rate of `0`. Set
`--mutation-rate 0` without `--vcf` to create a valid header-only VCF.gz.
Output is deterministic BGZF and is never appended to or overwritten.

## `build methdb`

Build one reusable prepared methylation profile as a MethDB file. This command
accepts the following options; unlike `run`, it has no master read-sampling
seed and does not accept `--methdb` as an input.

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| `-r` | `--reference` | FASTA path | required |
| `-o` | `--output` | new path | required |
| — | `--vcf` | one-sample diploid VCF path | none |
| — | `--cgmap` | CGmap path | none |
| — | `--bed-methyl` | bedMethyl path | none |
| — | `--asm` | ASM path | none; requires `--vcf` |
| — | `--asm-bed` | ASM BED path | none; requires `--vcf` |
| — | `--mutation-rate` | probability | `0.001`; mutually exclusive with `--vcf` |
| — | `--indel-fraction` | probability | `0.15` |
| — | `--indel-extension-probability` | probability | `0.15` |
| — | `--homozygous-only` | flag | off |
| — | `--seed-mut` | uint64 | `0` |
| — | `--seed-phase` | uint64 | `0` |
| — | `--seed-meth` | uint64 | `0` |
| — | `--beta-cg` | positive floats as `ALPHA,BETA` | `0.5,0.5` |
| — | `--beta-chg` | positive floats as `ALPHA,BETA` | `0.01,0.05` |
| — | `--beta-chh` | positive floats as `ALPHA,BETA` | `0.01,0.05` |
| — | `--cpg-only` | flag | off |
| — | `--cgmap-pool` | flag | off |
| — | `--methylation-model` | `bernoulli` or `bilstm` | `bernoulli` |
| — | `--no-update-variant-boundaries` | flag | off |
| — | `--core` | executable path | bundled `htsim-core` |

```bash
bsreadsim build methdb \
  -r reference.fa \
  -o sample.methdb \
  --mutation-rate 0 \
  --beta-cg 2,2 \
  --beta-chg 1,9 \
  --beta-chh 1,19 \
  --seed-meth 12
```

MethDB stores probabilities, not read-level Bernoulli realizations. Therefore
the state-model choice does not change the stored probability values.

## `build rrbs`

Export the exact RRBS candidate domain for external scoring. The final
`run rrbs` command must repeat every domain-defining option.

| Short | Long | Type | Default |
| --- | --- | --- | --- |
| `-r` | `--reference` | FASTA path | required |
| `-o` | `--output` | new candidate BED path | required |
| — | `--cut-site` | motif; repeatable | required |
| — | `--vcf` | one-sample diploid VCF path | none |
| — | `--single-end` | flag | paired-end |
| `-l` | `--read-length` | integer `1..10000` | `100` |
| — | `--insert-min` | positive integer | `100` |
| — | `--insert-mean` | positive integer | `400` |
| — | `--insert-max` | positive integer | `1000` |
| — | `--insert-sd` | non-negative float | `25` |
| — | `--max-ambiguous-fraction` | probability | `0.05` |
| — | `--mutation-rate` | probability | `0`; mutually exclusive with `--vcf` |
| — | `--indel-fraction` | probability | `0.15` |
| — | `--indel-extension-probability` | probability | `0.15` |
| — | `--homozygous-only` | flag | off |
| — | `--seed-mut` | uint64 | `0` |
| — | `--seed-phase` | uint64 | `0` |
| — | `--core` | executable path | bundled `htsim-core` |

```bash
bsreadsim build rrbs \
  -r reference.fa \
  -o candidates.bed \
  --cut-site 'C|CGG' \
  -l 100 \
  --insert-min 100 \
  --insert-mean 250 \
  --insert-max 500 \
  --mutation-rate 0
```

After an external program changes only the score column, consume the file with
`run rrbs --rrbs-candidates candidates.scored.bed --sampling score`.

## `export`

`export` has two kinds of source. `test-fasta` reads a release-bundled file,
so it needs only a destination. `methdb` decodes a user-supplied binary
methylation profile, so it requires both an input and an output.

| Target | Short | Long | Type | Default |
| --- | --- | --- | --- | --- |
| `test-fasta` | `-o` | `--output` | new FASTA path | required |
| `methdb` | `-i` | `--input` | existing MethDB path | required |
| `methdb` | `-o` | `--output` | new `.bed.gz` or `.bed` path | required |
| `methdb` | — | `--no-compression` | flag | off; write deterministic BGZF |
| `methdb` | — | `--core` | executable path | bundled `htsim-core` |

```bash
bsreadsim export test-fasta -o test.fa
bsreadsim export methdb -i sample.methdb -o sample.bed.gz
bsreadsim export methdb -i sample.methdb -o sample.bed --no-compression
```

`export test-fasta` verifies the bundled bytes against the release registry
and refuses to overwrite its destination. `export methdb` writes deterministic
BGZF by default and therefore requires a `.bed.gz` destination. With
`--no-compression`, it writes the same decoded bytes as plain text and requires
a `.bed` destination. Missing parent directories are created; existing output
files are never overwritten. The decoder streams both modes with bounded
buffers.

See [File formats](formats.md) for input contracts and
[Outputs](../outputs/index.md) for publication, manifest, and truth-file
details.
