# CLI reference

Use this page to look up the options accepted by each BSReadSim command. For
runnable commands, see [Tutorials](../simulation/tutorials.md). For guidance on
choosing settings, see [Customize](../simulation/customize.md). Input file
syntax is documented in [Input file formats](formats.md).

<div class="reference-jump" markdown>

**Jump to:** [commands](#commands) · [run simulations](#run-simulations) ·
[option combination rules](#option-combination-rules) · [build variants](#build-variants) ·
[build MethDB](#build-methdb) · [build RRBS candidates](#build-rrbs-candidates) ·
[export files](#export-files) · [core integration options](#core-integration-options)

</div>

## Commands

| Command form | Use it to | Result |
| --- | --- | --- |
| `bsreadsim run ASSAY` | Simulate `wgbs`, `rrbs`, `tbs`, `wgs`, `wes`, or `ts` reads | Reads and a run manifest |
| `bsreadsim build variants` | Generate variants or normalize and phase a VCF | Prepared VCF.gz |
| `bsreadsim build methdb` | Prepare a reusable methylation and variant snapshot | MethDB |
| `bsreadsim build rrbs` | Enumerate the RRBS candidate domain for external scoring | RRBS candidate BED |
| `bsreadsim export test-fasta` | Copy the bundled example reference | FASTA |
| `bsreadsim export methdb` | Decode a MethDB for inspection | Extended BED |

Only `run` creates a manifest. Use `bsreadsim --help`,
`bsreadsim COMMAND --help`, or `bsreadsim COMMAND SUBCOMMAND --help` to inspect
the interface installed in the current environment.

## Conventions

The **Default** column uses **Required** for mandatory options, **—** when no
optional input is selected, and **Disabled** for flags that are off unless
supplied.

- A probability is a finite number from `0` to `1`, inclusive.
- A uint64 seed is a decimal integer from `0` to `18446744073709551615`.
- Relative paths are resolved from the working directory.
- Existing destination artifacts are not overwritten.

### Help and version

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-h`, `--help` | Flag | Disabled | Shows help for the current command level and exits |
| `--version` | Flag | Disabled | Prints the version and exits; available only before a command |

## Run simulations

Options in this section apply to all six assays unless a narrower set is
stated. This table summarizes the assay-specific parts of each command:

| Command | Fragment sampling | Required assay input | Bisulfite options |
| --- | --- | --- | --- |
| `run wgbs` | `uniform` or `gc` | — | Yes |
| `run rrbs` | `uniform` or `score` | —; `--cut-site` defaults to <code>C&#124;CGG</code> | Yes |
| `run tbs` | `uniform` or `score` | `--targets` | Yes |
| `run wgs` | `uniform` or `gc` | — | No |
| `run wes` | `uniform` or `score` | `--targets` | No |
| `run ts` | `uniform` or `score` | `--targets` | No |

### Required inputs and read count { #required-inputs-and-dataset-size }

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-r`, `--reference` | FASTA path | Required | Loads the reference genome |
| `-o`, `--output` | Directory path | Required | Sets the output directory |
| `-n`, `--reads` | Positive integer | `1,000,000` | Sets the exact number of output read records |
| `-d`, `--depth` | Number greater than `0` | — | Derives the read count from mean depth, effective-region size, and read length |

`--reads` and `--depth` are mutually exclusive. If neither is supplied,
`--reads` defaults to `1,000,000` read records. In paired-end mode, `--reads`
must be even and counts both mates: `--reads 1000` produces 500 read pairs. In
single-end mode, it produces 1000 reads.

??? info "Details"

    `--reads` accepts up to `4294967295` single-end records or an even number
    up to `8589934590` in paired-end mode. Depth is resolved to complete
    fragments:

    ```text
    raw_reads = effective_reference_bases * D / read_length
    resolved_reads = emitted_mates * ceil(raw_reads / emitted_mates)
    resolved_fragments = resolved_reads / emitted_mates
    ```

    The effective region is eligible contig sequence for WGBS and WGS,
    eligible restriction-fragment envelopes for RRBS, and the union of target
    intervals for TBS, WES, and TS.

### Genome and variants { #genetic-variation }

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--vcf` | VCF path | — | Loads diploid variants and genotypes from VCF |
| `--mutation-rate` | Float in `[0, 1]` | `0.001` | Sets the probability of a de novo mutation event at each non-`N` reference position |
| `--indel-fraction` | Float in `[0, 1]` | `0.15` | Sets the proportion of generated mutation events that are indels |
| `--indel-extension-probability` | Float in `[0, 1]` | `0.15` | Sets the probability that an indel extends by each additional base, up to four bases |
| `--homozygous-only` | Flag | Disabled | Generates every de novo variant on both haplotypes |

??? info "Details"

    - By default, each generated event is homozygous with probability `1/3`.
      Otherwise, it is assigned to either haplotype with equal probability.
      `--homozygous-only` is still an active override in the current CLI: it
      forces every generated event onto both haplotypes.
    - When a VCF is provided, BSReadSim uses its diploid variant set instead
      of generating de novo variants. `--vcf` cannot be combined with
      `--mutation-rate`. Existing phasing is preserved; `--seed-phase`
      assigns unphased heterozygous variants to haplotypes.
    - VCF indels longer than four bases after normalization are rejected, as
      are MNPs and complex replacements.
    - The default mutation rate becomes `0` when `--gc-profile`,
      `--rrbs-candidates`, or `--methdb` is supplied. The indel controls and
      `--homozygous-only` affect only generated events; they do not modify
      variants loaded from VCF or MethDB.

### Methylation inputs

Available only for WGBS, RRBS, and TBS.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--methbg` | MethBG path | — | Loads methylation levels from MethBG |
| `--methbed` | MethBED path | — | Loads methylation levels from MethBED |
| `--bedmethyl` | bedMethyl path | — | Loads methylation levels from bedMethyl |
| `--cgmap` | CGmap path | — | Loads methylation levels from CGmap |
| `--methdb` | MethDB path | — | Reuses a methylation profile snapshot with embedded variants |
| `--asm` | CGmapTools ASS path | — | Adds `cgmaptools asm -m ass` rows called `TRUE` |
| `--asm-bed` | ASM BED path | — | Adds BED6+6 or BED6+10 ASM |

See [Input file formats](formats.md#methylation-profile-inputs) for row contracts,
[Customize](../simulation/customize.md#predefined-methylation) for profile
behavior, and [Option combination rules](#option-combination-rules) for allowed
combinations.

### Methylation profile and state model

Available only for WGBS, RRBS, and TBS.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--beta-cg` | Positive floats as `a,b` | `0.5,0.5` | Sets the CG-site Beta parameters |
| `--beta-chg` | Positive floats as `a,b` | `0.01,0.05` | Sets the CHG-site Beta parameters |
| `--beta-chh` | Positive floats as `a,b` | `0.01,0.05` | Sets the CHH-site Beta parameters |
| `--meth-model` | `bernoulli` | `bernoulli` | Selects the fragment-level methylation state model |
| `--cpg-only` | Flag | Disabled | Omits CHG and CHH sites from the prepared profile |
| `--pool-meth` | Flag | Disabled | Resamples input values within each contig and cytosine context |

`--pool-meth` requires one of the four text profile inputs; it cannot be used
with a generated profile or MethDB.

See [Customize](../simulation/customize.md#methylation) for probability
generation, profile fallback, pooling, and fragment-level state realization.

### Fragment length { #fragment-geometry }

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--insert-min` | Positive integer | `100` | Sets the minimum fragment length |
| `--insert-mean` | Positive integer | `400` | Sets the mean fragment length |
| `--insert-max` | Positive integer | `1000` | Sets the maximum fragment length |
| `--insert-sd` | Non-negative number | `25` | Sets the fragment-length standard deviation |

??? info "Details"

    - Fragment settings must satisfy
      `insert_min <= insert_mean <= insert_max`.
    - For whole-genome and targeted assays, `--insert-mean N --insert-sd 0`
      selects a fixed fragment length of `N`.
    - RRBS fragment lengths are determined by restriction sites; the minimum
      and maximum values define the retained length range. Its mean and SD do
      not reshape those restriction fragments.
    - Each read length must be no greater than `--insert-min`, except that a
      fixed-length whole-genome or targeted run compares it with
      `--insert-mean`.

### Read layout { #read-layout }

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-l`, `--read-length` | Integer from `1` to `10000` | `100` | Sets the number of bases in each read |
| `--max-ambiguous-fraction` | Float in `[0, 1]` | `0.05` | Sets the maximum allowed fraction of `N` bases per read |
| `--single-end` | Flag | Disabled | Selects single-end sequencing |

### Bisulfite conversion

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--conversion-rate` | Float in `[0, 1]` | `0.998` | Sets the probability that each unmethylated cytosine is converted |
| `--undirectional` | Flag | Disabled | Selects an undirectional bisulfite library |

Available only for WGBS, RRBS, and TBS. See
[Bisulfite conversion](../simulation/customize.md#bisulfite-conversion)
for molecule types and fragment-level conversion behavior.

### Base quality and sequencing errors

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-q`, `--phred` | Integer from `0` to `93` | `40` | Sets a fixed Phred score for every base |
| `--quality-model` | Quality-model JSON path | — | Samples each cycle's Phred score from a quality Markov model |
| `-e`, `--error-rate` | Float in `[0, 1]` | `0.005` | Sets a uniform substitution probability for every base |
| `--error-model` | Error-model JSON path | — | Samples each final base call from a Phred-specific base transition model |

The JSON inputs follow the [sequencing-model formats](formats.md#sequencing-model-inputs);
exclusive option pairs are listed under
[Option combination rules](#option-combination-rules).

### Read output

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-p`, `--prefix` | `[A-Za-z0-9._-]+`, up to 128 characters | `sim` | Sets the output filename prefix |
| `-f`, `--format` | `fastq`, `fastq.gz`, or `bam` | `fastq.gz` | Selects the read output format |
| `--gzip-level` | Integer from `0` to `9` | `6` | Sets the compression level for `fastq.gz` output |

FASTQ supports single- and paired-end reads; R2 is written only in paired-end
mode. BAM replaces FASTQ when `--format bam` is selected. `--gzip-level` has
no effect with uncompressed FASTQ or BAM. See
[Outputs](../outputs/index.md) for filenames, annotations, and completion
rules.

### BAM truth annotations

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--fragment-summary` | Flag | Disabled | Adds compact fragment metadata to BAM records |
| `--fragment-realization` | Flag | Disabled | Adds complete-fragment methylation and conversion states to BAM records |

Both options require BAM. Fragment realization is available only for
bisulfite assays and implies fragment summaries. See
[Annotated BAM](../outputs/index.md#annotated-bam) for tag contracts.

### Save and reuse simulation truth

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--save-methdb` | Flag | Disabled | Writes the methylation profile to MethDB with embedded variants |
| `--save-vcf` | Flag | Disabled | Writes the prepared, phased variant set to VCF |
| `--save-truth` | Flag | Disabled | Writes the variant set and, for bisulfite assays, the methylation profile to disk |

`--save-methdb` is available only for bisulfite assays. For non-bisulfite
assays, `--save-truth` exports only the prepared, phased variant set as VCF.
Reuse a saved VCF with `--vcf`, or reuse the prepared methylation profile and
its embedded variant set with `--methdb`. See
[Ground truth in simulation](../simulation/tutorials.md#ground-truth-in-simulation)
for runnable Save and Reuse examples.

### Random seeds

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-s`, `--seed` | uint64 | Generated and recorded | Sets the master seed for the simulation |
| `--seed-mut` | uint64 | Generated and recorded | Sets the seed for generating de novo variants |
| `--seed-phase` | uint64 | Generated and recorded | Sets the seed for assigning unphased variants to haplotypes |
| `--seed-meth` | uint64 | Generated and recorded | Sets the seed for preparing methylation probabilities |

`--seed-meth` is available only for WGBS, RRBS, and TBS. When omitted,
`--seed` is generated and each stage seed exposed by the assay is derived
independently from it. The manifest records every resolved seed. Explicit
stage seeds override derivation for their respective stages.

Fix all relevant seeds for byte-identical output. Keeping the biological seeds
fixed while changing only `--seed` draws different reads from the same
prepared snapshot.

### Execution

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-t`, `--threads` | Integer from `1` to `256` | `1` | Sets the total CPU thread budget |

Thread count changes resource use, not fixed-seed output bytes.

### Assay-specific fragment generation and sampling

#### Whole-genome assays

Available for whole-genome bisulfite sequencing (`run wgbs`) and whole-genome
sequencing (`run wgs`).

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--sampling` | `uniform` or `gc` | `uniform` | Selects uniform or GC-profile sampling |
| `--gc-profile` | Target-GC profile path | — | Supplies the distribution profile used by `--sampling gc` |

#### Reduced representation bisulfite sequencing (RRBS)

Available for `run rrbs`.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--sampling` | `uniform` or `score` | `uniform` | Selects uniform or score-weighted fragment sampling |
| `--cut-site` | Restriction-enzyme cut site | <code>C&#124;CGG</code> | Accepts a DNA motif with <code>&#124;</code> marking the cut position |
| `--rrbs-candidates` | RRBS candidate BED path | — | Supplies RRBS candidate scores for non-uniform fragment sampling |

`--cut-site` is case-insensitive and accepts one or more unique motifs in a
single comma-separated value, for example `--cut-site 'C|CGG,G|ANTC'`. Do not
repeat the option. Each motif must contain exactly one `|`, with only `A`, `C`,
`G`, `T`, or `N` on either side and at least one base after the cut.
`--sampling score` requires `--rrbs-candidates`; the candidate file must match
the run's candidate-defining settings.

#### Targeted assays

Available for targeted bisulfite sequencing (`run tbs`), whole-exome
sequencing (`run wes`), and targeted sequencing (`run ts`).

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--sampling` | `uniform` or `score` | `uniform` | Selects uniform or score-weighted target sampling |
| `--targets` | Capture target BED path | Required | Loads strand-aware capture targets from BED6 |
| `--center-sd` | Non-negative number | `50` | Sets the SD of fragment-center displacement, in bases |

## Option combination rules { #option-combination-rules }

**Genome and methylation inputs**

- Choose no more than one baseline methylation input and no more than one ASM
  representation.
- `--vcf` excludes `--mutation-rate`; ASM excludes nonzero de novo mutation
  generation but does not require a VCF.
- `--methdb` excludes VCF, de novo mutations, ASM, and methylation-value
  pooling.
- `--pool-meth` requires `--cgmap`, `--bedmethyl`, `--methbg`, or `--methbed`.
- With MethDB, beta parameters and `--cpg-only` do not alter the stored
  profile. `--meth-model` still controls read-level state realization.

**Read generation and truth output**

- `--quality-model` excludes `--phred`; `--error-model` excludes
  `--error-rate`.
- Fragment summaries and realization require BAM; realization is
  bisulfite-only.
- `--gzip-level` affects only `--format fastq.gz`.

**Fragment sampling**

- `--sampling gc` and `--gc-profile` must be supplied together.
- Variable-length GC-profile sampling does not support variants.
- RRBS score-weighted sampling requires `--rrbs-candidates`, and the candidate
  identities must match the current candidate domain.

**Assay availability**

- WGS, WES, and TS reject methylation inputs and bisulfite-only truth options.

Invalid combinations fail before final artifacts are written.

## Build reusable artifacts

All three `build` subcommands share the genome and variant controls below.
Unlike `run`, they have no master seed: their stage seeds default to `0`.
Build destinations must be new files in existing parent directories.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-r`, `--reference` | FASTA path | Required | Loads the reference genome |
| `--vcf` | VCF path | — | Loads variants instead of generating them |
| `--mutation-rate` | Float in `[0, 1]` | `0.001`; `0` for `build rrbs` | Sets the de novo event probability at each non-`N` reference position |
| `--indel-fraction` | Float in `[0, 1]` | `0.15` | Sets the proportion of generated events that are indels |
| `--indel-extension-probability` | Float in `[0, 1]` | `0.15` | Sets extension probability for generated indels, up to four bases |
| `--seed-mut` | uint64 | `0` | Sets the de novo mutation seed |
| `--seed-phase` | uint64 | `0` | Sets the seed for unphased VCF heterozygotes |
| `--homozygous-only` | Flag | Disabled | Forces every generated event onto both haplotypes |

`--vcf` and `--mutation-rate` are mutually exclusive. Indel controls,
`--seed-mut`, and `--homozygous-only` affect generated events only. The
[genome and variant rules](#genetic-variation) describe the generated genotype
model and accepted VCF variants.

### Build variants

`build variants` generates a reusable de novo variant set or normalizes and
phases variants from an existing VCF.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-o`, `--output` | New `.vcf.gz` path | Required | Sets the VCF output path |

Without a VCF, `--mutation-rate 0` creates a valid header-only VCF.gz. Output
is deterministic BGZF and is never appended to or overwritten.

### Build MethDB

`build methdb` prepares a reusable methylation profile and embedded prepared
variant set without generating reads. It accepts the shared build controls and
these additional options:

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-o`, `--output` | New MethDB path | Required | Sets the MethDB output path |
| `--methbg` | MethBG path | — | Loads methylation levels from MethBG |
| `--methbed` | MethBED path | — | Loads methylation levels from MethBED |
| `--bedmethyl` | bedMethyl path | — | Loads methylation levels from bedMethyl |
| `--cgmap` | CGmap path | — | Loads methylation levels from CGmap |
| `--asm` | CGmapTools ASS path | — | Adds `cgmaptools asm -m ass` rows called `TRUE` |
| `--asm-bed` | ASM BED path | — | Adds BED6+6 or BED6+10 ASM |
| `--beta-cg` | Positive floats as `a,b` | `0.5,0.5` | Sets the CG-site Beta parameters |
| `--beta-chg` | Positive floats as `a,b` | `0.01,0.05` | Sets the CHG-site Beta parameters |
| `--beta-chh` | Positive floats as `a,b` | `0.01,0.05` | Sets the CHH-site Beta parameters |
| `--meth-model` | `bernoulli` | `bernoulli` | Shares the build/run configuration field but does not change MethDB output |
| `--seed-meth` | uint64 | `0` | Sets the methylation-probability seed |
| `--cpg-only` | Flag | Disabled | Omits CHG and CHH sites from the prepared profile |
| `--pool-meth` | Flag | Disabled | Resamples text-profile values by contig and context |

Choose at most one baseline methylation input and at most one ASM input.
`--pool-meth` requires a text profile. `build methdb` does not accept MethDB as
an input. Because a snapshot stores probabilities rather than fragment-level
states, `--meth-model` has no effect on this command and should normally be
omitted. See [Option combination rules](#option-combination-rules) for input
combinations and
[Input file formats](formats.md#methdb) for the snapshot contract.

### Build RRBS candidates

`build rrbs` generates the exact candidate domain used by reduced
representation bisulfite sequencing (RRBS). It accepts the shared build
controls and these additional options:

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-o`, `--output` | New candidate BED path | Required | Sets the candidate BED output path |
| `--cut-site` | Comma-separated cut sites | <code>C&#124;CGG</code> | Sets one or more unique restriction motifs |
| `-l`, `--read-length` | Integer from `1` to `10000` | `100` | Sets the read length used to validate candidates |
| `--insert-min` | Positive integer | `100` | Sets the minimum retained restriction-fragment length |
| `--insert-mean` | Positive integer | `400` | Participates in shared geometry validation but does not change the candidate BED |
| `--insert-max` | Positive integer | `1000` | Sets the maximum retained restriction-fragment length |
| `--insert-sd` | Non-negative number | `25` | Is accepted for consistency but does not change the candidate BED |
| `--max-ambiguous-fraction` | Float in `[0, 1]` | `0.05` | Sets the maximum allowed `N` fraction in each candidate read |
| `--single-end` | Flag | Disabled | Validates one read per candidate instead of two |

The reference, variant source and focused seeds, cut sites, read layout and
length, insert bounds, and ambiguity threshold define candidate identities.
`--insert-mean` must still lie between the bounds because all fragment settings
share one validation contract; neither it nor `--insert-sd` needs to match the
later RRBS run. After external scoring, load the file with
`--sampling score --rrbs-candidates PATH` without changing any identity field.

## Export files

### Test FASTA

`export test-fasta` writes the bundled example reference and verifies its
bytes.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-o`, `--output` | New FASTA path | Required | Sets the new FASTA output path |

### MethDB inspection

`export methdb` decodes a MethDB snapshot as a human-readable extended BED.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-i`, `--input` | MethDB path | Required | Loads the source MethDB |
| `-o`, `--output` | New `.bed.gz` or `.bed` path | Required | Sets the new extended BED output path |
| `--no-compression` | Flag | Disabled | Writes plain extended BED |

??? info "Details"

    Compressed output must end in `.bed.gz` and is deterministic BGZF. With
    `--no-compression`, the output must end in `.bed`. Missing parent
    directories are created, but an existing destination is never
    overwritten. This export is intended for inspection; use the original
    MethDB to reuse the complete snapshot in another simulation.

## Core integration options { #core-integration-options }

These low-level options are intended for core integration and development.
Normal simulations and artifact builds should omit them.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--core` | Executable path | Bundled `htsim-core` | Overrides the core for every `run` and `build` command and for `export methdb` |
| `--no-update-variant-boundaries` | Flag | Disabled | Disables context updates; variant-aware bisulfite preparation rejects it |
