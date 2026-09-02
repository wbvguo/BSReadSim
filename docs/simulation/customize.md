# Customize

<div class="customize-page" markdown>

Use this page to customize genome and methylation setup, fragment sampling,
read generation, and output options. Each section explains the available options
and when to use them. See [Tutorials](tutorials.md) for runnable examples and
the [CLI reference](../reference/cli.md) for full syntax and option combination
rules.

In the **Default** column, `—` means the option is unused by default, while
**Off** means a flag must be included in the command to enable it.

## Genome { #genetic-variation }

BSReadSim resolves the reference genome into two haplotypes before it identifies
methylatable cytosines, methylation contexts, restriction sites, or eligible
fragment positions.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-r`,<br>`--reference` | FASTA path | Required | Loads the reference genome to construct haplotypes |

All coordinate-based inputs must use the same reference build. See the
[Reference FASTA contract](../reference/formats.md#fasta).

### Generate de novo variants { #generated-variants }

When preparing the haplotypes, BSReadSim can randomly introduce de novo SNVs and
indels at the rate set by `--mutation-rate`.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--mutation-rate` | Float in `[0, 1]` | `0.001` | Sets the probability of a de novo mutation event at each non-`N` reference position |
| `--indel-fraction` | Float in `[0, 1]` | `0.15` | Sets the proportion of generated mutation events that are indels |
| `--indel-extension-probability` | Float in `[0, 1]` | `0.15` | Sets the probability that an indel extends by each additional base, up to four bases |
| `--seed-mut` | uint64 | Randomly generated | Sets the seed for generating de novo variants |
| `--homozygous-only` | Flag | Off | Generates every de novo variant on both haplotypes |

By default, each generated event is placed on both haplotypes with probability
`1/3`; the remaining events are assigned to either haplotype with equal
probability. `--homozygous-only` overrides that default genotype model and
places every generated event on both haplotypes.

!!! tip "If genetic variants is not desired"
    Set `--mutation-rate 0` when fragments should be sampled directly from the reference genome without genetic variation.

??? info "When `--mutation-rate` defaults to `0`"

    `--mutation-rate` normally defaults to `0.001`. Some input modes change how
    this option behaves:

    - With `--sampling gc`, `--rrbs-candidates`, omitting
      `--mutation-rate` changes its default to `0`.
    - With `--vcf` or `--methdb`, BSReadSim uses the supplied variants and does not generate
      de novo variants.

    See [Genetic variation](../reference/cli.md#genetic-variation) for all
    option combination rules.

### Load variants from a VCF { #vcf-genome }

Use a one-sample diploid VCF when the simulation should incorporate a predefined
variant set.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--vcf` | VCF path | — | Loads diploid variants and genotypes from VCF |
| `--seed-phase` | uint64 | Randomly generated | Sets the seed for assigning unphased variants to haplotypes |

??? info "How VCF phasing is handled"

    - For a phased heterozygous genotype (`0|1` or `1|0`), BSReadSim preserves
      the existing allele assignment to haplotypes 1 and 2.
    - For an unphased heterozygous genotype (`0/1` or `1/0`), BSReadSim randomly
      assigns the ALT allele to haplotype 1 or 2. `--seed-phase` makes this
      assignment reproducible.

    See the [VCF contract](../reference/formats.md#vcf) for accepted
    records, normalization behavior, and boundary handling.

## Methylation { #methylation }

BSReadSim considers methylation at two levels:

1. A **methylation profile** assigns a methylation probability to each
   cytosine site included in the simulation.
2. A **fragment-level state** determines whether that cytosine is methylated
   on a particular sampled molecule.

The profile is prepared before fragments are sampled. Methylated or
unmethylated states are then drawn separately for each physical fragment.

### Generate methylation probabilities { #generated-methylation }

For a generated methylation profile, BSReadSim draws each eligible cytosine's methylation probability independently. The draw uses the Beta distribution specified for that cytosine's context (CG, CHG, or CHH).

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--beta-cg` | `a,b` | `0.5,0.5` | Sets the CG-site Beta parameters |
| `--beta-chg` | `a,b` | `0.01,0.05` | Sets the CHG-site Beta parameters |
| `--beta-chh` | `a,b` | `0.01,0.05` | Sets the CHH-site Beta parameters |
| `--seed-meth` | uint64 | Randomly generated | Sets the seed for generating methylation probabilities |
| `--cpg-only` | Flag | Off | Omits CHG and CHH sites from the prepared profile |

These defaults define a general synthetic profile rather than one specific to
a tissue or species.

??? info "Beta distribution parameters"

    Each option accepts the two positive shape parameters `a,b` of a
    [`Beta(a, b)` distribution](https://en.wikipedia.org/wiki/Beta_distribution).

### Load methylation profiles { #predefined-methylation }

Load site-level methylation values from an existing profile, or reuse a MethDB
profile snapshot from a previous BSReadSim simulation.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--cgmap` | CGmap path | — | Loads methylation levels from CGmap |
| `--bedmethyl` | bedMethyl path | — | Loads methylation levels from bedMethyl |
| `--methbed` | MethBED path | — | Loads methylation levels from MethBED |
| `--methbg` | MethBG path | — | Loads methylation levels from MethBG |
| `--methdb` | MethDB path | — | Reuses a methylation profile snapshot with embedded variants |
| `--seed-meth` | uint64 | Randomly generated | Sets the seed for fallback generation and pooled-value resampling |
| `--pool-meth` | Flag | Off | Resamples input values within each contig and cytosine context |

??? info "How methylation profiles are applied"

    - **Text profiles:** [CGmap](../reference/formats.md#cgmap),
      [bedMethyl](../reference/formats.md#bedmethyl),
      [MethBED](../reference/formats.md#methbed), and
      [MethBG](../reference/formats.md#methbg) contain site-level methylation
      values only. By default, each value applies to its listed genomic
      position, while unlisted eligible sites use generated context-specific
      methylation probabilities. These profiles may also be used with a VCF or
      an ASM input.
    - **Pooling:** `--pool-meth` is available for all four text formats. It
      groups values by contig and cytosine context, then samples with
      replacement across eligible sites instead of preserving the original
      positions. If the input contains no values for a contig-context group,
      sites in that group use the configured Beta distribution instead.
    - **MethDB snapshot:** [MethDB](../reference/formats.md#methdb) reuses the
      resolved methylation and variant snapshot from an earlier build or
      simulation. It loads the snapshot directly and ignores other variant and
      methylation settings.

### Add allele-specific methylation { #allele-specific-methylation }

BSReadSim also supports site-level allele-specific methylation (ASM) in ASM or
ASM BED format.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--asm` | ASM path | — | Loads an ASM profile |
| `--asm-bed` | ASM BED path | — | Loads an ASM BED profile |

??? info "How ASM is applied"

    Each ASM record links one methylation site to one heterozygous SNV and
    provides separate methylation probabilities for its REF and ALT alleles.

    - **Haplotypes:** With a VCF, both haplotypes are constructed from the
    reference and the complete VCF variant set. Without a VCF, they are
    constructed from the reference and the linked SNVs defined by the ASM
    records.
    - **Application:** At each listed site, each haplotype is assigned the methylation probability associated with its linked-SNV allele. Other sites keep their baseline methylation probabilities.
    - **Requirements:** The linked SNV must be heterozygous. The ASM target must
      remain present with the same CG, CHG, or CHH context on both haplotypes.
      See
      [ASM formats](../reference/formats.md#allele-specific-methylation-inputs)
      for complete input requirements.

### Realize methylation states { #methylation-states }

A methylation profile assigns a methylation probability to each site. On a
sampled fragment, each site is realized as either methylated or unmethylated.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--meth-model` | `bernoulli` or `bilstm` | `bernoulli` | Selects the requested fragment-level methylation state model |

??? info "How methylation states are realized"

    Bernoulli samples each site's methylation state independently. Before
    bisulfite sequencing, methylation states are sampled for each fragment.
    The same site may have different states across fragments, but overlapping
    mates from one fragment share the sampled state.

    `bilstm` is accepted as a requested model. The current runtime emits an
    explicit warning, records both requested and effective models in the
    manifest, and falls back to Bernoulli state realization.

## Fragment sampling { #supported-technologies }

Fragmentation defines the assay-specific geometry. Sampling then
determines which eligible start, restriction fragment, or target contributes
each output fragment.

### Generate fragments

<span id="wgbs"></span>**WGBS**

Random fragmentation places fragments across eligible whole-genome positions,
with lengths drawn from the configured distribution.

<span id="rrbs"></span>**RRBS**

Restriction-enzyme cut sites generate the candidate fragments;
`--insert-min` and `--insert-max` apply size selection.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--cut-site` | Restriction-enzyme cut site | <code>C&#124;CGG</code> | Accepts a DNA motif with <code>&#124;</code> marking the cut position |

??? info "Cut-site input"

    The default `C|CGG` represents the MspI recognition motif `CCGG`, cut
    after the first C. Input is case-insensitive and normalized to uppercase;
    only `A`, `C`, `G`, `T`, and `N` are accepted. Each motif must contain
    exactly one `|`.

    Specify multiple unique cut sites as one comma-separated value:

    ```bash
    --cut-site 'C|CGG,G|ANTC'
    ```

<span id="tbs"></span>**TBS**

Probe-enriched regions are supplied with `--targets`. For each selected
target, the fragment center is displaced according to `--center-sd`, and its
length is drawn from the configured distribution.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--targets` | Capture target BED path | Required | Loads strand-aware capture targets from BED6 |
| `--center-sd` | Non-negative number | `50` | Sets the SD of fragment-center displacement, in bases |

??? info "How TBS represents probe enrichment"

    TBS uses target intervals to represent probe-enriched regions and places
    fragments around each target. Individual probes and capture chemistry are
    not simulated. See the
    [target BED contract](../reference/formats.md#capture-target-bed) for the
    input format.

### Set fragment length { #fragment-geometry }

Configure the length of each physical fragment.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--insert-min` | uint32 | `100` | Sets the minimum fragment length |
| `--insert-mean` | uint32 | `400` | Sets the mean fragment length |
| `--insert-max` | uint32 | `1000` | Sets the maximum fragment length |
| `--insert-sd` | Non-negative number | `25` | Sets the fragment-length standard deviation |

??? info "How fragment length is determined"

    - **Length distribution:** For WGBS, TBS, WGS, WES, and TS, BSReadSim draws
      a normal deviate controlled by `--insert-mean` and `--insert-sd`, truncates
      its offset toward zero to obtain an integer length, and clamps values
      outside `--insert-min` and `--insert-max` to the nearest bound. It does not
      redraw out-of-range values, so the two bounds can carry extra probability
      mass.
    - **RRBS size selection:** RRBS lengths come from restriction sites;
      `--insert-min` and `--insert-max` define the retained size window.
      `--insert-sd 0` does not make RRBS fragments equal in length.
    - **Targeted assays:** For TBS, WES, and TS, fragment lengths follow the
      distribution above, while `--center-sd` controls their displacement
      from target centers.
    - **Fixed length:** To use one fragment length, set `--insert-mean N`
      together with `--insert-sd 0`.

    See [Fragment length](../reference/cli.md#fragment-geometry) for shared
    fragment-length constraints.

### Sample fragments

Fragments are sampled uniformly by default. Each assay provides options for
non-uniform sampling.

**WGBS**

`--sampling gc` enables non-uniform sampling to match the fragment-GC
distribution supplied by `--gc-profile`.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--sampling` | `uniform` or `gc` | `uniform` | Selects uniform or GC-profile sampling |
| `--gc-profile` | Target-GC profile path | — | Supplies the distribution profile used by `--sampling gc` |

??? info "How GC-profile sampling works"

    Rejection sampling accepts or rejects each proposed fragment with a
    probability calibrated to its GC content, producing the distribution
    specified in the
    [target-GC profile](../reference/formats.md#target-gc-profile).

**RRBS**

`--sampling score` uses the candidate BED `score` column to set relative
sampling probabilities. The scores can be generated by an external model.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--sampling` | `uniform` or `score` | `uniform` | Selects uniform or score-weighted candidate sampling |
| `--rrbs-candidates` | RRBS candidate BED path | — | Supplies RRBS candidate scores for non-uniform fragment sampling |

??? info "How score-weighted RRBS sampling works"

    `build rrbs` creates a candidate BED with every score initialized to `1`.
    Update only the `score` column with values from an external model, then
    load the scored BED with `--rrbs-candidates`. Sampling weights are
    calculated from the updated scores. See the
    [RRBS candidate contract](../reference/formats.md#rrbs-candidate-bed).

**TBS**

`--sampling score` samples targets in proportion to the values in BED column 5.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--sampling` | `uniform` or `score` | `uniform` | Selects uniform or score-weighted target sampling |

??? info "How score-weighted target sampling works"

    BED column 5 directly supplies each target's relative weight. The target
    is selected before fragment length and center displacement are drawn. See
    the [target BED contract](../reference/formats.md#capture-target-bed).

## Bisulfite conversion { #bisulfite-conversion }

Configure the library orientation and conversion rate for WGBS, RRBS, and TBS.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--conversion-rate` | Float in `[0, 1]` | `0.998` | Sets the probability that each unmethylated cytosine is converted |
| `--undirectional` | Flag | Off | Selects an undirectional bisulfite library |

??? info "How the bisulfite library and conversion are applied"

    The default directional library includes OT and OB molecules.
    `--undirectional` also includes CTOT and CTOB molecules.

    Bisulfite conversion is applied once to each complete physical fragment
    before read extraction. Unmethylated cytosines are converted according to
    `--conversion-rate`, whereas methylated cytosines remain unchanged.

## Read generation { #library-and-sequencing }

After fragment sampling and, where applicable, bisulfite conversion,
sequencing generates reads from the resulting physical fragments. The settings
below control read count, layout and length, base qualities, and sequencing
errors.

### Set number of reads { #dataset-size }

Set the output read count directly or derive it from sequencing depth:

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-n`,<br>`--reads` | Positive integer | `1,000,000` | Sets the exact number of output read records |
| `-d`,<br>`--depth` | Number greater than `0` | — | Derives the read count from mean depth, effective-region size, and read length |

??? info "How the number of reads is resolved"

    `--reads N` produces exactly `N` records. In paired-end mode, `N` must be
    even, yielding `N / 2` read pairs. When `--depth D` is specified, the
    number of reads is calculated automatically from the requested mean depth,
    the assay's effective-region size, and the configured read length.

    See [Required inputs and read count](../reference/cli.md#required-inputs-and-dataset-size)
    for the depth formula and assay-specific effective regions.

### Set read layout { #read-layout }

Choose single- or paired-end sequencing and set the read length.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-l`,<br>`--read-length` | Integer from `1` to `10000` | `100` | Sets the number of bases in each read |
| `--max-ambiguous-fraction` | Float in `[0, 1]` | `0.05` | Sets the maximum allowed fraction of `N` bases per read |
| `--single-end` | Flag | Off | Selects single-end sequencing |

??? info "How read layout is applied"

    Paired-end sequencing is the default and produces one read pair per
    fragment; `--single-end` produces one read per fragment. In either mode,
    every emitted read must stay within the `--max-ambiguous-fraction` limit.
    If any read exceeds the limit, the fragment is skipped and another
    fragment is sampled.

### Configure base quality and sequencing errors { #quality-and-error }

Use constant per-base settings or models estimated from sequencing data.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-q`,<br>`--phred` | Integer from `0` to `93` | `40` | Sets a fixed Phred score for every base |
| `--quality-model` | Quality-model JSON path | — | Samples each cycle's Phred score from a quality Markov model |
| `-e`,<br>`--error-rate` | Float in `[0, 1]` | `0.005` | Sets a uniform substitution probability for every base |
| `--error-model` | Error-model JSON path | — | Samples each final base call from a Phred-specific base transition model |

??? info "How estimated sequencing models are applied"

    Instead of using a fixed Phred score and a uniform sequencing-error rate,
    `--quality-model` and `--error-model` load models estimated from real
    sequencing data and apply them when generating synthetic reads:

    - `--quality-model` generates separate Phred-score sequences for R1 and
      R2. The first five cycles use individual initial distributions; each
      later score is drawn from the transition matrix row for the preceding
      score.
    - `--error-model` selects a base-transition matrix by mate and Phred score,
      then draws the final A, C, G, or T base from the row for the original
      base.

    See the
    [sequencing-model formats](../reference/formats.md#sequencing-model-inputs) for
    the required fields and dimensions.

## Output and reproducibility { #reproducibility }

Configure output, optional truth artifacts, and random seeds. A
manifest is written automatically for every successful run.

### Configure execution and output { #output-format }

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `-o`,<br>`--output` | Directory path | Required | Sets the output directory |
| `-p`,<br>`--prefix` | `[A-Za-z0-9._-]+`, up to 128 characters | `sim` | Sets the output filename prefix |
| `-f`,<br>`--format` | `fastq`, `fastq.gz`, or `bam` | `fastq.gz` | Selects the read output format |
| `--gzip-level` | Integer from `0` to `9` | `6` | Sets the compression level for `fastq.gz` output |
| `-t`,<br>`--threads` | Integer, `1`–`256` | `1` | Sets the number of threads |

??? info "FASTQ and BAM output"

    Read layout and output format are independent. `fastq` and `fastq.gz`
    support both single- and paired-end reads; R2 is written only for
    paired-end reads. `bam` also supports both layouts and is written instead
    of FASTQ.

??? info "Threads"

    BSReadSim uses the requested threads across fragment generation, read
    processing, output ordering, and BAM compression when applicable. The best
    value depends on available CPUs, memory, output format, and storage speed.
    Thread count changes resource use but not fixed-seed output bytes. See
    [Performance](../help/troubleshoot.md#performance) when tuning a run.

### Save simulation truth { #truth-artifacts }

Optionally add fragment-level truth to BAM or save reusable variant sets and
methylation profile snapshots.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--fragment-summary` | Flag | Off | Adds compact fragment metadata to BAM records |
| `--fragment-realization` | Flag | Off | Adds complete-fragment methylation and conversion states to BAM records |
| `--save-methdb` | Flag | Off | Writes the methylation profile to MethDB with embedded variants |
| `--save-vcf` | Flag | Off | Writes the prepared, phased variant set to VCF |
| `--save-truth` | Flag | Off | Writes the variant set and, for bisulfite assays, the methylation profile to disk |

??? info "BAM annotations and reusable truth"

    **BAM annotations:** Both fragment options require BAM.
    `--fragment-realization` is available only for bisulfite assays and implies
    `--fragment-summary`.

    **Reusable truth:** `--save-methdb` is available only for bisulfite assays;
    VCF truth is available for every assay.

    See [Outputs](../outputs/index.md) for filenames and file contracts.

### Control random variation { #random-seeds }

Separate seeds control preparation of the biological state and subsequent
fragment and read generation.

| Option | Value | Default | Description |
| --- | --- | --- | --- |
| `--seed-mut` | uint64 | Randomly generated | Sets the seed for generating de novo variants |
| `--seed-phase` | uint64 | Randomly generated | Sets the seed for assigning unphased variants to haplotypes |
| `--seed-meth` | uint64 | Randomly generated | Sets the seed for preparing methylation probabilities |
| `--seed` | uint64 | Randomly generated | Sets the master seed for the simulation |

??? info "How to reproduce or vary a run"

    When `--seed` is omitted, it is generated automatically. Each omitted
    biological seed is then derived independently from the resolved `--seed`.
    All resolved seeds are recorded in the run manifest and its fully expanded
    command. An explicitly supplied biological seed overrides derivation for
    that stage.

    Fix every relevant seed for an identical rerun. To draw different fragments
    and read-level states from the same prepared genome and methylation profile,
    keep `--seed-mut`, `--seed-phase`, and `--seed-meth` fixed and change only
    `--seed`.

</div>
