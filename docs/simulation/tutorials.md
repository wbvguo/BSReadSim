# Tutorials

These tutorials provide runnable commands for common bisulfite-sequencing
tasks. Replace the example paths before running them, and use inputs from the
same reference assembly. For WGS, WES, and TS examples, see
[Other assays](other-assays.md).

| Goal | Tutorial |
| --- | --- |
| Generate WGBS from a reference genome | [WGBS from a reference genome](#wgbs) |
| Incorporate existing variants and methylation | [Variant sets and methylation profiles](#variant-sets-and-methylation-profiles) |
| Add haplotype-specific methylation | [Allele-specific methylation](#allele-specific-methylation) |
| Simulate restriction or targeted sampling | [Enrichment-based assays](#enrichment-based-assays) |
| Export annotated BAM and reusable truth | [Save ground truth](#save) |
| Simulate new reads from a saved snapshot | [Reuse ground truth](#reuse) |

Use [Customize](customize.md) to choose models and review their default
behavior. Use the [CLI reference](../reference/cli.md) for accepted values and
complete option combination rules.

## WGBS from a reference genome { #wgbs }

This example generates 100,000 paired-end WGBS read records from a reference
FASTA, with de novo variants and a simulated methylation profile.

```bash
bsreadsim run wgbs \
  -r test.fa \
  -o test/wgbs-denovo \
  -n 100000 \
  --mutation-rate 0.001 \
  --seed-mut 7 \
  --seed-meth 11 \
  -s 42
```

With the default paired-end mode and `fastq.gz` format, the output directory
has this structure:

```text { .no-copy }
test/wgbs-denovo/
├── sim.R1.fastq.gz
├── sim.R2.fastq.gz
└── sim.manifest.json
```

??? info "How it works"

    **Generation:** With no VCF or methylation profile, BSReadSim generates
    [de novo variants](customize.md#generated-variants) and
    [context-specific methylation probabilities](customize.md#generated-methylation).
    [Whole-genome fragments](customize.md#generate-fragments) use the
    configured fragment-length distribution, and the default
    [Bernoulli model](customize.md#methylation-states) realizes methylation
    states before bisulfite conversion.

    **Sampling:** Eligible whole-genome starts are
    [sampled uniformly](customize.md#sample-fragments). The explicit
    [random seeds](customize.md#random-seeds) make genome preparation, fragment
    sampling, and read generation reproducible.

## Variant sets and methylation profiles

This recipe combines a one-sample diploid VCF with a site-level methylation
profile from the same reference assembly. It uses CGmap, but bedMethyl,
MethBED, and MethBG are also supported.

```bash
bsreadsim run wgbs \
  -r test.fa \
  -o test/wgbs-profile \
  -n 100000 \
  --vcf test.vcf.gz \
  --cgmap test.cgmap.gz \
  --seed-phase 7 \
  --seed-meth 11 \
  -s 42
```

??? info "How it works"

    When a VCF is provided, BSReadSim uses its diploid variant set instead of
    generating de novo variants. Existing phasing is preserved; unphased
    heterozygous variants are assigned to haplotypes using `--seed-phase`.

    The CGmap file supplies methylation levels. Unlisted eligible sites
    receive generated context-specific methylation probabilities rather than
    being treated as unmethylated.

    See [Load variants from a VCF](customize.md#vcf-genome),
    [Load methylation profiles](customize.md#predefined-methylation),
    and the [input file formats](../reference/formats.md).

## Allele-specific methylation

Use an ASM profile to assign different methylation probabilities to linked REF
and ALT alleles. A VCF can optionally be provided to preserve its complete
variant set and existing phasing.

```bash
bsreadsim run wgbs \
  -r test.fa \
  -o test/wgbs-asm \
  -n 100000 \
  --cgmap test.cgmap.gz \
  --asm test.asm.gz \
  --seed-phase 7 \
  --seed-meth 11 \
  -s 42
```

??? info "How it works"

    Without a VCF, the heterozygous SNVs linked by the ASM profile define the
    variant set. `--seed-phase` makes their haplotype assignment reproducible.

    At each ASM site, each haplotype receives the methylation probability
    associated with its linked-SNV allele. Other sites retain their baseline
    methylation probabilities.

    See [Add allele-specific methylation](customize.md#allele-specific-methylation)
    and the [ASM input formats](../reference/formats.md#allele-specific-methylation-inputs).

## Enrichment-based assays

### RRBS { #rrbs }

Reduced representation bisulfite sequencing (RRBS) samples restriction-enzyme
fragments inside a retained size window. This example uses the default MspI
cut site `C|CGG` and uniform sampling.

```bash
bsreadsim run rrbs \
  -r test.fa \
  -o test/rrbs \
  -n 100000 \
  --cut-site 'C|CGG' \
  --insert-min 100 \
  --insert-max 500 \
  --mutation-rate 0 \
  -s 42
```

??? info "How it works"

    **Generation:** `--cut-site 'C|CGG'` defines MspI digestion. RRBS retains
    restriction fragments between `--insert-min 100` and
    `--insert-max 500`.

    **Sampling:** Because no score model is selected, retained fragments are
    sampled uniformly. For non-uniform sampling, export the exact candidates
    with `build rrbs`, score them with an external model, then load the scored
    file with `--sampling score` and `--rrbs-candidates`.

    See [RRBS fragment sampling](customize.md#rrbs) and the
    [RRBS candidate contract](../reference/formats.md#rrbs-candidate-bed).

### TBS { #tbs }

Targeted bisulfite sequencing (TBS) uses target regions from a strand-aware
BED6 file to represent probe enrichment. This example samples eligible targets
uniformly.

```bash
bsreadsim run tbs \
  -r test.fa \
  -o test/tbs \
  -n 100000 \
  --targets targets.bed \
  --insert-mean 300 \
  --center-sd 50 \
  --mutation-rate 0 \
  -s 42
```

??? info "How it works"

    **Generation:** `--targets` supplies the enriched regions.
    `--center-sd 50` controls fragment-center displacement, and
    `--insert-mean 300` sets the mean fragment length.

    **Sampling:** Targets are sampled uniformly by default. `--sampling score`
    instead uses BED column 5 as each target's relative weight.

    This model represents enrichment through target-centered fragment
    placement; it does not simulate individual probes or capture chemistry.
    See [Targeted fragment generation](customize.md#tbs),
    [fragment sampling](customize.md#sample-fragments), and the
    [target BED contract](../reference/formats.md#capture-target-bed).

## Ground truth in simulation

BSReadSim provides ground truth for controlled simulations. The saved truth
can support method development and benchmarking or be reused in later
simulations.

### Save

This example saves reads as origin-annotated BAM and exports the prepared
variants and methylation profile as reusable VCF and MethDB artifacts.

```bash
bsreadsim run wgbs \
  -r test.fa \
  -o runs/benchmark \
  -n 100000 \
  --mutation-rate 0.001 \
  --seed-mut 7 \
  --seed-meth 11 \
  -s 42 \
  --format bam \
  --save-truth
```

The output directory has this structure:

```text { .no-copy }
runs/benchmark/
├── sim.bam
├── sim.manifest.json
└── truth/
    ├── sim.variants.vcf.gz
    └── sim.methdb
```

??? info "How it works"

    `--format bam` replaces FASTQ with origin-annotated BAM. `--save-truth`
    writes the prepared, phased variant set and methylation profile as
    reusable VCF and MethDB artifacts.

    Add [`--fragment-realization`](customize.md#truth-artifacts) when
    complete-fragment methylation and conversion states are needed in BAM.

    See [Save simulation truth](customize.md#truth-artifacts),
    [Outputs](../outputs/index.md), and
    [Annotated BAM](../outputs/index.md#annotated-bam) for the artifact and tag
    contracts.

### Reuse

Load the saved MethDB to reuse the same methylation profile and embedded
variants in another WGBS run:

```bash
bsreadsim run wgbs \
  -r test.fa \
  -o test/wgbs-reuse \
  -n 100000 \
  --methdb runs/benchmark/truth/sim.methdb \
  -s 43
```

`--methdb` loads the saved snapshot directly. Changing the master seed draws
new fragments and reads while retaining the prepared methylation profile and
variant set. Use the saved VCF with `--vcf` instead when only the prepared
variant set should be reused. See
[Load methylation profiles](customize.md#predefined-methylation) for input
combination rules.
