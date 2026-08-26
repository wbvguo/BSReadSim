# Customize

Customize follows the five stages introduced in the [Workflow](workflow.md).
BSReadSim provides complete defaults, so change only the components needed to represent your study.
This page covers bisulfite sequencing workflows, which traverse all five
stages. Standard sequencing workflows are documented separately under
[Other assays](other-assays.md).

## Genome { #genetic-variation }

Choose the genetic background from which fragments will be sampled.
Genetic variation is optional: keep the reference unchanged, introduce random variants, or incorporate variants from a diploid VCF.

| Goal | Strategy | Key controls |
| --- | --- | --- |
| Reference-only baseline | Keep the reference unchanged | `--mutation-rate 0` |
| Synthetic diploid genome | Introduce random SNVs and indels | `--mutation-rate`, `--indel-fraction`, `--indel-extension-probability`, `--seed-mut` |
| Sample-specific genome | Incorporate a one-sample diploid VCF | `--vcf`, `--seed-phase` |

Set the mutation policy explicitly for reference-only and synthetic simulations.
`--vcf` selects the predefined-variant path, automatically disables de novo
mutation generation, and is mutually exclusive with `--mutation-rate`.
Variants are resolved into two haplotypes before methylation contexts and assay-specific fragments are determined.

??? info "Reusable variant set"

    A normal run constructs the prepared variant set directly; add `--save-vcf` to publish that set as a phased VCF truth artifact.
    Build the VCF independently only when several simulations or external tools must share the same prepared variant set:

    ```bash
    bsreadsim build variants \
      -r reference.fa \
      -o variants.vcf.gz \
      --mutation-rate 0.001 \
      --seed-mut 7
    ```

    See the [VCF contract](../reference/formats/vcf-variants.md) for supported variants and genotypes.

## Methylome { #methylation }

Choose how methylation probabilities are assigned to eligible cytosines.
BSReadSim stores the prepared methylome as MethDB, then draws methylated or unmethylated states for each sampled fragment.

### Probability sources

| Goal | Strategy | Key controls |
| --- | --- | --- |
| Simulate without measured data | Generate context-specific probabilities | `--beta-cg`, `--beta-chg`, `--beta-chh`, `--cpg-only`, `--seed-meth` |
| Preserve measured loci | Incorporate a position-specific profile | `--cgmap` or `--bed-methyl` |
| Reproduce an empirical distribution | Sample context-specific pools from a profile | profile input plus `--cgmap-pool`, `--seed-meth` |
| Reuse a fixed methylome | Load a prepared MethDB | `--methdb` |

Choose only one baseline input: CGmap, bedMethyl, or MethDB.
Missing profile positions use context-specific generated probabilities.

### Methylation pattern

For a generated methylome, BSReadSim draws a probability for each eligible
cytosine from the Beta distribution assigned to its CG, CHG, or CHH context.
Together these site probabilities form MethDB. Position-specific profiles can
instead vary probability by locus, while context pooling preserves an
empirical distribution without preserving individual loci.

The current model generates site probabilities and realizes methylation states
independently. Probability assignment and state realization remain separate
model boundaries so sequence-dependent pattern models can be introduced
without changing the profile-input contracts.

### Allele-specific methylation

Use `--asm` or `--asm-bed` with a diploid VCF to overlay haplotype-specific
probabilities at linked sites. Allele-specific methylation overrides the
baseline at those sites, giving the effective precedence
`ASM > position-specific or pooled profile > generated fallback`.

??? info "Reusable methylation profile"

    A normal run constructs the prepared methylation profile in memory; add `--save-methdb` to publish it as a MethDB truth artifact, or use `--methdb` to load a saved profile.
    Build the profile independently when it must be inspected or distributed before read simulation:

    ```bash
    bsreadsim build methdb \
      -r reference.fa \
      -o sample.methdb \
      --mutation-rate 0 \
      --seed-meth 7
    ```

    See the [methylation input formats](../reference/formats.md#format-map) and [MethDB format](../reference/formats/methdb.md) for exact contracts.

## Fragmentation { #supported-technologies }

Choose the sequencing assay, then configure how many physical fragments are sampled and which fragment lengths are retained.

### WGBS { #wgbs }

Whole Genome Bisulfite Sequencing requires no assay-specific input beyond the reference FASTA.
Use `--gc-profile` when fragments should follow a target GC distribution rather than uniform opportunity sampling.

### RRBS { #rrbs }

Reduced Representation Bisulfite Sequencing requires one or more `--cut-site` motifs and retains enzyme-bounded fragments within the configured size range.
Use `--rrbs-candidates` with `--sampling score` when an external model supplies fragment weights.

### TBS { #tbs }

Targeted Bisulfite Sequencing requires strand-aware BED6 `--targets` and a fixed fragment length.
Use `--sampling score` for target-specific capture weights and `--fragment-center-stddev` to control placement around each target.

| Decision | Key controls |
| --- | --- |
| Number of source fragments | `--fragments` or `--depth` |
| Fragment length distribution | `--insert-min`, `--insert-mean`, `--insert-max`, `--insert-sd` |
| Ambiguous-base filtering | `--max-ambiguous-fraction` |

For fixed WGBS or TBS fragments, set
`--insert-mean N --insert-sd 0`.
RRBS fragment lengths are determined by restriction sites, while `--insert-min` and `--insert-max` define the retained size window.

When `--depth D` is used, the number of source fragments is:

```text
ceil(effective_reference_bases * D / (read_length * emitted_mates))
```

The effective region is eligible contig sequence for WGBS, the union of
eligible restriction-fragment envelopes for RRBS, and the union of target
intervals for TBS.

??? info "Reusable RRBS candidates"

    `bsreadsim build rrbs` exports the exact candidate set for external scoring or reuse.
    The build and run steps must use the same genome, variant, motif, and fragment-boundary settings.
    See the [RRBS tutorial](tutorials.md#rrbs) for the scored-candidate workflow.

## Sequencing { #library-and-sequencing }

Configure how each sampled fragment becomes one or two sequencing reads.
BSReadSim draws fragment-level methylation states and applies bisulfite
conversion before extracting mates. Every mode assigns base qualities and
introduces sequencing substitutions.

| Component | Choices | Key controls |
| --- | --- | --- |
| Read layout | Paired-end or single-end, with configurable read length | `--single-end`, `--read-length` |
| Bisulfite conversion | Directional or undirectional libraries in WGBS, RRBS, and TBS | `--conversion-rate`, `--undirectional` |
| Base quality | Fixed Phred score or empirical quality model | `--phred` or `--quality-model` |
| Sequencing error | Uniform substitution rate or quality-dependent confusion model | `--error-rate` or `--error-model` |

Methylation and conversion are realized on the complete physical fragment
before mate extraction, so overlapping mates share the same underlying event.
An empirical quality model replaces `--phred`, and an empirical error model replaces `--error-rate`.
See the [sequencing-model contracts](../reference/formats/sequencing-models.md) for the required files.

## Output { #reproducibility }

Choose the read representation and any reusable simulation truth artifacts.
Every successful run writes a manifest automatically; no additional option is required.

| Goal | Published files | Key controls |
| --- | --- | --- |
| Run an existing analysis pipeline | Compressed or uncompressed FASTQ | `--format fastq.gz` or `--format fastq` |
| Benchmark against read-level truth | Annotated BAM | `--format bam`, `--fragment-summary`, `--fragment-realization` |
| Preserve simulation truth | Saved variant set and, for bisulfite modes, saved methylation profile beneath `truth/` | `--save-vcf`, `--save-methdb`, `--save-truth` |
| Control output names | One read representation plus the manifest | `--output`, `--prefix` |

BAM replaces FASTQ rather than accompanying it.
The manifest is published last and records both the received command and a
full command expanded from the effective configuration, together
with input identities, seeds, versions, a data summary, output paths, and checksums.

| Seed | Controls |
| --- | --- |
| `--seed-mut` | Random variant generation |
| `--seed-phase` | Assignment of unphased VCF heterozygotes |
| `--seed-meth` | Generated methylation probabilities and empirical-pool draws in bisulfite modes |
| `-s`, `--seed` | Fragment selection, quality, and error; also methylation and conversion in bisulfite modes |

Fix every relevant seed for an identical rerun.
Keep the focused biological seeds fixed while changing only `--seed` to draw
different reads from the same prepared genome and, when applicable, methylome.

??? note "Performance controls"

    `--core-workers` controls fragment generation, `--workers` controls read processing, and `--chunk-size` plus `--max-in-flight-fragments` bound work batches and queued fragments.
    These options change resource use, not fixed-input, fixed-seed simulation semantics.

Use the [CLI reference](../reference/cli.md) for complete option types, defaults, and compatibility rules.
See [Tutorials](tutorials.md) for complete simulation recipes and [Outputs](../outputs/index.md) for the published file contracts.
