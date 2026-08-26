# Tutorials

These tutorials are organized around common simulation goals. Start with
Synthetic WGBS for a generated baseline, then progressively incorporate a
variant set, a methylation profile, or allele-specific methylation. The
tutorials on RRBS/TBS and ground truth can be followed independently.

| Goal | Tutorial |
| --- | --- |
| Generate WGBS data without measured profiles | [Synthetic WGBS](#wgbs) |
| Incorporate measured variants and methylation | [Variant sets and methylation profiles](#variant-sets-and-methylation-profiles) |
| Add haplotype-specific methylation | [Allele-specific methylation](#allele-specific-methylation) |
| Simulate restriction or capture enrichment | [RRBS and TBS](#rrbs-and-tbs) |
| Generate reads with detailed truth | [Ground-truth benchmarking](#ground-truth-benchmarking) |

Each section provides a complete command template, its required inputs, and
the expected result. Replace example paths with inputs from the same reference
assembly, and use [Customize](customize.md) for the full set of model and
parameter choices.

## Synthetic WGBS { #wgbs }

**When to use.** Generate a WGBS baseline when no measured genetic or
methylation profile is available. BSReadSim introduces variants at the
requested mutation rate and generates context-specific methylation
probabilities from its default Beta distributions.

**Inputs.** A reference FASTA. No VCF or methylation profile is required.

**Methylation pattern.** Without a methylation input, BSReadSim draws a
probability for every eligible cytosine from the Beta distribution assigned to
its CG, CHG, or CHH context. These probabilities form the MethDB used to
realize methylation on sampled fragments. `--seed-meth` makes the generated
MethDB reproducible.

Adjust the three distributions to create a different synthetic pattern:

| Pattern | `--beta-cg` | `--beta-chg` | `--beta-chh` |
| --- | --- | --- | --- |
| High CpG, low non-CpG | `8,2` | `1,19` | `1,19` |
| Low CpG, low non-CpG | `2,8` | `1,19` | `1,19` |

The mean of `Beta(a, b)` is `a / (a + b)`, so the two CpG distributions are
centered at `0.8` and `0.2`. Probabilities are generated independently at
individual sites in the current model.

**Run.**

```bash
bsreadsim run wgbs \
  -r reference.fa \
  -o runs/synthetic-wgbs \
  -n 100000 \
  --mutation-rate 0.001 \
  --seed-mut 7 \
  --seed-meth 11 \
  -s 42
```

**Expected result.** `runs/synthetic-wgbs/` contains paired compressed FASTQ
and `sim.manifest.json`. The manifest records the generated variant and
methylation settings plus the effective seeds. Add `--save-methdb` when the
generated MethDB should also be published beneath `truth/`.

**Customize further.** Adjust `--mutation-rate` to control synthetic genetic
variation, or set it to `0` when no variants are desired. Keep `--seed-mut`
and `--seed-meth` fixed to reuse the same genome and methylome while changing
`-s` to draw different reads. See [Genome](customize.md#genetic-variation) and
[Methylome](customize.md#methylation) for the available models.

## Variant sets and methylation profiles

**When to use.** Incorporate a measured variant set and methylation
probabilities into an input-informed WGBS simulation.

**Inputs.** A reference FASTA, a one-sample diploid VCF variant set, and either
a bedMethyl or CGmap methylation profile from the same assembly. A VCF
automatically disables de novo mutation generation.

**Run.**

```bash
bsreadsim run wgbs \
  -r reference.fa \
  -o runs/profile-based \
  -n 100000 \
  --vcf sample.vcf.gz \
  --bed-methyl sample.bed.gz \
  --seed-phase 7 \
  --seed-meth 11 \
  -s 42
```

**Expected result.** `runs/profile-based/` contains paired compressed FASTQ
and a manifest identifying both measured inputs. Profile positions use their
supplied probabilities; missing positions use the context-specific generated
fallback.

**Customize further.** Use `--cgmap` instead of `--bed-methyl`, or enable
context pooling when an empirical distribution matters more than individual
loci. See the [VCF](../reference/formats/vcf-variants.md),
[bedMethyl](../reference/formats/bed-methyl-profile.md), and
[CGmap](../reference/formats/cgmap-profile.md) contracts.

## Allele-specific methylation

**When to use.** Simulate selected sites with different methylation
probabilities on the REF and ALT haplotypes.

**Inputs.** The reference, diploid VCF, and baseline methylation profile from
the profile-based scenario, plus an ASM or ASM BED file. Every ASM row must
link to a heterozygous SNV in the VCF.

**Run.**

```bash
bsreadsim run wgbs \
  -r reference.fa \
  -o runs/asm \
  -n 100000 \
  --vcf sample.vcf.gz \
  --bed-methyl sample.bed.gz \
  --asm sample.asm.gz \
  --seed-phase 7 \
  --seed-meth 11 \
  -s 42
```

**Expected result.** `runs/asm/` contains reads generated from the
profile-informed methylome with haplotype-specific probabilities overlaid at
linked sites. The manifest records the VCF, baseline profile, and ASM input.

**Customize further.** Use `--asm-bed` for the BED representation. ASM
overrides the baseline profile at linked sites. See the
[ASM](../reference/formats/asm-profile.md) and
[ASM BED](../reference/formats/asm-bed-profile.md) contracts.

## RRBS and TBS

### RRBS { #rrbs }

**When to use.** Simulate a library enriched for restriction-enzyme fragments
within a selected size range.

**Inputs.** A reference FASTA and one or more cut motifs. This baseline uses
MspI `C|CGG` and uniform sampling over eligible fragments.

**Run.**

```bash
bsreadsim run rrbs \
  -r reference.fa \
  -o runs/rrbs-uniform \
  -n 100000 \
  --cut-site 'C|CGG' \
  --insert-min 100 \
  --insert-mean 250 \
  --insert-max 500 \
  --mutation-rate 0 \
  -s 42
```

**Expected result.** `runs/rrbs-uniform/` contains reads sampled from
MspI-bounded fragments between 100 and 500 bases, plus the run manifest.

**Customize further.** Add motifs for a multi-enzyme library or use externally
scored candidates. Score-based runs must share genome, variant, motif, and
fragment-boundary settings with the candidate build. See
[Customize](customize.md#rrbs) and the
[RRBS candidate format](../reference/formats/rrbs-candidate-bed.md).

### TBS { #tbs }

**When to use.** Simulate a bisulfite library enriched for predefined genomic
targets.

**Inputs.** A reference FASTA and a strand-aware BED6 target file. This
baseline assigns equal sampling mass to every eligible target.

**Run.**

```bash
bsreadsim run tbs \
  -r reference.fa \
  -o runs/tbs-uniform \
  -n 100000 \
  --targets targets.bed \
  --insert-mean 300 \
  --insert-sd 0 \
  --mutation-rate 0 \
  -s 42
```

**Expected result.** `runs/tbs-uniform/` contains reads sampled around the BED
targets with a fixed 300-base fragment length, plus the run manifest.

**Customize further.** Use score-based sampling when BED column 5 should
supply relative capture weights. Scores are relative and do not need to sum to
one. See [Customize](customize.md#tbs), the
[TBS target format](../reference/formats/tbs-catalog.md), and the
[target-score sampling guide](../reference/advanced/target-score-sampling.md).

## Ground-truth benchmarking

**When to use.** Generate reads with origin, methylation, conversion, variant,
and sequencing-error truth for method development or accuracy evaluation.

**Inputs.** A reference FASTA. This example generates variants and a
methylome internally; measured inputs can be substituted from the earlier
profile-based scenarios.

**Run.**

```bash
bsreadsim run wgbs \
  -r reference.fa \
  -o runs/benchmark \
  -n 100000 \
  --mutation-rate 0.001 \
  --seed-mut 7 \
  --seed-meth 11 \
  -s 42 \
  --format bam \
  --save-truth
```

**Expected result.** `runs/benchmark/` contains `sim.bam`,
`sim.manifest.json`, `truth/sim.variants.vcf.gz`, and `truth/sim.methdb`. BAM
annotations retain read-level truth, while the VCF and MethDB preserve the
biological models used by the run.

**Customize further.** Add fragment summaries or complete-fragment
realizations for more detailed BAM truth. See
[Output](customize.md#reproducibility), [Outputs](../outputs/index.md), and
[Annotated BAM](../reference/formats/bam.md) for the exact artifacts and tags.
