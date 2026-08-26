# Workflow

BSReadSim generates both bisulfite and ordinary sequencing reads through a workflow that mirrors the key biological and technical stages of an experiment. It first prepares the diploid genome, then applies assay-specific fragment sampling. Bisulfite modes additionally prepare a methylation landscape, realize fragment methylation, and apply conversion; standard modes skip that branch. Every mode finishes with modeled base qualities and sequencing errors.
The resulting reads are exported as FASTQ or annotated BAM, together with a manifest that records the inputs, settings, seeds, and checksums required for reproducibility.

## From reference genome to sequencing reads

<figure class="workflow-figure">
  <a href="../../img/BSReadSim_workflow.png" title="Open the full-resolution BSReadSim workflow">
    <img src="../../img/BSReadSim_workflow.png"
         alt="BSReadSim workflow from a reference genome through variants, haplotypes, optional methylation and bisulfite conversion, assay-specific sampling, sequencing errors and read output">
  </a>
  <figcaption>BSReadSim resolves the genome and methylome before sampling physical fragments, realizing reads, and publishing FASTQ or annotated BAM with a run manifest.</figcaption>
</figure>

```text
reference FASTA
  -> diploid haplotypes
  -> choose an assay
     -> WGBS/RRBS/TBS: methylome + bisulfite chemistry
     -> WGS/WES/TS: ordinary fragment sampling
  -> base quality and sequencing error
  -> FASTQ or annotated BAM + run manifest
```

## Five stages of read simulation

| Stage            | Description                                                                                                                               |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Genome           | Construct a prepared variant set, then resolve two haplotypes from the reference.                                                    |
| Methylome        | For bisulfite modes, construct a prepared methylation profile from generated, measured, and allele-specific probabilities.             |
| Fragmentation    | Sample physical DNA fragments using whole-genome, restriction-fragment, or BED-targeted capture models.                                   |
| Sequencing       | In bisulfite modes, realize methylation and conversion; in every mode, add base qualities and sequencing errors.                          |
| Output           | Write reads as FASTQ or annotated BAM and record the simulation in a provenance manifest.                                                 |

This order reflects the underlying data-generating process: the genome defines the fragments available to the assay and, for bisulfite runs, the methylation contexts; those fragments determine the reads that are ultimately observed.
It also preserves variant-dependent edge cases: genetic variants can create or disrupt cytosine contexts, introduce methylatable cytosines, alter RRBS cut sites, and shift fragment sampling probabilities.
BSReadSim therefore resolves each upstream layer before simulating the next.

## Six supported sequencing assays

BSReadSim supports three bisulfite assays and their ordinary whole-genome or
targeted counterparts. WES and TS share the BED-targeted capture algorithm;
their distinct identities make the intended experiment explicit in the
protocol and manifest.

<figure class="workflow-figure">
  <a href="../../img/BS_seqtech.png" title="Open the full-resolution sequencing technology comparison">
    <img src="../../img/BS_seqtech.png"
         alt="Comparison of WGBS random fragmentation, RRBS restriction-enzyme digestion, and TBS probe-based enrichment">
  </a>
  <figcaption>
    WGBS samples broadly across the genome, RRBS selects restriction-enzyme fragments, and TBS enriches probe-captured regions.
  </figcaption>
</figure>

| Technology                                         | Key characteristics                                                  |
| -------------------------------------------------- | -------------------------------------------------------------------- |
| Whole Genome Bisulfite Sequencing (WGBS)           | Provides genome-wide methylation profiling at single-base resolution |
| Reduced Representation Bisulfite Sequencing (RRBS) | Enriches CpG-rich regions through enzyme digestion and size selection |
| Targeted Bisulfite Sequencing (TBS)                | Profiles predefined genomic regions through probe-based capture      |
| Whole Genome Sequencing (WGS)                      | Samples ordinary reads across the genome without bisulfite chemistry |
| Whole Exome Sequencing (WES)                       | Captures exon intervals supplied as a target BED                      |
| Targeted Sequencing (TS)                           | Captures a general panel supplied as a target BED                     |

See [Other assays](other-assays.md) for the WGS, WES, and TS
simulation paths and complete starter commands.
Continue with a study-oriented bisulfite example from [Tutorials](tutorials.md), then use [Customize](customize.md) to control each stage.
See [Outputs](../outputs/index.md) for the published files and simulation truth artifacts.
