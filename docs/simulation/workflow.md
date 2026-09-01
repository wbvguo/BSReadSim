# Workflow

BSReadSim generates sequencing reads through a workflow that mirrors the key biological and technical stages of an experiment. It first prepares a diploid genome, then applies assay-specific fragment sampling. Bisulfite modes additionally construct site-specific methylation probabilities, draw methylation states, and apply bisulfite conversion. Finally, it adds base qualities and sequencing errors before exporting reads as FASTQ or origin-annotated BAM.


## From reference genome to sequencing reads

<figure class="workflow-figure">
  <a href="../../img/BSReadSim_workflow.png" title="Open the full-resolution BSReadSim workflow">
    <img src="../../img/BSReadSim_workflow.png"
         alt="BSReadSim workflow from a reference genome through variants, haplotypes, optional methylation and bisulfite conversion, assay-specific sampling, sequencing errors and read output">
  </a>
  <figcaption>BSReadSim resolves the genome and methylome before sampling physical fragments, generating reads, and writing FASTQ or origin-annotated BAM with a run manifest.</figcaption>
</figure>

## Five stages of read simulation

| Stage            | Description                                                                                                                               |
| ---------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| Genome           | Resolve two haplotypes from the reference genome, optionally incorporating VCF variants or de novo mutations. |
| Methylome        | For bisulfite modes, prepare methylation probabilities using beta distributions and any provided measured or allele-specific profiles. |
| Fragmentation    | Generate physical DNA fragments across whole-genome, restriction-fragment, or targeted domains using uniform or profile-based sampling. |
| Sequencing       | Sequence fragments according to the selected assay protocol, then assign base qualities and introduce sequencing errors. |
| Output           | Output reads as FASTQ or annotated BAM, with a provenance manifest. |

The order matters because variants can change cytosine contexts, methylatable
sites, RRBS cut sites, and fragment opportunities. Each upstream layer is
therefore resolved before the next one.

<!-- !!! note "Current model scope"
    See [Customize](customize.md) for model assumptions and assay-specific limitations. -->

## Six supported sequencing assays

BSReadSim supports three bisulfite assays (WGBS, RRBS, and TBS) and three non-bisulfite assays (WGS, WES, and TS), spanning whole-genome and region-enrichment protocols. The figure below summarizes the experimental strategies used in the three bisulfite assays.

<figure class="workflow-figure">
  <a href="../../img/BS_seqtech.png" title="Open the full-resolution sequencing technology comparison">
    <img src="../../img/BS_seqtech.png"
         alt="Comparison of WGBS genome-wide sampling, RRBS restriction-enzyme digestion, and TBS BED-targeted sampling">
  </a>
  <figcaption>
    WGBS uses random fragmentation for genome-wide coverage; RRBS enriches CpG-rich regions through restriction-enzyme digestion and size selection; and TBS enriches predefined targets through probe enrichment.
  </figcaption>
</figure>

See [Tutorials](tutorials.md) for runnable examples, [Customize](customize.md) for model options, and [Other assays](other-assays.md) for non-bisulfite workflows.