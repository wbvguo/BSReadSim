<div class="docs-hero" markdown>

<img class="docs-hero__logo" src="img/BSReadSim.png" alt="BSReadSim logo">

# BSReadSim

BSReadSim is a versatile simulator built for bisulfite sequencing. It combines
genetic variation, DNA methylation, assay-specific sampling, bisulfite
chemistry, and sequencing errors to produce realistic reads with traceable
ground truth. The same simulation engine also provides complementary support
for ordinary genomics assays.
Learn more in the [BSReadSim preprint](https://doi.org/10.1101/2024.12.24.627620).

<div class="docs-actions" markdown>

[Quick start](getting-started/quickstart.md){ .md-button .md-button--primary }
[Workflow](simulation/workflow.md){ .md-button }

</div>

</div>

## Supported technology assays

BSReadSim focuses on WGBS, RRBS, and TBS with bisulfite chemistry. The same
simulation engine also supports WGS, WES, and targeted sequencing when a study
needs ordinary sequencing reads or a non-bisulfite control.

<div class="technology-grid" markdown>

<div class="technology-card" markdown>

### WGBS

Whole-genome bisulfite sequencing profiles DNA methylation at single-base resolution across the genome.

[Configure WGBS](simulation/customize.md#wgbs)

</div>

<div class="technology-card" markdown>

### RRBS

Reduced representation bisulfite sequencing enriches CpG-rich fragments through restriction-enzyme digestion and size selection.

[Configure RRBS](simulation/customize.md#rrbs)

</div>

<div class="technology-card" markdown>

### TBS

Targeted bisulfite sequencing enriches predefined genomic regions through probe-based capture.

[Configure TBS](simulation/customize.md#tbs)

</div>

<div class="technology-card technology-card--secondary" markdown>

### Other assays

Additional support for ordinary whole-genome, whole-exome, and targeted sequencing.

[Learn more](simulation/other-assays.md)

</div>

</div>

## Installation

Install the current release from source on Linux or WSL2.
A Bioconda package is coming soon.

[Installation guide](getting-started/installation.md){ .md-button .md-button--primary }
[Supported platforms](getting-started/platforms.md){ .md-button }

## Customize the simulation

Configure each biological and technical layer to match your study.

<ol class="simulation-flow">
  <li><strong><a href="simulation/customize/#genetic-variation">Genetic variation</a></strong>
  Generate SNVs and indels or incorporate sample variants from a diploid VCF.</li>
  <li><strong><a href="simulation/customize/#methylation">DNA methylation</a></strong>
  For bisulfite runs, generate context-specific methylation or incorporate measured and allele-specific profiles.</li>
  <li><strong><a href="simulation/customize/#supported-technologies">Sequencing assay</a></strong>
  Choose WGBS, RRBS, or TBS and configure its fragment-sampling model. Standard sequencing is covered under <a href="simulation/other-assays/">Other assays</a>.</li>
  <li><strong><a href="simulation/customize/#library-and-sequencing">Library and sequencing</a></strong>
  Configure fragment length, optional bisulfite chemistry, base quality, and sequencing errors.</li>
</ol>

[Customize](simulation/customize.md){ .md-button }
[View tutorials](simulation/tutorials.md){ .md-button }

## Understand the output

Every successful run produces reads as either FASTQ or annotated BAM, together with a run manifest.

- **FASTQ** provides standard reads and quality scores for existing analysis pipelines.
- **Annotated BAM** maps reads to their simulated origins and records applicable methylation/conversion plus variant and sequencing-error ground truth.
- **Simulation truth artifacts** preserve the prepared variant set as phased VCF and, for bisulfite modes, the prepared methylation profile as MethDB.
- **Run manifest** is pretty JSON recording the received command, an expanded full command, inputs, effective settings, random seeds, software versions, counts, and checksums.

[Explore outputs](outputs/index.md){ .md-button }
[View format reference](reference/formats.md){ .md-button }
