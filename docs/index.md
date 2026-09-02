<div class="docs-hero" markdown>

<img class="docs-hero__logo" src="img/BSReadSim.png" alt="BSReadSim logo">

# BSReadSim

BSReadSim is a versatile and efficient read simulator for genomic sequencing, supporting both conventional and bisulfite-based assays.
It combines configurable biological and technical models to produce realistic reads with traceable ground truth.
The resulting datasets can be used to guide experimental design, develop bioinformatics tools, and benchmark their performance under controlled conditions.
Learn more in the [BSReadSim preprint](https://doi.org/10.1101/2024.12.24.627620).

<div class="docs-actions" markdown>

[Quick start](getting-started/quickstart.md){ .md-button .md-button--primary }
[Workflow](simulation/workflow.md){ .md-button }

</div>

</div>

## Supported assays

<div class="technology-grid" markdown>

<div class="technology-card" markdown>

<p class="technology-card__title">WGBS</p>

Whole-genome bisulfite sequencing profiles DNA methylation at single-base resolution across the genome.

[Configure WGBS](simulation/customize.md#wgbs)

</div>

<div class="technology-card" markdown>

<p class="technology-card__title">RRBS</p>

Reduced representation bisulfite sequencing enriches CpG-rich regions through restriction-enzyme digestion and size selection.

[Configure RRBS](simulation/customize.md#rrbs)

</div>

<div class="technology-card" markdown>

<p class="technology-card__title">TBS</p>

Targeted bisulfite sequencing enriches predefined genomic regions through probe-based capture.

[Configure TBS](simulation/customize.md#tbs)

</div>

<div class="technology-card technology-card--secondary" markdown>

<p class="technology-card__title">Other assays</p>

Additional support includes conventional whole-genome, whole-exome, and targeted sequencing.

[Learn more](simulation/other-assays.md)

</div>

</div>

## Installation

Install the current release on Linux or WSL2.

[Installation guide](getting-started/installation.md){ .md-button .md-button--primary }

## Customization

Tailor models and parameters to your study in **Customize**. For complete commands organized by simulation goal, see **Tutorials**.

[Customize](simulation/customize.md){ .md-button }
[View tutorials](simulation/tutorials.md){ .md-button }

## Output

Each run produces reads as FASTQ files or an origin-annotated BAM, together with a manifest. The underlying variant set and methylation profile can optionally be saved as VCF and MethDB files for reuse.

[Explore outputs](outputs/index.md){ .md-button }
[View file formats](reference/formats.md){ .md-button }

## Citation

If you use BSReadSim in your research, please cite the [BSReadSim preprint](https://doi.org/10.1101/2024.12.24.627620):

```bibtex
@article{guo2024bsreadsim,
  title = {BSReadSim: a versatile and efficient simulator to generate realistic bisulfite sequencing reads},
  author = {Guo, Wenbin and Pellegrini, Matteo},
  journal = {bioRxiv},
  year = {2024},
  publisher = {Cold Spring Harbor Laboratory},
  doi = {10.1101/2024.12.24.627620}
}
```
