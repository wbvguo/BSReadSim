<img src="docs/img/BSReadSim.png" alt="BSReadSim logo" width="150" align="right">

# BSReadSim

BSReadSim is a versatile and efficient read simulator for genomic sequencing,
supporting both conventional and bisulfite-based assays. It combines
configurable biological and technical models to produce realistic reads with
traceable ground truth. The resulting datasets can be used to guide
experimental design, develop bioinformatics tools, and benchmark their
performance under controlled conditions. Learn more in the
[BSReadSim preprint](https://doi.org/10.1101/2024.12.24.627620).

## Highlights

- Simulates WGBS, RRBS, and TBS together with matched WGS, WES, and TS controls
  in one framework.
- Maintains a coherent, haplotype-aware diploid genome and methylome, including
  variant phasing, short indels, and variant-induced methylation-context
  changes.
- Supports profile-based simulation that can incorporate user-supplied genetic
  variants, site-level methylation and allele-specific methylation,
  fragment-selection profiles, and empirical quality and error models.
- Models assay-specific library selection, including restriction-aware RRBS,
  capture-aware targeted assays, and GC- or score-weighted fragment sampling.
- Produces reproducible FASTQ or origin-annotated BAM, with reusable ground truth
  files and a manifest recording the effective configuration and output
  identities.
- High computational efficiency

## Installation and usage

For installation instructions, the quick start, simulation tutorials, CLI
options, file formats, outputs, and troubleshooting, see the complete user
guide:

[**wbvguo.github.io/BSReadSim/**](https://wbvguo.github.io/BSReadSim/)

## Citation

If you use BSReadSim in your research, please cite:

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

## License

BSReadSim is available under the MIT license.
