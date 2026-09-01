# Other assays

Besides bisulfite sequencing, BSReadSim simulates whole-genome sequencing
(WGS), whole-exome sequencing (WES), and targeted sequencing (TS).

## Non-bisulfite workflow { #how-standard-sequencing-is-simulated }

The non-bisulfite workflow shares the genome, fragment, read, output, and seed
controls without methylation modeling or bisulfite conversion.

```text
reference FASTA
  -> prepared diploid genome
  -> fragment generation and sampling
  -> read extraction, base quality, and substitution errors
  -> FASTQ or annotated BAM + manifest
```

- **WGS:** Samples fragments from eligible positions across the genome and
  requires no target file.
- **WES:** Samples fragments around exome intervals supplied with `--targets`.
- **TS:** Samples fragments around custom panel intervals supplied with
  `--targets`.

WES, TS, and TBS use the same
[target-centered fragment model](customize.md#tbs).

## Choose an assay

### WGS { #wgs }

Use `wgs` for whole-genome sampling. This example uses the default uniform
sampling.

```bash
bsreadsim run wgs \
  -r reference.fa \
  -o runs/wgs \
  -n 100000 \
  --mutation-rate 0 \
  -s 42
```

See [Generate fragments](customize.md#generate-fragments) and
[Sample fragments](customize.md#sample-fragments) for whole-genome fragment
geometry and sampling controls.

### WES { #wes }

Use `wes` with exome intervals in a strand-aware BED6 file. This example uses
the default uniform target sampling.

```bash
bsreadsim run wes \
  -r reference.fa \
  -o runs/wes \
  -n 100000 \
  --targets exome.bed \
  --insert-mean 300 \
  --center-sd 50 \
  --mutation-rate 0 \
  -s 42
```

### TS { #ts }

Use `ts` with a custom target panel in a strand-aware BED6 file. This example
weights targets by the scores in BED column 5.

```bash
bsreadsim run ts \
  -r reference.fa \
  -o runs/ts \
  -n 100000 \
  --sampling score \
  --targets panel.bed \
  --insert-mean 300 \
  --center-sd 50 \
  --mutation-rate 0 \
  -s 42
```

## Shared simulation controls

The following controls also apply to non-bisulfite assays:

| Control | Use in non-bisulfite assays |
| --- | --- |
| [Genome](customize.md#genetic-variation) | Keep the reference unchanged, generate de novo variants, or load a diploid VCF |
| [Fragment generation](customize.md#generate-fragments) | Use whole-genome geometry for WGS and target-centered geometry for WES and TS |
| [Fragment length](customize.md#fragment-geometry) | Set the minimum, mean, maximum, and standard deviation |
| [Fragment sampling](customize.md#sample-fragments) | Use `uniform` by default, `gc` for WGS, or `score` for WES and TS |
| [Read count](customize.md#dataset-size) | Set an exact read count or assay-aware depth |
| [Read layout](customize.md#read-layout) | Choose single- or paired-end reads and set the read length |
| [Quality and errors](customize.md#quality-and-error) | Use fixed values or data-derived models for Phred scores and substitutions |
| [Output and truth](customize.md#reproducibility) | Write FASTQ or annotated BAM, a manifest, and optional variant truth |
| [Random seeds](customize.md#random-seeds) | Control genome preparation and fragment and read generation; the methylation seed is unused |

??? info "How depth is calculated"

    WGS depth uses eligible whole-genome sequence. WES and TS depth use the
    union of eligible target intervals. Use `--reads` instead when an exact
    output record count is required.

## Outputs and simulation truth

Each run writes FASTQ or an origin-annotated BAM, together with a manifest.
With the default prefix `sim`, a paired-end FASTQ run using `--save-truth` has
the following structure:

```text { .no-copy }
OUTPUT/
├── sim.R1.fastq.gz
├── sim.R2.fastq.gz
├── sim.manifest.json
└── truth/
    └── sim.variants.vcf.gz
```

Single-end FASTQ omits `sim.R2.fastq.gz`. BAM output replaces the FASTQ files
with `sim.bam` and omits methylation and conversion annotations. For
non-bisulfite assays, `--save-truth` exports only the prepared, phased variant
set as VCF.

See [Choose read output](customize.md#output-format),
[Save simulation truth](customize.md#truth-artifacts), and
[Outputs](../outputs/index.md) for the available representations, artifact
names, and completion rules.
