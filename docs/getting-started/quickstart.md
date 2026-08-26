# Quick start

## 1. Prepare a reference genome

For this demo, use the toy reference bundled with BSReadSim:

```bash
bsreadsim export test-fasta -o test.fa
```

For real applications, prepare a reference genome FASTA for the species and assembly used by your study; genome sequences are available from [Ensembl](https://www.ensembl.org/info/data/ftp/index.html?redirect=no) and [GENCODE](https://www.gencodegenes.org/).

## 2. Run a WGBS simulation

```bash
bsreadsim run wgbs -r test.fa -o test/ -n 1000 -s 42
```

- `run wgbs` simulates Whole Genome Bisulfite Sequencing reads.
- `-r` specifies the reference genome in FASTA format.
- `-o` specifies the output directory for the simulated dataset.
- `-n` specifies the number of source DNA fragments to sample.
- `-s` fixes the master seed for fragment selection and read realization.

All unspecified fragment, read, sequencing, and output options use the WGBS defaults.

## 3. Inspect the output

The output directory contains paired FASTQ files and a run manifest:

```text
test/
├── sim.R1.fastq.gz
├── sim.R2.fastq.gz
└── sim.manifest.json
```

Each FASTQ record contains a read name, simulated sequence, and quality scores. Preview the first R1 record:

```bash
zcat test/sim.R1.fastq.gz | head -n 4
```

The manifest records the effective configuration, seeds, inputs, counts, and checksums. Preview its first 40 lines:

```bash
head -n 40 test/sim.manifest.json
```

## Continue with customized simulation

Choose the technology that matches your experiment:

- [Whole Genome Bisulfite Sequencing (WGBS)](../simulation/customize.md#wgbs): genome-wide methylation profiling
- [Reduced Representation Bisulfite Sequencing (RRBS)](../simulation/customize.md#rrbs): restriction-enzyme enrichment of CpG-rich fragments
- [Targeted Bisulfite Sequencing (TBS)](../simulation/customize.md#tbs): probe-based enrichment of selected regions
- [Other genomic assays](../simulation/other-assays.md): ordinary whole-genome, whole-exome, or panel-enriched targeted sequencing

Open [Customize](../simulation/customize.md) to control genetic variation, methylation, fragment sampling, sequencing behavior, and output format.
See [Tutorials](../simulation/tutorials.md) for bisulfite simulation recipes.
See [Outputs](../outputs/index.md) for descriptions of the generated files.
