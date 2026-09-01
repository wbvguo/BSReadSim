# Quick start

## 1. Prepare a reference genome

Obtain a reference genome in FASTA format. This demo uses the repository's
bundled [test.fa](https://github.com/wbvguo/BSReadSim/blob/main/data/example/test.fa).
For real analyses, prepare or download a reference genome from
[Ensembl](https://www.ensembl.org/info/data/ftp/index.html?redirect=no) or
[GENCODE](https://www.gencodegenes.org/) that matches the species and build in your study.

## 2. Run a WGBS simulation

```bash
bsreadsim run wgbs -r test.fa -o test/ -n 1000 -s 42
```

??? info "Command options"
    `run wgbs` simulates WGBS reads. In this command:

    - `-r` specifies the reference FASTA;
    - `-o` specifies the output directory;
    - `-n` specifies the total number of reads to generate (n/2 read pairs);
    - `-s` specifies the master seed;

    All other parameters use the WGBS defaults. See
    [Customize](../simulation/customize.md) for more configuration options.

## 3. Inspect the output

This command uses the default paired-end mode and `fastq.gz` format, so the
output directory contains two FASTQ files and a run manifest:

```text
test/
├── sim.R1.fastq.gz
├── sim.R2.fastq.gz
└── sim.manifest.json
```

Preview the first R1 record:

```bash
zcat test/sim.R1.fastq.gz | head -n 4
```

The manifest records the effective configuration, provenance, and run summary:

```bash
head -n 40 test/sim.manifest.json
```

## Continue with customized simulation

Choose [WGBS](../simulation/customize.md#wgbs),
[RRBS](../simulation/customize.md#rrbs),
[TBS](../simulation/customize.md#tbs), or a
[non-bisulfite assay](../simulation/other-assays.md). Then check the following

- [Tutorials](../simulation/tutorials.md) for complete task-oriented commands;
- [Customize](../simulation/customize.md) to choose models and parameters;
- [CLI parameters](../reference/cli.md) for default and optional configurations.
