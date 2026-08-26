# Inspect outputs

Every successful run publishes its manifest last. Treat that file as the
run-level record of effective settings, seeds, input identities, summary, and
output checksums.

## Read the manifest

```bash
python -m json.tool runs/wgbs/sim.manifest.json | less
```

Keep the manifest with the reads. It records generated seeds in the
`command.full_command` and `details.randomness` even when no seed
was supplied on the command line.

## Preview FASTQ

```bash
gzip -dc runs/wgbs/sim.R1.fastq.gz | head -n 8
gzip -dc runs/wgbs/sim.R2.fastq.gz | head -n 8
```

Two FASTQ records should share the same fragment identity in their read names.
The exact naming format is described in [Read names](../reference/formats/read-name.md).

## Inspect annotated BAM

Request BAM with `--format bam`. BSReadSim then emits `sim.bam` instead of the FASTQ
pair. The stream is deliberately unsorted:

```bash
samtools view -H runs/example-bam/sim.bam
samtools view runs/example-bam/sim.bam | head
samtools sort -o runs/example-bam/sim.sorted.bam runs/example-bam/sim.bam
samtools index runs/example-bam/sim.sorted.bam
```

Every record contains `zt:Z` base-level state and `zr:B:S` read summaries.
`--fragment-summary` adds `zf:B:S`; `--fragment-realization` adds `zx:Z` and
implies fragment summaries. See [Annotated BAM](../reference/formats/bam.md) for the bit
and field definitions.

FASTQ can be recovered from paired BAM records:

```bash
samtools fastq \
  -1 recovered.R1.fastq \
  -2 recovered.R2.fastq \
  -0 /dev/null \
  -s /dev/null \
  -n runs/example-bam/sim.bam
```

Use the original manifest, rather than regenerated FASTQ metadata, for the
authoritative configuration and artifact checksums.
