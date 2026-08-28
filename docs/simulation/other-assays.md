# Other assays

BSReadSim focuses on bisulfite sequencing, but the same genome, fragment, and
sequencing models can also generate standard sequencing reads. This page
covers the three standard assays available alongside the bisulfite workflows:
Whole Genome Sequencing (WGS), Whole Exome Sequencing (WES), and Targeted
Sequencing (TS).

“Other assays” is a documentation grouping, not a fourth command. Choose
`wgs`, `wes`, or `ts` explicitly so the selected technology is preserved in
the run manifest and annotated BAM.

## How standard sequencing is simulated

All three assays retain the parts of the BSReadSim workflow that do not depend
on bisulfite chemistry:

```text
reference FASTA
  -> reference-only, de novo, or VCF-derived diploid genome
  -> WGS whole-genome sampling or WES/TS BED-targeted capture
  -> single- or paired-end read extraction
  -> base quality and sequencing error
  -> FASTQ or annotated BAM + run manifest
```

Unlike WGBS, RRBS, and TBS, these assays do not construct a methylome, realize
cytosine methylation, or apply bisulfite conversion. They omit those stages
entirely rather than approximating standard sequencing with a zero conversion
rate. Consequently, methylation profiles, MethDB, ASM, and bisulfite-specific
controls are not accepted by `wgs`, `wes`, or `ts`.

| Assay | Fragment domain | Assay input | Related bisulfite mode |
| --- | --- | --- | --- |
| **WGS** | eligible starts across the genome | reference FASTA; optional target-GC profile | WGBS uses the same whole-genome sampler |
| **WES** | exon or other exome intervals | strand-aware BED6 target file | TBS uses the same capture sampler with bisulfite chemistry |
| **TS** | intervals in a general capture panel | strand-aware BED6 target file | TBS uses the same capture sampler with bisulfite chemistry |

## Choose an assay

### WGS { #wgs }

Use WGS when reads should originate across the eligible genome without
bisulfite chemistry. Fragment starts are sampled uniformly by default. An
optional `--gc-profile` requests a target fragment-GC distribution, and WGS
supports either fixed or clamped-normal variable insert lengths.
When a GC profile is combined with variable inserts, the run must currently be
reference-only; fixed-insert GC-profile runs also support VCF and de novo
variants.

```bash
bsreadsim run wgs \
  -r reference.fa \
  -o runs/wgs \
  -n 100000 \
  --mutation-rate 0 \
  -s 42
```

### WES { #wes }

Use WES for exome capture. Supply exon or other exome intervals as a
strand-aware BED6 file. The simulator selects an eligible target and places a
fixed-length fragment around its center; `--fragment-center-stddev` controls
the placement spread. Uniform sampling gives every eligible BED row equal
mass, while `--sampling score` uses BED column 5 as a relative capture weight.

```bash
bsreadsim run wes \
  -r reference.fa \
  -o runs/wes \
  -n 100000 \
  --targets exome.bed \
  --sampling uniform \
  --insert-mean 300 \
  --insert-sd 0 \
  --mutation-rate 0 \
  -s 42
```

### TS { #ts }

Use TS for a general targeted panel, such as a disease, validation, or custom
capture panel. It uses the same target eligibility, scoring, and placement
model as WES. The separate `ts` command records the experiment as Targeted
Sequencing rather than Whole Exome Sequencing.

```bash
bsreadsim run ts \
  -r reference.fa \
  -o runs/ts \
  -n 100000 \
  --targets panel.bed \
  --sampling score \
  --insert-mean 300 \
  --insert-sd 0 \
  --fragment-center-stddev 50 \
  --mutation-rate 0 \
  -s 42
```

With score-based sampling, BED scores are relative weights and do not need to
sum to one.

## Shared simulation controls

WGS, WES, and TS share the standard parts of the simulation model. This makes
them useful as non-bisulfite controls and for benchmarking tools that operate
on ordinary sequencing reads.

| Layer | Available choices |
| --- | --- |
| Genetic background | unchanged reference, de novo SNVs and indels, or a one-sample diploid VCF |
| Dataset size | exact total read count with `--reads`, or assay-aware mean depth with `--depth` |
| Read layout | single-end or paired-end reads with a configurable read length |
| Sequencing | fixed or empirical base quality and uniform or model-based substitution errors |
| Reproducibility | independent variant, phasing, and master seeds plus an automatic run manifest |

For WGS, depth is calculated over contigs with eligible whole-genome fragment
starts. For WES and TS, it is calculated over the union of eligible target
intervals. Use `--reads` when the exact output record count is more important
than assay-aware depth.

## Outputs and ground truth

Every run emits either FASTQ or annotated BAM plus a manifest. Annotated BAM
preserves fragment origin, variants, and sequencing-error truth where
applicable, but it contains no methylation or conversion events for these
assays. `--save-truth` publishes the prepared variant set as a normalized,
phased VCF truth artifact; it does not create a methylation profile or MethDB.

See [Customize](customize.md) for all shared controls, the
[target BED contract](../reference/formats/tbs-catalog.md) for WES and TS input rules,
and [Outputs](../outputs/index.md) for the published file contracts.

<div class="next-step" markdown>

**Next:** [Customize the genome, fragmentation, sequencing, and output stages](customize.md).

</div>
