# Troubleshoot

Start with the first `bsreadsim: error:` message. Later shell messages or
downstream-tool failures are often consequences of that first error. Use the
matching symptom below, then follow the linked reference only when the short
fix does not resolve it.

<div class="reference-jump" markdown>

**Jump to:** [commands and setup](#commands-and-setup) ·
[input files](#input-files) · [read generation and sampling](#read-generation-and-sampling) ·
[outputs](#outputs) · [performance](#performance)

</div>

## Commands and setup

### A command or option is not recognized

Run `--help` at the same level as the failing command:

```bash
bsreadsim --help
bsreadsim run --help
bsreadsim run wgbs --help
bsreadsim build methdb --help
```

The current interface expects an action and target before its options, for
example `bsreadsim run wgbs [OPTIONS]`. Old flat commands and JSON
run-configuration options are not accepted. See the
[CLI reference](../reference/cli.md#commands) for the current hierarchy.

### Installation or a source build fails

First confirm that the operating system and toolchain meet the
[installation requirements](../getting-started/installation.md). Native
Windows, macOS, and Linux ARM64 are not currently supported.

If the error reports missing HTSlib or htscodecs sources, initialize the
pinned submodules from the repository root, then rerun the build:

```bash
git submodule update --init --recursive --depth 1
```

BSReadSim builds against the repository-pinned HTSlib rather than a system
HTSlib installation.

## Input files

### A contig or coordinate does not match

Use the first whitespace-delimited token in each FASTA header as the contig
name. Names are case-sensitive: `1`, `chr1`, and `chr1 description` are not
interchangeable.

Confirm that every coordinate-bearing input was prepared from the exact
reference assembly used by the command. A matching contig name does not make
coordinates from another assembly valid.

### An input row is rejected

Confirm that the option selects the file's actual serialization; a filename
suffix does not select the parser. Then compare the first rejected row with
the relevant contract in [Input file formats](../reference/formats.md).

Check these points first:

- contig names and coordinates agree with the reference;
- VCF and CGmap positions are one-based, while BED-derived intervals are
  zero-based and half-open;
- the field count is consistent with the selected format; and
- row ordering, numeric ranges, headers, and compression match that format's
  contract.

BSReadSim rejects malformed or ambiguous rows rather than inferring their
intended meaning.

### Input options conflict

- **Baseline methylation:** Choose at most one of MethBG, MethBED, bedMethyl,
  CGmap, and MethDB.
- **VCF:** Remove `--mutation-rate` whenever `--vcf` is present, including an
  explicit value of zero.
- **ASM:** A VCF is optional. Without one, retained ASM links define the
  required heterozygous SNVs; with one, every link must match the VCF. Omit a
  nonzero `--mutation-rate` when using ASM.
- **MethDB:** Use the matching reference and remove explicit VCF, baseline
  profile, ASM, and pooling inputs. MethDB also rejects nonzero de novo
  mutation.
- **Pooling:** `--pool-meth` requires CGmap, bedMethyl, MethBG, or MethBED; it
  cannot operate on a generated profile or MethDB.
- **Standard sequencing:** WGS, WES, and TS reject methylation inputs and
  bisulfite-only truth options.

See [CLI option combination rules](../reference/cli.md#option-combination-rules)
for the complete set of exclusions.

## Read generation and sampling

### The read count is rejected

Use either `--reads` or `--depth`, not both. In paired-end mode, `--reads`
counts individual records and must be even; `--reads 1000` produces 500 read
pairs. See
[Required inputs and read count](../reference/cli.md#required-inputs-and-dataset-size)
for limits and depth rounding.

### Read or fragment lengths conflict { #read-and-insert-lengths-conflict }

First ensure that
`insert_min <= insert_mean <= insert_max`.

- For WGBS, WGS, TBS, WES, and TS, `--insert-sd 0` selects a fixed fragment
  length, so `--read-length` must not exceed `--insert-mean`. With a positive
  SD, the read length must not exceed `--insert-min`.
- For RRBS, `--read-length` must not exceed `--insert-min`. The minimum and
  maximum retain restriction fragments; the mean must still lie between them
  even though it does not reshape the candidate fragments.

Paired reads may overlap, but both mates must fit the physical fragment. See
[Fragment length](../reference/cli.md#fragment-geometry) for assay-specific
behavior.

### GC-profile sampling is rejected

Whole-genome sampling requires `--sampling gc` and `--gc-profile PATH`
together. With variable fragment lengths, GC-profile sampling does not support
VCF, ASM, or de novo variants. Either remove those variant sources and use
mutation rate zero, or select a fixed fragment length with `--insert-sd 0`.

If the error mentions a positive-probability GC bin with no eligible fragment,
adjust the profile or fragment geometry. See the
[target-GC profile contract](../reference/formats.md#target-gc-profile).

### An RRBS candidate file is rejected

Regenerate the candidate BED with the same reference, VCF or mutation policy,
`--seed-mut`, `--seed-phase`, cut sites, read layout and length, insert minimum
and maximum, and ambiguity threshold used by the run. `--insert-mean` and
`--insert-sd` do not need to match the later run, although the mean must remain
within its bounds.

External scoring may change only the score column. Missing, extra, duplicated,
or modified candidate identity fields are rejected; row order may change. See
the [RRBS candidate BED contract](../reference/formats.md#rrbs-candidate-bed).

### Targets produce no fragments

Start with `--sampling uniform`, then check targets near contig boundaries,
fragment and read geometry, `--center-sd`, and N-rich reference regions. A
valid BED interval may still be ineligible when no complete fragment and read
layout can be placed around it.

With `--sampling score`, at least one eligible target must have a positive
weight. See the
[capture-target BED contract](../reference/formats.md#capture-target-bed).

### A quality or error model is rejected

- Use either `--quality-model FILE` or an explicit `--phred`.
- Use either `--error-model FILE` or an explicit `--error-rate`.
- Ensure that the error model contains every Phred value that the selected
  quality policy can emit.

When both empirical models are used, omit both uniform settings. Validate the
JSON structure against the
[sequencing-model contracts](../reference/formats.md#sequencing-model-inputs).

## Outputs

### Destination files already exist

For `run`, the output directory may already exist, but none of the selected
prefix's destination artifacts may exist. Choose another directory, change
`--prefix`, or move the old artifacts before rerunning.

For a standalone `build` command, the parent directory must already exist and
the destination file must be new. BSReadSim never replaces an existing
artifact.

### The manifest is missing

Treat the run as incomplete. The manifest is published last, and a completed
dataset has `<prefix>.manifest.json` with `"status": "complete"`.

An abruptly interrupted process may leave `.PREFIX.staging-*` data. Remove it
only after confirming that no BSReadSim process still owns the directory. See
the [run-manifest contract](../outputs/index.md#run-manifest).

### BAM tools report unsorted input

BSReadSim intentionally writes unsorted BAM. Sort and index it before using a
coordinate-based tool:

```bash
samtools sort -o sim.sorted.bam sim.bam
samtools index sim.sorted.bam
```

See [Annotated BAM](../outputs/index.md#annotated-bam) for the BAM contract and
[Read files](../outputs/index.md#read-files) for inspection commands.

## Performance

### A run is slow or uses too much memory

Start with `--threads 1`, then increase the thread budget while measuring CPU,
memory, and disk throughput. BSReadSim distributes this single budget across
the pipeline.

Omit `--fragment-realization` when complete-fragment BAM truth is unnecessary.
For compressed FASTQ, a lower `--gzip-level` reduces CPU work at the cost of
larger files. Output uses disk-backed staging, so slow storage can also limit
throughput.
