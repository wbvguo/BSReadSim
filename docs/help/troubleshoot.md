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

### `bsreadsim` is not found

Confirm that the environment containing BSReadSim is active and that its
executable is on `PATH`:

```bash
command -v bsreadsim
python -m pip show bsreadsim
bsreadsim --version
```

If `pip show` succeeds but `command -v` does not, the environment's script
directory may be absent from `PATH`, or the shell may be using a different
environment. Activate the intended Conda or virtual environment, then retry.
Do not mix a system `pip` with the Python executable from another environment.

### A command or option is not recognized

Check the installed version, then run `--help` at the same level as the
failing command:

```bash
bsreadsim --version
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

First confirm that the platform, Python version, and toolchain meet the
[installation requirements](../getting-started/installation.md).

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
name. For example, the input name for `>chr1 description` is `chr1`, not the
complete header. Names are case-sensitive, so `1` and `chr1` are different.

Confirm that every coordinate-bearing input was prepared from the exact
reference assembly used by the command. A matching contig name does not make
coordinates from another assembly valid.

An input may cover only a subset of FASTA contigs. Omitted contigs are not an
error, but every present contig block must follow the relative FASTA order.
Do not use lexical order when it would place `chr10` before `chr2`. Within a
contig, VCF positions are nondecreasing; methylation and ASM target positions
are strictly increasing and unique.

### An input file or row is rejected

When supported, isolate the problem with `validate` before running a
simulation:

```bash
bsreadsim validate \
  -r reference.fa.gz \
  --vcf variants.vcf.gz \
  --strict
```

Supply the same VCF, text methylation-profile, and ASM options that need to be
checked. `validate` does not directly accept MethDB, target-GC profiles, RRBS
candidate BED, capture-target BED, or sequencing-model JSON; those files are
checked by the command that consumes them. See
[Validate inputs](../reference/cli.md#validate-inputs) for its complete scope.

Confirm that each option selects the file's actual serialization; a filename
suffix does not select the parser. Then compare the first rejected row with
the relevant contract in [File formats](../reference/formats.md).

Check these points first:

- contig names and coordinates agree with the reference;
- VCF and CGmap positions are one-based, while BED-derived intervals are
  zero-based and half-open;
- the field count is consistent with the selected format; and
- row ordering, numeric ranges, headers, and compression match that format's
  contract.

BSReadSim rejects malformed or ambiguous rows rather than inferring their
intended meaning.

### VCF records are reported as skipped

Normal validation accepts but skips MNPs, complex replacements, and indels
longer than four bases after normalization. The validation summary reports
their counts, and they do not appear in saved variant truth. Use `--strict`
when any skipped VCF record should make validation fail.

Malformed rows, unknown contigs, invalid genotypes, retained REF mismatches,
and overlapping retained events are errors in both modes; `--strict` does not
change those checks.

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
VCF, ASM, nonzero de novo mutation, or a MethDB containing embedded variants.
Use a variant-free input world, or select a fixed fragment length with
`--insert-sd 0`.

With a fixed insert length, every positive-probability bin must have an
eligible fragment. A reference-only variable-insert run drops unreachable
positive mass and renormalizes the remaining bins, but still fails when no
positive-probability bin is reachable. Adjust the profile or fragment
geometry when either condition fails. See the
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

### Output cannot be created

For `run`, the output directory may already exist, but none of the selected
prefix's destination artifacts may exist. Choose another directory, change
`--prefix`, or move the old artifacts before rerunning.

For a standalone `build` command, the parent directory must already exist and
the destination file must be new. BSReadSim never replaces an existing
artifact. Also confirm that the destination is writable and has enough free
space for the selected read and truth outputs.

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

Compare `-t 1`, `-t 2`, and `-t 4` while measuring CPU, memory, and disk
throughput. `-t` sets the number of threads used by the pipeline; increasing
it may also increase memory use.

Omit `--fragment-realization` when complete-fragment BAM truth is unnecessary.
For compressed FASTQ, a lower `--gzip-level` reduces CPU work at the cost of
larger files. Output uses disk-backed staging, so slow storage can also limit
throughput.
