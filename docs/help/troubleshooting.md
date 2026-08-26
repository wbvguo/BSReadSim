# Troubleshooting

Start with the first BSReadSim error, not a later shell or downstream-tool
message. Configuration, input formats, and destination collisions are checked
before final artifacts are published.

| Problem | Jump to |
| --- | --- |
| Command or option is unknown | [Current command hierarchy](#a-command-or-option-is-not-recognized) |
| Output collision or missing manifest | [Destinations](#the-output-directory-already-exists) · [Completion](#the-manifest-is-missing) |
| FASTA, BED, VCF, or methylation input | [Contigs](#a-contig-name-does-not-match) · [Formats](#an-input-format-is-rejected) |
| Variant or MethDB mismatch | [VCF](#vcf-and-mutation-options-conflict) · [MethDB](#a-methdb-is-rejected) |
| Fragment geometry | [Read and insert lengths](#read-and-insert-lengths-conflict) · [TBS fixed insert](#tbs-rejects-the-default-insert-settings) |
| RRBS or TBS selection | [RRBS candidates](#an-rrbs-candidate-file-is-rejected) · [TBS targets](#a-tbs-target-is-valid-bed-but-produces-no-fragments) |
| Sequencing models | [Quality and error](#quality-or-error-model-options-conflict) |
| Performance or BAM tooling | [Resources](#a-run-is-slow-or-uses-too-much-memory) · [Sorted BAM](#bam-tools-expect-sorted-input) |
| Installation | [HTSlib](#a-source-build-cannot-find-htslib) · [Platforms](#installation-fails-on-windows-or-macos) |

## A command or option is not recognized

The current hierarchy always puts an action and usually a target before its
options:

```text
bsreadsim run wgbs [OPTIONS]
bsreadsim run rrbs [OPTIONS]
bsreadsim run tbs [OPTIONS]
bsreadsim build variants [OPTIONS]
bsreadsim build methdb [OPTIONS]
bsreadsim build rrbs [OPTIONS]
```

Use help at the same level as the problem:

```bash
bsreadsim --help
bsreadsim run --help
bsreadsim run wgbs --help
bsreadsim build methdb --help
```

Old flat commands and JSON run-configuration options are not part of the
current user CLI. Convert them to direct parameters under the appropriate
`run` or `build` subcommand. The [CLI reference](../reference/cli.md) contains
the complete short/long/type/default tables.

## The output directory already exists

An existing directory is allowed; an existing destination artifact is not.
For prefix `sim`, a second paired FASTQ run collides with files such as
`sim.R1.fastq.gz` and `sim.manifest.json`. Choose a new `--output` directory,
change `--prefix`, or deliberately move the old artifacts before rerunning.
BSReadSim never replaces them.

`run` creates missing parent directories. By contrast, each standalone
`build` command expects the output parent to exist and the destination file to
be absent.

## The manifest is missing

The manifest is published last. If it is absent, do not treat visible read or
truth files as a committed BSReadSim dataset. A successful run has
`PREFIX.manifest.json` with `"status": "complete"`.

Private `.PREFIX.staging-*` directories are implementation staging areas, not
completed output. A normally handled error removes its own staging data. If a
process was killed abruptly, confirm that no BSReadSim process still owns a
staging directory before cleaning it up.

## A contig name does not match

Every biological input must use the first whitespace-delimited token from the
corresponding FASTA header. Names are case-sensitive; `1`, `chr1`, and
`chr1 description` are not interchangeable aliases.

For multiple inputs, compare their contig sets against the exact FASTA used by
the run. A coordinate that exists in a different genome build is still an
error even when the contig name happens to match.

## An input format is rejected

Select the option naming the serialization. For example, use `--cgmap` for an
eight-column CGmap file and `--bed-methyl` for BED9+2 or BED9+9. Compression
and filename suffixes do not select a parser.

Check structural requirements before changing biological parameters:

- FASTA sequence accepts only A, C, G, T, and N;
- VCF and CGmap positions are one-based, while BED intervals are zero-based
  and half-open;
- all coordinates must lie within the exact named contig;
- TBS requires exactly BED6 with `start < end`, a nonempty name, a finite
  non-negative score, and strand `+`, `-`, or `.`; and
- compressed input support is defined separately by each format contract.

Follow [File formats](../reference/formats.md). BSReadSim fails closed instead
of guessing how a malformed row was intended.

## VCF and mutation options conflict

A predefined VCF and de novo mutation generation cannot both be active. Use
`--vcf sample.vcf.gz` by itself; it automatically selects an internal mutation
rate of zero. Remove `--mutation-rate` entirely when `--vcf` is present,
including an explicit zero.

`--asm` and `--asm-bed` require an input VCF containing the linked
heterozygous SNVs. They cannot link to de novo variants. Also remove
`--no-update-variant-boundaries` whenever a VCF or positive mutation rate is
active.

## A MethDB is rejected

MethDB stores a prepared methylation profile bound to the reference, prepared
variant set, and methylation settings used when it was built. Repeat the exact
reference, mutation or VCF choice, and relevant `--seed-mut`, `--seed-phase`,
and `--seed-meth` values.

Do not combine `--methdb` with `--cgmap`, `--bed-methyl`, `--asm`,
`--asm-bed`, or `--cgmap-pool`. A MethDB is already the complete normalized
methylation profile. Historical incompatible representations must be
regenerated with the current `build methdb` or `--save-methdb` implementation.

## Read and insert lengths conflict

The current shared geometry contract requires:

```text
insert_min <= insert_mean <= insert_max
```

For variable WGBS and RRBS, each read must fit `insert_min`. For fixed WGBS
and TBS (`--insert-sd 0`), each read must fit `insert_mean`. With paired-end
data, R1 and R2 currently have the same `--read-length`; overlap is allowed,
but each mate must fit the physical fragment.

For RRBS, min/max filter restriction fragments even though mean/SD do not draw
their lengths. If you lower `--insert-max` below the default mean 400, set an
in-range `--insert-mean` as well.

## TBS rejects the default insert settings

TBS requires one fixed insert. The shared defaults use a positive SD, so set
both the length and fixed model, for example
`--insert-mean 300 --insert-sd 0`. Ensure that `--read-length` does not exceed
the mean and that the fixed fragment can be placed around at least one target.

## An RRBS candidate file is rejected

`--sampling score` requires `--rrbs-candidates`. The candidate BED is verified
against a freshly regenerated domain. Repeat the same reference, cut sites,
read and insert geometry, ambiguity threshold, VCF/mutation settings, and
focused seeds used for `bsreadsim build rrbs`.

External tools may modify only the score column. Missing, extra, reordered
identity fields, duplicate, or otherwise modified candidates are rejected.
Input row order itself may change because matching is by deterministic
identity. The run's master `--seed` controls later sampling and is not part of
candidate construction.

## A TBS target is valid BED but produces no fragments

A syntactically valid target can still be ineligible when the fixed insert
cannot be placed around its center, projection crosses an unsupported
haplotype boundary, or a mate exceeds `--max-ambiguous-fraction`.

Start with `--sampling uniform`, inspect edge-of-contig targets, verify the
fixed insert and read length, and examine N-rich reference regions. With
`--sampling score`, at least one eligible target must retain positive mass.

## Quality or error model options conflict

An empirical model replaces its uniform counterpart:

- use either `--quality-model FILE` or an explicit `--phred`, not both;
- use either `--error-model FILE` or an explicit `--error-rate`, not both; and
- ensure the error model contains every Phred value the quality policy can
  emit.

Omit `--phred` and `--error-rate` when selecting both empirical files. The
[sequencing model contract](../reference/formats/sequencing-models.md) explains strict
JSON fields and dimensions.

## A run is slow or uses too much memory

Return to one worker in each stage, then reduce
`--max-in-flight-fragments`. If core generation chunks are too large, reduce
`--chunk-size`. Increase `--core-workers` or `--workers` only after a baseline
run succeeds and measure whether CPU or storage is the bottleneck.

Output staging is disk-backed and streaming; it does not hold a complete BAM
or FASTQ in memory. BAM fragment realization, empirical models, queued
fragments, and large batches can still increase the active working set. Lower
gzip levels can reduce CPU work at the cost of larger files.

## BAM tools expect sorted input

BSReadSim emits unsorted BAM. Sort and index it for coordinate-based tools:

```bash
samtools sort -o sim.sorted.bam sim.bam
samtools index sim.sorted.bam
```

The BAM MAPQ describes simulated-origin annotation and is not calibrated
mapping confidence from an aligner.

## A source build cannot find HTSlib

Initialize both submodule levels and retry:

```bash
git submodule update --init --recursive --depth 1
```

The build uses the pinned repository HTSlib rather than a system provider.

## Installation fails on Windows or macOS

Native Windows and macOS are outside the current production support matrix.
On Windows, build and run inside WSL2 with a Linux Python interpreter. Review
the [supported platforms](../getting-started/platforms.md) before debugging the
compiler toolchain.
