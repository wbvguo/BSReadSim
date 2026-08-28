# Outputs

A successful run publishes one read representation—FASTQ or annotated BAM—and
then publishes `<prefix>.manifest.json` as the completion record. The command
prints the manifest path to standard output, for example:

```text
/work/project/runs/example-fastq/sim.manifest.json
```

The default prefix is `sim`. Change it with `--prefix`; use `--format fastq`
when uncompressed `.fastq` files are more convenient.

The manifest is UTF-8, two-space-indented JSON with one object key per line.
Its `command` section records canonical, shell-quoted forms of the CLI that was
received and the fully expanded reproducible command. Shell-only syntax such as
environment-variable expressions cannot be recovered after the shell expands
it. Runs started through the Python API record
`"interface": "python-api"` instead.

## What can be produced?

| Output choice | Select with | Files for the default `sim` prefix |
| --- | --- | --- |
| Paired compressed FASTQ | Default (`--format fastq.gz`) | `sim.R1.fastq.gz`, `sim.R2.fastq.gz`, `sim.manifest.json` |
| Paired uncompressed FASTQ | `--format fastq` | `sim.R1.fastq`, `sim.R2.fastq`, `sim.manifest.json` |
| Single-end FASTQ | `--single-end` | `sim.R1.fastq.gz`, `sim.manifest.json` |
| Annotated BAM | `--format bam` | `sim.bam`, `sim.manifest.json` |
| Saved methylation profile | Add `--save-methdb` | `truth/sim.methdb`, plus one of the sets above |
| Saved variant set | Add `--save-vcf` | `truth/sim.variants.vcf.gz`, plus one of the sets above |
| All simulation truth artifacts | Add `--save-truth` | both files beneath `truth/`, plus one read set above |

BAM replaces the FASTQ files; it is not an additional sidecar. A saved
methylation profile is different: its MethDB file is an optional reusable
simulation truth artifact and does not replace the selected read output.

Standalone `bsreadsim build variants`, `build methdb`, and `build rrbs`
commands write only their requested artifact and print its path. The
transactional file set and run manifest described on this page belong to
`bsreadsim run`.

## Example A: paired FASTQ

```bash
bsreadsim run wgbs \
  --reference reference.fa \
  --output runs/example-fastq \
  --fragments 2 \
  --read-length 5 \
  --insert-mean 8 \
  --insert-sd 0 \
  --mutation-rate 0 \
  --phred 30 \
  --error-rate 0 \
  --seed 42
```

The output directory has one file per mate and one run manifest:

```text
runs/example-fastq/
├── sim.R1.fastq.gz
├── sim.R2.fastq.gz
└── sim.manifest.json
```

Exact loci and bases depend on the reference and seed. The short records below
show the emitted format.

### `sim.R1.fastq.gz`

```bash
gzip -dc runs/example-fastq/sim.R1.fastq.gz | head -n 4
```

```text
@chr1:101-108:0/1
ACTGT
+
?????
```

This is mate 1. The sequence and quality strings have equal length; `?` is
Phred 30 in Phred+33 encoding.

### `sim.R2.fastq.gz`

```bash
gzip -dc runs/example-fastq/sim.R2.fastq.gz | head -n 4
```

```text
@chr1:101-108:0/2
CGACA
+
?????
```

The two names share `chr1:101-108:0`, so they came from the same simulated
fragment. `/1` and `/2` identify the mates. Coordinates are one-based and
inclusive; the final field is a zero-based hexadecimal fragment ordinal. See
[read names](../reference/formats/read-name.md) for the complete contract.

### `sim.manifest.json`

The manifest is already pretty-printed. Inspect it directly or use
`python -m json.tool` for validation:

```bash
python -m json.tool runs/example-fastq/sim.manifest.json | less
```

Selected fields from a two-fragment paired run have this shape:

```json
{
  "command": {
    "full_command": "bsreadsim run wgbs --reference /work/project/reference.fa --output /work/project/runs/example-fastq --fragments 2 --seed 42 ...",
    "interface": "cli",
    "user_command": "bsreadsim run wgbs --reference reference.fa --output runs/example-fastq --fragments 2 --seed 42"
  },
  "details": {
    "configuration": {
      "coverage": {
        "type": "uniform"
      },
      "methylation": {
        "beta": {
          "CG": "0.5,0.5",
          "CHG": "0.01,0.05",
          "CHH": "0.01,0.05"
        }
      },
      "technology": "WGBS"
    },
    "randomness": {
      "master_seed": "42",
      "methylation_seed": "0",
      "mutation_seed": "0",
      "phasing_seed": "0"
    }
  },
  "inputs": [
    {
      "format": "fasta",
      "path": "/work/project/reference.fa",
      "role": "reference",
      "sha256": "<64 hexadecimal characters>",
      "size_bytes": 123
    }
  ],
  "outputs": [
    {
      "path": "/work/project/runs/example-fastq/sim.R1.fastq.gz",
      "record_count": 2,
      "role": "read1",
      "sha256": "<64 hexadecimal characters>",
      "size_bytes": 123
    },
    {
      "path": "/work/project/runs/example-fastq/sim.R2.fastq.gz",
      "record_count": 2,
      "role": "read2",
      "sha256": "<64 hexadecimal characters>",
      "size_bytes": 124
    }
  ],
  "run_id": "<UUID>",
  "status": "complete",
  "summary": {
    "fragment_count": 2,
    "methylation_site_count": 24,
    "output_file_count": 2,
    "output_format": "fastq.gz",
    "output_size_bytes": 247,
    "paired_end": true,
    "read_base_count": 20,
    "read_count": 4,
    "skipped_fragment_count": 0,
    "technology": "WGBS",
    "template_base_count": 16
  },
  "version": 2
}
```

This is a selected-field view, not the entire document. The five user-facing
blocks are `summary`, `inputs`, `outputs`, `command`, and `details`. `summary`
contains one final set of run totals; per-file record counts remain with the
corresponding `outputs` entries. `details` contains the effective configuration,
contig identities, model choices, seeds, software and protocol versions,
contracts, stream identity, and the reproducibility digest.

Within `command`, `user_command` preserves a canonical rendering of the
arguments received from the user. `full_command` is generated from the
effective configuration and contains every effective value-bearing run option,
including normalized absolute paths, all deterministic defaults, focused
seeds, execution settings, output policy, and exactly one effective `--seed`.
Boolean flags are present when enabled and absent when disabled. Re-parsing
`full_command` produces the same canonical configuration represented by
`details.configuration` and `details.randomness`.

`details.configuration` deliberately omits `seed` and `seeds`; their structured
record is `details.randomness`. Model declarations use `type` (for example,
`"type": "uniform"`) rather than the less specific `kind`. Beta pairs use the
same compact comma-separated representation as the CLI.

Output destinations are write-once. Before executing `command.full_command`,
move the original output directory aside, or change only `--output` to a new
empty destination when reproducing the simulated reads.

| Output role | `record_count` means |
| --- | --- |
| `read1`, `read2` | accepted source fragments |
| `bam` | read alignment records |
| `truth.methdb` | stored methylation sites |
| `truth.vcf` | non-header VCF records; zero is valid |

## Example B: annotated BAM

```bash
bsreadsim run wgbs \
  --reference reference.fa \
  --output runs/example-bam \
  --fragments 2 \
  --read-length 5 \
  --insert-mean 8 \
  --insert-sd 0 \
  --mutation-rate 0 \
  --seed 42 \
  --format bam
```

```text
runs/example-bam/
├── sim.bam
└── sim.manifest.json
```

### `sim.bam`

BAM is binary, so inspect it with samtools rather than `head`:

```bash
samtools view -H runs/example-bam/sim.bam
samtools view runs/example-bam/sim.bam | head -n 1
```

A representative record rendered as SAM looks like this (fields are
tab-separated):

```text
chrMini:11-16:7	0	chrMini	11	60	2M1I1M2D1M	*	0	0	ACGTN	ABCDE	RG:Z:00000000-0000-4000-8000-000000000002	AS:i:5	XG:Z:CT	XR:Z:CT	YS:Z:OT	zt:Z:AAAAA	zr:B:S,5,0,0,0,0,0,0,0,0,0,0,0
```

The standard SAM columns retain the simulated origin, CIGAR, sequence, and
quality. BSReadSim adds:

| Tag | Present when | Meaning |
| --- | --- | --- |
| `XG:Z` | Bisulfite BAM records | Bismark-compatible genome conversion: `CT` or `GA` |
| `XR:Z` | Bisulfite BAM records | Bismark-compatible conversion of the current read: `CT` or `GA` |
| `YS:Z` | Bisulfite BAM records | Complete library strand: `OT`, `OB`, `CTOT`, or `CTOB` |
| `zt:Z` | Every BAM record | One base-state character per BAM `SEQ` base |
| `zr:B:S` | Every BAM record | Twelve-value summary for that read |
| `zf:B:S` | `--format bam --fragment-summary` | Twelve-value summary for the complete physical fragment, copied to both mates |
| `zx:Z` | `--fragment-realization` | Packed complete-fragment methylation and conversion realization, copied to both mates |

`--fragment-realization` requires BAM and implies fragment summaries. The BAM stream is
unsorted; sort and index it before using coordinate-indexed tools:

```bash
samtools sort -o runs/example-bam/sim.sorted.bam runs/example-bam/sim.bam
samtools index runs/example-bam/sim.sorted.bam
```

Reads and qualities remain recoverable with `samtools fastq`. See
[annotated BAM](../reference/formats/bam.md) for every tag bit and field.
FASTQ and QNAME do not duplicate the `XG`, `XR`, or `YS` truth annotations.

## Example C: saved methylation profile

Save the prepared methylation profile as a reusable MethDB truth artifact
without changing the chosen read container:

```bash
bsreadsim run wgbs \
  --reference reference.fa \
  --output runs/with-methdb \
  --fragments 2 \
  --mutation-rate 0 \
  --seed-meth 7 \
  --seed 42 \
  --save-methdb
```

```text
runs/with-methdb/
├── sim.R1.fastq.gz
├── sim.R2.fastq.gz
├── truth/
│   └── sim.methdb
└── sim.manifest.json
```

### `truth/sim.methdb`

MethDB is an opaque binary file, not a table. Its first six bytes identify
the format and the seventh byte is the internal format version:

```bash
xxd -l 7 runs/with-methdb/truth/sim.methdb
```

```text
00000000: 6d65 7468 6462 01                        methdb.
```

Do not edit or concatenate it. Load the exact saved profile in a later run with
`--methdb runs/with-methdb/truth/sim.methdb`. The creating run records it in
the manifest output list with role `truth.methdb`, its record count, size, and
SHA-256. Decode a human-readable, BGZF-compressed extended BED with:

```bash
bsreadsim export methdb \
  -i runs/with-methdb/truth/sim.methdb \
  -o sim.methdb.bed.gz
```

See [fixed MethDB](../reference/formats/methdb.md) for the BED columns and the plain-text
`--no-compression` mode.

## Example D: saved variant set

`--save-vcf` publishes the exact prepared variant set used by the run as a
normalized and phased VCF truth artifact:

```bash
bsreadsim run wgbs \
  --reference reference.fa \
  --output runs/with-variants \
  --fragments 1000 \
  --mutation-rate 0.001 \
  --seed-mut 7 \
  --seed 42 \
  --save-vcf
```

```text
runs/with-variants/
├── sim.R1.fastq.gz
├── sim.R2.fastq.gz
├── truth/
│   └── sim.variants.vcf.gz
└── sim.manifest.json
```

Inspect the compressed VCF as text:

```bash
gzip -dc runs/with-variants/truth/sim.variants.vcf.gz | head
```

For a run using `--vcf INPUT`, the saved variant set is validated, normalized,
and deterministically phased; it is not necessarily a byte copy of `INPUT`.
For a run with `--mutation-rate 0`, it is a valid header-only VCF.gz and its
manifest record count is zero. `--save-truth` combines this artifact with the
saved methylation profile from Example C.

## When is a run complete?

BSReadSim stages data files, publishes them, and publishes the manifest last.
Treat a present manifest with `"status": "complete"` as the completion marker.
Existing destination artifacts are never overwritten; a conflicting prefix is
reported as an error. Keep the manifest beside the reads so configuration,
seeds, summary, and checksums stay attached to the dataset.

Staging is disk-backed and streaming: each FASTQ/BAM batch is written
immediately to a private file in the destination filesystem, with its digest
updated as bytes pass. Successful completion publishes those files and then
the manifest. Transactional publication therefore does not hold the complete
BAM or FASTQ in RAM; memory remains bounded by active batches and in-flight
work. Writing directly to the final filename would expose a plausible-looking
partial file after a crash without materially reducing memory use.

<div class="next-step" markdown>

**Next:** [Inspect the generated files](inspect.md), or read
the [reproducibility section](../simulation/customize.md#reproducibility).

</div>
