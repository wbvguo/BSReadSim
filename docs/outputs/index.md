# Outputs

Every successful `run` writes simulated reads and a manifest. The variant set
and, when applicable, methylation profile may also be exported as truth
artifacts. See
[Output and reproducibility](../simulation/customize.md#reproducibility) to
configure the output format and optional artifacts; this page explains the
resulting files.

## Output directory

With the default prefix `sim`, a paired-end bisulfite run using compressed
FASTQ and exporting both truth artifacts has this structure:

```text { .no-copy }
OUTPUT/
├── sim.R1.fastq.gz
├── sim.R2.fastq.gz
├── sim.manifest.json
└── truth/
    ├── sim.variants.vcf.gz
    └── sim.methdb
```

Single-end FASTQ omits `sim.R2.fastq.gz`. BAM output replaces the FASTQ files
with `sim.bam`. The manifest is written for every successful run, while
exported truth artifacts are stored under `truth/`.

??? info "Related options"

    - [`--prefix`](../simulation/customize.md#output-format) changes the
      filename prefix.
    - [`--single-end`](../simulation/customize.md#read-layout) selects
      single-end output. Single-end FASTQ contains R1 only; paired-end FASTQ
      contains matching R1 and R2 files.
    - [`--format`](../simulation/customize.md#output-format) selects plain-text
      `.fastq`, gzip-compressed `.fastq.gz`, or BAM.
    - [`--save-vcf`](../simulation/customize.md#truth-artifacts) exports the
      VCF; [`--save-methdb`](../simulation/customize.md#truth-artifacts)
      exports MethDB; and
      [`--save-truth`](../simulation/customize.md#truth-artifacts) exports all
      applicable truth artifacts.

??? info "Inspect output files"

    Preview a manifest and compressed FASTQ with standard command-line tools:

    ```bash
    python -m json.tool runs/wgbs/sim.manifest.json
    gzip -dc runs/wgbs/sim.R1.fastq.gz | head -n 8
    ```

    BAM inspection, sorting, and indexing require
    [samtools](https://www.htslib.org/), which is not installed by BSReadSim:

    ```bash
    samtools view -H runs/example-bam/sim.bam
    samtools view runs/example-bam/sim.bam | head
    samtools sort -o runs/example-bam/sim.sorted.bam runs/example-bam/sim.bam
    samtools index runs/example-bam/sim.sorted.bam
    ```

## Read files

### FASTQ

Each read uses the standard four-line FASTQ format:

```text { .no-copy }
@chr1:101-108:2a/1
ACGTCGTA
+
IIIIIIII
```

The four lines contain the read identifier, sequence, `+` separator, and Phred
quality string. Read identifiers have this form:

```text { .no-copy }
@<contig>:<start>-<end>:<ordinal-hex>/<mate>
```

`start` and `end` are the one-based inclusive fragment envelope. The fragment
ordinal is zero-based lowercase hexadecimal, and `/1` or `/2` identifies the
mate. Both mates share the contig, envelope, and ordinal.

### Annotated BAM

Annotated BAM stores reads with truth annotations in one unsorted file.
Standard SAM fields record each read's sequence, quality, and simulated origin;
additional tags retain
applicable strand, methylation, conversion, variant, and sequencing-error
truth. Optional fragment annotations retain a summary or the complete
methylation and conversion realization. The file conforms to SAM/BAM 1.6 and
is readable by HTSlib and samtools.

??? info "BAM header and standard fields"

    Ordinary `@HD`, `@SQ`, `@RG`, and `@PG` records are followed by BSReadSim
    format comments:

    ```text
    @CO    BSREADSIM_BAM_CONTRACT=bsreadsim-bam
    @CO    BSREADSIM_ZT=state64;ALPHABET=ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_
    @CO    BSREADSIM_XG=bismark-genome-conversion;ENABLED=0|1;VALUES=CT|GA
    @CO    BSREADSIM_XR=bismark-read-conversion;ENABLED=0|1;VALUES=CT|GA
    @CO    BSREADSIM_YS=bismark-strand-id;ENABLED=0|1;VALUES=OT|OB|CTOT|CTOB
    @CO    BSREADSIM_ZR=u16x12;REQUIRED=1
    @CO    BSREADSIM_ZF=u16x12;ENABLED=0|1
    @CO    BSREADSIM_ZX=packed-b64url;ENABLED=0|1;BIT_ORDER=LSB0
    ```

    The `@RG ID` is the run UUID. The manifest records the same UUID and
    normalized configuration.

    Viewed as SAM, a representative single-end record has this shape (the run
    UUID is abbreviated):

    ```text
    chr1:101-104:0	0	chr1	101	60	4M	*	0	0	ACGT	IIII	AS:i:4	RG:Z:RUN_UUID	XG:Z:CT	XR:Z:CT	YS:Z:OT	zt:Z:AAAA	zr:B:S,0,0,0,0,0,0,0,0,0,0,0,0
    ```

    The first eleven fields are standard SAM. The remaining fields are typed
    tags; paired records additionally carry mate fields and tags.

    Both mates share
    `<contig>:<one-based-inclusive-start>-<one-based-inclusive-end>:<ordinal-hex>`.
    The ordinal is variable-width lowercase hexadecimal without `0x`; SAM
    flags identify mate 1 and mate 2. Records are mapped to the simulated
    origin with an indel-aware, query-complete CIGAR and reference-forward
    SEQ/QUAL.

    - `MAPQ=60` denotes a known simulated origin, not empirical confidence.
    - `AS:i` equals query length, the maximum simulation-origin score.
    - `RG:Z` links the record to the run.
    - Paired records use standard flags, RNEXT/PNEXT/TLEN, `MQ:i`, and `MC:Z`.

    FASTQ records and QNAMEs do not carry strand-origin truth. The readable
    truth annotation is available only in annotated BAM.

??? info "BAM truth tags"

    **Bisulfite `XG`, `XR`, and `YS` tags**

    WGBS, RRBS, and TBS use the Bismark conventions: `XG` records genome
    conversion, `XR` records read conversion, and `YS` identifies the library
    strand. Mates share `XG` and `YS`, while `XR` is mate-specific.

    | `YS` library strand | `XG` | R1 `XR` | R2 `XR` | Enabled by |
    | --- | --- | --- | --- | --- |
    | `OT` | `CT` | `CT` | `GA` | directional and undirectional |
    | `OB` | `GA` | `CT` | `GA` | directional and undirectional |
    | `CTOT` | `CT` | `GA` | `CT` | undirectional only |
    | `CTOB` | `GA` | `GA` | `CT` | undirectional only |

    **Required `zt:Z` per-read state**

    `zt` is a per-base truth string with the same length and order as the BAM
    `SEQ` field. Each character describes the base at the corresponding
    position, including for reverse-strand reads.

    | Bits | Meaning |
    | --- | --- |
    | 0-1 | context: 0 none, 1 CG, 2 CHG, 3 CHH |
    | 2 | methylated |
    | 3 | bisulfite conversion succeeded |
    | 4 | base is produced or affected by a variant event |
    | 5 | sequencing error |

    **Read and fragment summaries (`zr:B:S` and `zf:B:S`)**

    `zr` and `zf` are numeric-array BAM auxiliary tags. In `B:S`, `B` denotes
    an array and `S` selects unsigned 16-bit elements. BSReadSim uses a fixed
    12-value layout: one flags value and 11 counts, all integers in
    `[0, 65535]`. In SAM text, the tag has this form:

    ```text
    zr:B:S,VALUE_0,VALUE_1,...,VALUE_11
    ```

    `zr` is required and summarizes the current read. `zf` is optional and
    summarizes the complete physical fragment; both mates receive the same
    value.

    The first element packs the following flags:

    - bits 0-1: haplotype
    - bits 2-3: informative strand (`0` none, `1` Watson, `2` Crick)
    - bits 4-6: conversion mode (`0` C-to-T, `1` G-to-A, `2` none)
    - bits 7-11: ASM, SNV, insertion, deletion, and methylation conversion
    - bit 13: count overflow

    The array elements are ordered as follows; `_m` denotes the methylated
    count:

    ```text
    0       flags
    1-2     n_cg, n_cg_m
    3-4     n_chg, n_chg_m
    5-6     n_chh, n_chh_m
    7-8     n_conversion_success, n_conversion_failure
    9       n_variants
    10      n_sequencing_errors
    11      n_false_methylation_errors
    ```

    **Complete-fragment realization (`zx:Z`)**

    `zx:Z` stores realized methylation states and bisulfite-conversion outcomes
    across the complete physical fragment, including positions outside the
    emitted reads. Both mates receive the same value.
    [`--fragment-realization`](../simulation/customize.md#truth-artifacts)
    enables the tag and implies
    [`--fragment-summary`](../simulation/customize.md#truth-artifacts).

    The value has four dot-separated fields:

    ```text
    <site-count-hex>.<convertible-count-hex>.<methylation-bits>.<conversion-bits>
    ```

    `site-count-hex` counts methylation sites. `convertible-count-hex` counts C
    bases in C-to-T mode or G bases in G-to-A mode. These counts may differ
    because they describe different sets of positions.

    A `1` means methylated in `methylation-bits` and successfully converted in
    `conversion-bits`. Counts use lowercase hexadecimal; bit vectors are
    packed in LSB0 order and encoded as unpadded base64url. Here, `64` names the
    encoding alphabet, not a site limit.

    Take a three-site example:

    ```text
    site-count          3
    convertible-count   3
    methylation         1,0,1 -> 00000101 -> 000001(B) 010000(Q) -> BQ
    conversion          0,1,0 -> 00000010 -> 000000(A) 100000(g) -> Ag
    zx                  3.3.BQ.Ag
    ```

## Run manifest

The manifest is written last and records the effective configuration and
provenance needed to inspect and reproduce the run. Treat the dataset as
complete only when the manifest contains `"status": "complete"`.

??? info "Manifest fields"

    The manifest is UTF-8 JSON with the filename
    `<prefix>.manifest.json`.

    | Block | Contents |
    | --- | --- |
    | `summary` | Assay, read format, output counts, and sizes |
    | `inputs` | Input roles, resolved paths, byte sizes, and SHA-256 digests |
    | `outputs` | Output roles, paths, record counts, byte sizes, and SHA-256 digests |
    | `command` | Received command and fully expanded effective command |
    | `details` | Normalized settings, resolved seeds, contigs, versions, and contracts |

    The manifest and BAM header share a run UUID.
    `command.full_command` records resolved paths and defaults, while
    `details.randomness` records the master seed and resolved stage seeds.

    **Manifest example**

    This abridged example shows the main structure. Actual manifests include
    additional counts, file identities, configuration fields, and contracts.

    ```json
    {
      "version": 2,
      "status": "complete",
      "run_id": "12345678-1234-4234-8234-123456789abc",
      "summary": {
        "technology": "WGBS",
        "output_format": "fastq.gz",
        "paired_end": true,
        "fragment_count": 500000,
        "read_count": 1000000
      },
      "inputs": [
        {
          "role": "reference",
          "format": "fasta",
          "path": "/data/reference.fa",
          "sha256": "..."
        }
      ],
      "outputs": [
        {
          "role": "read1",
          "path": "/data/run/sim.R1.fastq.gz",
          "record_count": 500000,
          "sha256": "..."
        },
        {
          "role": "read2",
          "path": "/data/run/sim.R2.fastq.gz",
          "record_count": 500000,
          "sha256": "..."
        }
      ],
      "command": {
        "interface": "cli",
        "user_command": "bsreadsim run wgbs ...",
        "full_command": "bsreadsim run wgbs ..."
      },
      "details": {
        "randomness": {
          "master_seed": "104729",
          "mutation_seed": "...",
          "phasing_seed": "...",
          "methylation_seed": "..."
        }
      }
    }
    ```

## Truth artifacts { #saved-truth-artifacts }

The exported VCF contains the normalized and phased variant set used by the
run, so it may not be a direct copy of the input VCF. A run without variants
produces a valid header-only file.

For WGBS, RRBS, and TBS, the exported MethDB contains the methylation profile
and embedded variants. Load it with `--methdb` to reuse the same snapshot in a
later run.

See the [VCF](../reference/formats.md#vcf) and
[MethDB](../reference/formats.md#methdb) file specifications.
