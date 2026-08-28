# Annotated BAM

Annotated BAM replaces FASTQ output when `--format bam` is selected. It is an
unsorted SAM/BAM 1.6 stream readable by HTSlib and samtools. Reads and
qualities are recoverable with `samtools fastq`; the manifest remains the
authoritative run-level audit record.

## Artifact set

With `--format fastq` or `--format fastq.gz`, a run publishes R1, optional R2,
and the manifest. With `--format bam`, it publishes `<prefix>.bam` and
`<prefix>.manifest.json`. The `--save-methdb`, `--save-vcf`, and `--save-truth`
flags add automatically named artifacts beneath `OUTPUT/truth/`.

## Header

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

The `@RG ID` is the run UUID. The manifest records the same UUID and complete
normalized configuration.

## QNAME and standard SAM fields

Both mates share
`<contig>:<one-based-inclusive-start>-<one-based-inclusive-end>:<ordinal-hex>`.
The ordinal is variable-width lowercase hexadecimal without `0x`; flags
identify mate 1 and mate 2. Records are mapped to the simulated origin with an
indel-aware, query-complete CIGAR and reference-forward SEQ/QUAL.

- `MAPQ=60` denotes a known simulated origin, not empirical confidence.
- `AS:i` equals query length, the maximum simulation-origin score.
- `RG:Z` links the record to the run.
- Paired records use standard flags, RNEXT/PNEXT/TLEN, `MQ:i`, and `MC:Z`.

FASTQ records and QNAMEs deliberately do not carry strand-origin truth. The
readable truth annotation exists only in annotated BAM.

## Bisulfite `XG`, `XR`, and `YS` tags

WGBS, RRBS, and TBS records use the Bismark conversion conventions:

- `XG:Z:CT|GA` is the genome conversion and is identical on both mates.
- `XR:Z:CT|GA` is the current read conversion and normally differs between
  paired mates.
- `YS:Z:OT|OB|CTOT|CTOB` is the complete library-strand identity and is
  identical on both mates. Bismark calls this its strand-ID tag.

The paired-end mapping is:

| `YS` library strand | `XG` | R1 `XR` | R2 `XR` | Enabled by |
| --- | --- | --- | --- | --- |
| `OT` | `CT` | `CT` | `GA` | directional and undirectional |
| `OB` | `GA` | `CT` | `GA` | directional and undirectional |
| `CTOT` | `CT` | `GA` | `CT` | undirectional only |
| `CTOB` | `GA` | `GA` | `CT` | undirectional only |

The three tags are always present on bisulfite BAM records. WGS, WES, and TS
records omit them rather than inventing a non-bisulfite conversion value.

## Required `zt:Z` per-read state

`zt` has exactly one state64 character per BAM SEQ base. Character index packs:

| Bits | Meaning |
| --- | --- |
| 0-1 | context: 0 none, 1 CG, 2 CHG, 3 CHH |
| 2 | methylated |
| 3 | bisulfite conversion succeeded |
| 4 | base is produced or affected by a variant event |
| 5 | sequencing error |

For reverse-strand records, `zt` is reversed with BAM SEQ so offsets remain
aligned. No separators or inner wrapper are used.

## Required `zr:B:S` and optional `zf:B:S`

Both use twelve unsigned 16-bit values. `zr` summarizes the current read;
`zf`, when selected, summarizes the complete physical fragment and is copied
to both mate records.

| Index | Meaning |
| ---: | --- |
| 0 | flags |
| 1 | `n_cg` |
| 2 | `n_cg_m` |
| 3 | `n_chg` |
| 4 | `n_chg_m` |
| 5 | `n_chh` |
| 6 | `n_chh_m` |
| 7 | `n_conversion_success` |
| 8 | `n_conversion_failure` |
| 9 | `n_variants` |
| 10 | `n_sequencing_errors` |
| 11 | `n_false_methylation_errors` |

Summary flag bits are: haplotype 0-1, resolved informative strand 2-3
(0 none, 1 Watson, 2 Crick), conversion mode 4-6, ASM 7, SNV 8, insertion 9,
deletion 10, any methylation conversion 11, saturation 13; all other bits are
zero. An explicit simulator capture strand is retained; otherwise the
informative strand is resolved from the fragment chemistry.
Conversion modes are 0 for C-to-T, 1 for G-to-A, and 2 for no bisulfite
chemistry. WGS, WES, and TS records use mode 2 and zero methylation/conversion
counts.

## Optional `zx:Z` complete-fragment realization

`--fragment-realization` implies `--fragment-summary` and BAM output. The same
`zx` value is attached to both mates so either record is independently useful.
Its compact ASCII payload is:

```text
<site-count-hex>.<convertible-count-hex>.<methylation-bits>.<conversion-bits>
```

The two bit vectors are unpadded base64url strings with LSB0 bit order. The
first vector follows all MethDB sites on the physical fragment in fragment
order and records sampled methylation state. The second follows all canonical
convertible template bases and records whether conversion succeeded. Counts
define the exact useful bit lengths and make padding unambiguous.

Sequencing errors are intentionally absent from `zx`: they are read-specific,
already represented by BAM SEQ/QUAL and `zt`, and do not exist in an uncovered
fragment interior. Variant definitions and fixed methylation probabilities are
part of the prepared variant set and methylation profile supplied by VCF and
MethDB; they are not repeated per fragment.

## Interoperability

The stream is unsorted. Coordinate-indexed consumers should run:

```sh
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```
