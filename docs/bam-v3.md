# annotated BAM v3

Status: normative contract for the detailed BSReadSim details artifact.

Annotated BAM replaces FASTQ sidecars when `--bam` is selected. It is an
unsorted SAM/BAM 1.6 stream readable by HTSlib and samtools. Reads and
qualities are recoverable with `samtools fastq`; the manifest remains the
authoritative run-level audit record.

## Artifact set

Without `--bam`, a run publishes R1, optional R2, and the manifest. With
`--bam`, it publishes `<prefix>.bam` and `<prefix>.manifest.json`. A fixed
MethDB snapshot is an additional, explicit artifact only when
`--save-methdb PATH` is selected.

## Header

Ordinary `@HD`, `@SQ`, `@RG`, and `@PG` records are followed by deterministic
contract comments:

```text
@CO    BSREADSIM_BAM=bsreadsim-bam-v3
@CO    BSREADSIM_ZT=state64-v1;ALPHABET=ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_
@CO    BSREADSIM_ZR=u16x12-v1;REQUIRED=1
@CO    BSREADSIM_ZF=u16x12-v1;ENABLED=0|1
@CO    BSREADSIM_ZX=packed-b64url-v1;ENABLED=0|1;BIT_ORDER=LSB0
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

Summary flag bits are: haplotype 0-1, capture strand 2-3, conversion mode
4-6, ASM 7, SNV 8, insertion 9, deletion 10, any methylation conversion 11,
saturation 13; all other bits are zero.

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
catalog information supplied by VCF/MethDB, not repeated per fragment.

## Transaction and interoperability

BAM is streamed through the HTSlib-backed `htsim-core --sam-to-bam` helper.
Finalization, digest/count reconciliation, and manifest-last publication form
one transaction. The stream is unsorted; coordinate-indexed consumers run:

```sh
samtools sort -o sample.sorted.bam sample.bam
samtools index sample.sorted.bam
```
