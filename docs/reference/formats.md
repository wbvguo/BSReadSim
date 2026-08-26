# File formats

Choose an input parser with the option name, not only the filename suffix.
For example, `.bed` does not distinguish bedMethyl, ASM BED, RRBS candidates,
and TBS targets. BSReadSim validates the selected format strictly and records
the resolved file identity in each run manifest.

In this documentation, a **variant set** is the normalized, phased collection
of variants used to prepare the simulated haplotypes. A **methylation profile**
is the resolved collection of site-level methylation probabilities. When saved
as VCF and MethDB, respectively, they are reusable **simulation truth
artifacts**.

## Format map

| Data | Consume with | Produce with | Exact contract |
| --- | --- | --- | --- |
| Reference FASTA | `--reference` | external | A/C/G/T/N sequence; first header token is the contig name |
| Prepared variant set | `--vcf` | `build variants`, `--save-vcf`, `--save-truth` | [VCF variant sets](formats/vcf-variants.md) |
| CGmap methylation | `--cgmap` | external | [CGmap profile](formats/cgmap-profile.md) |
| bedMethyl methylation | `--bed-methyl` | external | [bedMethyl profile](formats/bed-methyl-profile.md) |
| Prepared methylation profile | `--methdb` | `build methdb`, `--save-methdb`, `--save-truth`; inspect with `export methdb` | [MethDB profile and exported BED](formats/methdb.md) |
| Allele-specific methylation | `--asm` | external | [ASM profile](formats/asm-profile.md) |
| Allele-specific methylation BED | `--asm-bed` | external | [ASM BED profile](formats/asm-bed-profile.md) |
| WGBS/WGS target-GC weights | `--gc-profile` | external | [Target-GC profile](formats/coverage-profile-target.md) |
| RRBS candidate weights | `--rrbs-candidates` | `build rrbs`, followed by external score editing | [Candidate BED](formats/rrbs-candidate-bed.md) |
| TBS/WES/TS capture targets | `--targets` | external | [Target BED](formats/tbs-catalog.md) |
| Quality Markov model | `--quality-model` | external | [Sequencing models](formats/sequencing-models.md) |
| Quality-conditioned error | `--error-model` | external | [Sequencing models](formats/sequencing-models.md) |

## Coordinate conventions

Do not convert every input to one coordinate convention blindly:

| Format family | Coordinates |
| --- | --- |
| FASTA | sequence indexed internally from zero; no coordinate column |
| VCF and CGmap | one-based positions |
| BED-derived inputs | zero-based, half-open intervals |
| FASTQ read names | one-based, inclusive displayed fragment coordinates |
| BAM/SAM | standard SAM coordinate and CIGAR conventions |

Each linked contract is authoritative for edge cases. Contig names are always
case-sensitive and must match the first whitespace-delimited token in the
FASTA header.

## Compression and file identity

Compression support belongs to each contract; a `.gz` suffix alone does not
select or guarantee a parser. In particular, `build variants` always requires
a new `.vcf.gz` destination and writes deterministic BGZF. MethDB is already
an internally compressed opaque binary format and should not be wrapped or
edited as text. `export methdb` decodes it to deterministic BGZF `.bed.gz` by
default, or plain `.bed` with `--no-compression`.

Users provide paths only. BSReadSim computes SHA-256 identities internally,
checks model metadata where applicable, and records the identities in the run
manifest. There are no checksum CLI parameters.

## Read and evidence formats

| Output | Selected with | Contract |
| --- | --- | --- |
| FASTQ or FASTQ.gz | `--format fastq` or `--format fastq.gz` | [Read names](formats/read-name.md) and standard four-line records |
| Annotated BAM | `--format bam` | [BAM fields and tags](formats/bam.md) |
| Run file set | every `run` | [Publication, truth, and manifest](../outputs/index.md) |

Invalid or ambiguous inputs fail before final output publication. When an
error names a row or field, fix the selected format rather than changing the
filename extension.
