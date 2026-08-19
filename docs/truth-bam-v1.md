# Truth-aligned BAM v1

Status: released opt-in output contract.

`bsreadsim run --truth-bam` adds `<prefix>.truth.bam` to the normal FASTQ and
manifest outputs. It is a SAM/BAM 1.6 alignment file serialized directly by
the pinned HTSlib; no intermediate SAM file is created. The manifest role is
`truth_bam`, and its `record_count` is the number of mates rather than the
number of biological fragments.

This is a truth-origin alignment, not a realignment of the simulated reads.
Coordinates and CIGAR come from HTSIM's haplotype-to-reference projection;
bisulfite conversion and sequencing substitutions remain in `SEQ` as observed
read bases. The contract follows the official
[SAM/BAM specification](https://github.com/samtools/hts-specs).

## CLI and output modes

```sh
bsreadsim run \
  --reference GRCh38.fa \
  --output runs/wgbs-bam \
  --read-pairs 1000000 \
  --read-length 150 \
  --insert-size 300 \
  --seed 42 \
  --truth-bam
```

`--truth-bam` is independent of `--mode`:

- production plus `--truth-bam` emits FASTQ, truth BAM, and the manifest, but
  does not emit Truth JSONL;
- debug plus `--truth-bam` emits FASTQ, Truth JSONL, truth BAM, and the
  manifest;
- without `--truth-bam`, the existing production/debug artifact sets are
  unchanged.

Truth BAM requires the Full Truth coordinate projection internally. It
therefore does not use the FASTQ-only common-column fast path and can cost
substantially more CPU and shared-memory bandwidth than ordinary production
output. `--gzip-level 0..9` also selects the BAM BGZF compression level.

## Header and ordering

The header contains:

- `@HD VN:1.6 SO:unsorted`;
- one `@SQ` per FASTA contig, in FASTA order, with exact name and length;
- one `@RG`, whose ID is the run UUID and whose sample name is the output
  prefix;
- one `@PG ID:bsreadsim` with the Python package version;
- an `@CO` declaration that MAPQ 60 is a simulated-origin marker, not a
  calibrated alignment confidence.

Because contig names become SAM `SN`/`RNAME` values, `--truth-bam` also
requires the SAM 1.6 reference-name character contract. A FASTA name containing
quotes, commas, backslashes, or brackets, or beginning with `*` or `=`, fails
before output staging rather than creating a BAM that other tools parse
inconsistently.

Records are emitted by fragment ordinal and paired records are adjacent. The
file is deliberately not coordinate-sorted and no BAI/CSI index is emitted.
Tools that require indexed genomic access should create a derived file:

```sh
samtools quickcheck runs/wgbs-bam/sim.truth.bam
samtools sort -o runs/wgbs-bam/sim.truth.sorted.bam \
  runs/wgbs-bam/sim.truth.bam
samtools index runs/wgbs-bam/sim.truth.sorted.bam
```

Do not index the original unsorted BAM. Sorting changes the BAM byte identity,
so the derived sorted file is intentionally outside the original run
manifest.

## Alignment records

For a fragment on `chr_id` spanning the one-based inclusive interval
`c1-c4` with ordinal `N`, both mates use QNAME `chr_id:c1-c4:N`. FASTQ appends
`/1` or `/2` to that core identifier. This shared BAM QNAME is required for
standard paired-read consumers and is specified in
[read-name v1](read-name-v1.md).

The record contract is:

- all records are primary and mapped to their simulated source contig;
- paired records set `0x1`, `0x2`, `0x40`/`0x80`, `0x10`, and `0x20` according
  to their generated mate geometry; single-end records use only `0x10` when
  reverse-strand;
- `POS` is the one-based first projected reference position;
- CIGAR uses `M` for every reference-projected read base, `I` for a haplotype
  insertion base, and `D` for gaps between consecutive projected positions;
- reverse-strand BAM `SEQ` is reverse-complemented and `QUAL` is reversed, as
  required by SAM; converting the BAM back to FASTQ therefore restores the
  exact emitted FASTQ read;
- paired `RNEXT`, `PNEXT`, and signed `TLEN` describe the other generated mate;
- MAPQ is fixed at 60 to pass conventional mapped-read filters, but it must not
  be interpreted as an empirical or aligner-calibrated probability;
- every record has standard `RG` and `PG` tags; paired records also carry the
  standard `MC` mate-CIGAR tag.

A mate consisting entirely of inserted haplotype bases is represented by a
pure-`I` CIGAR at its insertion anchor. HTSlib writes the BAM bin using the SAM
rule that an alignment with no reference-consuming CIGAR operation has an
effective reference span of one base.

## Downstream boundary

The file is consumable by general SAM/BAM readers and workflows that accept
an unsorted mapped BAM. It intentionally omits aligner-derived `AS`, `NM`, and
`MD` tags and bisulfite-aligner-specific tags such as `XM`, `XR`, or `XG`.
Downstream tools that require those tags must calculate their own reference-
aware annotations or run their documented preparation step. Variant callers
or methylation callers should also be configured with the same reference and
must not treat synthetic MAPQ as mapping uncertainty.

Publication remains transactional: BAM finalization must write a valid BGZF
EOF block, all output counts and SHA-256 identities must reconcile, and the
complete manifest is linked last. A malformed SAM record, HTSlib failure, or
output collision aborts the entire staged artifact set.
