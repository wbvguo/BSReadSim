# Read-name contract v1

Status: released output contract.

BSReadSim assigns the final read name only when FASTQ or truth BAM records are
formatted. HTSIM supplies the contig, fragment reference envelope, and
fragment ordinal; it does not serialize a presentation-specific identifier in
the binary protocol.

The contract follows slide 4 of `docs/design/bsreadsim.pptx`:

```text
FASTQ R1: @<chr_id>:<c1>-<c4>:<read_num>/1
FASTQ R2: @<chr_id>:<c1>-<c4>:<read_num>/2
BAM QNAME: <chr_id>:<c1>-<c4>:<read_num>
```

The fields are:

- `chr_id`: the exact FASTA contig name;
- `c1`: the fragment's leftmost reference coordinate, one-based and
  inclusive;
- `c4`: the fragment's rightmost reference coordinate, one-based and
  inclusive;
- `read_num`: the zero-based fragment ordinal within the current simulation;
- the optional `/1` or `/2`: the FASTQ mate number.

Internally, HTSIM and the protocol use a zero-based half-open fragment envelope
`[reference_start, reference_end)`. A non-empty envelope is rendered as
`c1 = reference_start + 1` and `c4 = reference_end`. A physical fragment made
only of inserted bases has a zero-width envelope; both fields are its one-based
insertion anchor, `reference_start + 1`.

For example, internal fragment `[100,108)` with ordinal 7 produces
`@chr1:101-108:7/1` and, when paired, `@chr1:101-108:7/2`. Both corresponding
BAM records use the shared QNAME `chr1:101-108:7`, as required by paired-read
consumers.

`output.prefix` remains the output filename prefix and the truth-BAM sample
name. It is deliberately not part of the read identifier. Contig names used in
read names must be non-empty printable text without whitespace;
truth BAM additionally applies the SAM 1.6 reference-name and QNAME limits.

Consumers should treat the identifier as opaque unless they implement this
versioned contract. If a contig name itself contains `:`, parse the fixed
coordinate and ordinal fields from the right.
