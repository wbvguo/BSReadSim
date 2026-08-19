# CGmap methylation profile v1

Status: normative predefined methylation-level input contract.

For the equivalent BED representation, use
[bedMethyl profile v1](bed-methyl-profile-v1.md). CGmap and bedMethyl inputs are
mutually exclusive and normalize to the same position-specific overlay.

The C++ core is the sole owner of CGmap parsing and MethDB overlay. Python
hashes the input during preparation, projects the path and SHA-256 to the core,
and consumes only typed protocol methylation sites. It MUST NOT parse CGmap or
construct a second methylation database.

## Accepted v1 text contract

Input is a plain or gzip-compressed regular file whose raw bytes match the
declared SHA-256. Empty lines and lines beginning with `#` are ignored. Every
data row has exactly the standard eight tab-separated CGmap fields:

```text
CHR  NUC  POS  CONTEXT  DINUC  METH  MC  NC
```

- `CHR` must name a FASTA contig. Rows follow FASTA contig order and strictly
  increasing positions within a contig; duplicate or returning contigs fail.
- `NUC` is `C` or `G`; `POS` is a positive one-based `uint32` coordinate.
- `CONTEXT` is `CG`, `CHG`, or `CHH`; strand is derived from `NUC` and encoded
  as the corresponding protocol context.
- `DINUC` is one of `CA`, `CC`, `CG`, or `CT`. Its cytosine-oriented second
  base must agree with the verified FASTA (reverse-complementing the upstream
  reference base for a `G` row), and `CG` must agree with the `CG` context
  class.
- `METH` is a finite decimal in `[0,1]` or lowercase `na`. `na` means that no
  predefined level is available and the site uses the Beta fallback.
- `MC` and `NC` are normalized as `uint32` decimal counts and must satisfy
  `MC <= NC`; this deliberately permits deep-coverage counts beyond the
  historical format table's 12-bit recommendation.

Before the protocol header, every row is checked against the full reference
contig: its coordinate, C/G strand, and three-class context must agree exactly.
Unknown contigs, malformed values, unresolved boundary/N contexts, and any
reference mismatch fail with zero protocol bytes. Rows for CHG/CHH are still
validated when `collect_non_cpg=false`, but they are not emitted.

## MethDB overlay

The base catalog is first populated with deterministic context-specific Beta
levels. A defined CGmap value then replaces the level at the same reference
position and exact reference context and changes its protocol source from
`BETA` to `CGMAP`. Missing positions and `na` retain the Beta value.

When `cgmap_pool=true`, the direct positional overlay above is replaced by the
addressed, per-contig context-class sampling contract in
[cgmap-pool-v1.md](cgmap-pool-v1.md). This mode is explicit because it no
longer preserves the input level at the same genomic position.

With a diploid VCF projection, the CGmap value applies only to a site that is
reference-equivalent on that haplotype. A variant-created/changed context or an
inserted site retains its independently addressed Beta value. Therefore a
heterozygous context change can yield one `CGMAP` reference-haplotype site and
one `BETA` alternate-haplotype site. The reference MethDB is filled before
variant-boundary updates; no field is overloaded.

## Resource and width boundary

Rows are normalized into an unlinked fixed-record temporary spool. Only one
contig's records are materialized while its MethDB is built; a whole-genome
CGmap is not retained in RAM. Positions and per-contig row counts are `uint32`.
Total rows, spool offsets, and byte-size arithmetic are checked `uint64`; no
u96 representation is used.
