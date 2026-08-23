# htsim core

`htsim-core` owns reference loading and biological fragment generation. It
writes the sole supported binary protocol to stdout; Python owns sequencing
errors, FASTQ and optional details output, alignment formatting, and publication.
For BAM, Python streams SAM records back through the executable's narrow
HTSlib serializer; htsim still does not choose output paths.

## Source layout

`htsim/src` is intentionally flat. Each implementation domain owns one public
header and one source file:

- `core.{h,cpp}`: command-line parsing and top-level orchestration
- `bam.{h,cpp}`: strict streaming SAM-to-BAM serialization through HTSlib
- `variant.{h,cpp}`: variants, mutations, and haplotype projection
- `reference.{h,cpp}`: FASTA loading and immutable contig storage
- `methdb.{h,cpp}`: methylation sites, CGmap/bedMethyl, ASM/ASM BED, and diploid MethDB catalogs
- `fragment.{h,cpp}`: allocation, insert lengths, and fragment construction
- `wgbs.{h,cpp}`: reference and haplotype-aware WGBS selection
- `rrbs.{h,cpp}`: restriction-site discovery, catalog, and sampling
- `tbs.{h,cpp}`: target BED projection, catalog, and sampling
- `protocol.{h,cpp}`: column projection, ordered emission, and wire I/O
- `utilities.{h,cpp}`: RNG, distributions, hashes, and HTSlib-backed text snapshots
- `types.h`: shared biological value types
- `htsim.cpp`: the minimal executable entry point

This is the smallest conventional C++17 layout that still separates public
declarations from implementation. Named section markers inside each domain
file preserve the smaller component boundaries without nested source folders.

The root build exposes the domains as one private `htsim` target, so
the directory layout documents responsibilities without duplicating the link
dependency graph.

Component-owned C++ tests live under `htsim/tests` and mirror the source
domains. Cross-language and end-to-end tests remain under the repository-level
`tests/integration` suite.

Build from the repository root:

```sh
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The build uses the commit-pinned, top-level `htslib/` submodule. HTSlib
owns format detection, gzip/BGZF decoding, and VCF syntax; htsim retains the
raw-byte digest, stable-descriptor, and simulator-specific semantic checks.
