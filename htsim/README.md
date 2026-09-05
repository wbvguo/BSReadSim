# htsim core

`htsim` owns the reference, the variants, the methylation profile, and
where each fragment comes from. It writes a fragment stream to stdout (format
in `src/stream.h`); Python owns methylation states, bisulfite conversion,
qualities, sequencing errors, FASTQ and BAM formatting, and publication.
For BAM, Python streams SAM records back through `htsim --sam-to-bam`.

The whole-genome sampler backs WGBS and WGS, the restriction-fragment sampler
backs RRBS, and the target sampler backs TBS, WES, and TS. Standard
technologies emit no methylation sites, so Python skips methylation and
bisulfite conversion for them.

## Running htsim on its own

Installing the `bsreadsim` package also installs the `htsim` command (the
bundled executable); a source build leaves it at `build/bin/htsim`. Without
bisulfite (WGS, WES, TS), `htsim` can also write reads itself, to files named
after `--output PREFIX`:

```sh
htsim --reference genome.fa --technology WGS --fragments 1000000 \
    --output sim --format fastq.gz,bam --threads 4
```

writes `sim.R1.fastq.gz`, `sim.R2.fastq.gz`, and `sim.bam`. `--format` takes one
or more of `fastq.gz` (the default), `fastq` (`sim.R1.fastq`, `sim.R2.fastq`),
`bam`, and `sam`. Only `--reference` is required; every other option has the default
shown by `htsim --help` (those of the parameter structs in `src/struct.h`).
The reads are cut, named, and aligned as the Python package cuts, names, and
aligns them (`src/reads.h`), with one base quality (`--phred`) and uniform
substitutions (`--error-rate`). BAM and SAM are unsorted, with mates adjacent, and
has the package's alignment fields and AS/MQ/MC tags but not its truth tags
(zt, zr, zf, zx); those, quality and error models, and bisulfite assays need the
`bsreadsim` package.

## Source layout

The core is derived from WGSIM and grew out of the 2024 htsim. One header and
one source file per topic, in reading order:

- `struct`: constants, lookup tables, and the shared structs, each with the
  small functions that read or encode its fields
- `utility`: general tools that know nothing of BSReadSim's data: errors,
  number parsing, the random generator, threads, text input and output, and
  pipes
- `haplo`: VCF input, de novo variants, the two haplotypes of a contig, the
  walk that reads a fragment from a haplotype, and VCF output
- `methdb`: the methylation level of every site (Beta or pooled draws, sites
  changed by variants), the MethDB snapshot, and its BED export
- `profile`: methylation inputs (CGmap, bedMethyl, MethBG, MethBED, and ASM),
  checked against the reference and applied to the methdb
- `mode`: capture targets, the GC profile, how many fragments each contig
  gets, and where they are placed
- `rrcut`: RRBS cut sites and restriction-fragment candidates
- `stream`: the fragment stream read by Python
- `reads`: reads written by htsim itself (`--output`, `--format`; no bisulfite)
- `option`: command-line options and their allowed combinations
- `htsim.cpp`: reading the inputs, one contig at a time, and the commands

A contig is processed at a time, and `run` prepares every contig the way
`methdb-build` does (`prepare_contig`: variants, haplotypes, methylome, with
every input checked) before sampling its fragments, whatever its share of them.
Each haplotype is a `mut_t` array with one
16-bit word per reference base (the base, or a substitution, a deletion, or up
to four bases inserted after it). `meth_vec` holds one record per reference
site in position order (and per C/G near a variant, which it may turn into a
site), with its context, level, and source on both haplotypes;
`snpmeth_map` holds the variants and the sites of the bases they add. Sites
are found by binary search (`find_site`); fragment generation, the hot path,
starts from a coarse index with one entry per 64 bases (`site_index`) and walks
forward. On GRCh38 chr1 this is as fast as the former per-base index
(`posidx_arr`, 4 bytes per base) and needs 27% less memory.

RRBS and targeted runs build the methylome only where their fragments fall
(`frag_domain::reach`: the RRBS candidates, and each target's centre plus or
minus the longest insert and 6 centre SDs, a fragment reaching beyond being
drawn again), and at the ASM targets, which are checked wherever they are; a
MethDB output still holds every site. On GRCh38 chr1 at 4 threads, the setup
takes 1.9 s and 1.3 GB instead of 7.1 s and 2.4 GB for TBS, 7.4 s and 1.8 GB
instead of 10.2 s and 2.4 GB for RRBS, and 1.3 GB instead of 2.4 GB for
`--cpg-only` WGBS. Most of what remains is the two haplotypes (2 bytes per base
each).

## Reproducibility

Every random decision draws from a generator seeded by hashing a stage seed
with what is being decided: a contig and a fragment's ordinal within it, or a
contig and a methylation site. Output therefore depends only on the version,
the platform, the options, and the seeds, never on the thread count.

## Building

From the repository root:

```sh
git submodule update --init --recursive
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The build uses the pinned `htslib/` submodule for gzip/BGZF input and output,
MD5 digests, and SAM-to-BAM conversion. The semantic suite under
`tests/semantic` checks the core's behavior through the public interface.
