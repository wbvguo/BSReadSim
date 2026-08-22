# BSReadSim

BSReadSim defaults to a WGBS generator with independent Bernoulli methylation
state, using one supported architecture:

- `htsim-core` (C++17) generates haplotypes, methylation metadata, and raw
  fragment templates;
- `bsreadsim` (Python) validates the core stream, applies bisulfite conversion,
  quality and sequencing-error models, and publishes FASTQ, optional details
  artifacts, and a reproducibility manifest atomically.

The two components communicate only through the current binary protocol,
documented in [docs/protocol.md](docs/protocol.md). Protocol v1 compatibility
and the 0.1.0 implementation remain in Git history, not in the 0.4.0 runtime.
The complete ownership and stage order is in
[docs/architecture.md](docs/architecture.md).

RRBS and TBS are explicit, isolated C++ generation components rather than
alternate Python pipelines. All supported technologies use independent
Bernoulli methylation site states; speculative model-plugin machinery is not
part of the runtime.

| Technology | Fragment source | Released coverage | Count mode | Current gate |
| --- | --- | --- | --- | --- |
| WGBS (default) | reference/haplotype spans | uniform or target GC distribution | `read_pairs` or depth | uniform supports fixed/variable insert and variants; target GC is reference-only, exact for fixed inserts and approximate for variable inserts |
| RRBS | restriction motifs scanned on each haplotype | uniform or external fragment scores | `read_pairs` | physical enzyme-bounded `[insert_min, insert_max]`; VCF and de novo variants supported |
| TBS | verified BED6 targets, then conditional haplotype projection | uniform or exact target output weights | `read_pairs` | fixed physical insert; capture strand, VCF, and de novo variants supported |

All three feed the same protocol and Python post-processing stages. Independent
Bernoulli site state is the default for each technology.

## Build and test

The default CMake build contains only the supported C++ core and its component
tests. It requires a C++17 compiler, CMake 3.20+, zlib, Autoconf, and Make.
HTSlib 1.24 is a commit-pinned recursive Git submodule and is built statically,
so the installed core has no checkout-specific HTSlib `RUNPATH`. For a smaller
fresh checkout, clone both submodule levels shallowly. In an existing checkout,
initialize the same pinned revisions with `--depth 1`:

```sh
git clone --recurse-submodules --shallow-submodules \
  https://github.com/wbvguo/BSReadSim.git
cd BSReadSim

# Existing checkout:
git submodule update --init --recursive --depth 1
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

Build dependencies for a source install are CMake 3.20+, a C++17 compiler,
zlib development files, Autoconf, Make, and a Python C-extension toolchain.
Official source distributions carry the pinned submodule sources. A regular
Python install builds and bundles both the matching private core, its HTSlib
code and license notices, and the required native hot loops:

```sh
python3 -m pip install .
python3 -m unittest discover -s tests -v
```

Institutional builds may explicitly use an installed HTSlib 1.24+ instead:

```sh
cmake -S . -B build -DBSREADSIM_HTSLIB_PROVIDER=SYSTEM
```

## Input guide

Choose the option that names the file's actual serialization. BSReadSim never
infers a parser from `.bed`, `.gz`, or another filename suffix. Contig names in
every biological input must exactly match the first whitespace-delimited name
in the corresponding FASTA header.

| Input role | CLI option | Accepted user-facing contract |
| --- | --- | --- |
| Reference | `-r`, `--reference` | Required FASTA. Sequence is normalized to A/C/G/T/N; its contig order defines the required order for VCF and methylation rows. |
| Predefined variants | `--vcf` | Strict one-sample VCF 4.2/4.3, plain or gzip. `--mutation-rate` must be zero (selected automatically when omitted); REF alleles are checked against FASTA. Use an explicit `--seed` when a run must be repeated or matched to an RRBS catalog. |
| Site methylation | `--cgmap` | Eight-column CGmap: `CHR NUC POS CONTEXT DINUC METH MC NC`. `POS` is one-based and `METH` is `[0,1]` or `na`. |
| Site methylation in BED | `--bed-methyl` | UCSC/ENCODE bedMethyl BED9+2 or BED9+9, plain or gzip. Intervals are zero-based, half-open, and exactly one base; `percentModified` is `[0,100]`. |
| Fixed MethDB | `--methdb` | BSReadSim binary snapshot containing exact normalized site probabilities, source, context, and allele. It is SHA-256 and reference/catalog-identity checked. |
| Allele-specific methylation | `--asm` | Fourteen-column HTSIM ASM profile with one-based target and linked-SNV positions. This is not raw CGmapTools `asm` output. |
| Allele-specific methylation in BED | `--asm-bed` | Project-defined BED6+6, plain or gzip. It carries a one-base target, a one-base linked VCF SNV, REF/ALT, `REF_METH`, and `ALT_METH`. It is not a generic BED3 file. |
| TBS targets | `--targets` | Strict BED6 used with `--technology TBS`; zero-based half-open intervals and strand `+`, `-`, or `.`. |
| WGBS target GC distribution | `--coverage-profile` | Strict TSV containing one probability per line. Fixed inserts are calibrated on post-haplotype opportunities and support VCF, de novo mutation, and ASM; variable insert plus variants remains fail-closed. |
| RRBS candidate scores | `--rrbs-candidates` | Ten-column BED generated first by `bsreadsim catalog rrbs`; it is validated against a regenerated native catalog and is not a general fragment BED. |

`--cgmap` and `--bed-methyl` are alternatives and cannot be supplied together.
With neither option, MethDB uses the configured context-specific Beta
distributions. `--cgmap-pool` accepts either representation and changes the
input from a position-specific overlay into a context-class sampling pool.
`--methdb-seed` fixes probability-catalog and variant-phasing randomness;
`--seed` independently controls fragment selection and downstream realization.
Use `--save-methdb fixed.methdb` to export and immediately run against that
snapshot, or `--methdb fixed.methdb` to load it later.

`--asm` and `--asm-bed` are also alternatives. Either form requires `--vcf`,
and every ASM row must resolve to one exact heterozygous VCF SNV. A site-level
CGmap or bedMethyl profile may be used at the same time; precedence is
`ASM > CGmap/bedMethyl > Beta`.

The following schematic rows are tab-separated. The CGmap and bedMethyl rows
describe the same one-based position 100 / zero-based interval `[99,100)`:

```text
# CGmap: CHR NUC POS CONTEXT DINUC METH MC NC
chr1	C	100	CG	CG	0.35	7	20

# bedMethyl BED9+2: chrom start end name score strand thickStart thickEnd RGB coverage percent
chr1	99	100	m	20	+	99	100	0	20	35

# ASM BED6+6: target BED6, linked SNV interval, REF, ALT, REF_METH, ALT_METH
chr1	99	100	target-100	0	+	499	500	A	G	0.2	0.8
```

The ASM BED example targets the same cytosine and links it to VCF position 500
(`A>G`), represented in BED coordinates as `[499,500)`.

Rows are validated against the complete FASTA before protocol output. Keep VCF,
CGmap, bedMethyl, and ASM rows in FASTA contig order and increasing coordinate
order; duplicate targets and context/reference mismatches fail closed. The
`bsreadsim` CLI computes and records SHA-256 identities automatically. Only a
developer invoking `htsim-core` directly must also provide each matching
`--*-sha256` option.

Full column and validation rules are in
[VCF catalog v1](docs/vcf-catalog-v1.md),
[CGmap profile v1](docs/cgmap-profile-v1.md),
[bedMethyl profile v1](docs/bed-methyl-profile-v1.md),
[HTSIM ASM profile v1](docs/asm-profile-v1.md),
[ASM BED profile v1](docs/asm-bed-profile-v1.md),
[TBS target catalog v1](docs/tbs-catalog-v1.md), and
[RRBS candidate BED v1](docs/rrbs-candidate-bed-v1.md). The WGBS profile is
specified separately in
[WGBS target GC distribution v1](docs/coverage-profile-target-v2.md).

## Run

Run directly from command-line parameters; no input JSON file is required.
Without `--bam`, a run emits FASTQ plus the reproducibility manifest
without constructing Full Details. `-n` counts biological fragments: in
paired-end mode each fragment produces one R1 and one R2 record. FASTQ
identifiers use a variable-width lowercase hexadecimal fragment ordinal; see
the [read-name v2 contract](docs/read-name-v2.md). The output prefix controls
file names, not read identifiers.

```sh
bsreadsim run \
  --reference GRCh38.fa \
  --output runs/wgbs \
  --read-pairs 1000000 \
  --read-length 150 \
  --insert-size 300 \
  --seed 42

# Target-GC WGBS; the TSV SHA-256 is computed and recorded automatically:
bsreadsim run \
  --reference GRCh38.fa \
  --output runs/wgbs-profile \
  --read-pairs 1000000 \
  --read-length 150 \
  --insert-min 150 --insert-mean 400 --insert-max 1000 --insert-stddev 25 \
  --coverage-profile data/experiments/wgbs-gc-target-mock.tsv \
  --seed 42

# Downstream-compatible annotated BAM (replaces FASTQ sidecars):
bsreadsim run -r GRCh38.fa -o runs/bam -n 1000 \
  --seed 42 --bam

# An eight-column CGmap profile supplies position-specific levels:
bsreadsim run -r GRCh38.fa -o runs/cgmap -n 100000 \
  --cgmap sample.cgmap.gz --seed 42

# A standard bedMethyl profile is an explicit alternative to --cgmap:
bsreadsim run -r GRCh38.fa -o runs/bedmethyl -n 100000 \
  --bed-methyl sample.bedmethyl.gz --seed 42

# Freeze the normalized probability catalog, then reuse it with new run seeds:
bsreadsim run -r GRCh38.fa -o runs/freeze -n 100000 \
  --methdb-seed 7 --save-methdb fixed.methdb --seed 42
bsreadsim run -r GRCh38.fa -o runs/replay -n 100000 \
  --methdb-seed 7 --methdb fixed.methdb --seed 99

# Allele-specific BED6+6 retains an exact typed VCF SNV link:
bsreadsim run -r GRCh38.fa -o runs/asm-bed -n 100000 \
  --vcf sample.vcf.gz --mutation-rate 0 \
  --bed-methyl sample.bedmethyl.gz --asm-bed sample.asm.bed.gz --seed 42

# Export exact RRBS fragments for an external GAM. The BED has no hash field;
# no JSON configuration file is involved:
bsreadsim catalog rrbs -r GRCh38.fa -o rrbs-candidates.bed \
  --cut-site 'C|CGG'

# After the predictor replaces only BED column 5 (score), sample directly:
bsreadsim run -r GRCh38.fa -o runs/rrbs-profile -n 100000 \
  --technology RRBS --cut-site 'C|CGG' \
  --rrbs-candidates rrbs-candidates.bed --rrbs-score

# For uniform sampling, omit --rrbs-score; BED score is ignored:
bsreadsim run -r GRCh38.fa -o runs/rrbs-uniform -n 100000 \
  --technology RRBS --cut-site 'C|CGG' \
  --rrbs-candidates rrbs-candidates.bed

# Developer override:
bsreadsim run -r GRCh38.fa -o runs/dev -n 1000 \
  --core build/bin/htsim-core
```

Installed wheels discover their bundled core automatically. Use
`bsreadsim run --help` for all WGBS, RRBS, TBS, model, and execution flags.
The RRBS candidate export and sampling workflow above is CLI-only and needs no
saved JSON configuration. If catalog-defining fragment or variant flags are
changed from their defaults, repeat the same flags on `catalog rrbs` and
`run`; VCF or de novo workflows also repeat the same explicit `--seed`.
`bsreadsim run` is the sole simulation entry point; there is no JSON
run-configuration command or file loader. CLI arguments enter the internal
[normalization schema](src/bsreadsim/run-config.schema.json), which provides
shared defaults and scientific validation without becoming a user input
format. Immutable biological inputs and model artifacts are SHA-256 verified
before launch; the deliberately mutable RRBS candidate-score exchange is
instead checked by exact regenerated-row matching. Successful data files become
visible only when a complete manifest can be committed.

The output contract is documented in
[docs/output-policy.md](docs/output-policy.md). A normal run emits R1, optional
R2, and the reproducibility manifest. There is no production/debug selector
and no per-fragment Details JSONL artifact.

`--bam` switches the data product to `<prefix>.bam` plus the
manifest; FASTQ sidecars are not emitted because `samtools fastq` can recover
the reads and qualities. The unsorted SAM/BAM 1.6 records retain source
coordinates, indel-aware CIGAR, paired flags/mate fields, qualities, fixed
`zt` and `zr` tags, optional `zf` summaries selected by
`--fragment-summary`, and optional complete-fragment `zx` realizations selected
by `--fragment-realization`. The latter implies BAM and `zf`.
Coordinate-indexed consumers should run
`samtools sort` and then `samtools index`. The exact contract is in
[docs/bam-v3.md](docs/bam-v3.md).

`--methylation-model bilstm` currently emits a warning and records an effective
fallback to the built-in Bernoulli model. The stable model interface is ready
for a future correlated plugin; this release does not claim BiLSTM behavior.

For storage/speed tuning, `output.gzip_level` accepts `0` through `9` and
defaults to `6`. A lower value such as `1` is an explicit, manifest-visible
choice: it preserves decompressed FASTQ content but produces larger files.

The two bounded compute stages are configured independently. For example, the
measured FASTQ-only WGBS setting on the benchmark host was:

```sh
bsreadsim run -r GRCh38.fa -o runs/wgbs -n 1000000 \
  --core-workers 2 \
  --workers 2 \
  --chunk-size 64 \
  --max-in-flight-fragments 256
```

`core_workers` prepares C++ protocol frames; `workers` post-processes fragments
in Python. With `workers = 1`, all policies run through the raw-payload batch
core inline from one reusable local buffer without a process pool. Eligible
uniform-policy FASTQ-only blocks use its native NumPy branch; other policies use
the typed-object branch inside that same core. Larger worker counts use
spawned processes and shared-memory slots. There is no one-to-one assignment
between C++ and Python workers, and the defaults remain one each for
conservative, low-process-count runs. The example is the measured optimum
for the native, uniform-policy, FASTQ-only benchmark on one host; model cost,
Details output, compression, CPU topology, and memory budget can change the best
ratio.

Experiment drivers, benchmark harnesses, and all generated research data belong
under [data/](data/README.md). Large inputs and outputs remain outside Git.
Runnable, deterministic RRBS and TBS mock configurations are documented in
[data/experiments/](data/experiments/README.md).

## History

Retired implementations are kept in Git rather than in the active source
tree. Protocol v1 and the historical `htsim`/`rrcut` programs remain
recoverable from the `v0.1.0` tag.

BSReadSim is available under the MIT license. Local design-paper assets, when
present, live under the ignored `docs/design/` directory.
