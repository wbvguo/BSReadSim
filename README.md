# BSReadSim

BSReadSim generates reproducible synthetic reads for bisulfite sequencing
(WGBS, RRBS, and TBS) and non-bisulfite sequencing (WGS, WES, and targeted
sequencing). It combines configurable variation, methylation, fragment,
chemistry, quality, and sequencing-error models, then writes FASTQ or
origin-annotated BAM. Every completed run includes a manifest describing the
effective settings, inputs, seeds, counts, and output checksums.

## Highlights

- WGBS, RRBS, TBS, WGS, WES, and TS through one command-line interface
- CGmap, bedMethyl, MethBED, MethDB, CGmapTools ASS, ASM BED, diploid VCF,
  target BED, and scored RRBS inputs
- fixed or variable insert lengths and single- or paired-end reads
- deterministic output across supported thread counts for a fixed configuration
- manifests containing effective settings, seeds, counts, and SHA-256 identities
- strict input validation with no parser guessing from filename suffixes

## Documentation

The searchable user guide is available at
[wbvguo.github.io/BSReadSim](https://wbvguo.github.io/BSReadSim/).

Repository versions of the main entry points are also available directly:

- [Installation](docs/getting-started/installation.md)
- [Quick start](docs/getting-started/quickstart.md)
- [Simulation](docs/simulation/workflow.md)
- [Other assays](docs/simulation/other-assays.md)
- [Outputs](docs/outputs/index.md)
- [CLI parameters and defaults](docs/reference/cli.md)
- [Input file formats](docs/reference/formats.md)

The complete site starts at [docs/index.md](docs/index.md).

<details>
<summary>Preview the documentation locally</summary>

Run the following commands from the repository root:

```bash
python3 -m pip install --upgrade "pip>=25.1"
python3 -m pip install --group docs
python3 -m mkdocs serve
```

Then open <http://127.0.0.1:8700/> in a browser. MkDocs automatically reloads
the site when documentation files change; press `Ctrl+C` to stop the server.

</details>

## Install

BSReadSim supports Linux x86-64 with Python 3.10 or newer. Windows users
should build and run it inside WSL2. A source installation requires CMake 3.20
or newer, a C++17 toolchain, zlib development files, Autoconf, and Make.

```bash
git clone --recurse-submodules --shallow-submodules \
  https://github.com/wbvguo/BSReadSim.git
cd BSReadSim
python3 -m pip install .
```

This installs into the Python environment of your choice; `venv` is optional.
A Bioconda package is planned as the second supported installation path, but
has not been published yet. BSReadSim will not be distributed through PyPI.

## Quick example

```bash
bsreadsim run wgbs \
  --reference GRCh38.fa \
  --output runs/example \
  --reads 100000 \
  --read-length 150 \
  --insert-mean 300 \
  --insert-sd 0 \
  --mutation-rate 0 \
  --seed 42
```

This writes 100,000 read records (50,000 paired-end read pairs) and a
reproducibility manifest without introducing variants. Add `--format bam` to
write an unsorted, origin-annotated BAM instead of FASTQ.

Use `bsreadsim run wgbs --help` for the installed command reference.

For non-bisulfite whole-genome reads, use the same sampling and sequencing models
without constructing a methylome or applying bisulfite conversion:

```bash
bsreadsim run wgs \
  --reference GRCh38.fa \
  --output runs/wgs \
  --reads 100000 \
  --mutation-rate 0 \
  --seed 42
```

## Citation

If BSReadSim contributes to research, cite the
[BSReadSim preprint](https://doi.org/10.1101/2024.12.24.627620) and report the
software version used. Machine-readable metadata is provided in
[`CITATION.cff`](CITATION.cff).

BSReadSim is available under the MIT license.
