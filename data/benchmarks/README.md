# Benchmarking

[`benchmark.py`](benchmark.py) measures the only current end-to-end pipeline.
It accepts the same parameter-style arguments as `bsreadsim run` after `--`,
makes every referenced path absolute, materializes an omitted seed once, and
writes each repetition to an isolated directory. Every repetition must produce
identical scientific counts and byte-identical artifact hashes before timing
is reported. No run-configuration file is read or generated.

For a one-CPU FASTQ measurement:

```sh
taskset -c 0 python3 data/benchmarks/benchmark.py \
  --core build/bin/htsim-core \
  --warmup --repetitions 5 \
  --json-out /tmp/bsreadsim-production.json -- \
  -r /path/to/reference.fa -n 500000 --seed 17 \
  --read-length 150 --insert-size 300 --workers 1 --core-workers 1
```

Use `--bam` after `--` to measure annotated BAM instead of FASTQ.
Worker counts, C++ worker counts, batch sizes, compression, and fragment counts
come from the arguments after `--`, so the exact workload is preserved in
every run manifest. The benchmark owns `-o/--output` and `--core`; specify
those before `--`.
`--memory-sample-ms 25` (the default) samples aggregate Linux process-tree RSS
and PSS; set it to `0` to disable sampling. The report also records raw wall
times, CPU accounting, process I/O counters, exact output sizes and hashes,
the core binary digest, Git state, and unaggregated measurements.

For comparisons, benchmark separate commits or installed wheels in alternating
order on the same host and input. Set `BSREADSIM_BENCHMARK_USE_INSTALLED=1` to
prevent the harness from importing the working tree. Do not compare only the
reported medians: retain all paired samples and reject any output-hash or count
difference first. Historical protocol-promotion harnesses and reports remain
available in Git history; they are not part of the current source tree.

[`compare_wgbs_coverage.py`](compare_wgbs_coverage.py) is the narrower paired
uniform-versus-target-GC harness. It refuses to run unless inherited affinity
contains exactly one logical CPU, warms both modes, and measures them in
`U-P / P-U / U-P` adjacent pairs. It reports two lanes: raw C++ core production
to `/dev/null`, and the complete production FASTQ path with gzip level 1.
Within each end-to-end mode, all repetitions must retain identical counts and
FASTQ hashes. Formal runs also require the Python native extension and record
its resolved path and SHA-256; this prevents a pure-Python fallback from being
misreported as an algorithmic throughput regression.

```sh
taskset -c 0 python3 data/benchmarks/compare_wgbs_coverage.py \
  --reference /path/to/chr21.fa.gz \
  --profile data/experiments/wgbs-gc-target-mock.tsv \
  --core build/bin/htsim-core \
  --fragments 500000 --warmup-fragments 50000 \
  --output-json /tmp/wgbs-coverage-paired.json
```

## Retained chr21 one-CPU comparison

[`chr21-wgbs-uniform-vs-target-1cpu.json`](chr21-wgbs-uniform-vs-target-1cpu.json)
retains the interleaved measurements from commit `4b21577`. The workload was
500,000 PE150 fragments (1,000,000 reads), fixed insert 300, no mutation or
sequencing error, gzip level 1, and one C++ plus one Python worker, with every
process and thread pinned to logical CPU 0. Both modes were warmed before the
`U-P / P-U / U-P` formal order.

| Lane | Uniform fragments/s | Target fragments/s | Target throughput change | Target wall overhead |
| --- | ---: | ---: | ---: | ---: |
| C++ producer to `/dev/null` | 40,481.99 | 33,805.73 | -16.49% | +19.75% |
| production FASTQ end to end | 20,352.30 | 18,442.23 | -9.39% | +10.36% |

The target run rejected 405,921 proposals for 500,000 outputs, an observed
acceptance of 55.19%. Core pairwise profile/uniform throughput ratios were
`0.8406`, `0.8350`, and `0.8375`. End-to-end ratios were `0.9061`, `0.8988`,
and `1.0126`; the third pair had elevated system/I/O time for both modes and a
particularly slow uniform sample. It remains in the report. The median result
supports an approximately 9% to 10% end-to-end target-profile penalty on this
host and workload, but it is not a universal constant: target acceptance,
reference length, fragment count, compression, and host contention all matter.

## Current 0.3 one-million-read evidence

The reports dated 2026-08-19 bind the reconciled 0.3 runtime to its core,
native-extension, input, output, and pre-squash source identities:

- `bsreadsim-0.3.0-wgbs-components-1m-reads-2026-08-19.json` separates the
  C++ protocol producer from FASTQ end to end for uniform and target-GC
  variable-insert WGBS. Each cell has three interleaved 1,000,000-read runs.
- `bsreadsim-0.3.0-post-refactor-bam-1m-reads-2026-08-20.json` records
  the reorganized annotated BAM-only path: an uninstrumented one-CPU result,
  a matching Full Truth core-producer lane, process-boundary timings, process-
  tree RSS/PSS, native artifact identities, and byte-identical BAM records
  across the typed-safe and columnar protocol batch bounds.
- `archive/bsreadsim-0.3.0-output-modes-1m-reads-2026-08-19.json` is retained
  as pre-cleanup historical evidence. It describes the retired production/debug
  and Truth JSONL surface and is not a current performance gate.
- `bsreadsim-0.3.0-technologies-1m-reads-2026-08-19.json` compares WGBS, RRBS,
  and TBS in a balanced three-round order on one deterministic synthetic
  corpus. Core fragment, mate, template-base, methylation-site, and rejection
  counts are equal across technologies.

Generate the exact 10 Mb technology corpus outside Git with:

```sh
python3 data/benchmarks/generate_technology_corpus.py \
  --output-directory /tmp/bsreadsim-technology-corpus
```

The generator freezes and verifies both reference and target SHA-256 digests.
The corpus is a throughput control, not an empirical biological model.

The generic driver accepts `--warmup-fragments N`, which implies `--warmup`
and replaces only an explicit `--read-pairs` count for the warmup run. This
avoids duplicating a full 1,000,000-read truth-BAM output merely to
warm code and filesystem caches. Truth-BAM output is an accepted benchmark
role and remains covered by the same count/hash repetition gate.
