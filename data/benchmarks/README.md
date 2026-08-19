# Benchmarking

[`benchmark.py`](benchmark.py) measures the only current end-to-end pipeline.
It accepts the same parameter-style arguments as `bsreadsim run` after `--`,
makes every referenced path absolute, materializes an omitted seed once, and
writes each repetition to an isolated directory. Every repetition must produce
identical scientific counts and byte-identical FASTQ/Truth hashes before timing
is reported. No run-configuration file is read or generated.

For a one-CPU production measurement:

```sh
taskset -c 0 python3 data/benchmarks/benchmark.py \
  --core build/bin/htsim-core \
  --mode production --warmup --repetitions 5 \
  --json-out /tmp/bsreadsim-production.json -- \
  -r /path/to/reference.fa -n 500000 --seed 17 \
  --read-length 150 --insert-size 300 --workers 1 --core-workers 1
```

Use `--mode debug` only when measuring Full Truth. Worker counts, C++ worker
counts, batch sizes, compression, and fragment counts come from the arguments
after `--`, so the exact workload is preserved in every run manifest. The
benchmark owns `-o/--output`, `--core`, and `--mode`; specify those before `--`
where applicable rather than among the run arguments.
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
FASTQ hashes.

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
