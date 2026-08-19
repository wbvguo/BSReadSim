# Experiments

Place scientific experiments here, one generated run directory per invocation.
Small top-level drivers and intentionally tiny `*.fa`, `*.bed`, and
`*.vcf` mock fixtures may be tracked; nested inputs, model artifacts, and
results are ignored.

An experiment is not a regression test. Promote a minimal deterministic case
to `tests/` only when it becomes a correctness contract.

## Runnable RRBS/TBS mocks

The checked-in mock is intentionally small but exercises real biological input
boundaries:

- `mock-reference.fa` contains seven MspI `C|CGG` cut sites on `chrMock` and a
  control contig with no RRBS cuts or TBS targets;
- the RRBS CLI example samples enzyme-bounded fragments of length 8, 16, or 24;
- `mock-targets.bed` contains forward, reverse, and unknown-strand BED6 targets;
- the TBS CLI example uses exact BED score weights `3:1:0`, so the zero-weight
  control target cannot be selected;
- `mock-variants.vcf` contains phased insertions used to verify that haplotypes
  are constructed before motif discovery or target-centered fragmentation.

From the repository root after installing the package, the same mock components
can be run entirely with command-line arguments:

```sh
bsreadsim run \
  -r data/experiments/mock-reference.fa \
  -o data/experiments/runs/rrbs-cli \
  -n 24 --seed 20260813 --mutation-rate 0 \
  --technology RRBS --cut-site 'C|CGG' \
  --read-length 4 \
  --insert-min 8 --insert-mean 16 --insert-max 24 --insert-stddev 0 \
  --max-ambiguous-fraction 0 --error-rate 0 --compression none

bsreadsim run \
  -r data/experiments/mock-reference.fa \
  -o data/experiments/runs/tbs-cli \
  -n 24 --seed 20260813 --mutation-rate 0 \
  --technology TBS --targets data/experiments/mock-targets.bed \
  --target-score --fragment-center-stddev 0 \
  --read-length 4 --insert-size 12 \
  --max-ambiguous-fraction 0 --error-rate 0 --compression none
```

Results are written below `data/experiments/runs/` and remain ignored. Both
commands use independent Bernoulli state, the normalized default. Add
`--vcf data/experiments/mock-variants.vcf` to either command to exercise the
phased-variant path; the CLI automatically selects de novo mutation rate zero
when `--vcf` is present and no explicit mutation rate is given.

These commands are the complete runnable interface; no parallel JSON
run-configuration fixtures are maintained.

## WGBS target-GC mock

`wgbs-gc-target-mock.tsv` is an illustrative chr21-shaped 100-bin target
distribution. It is not fitted from ERR2359938 or another empirical dataset
and must not be presented as the paper profile. Each physical line is one
output probability; the zero-based line number is the bin, and all 100 values
sum to one. The profile has zero mass on GC bins that the PE150/insert-300
chr21 domain cannot reach.

The target was constructed from the old asymmetric acceptance mock and the
eligible hg38 chr21 start histogram: `p_i` is proportional to `N_i*a_i`.
This preserves that mock's intended chr21 bias while changing the public
contract to the directly interpretable requested output distribution. At run
time the core independently calibrates `p_i/N_i`; it does not trust or reuse
the construction-time acceptance values.

The CLI calculates and records the profile SHA-256. The normalizer,
preparation stage, and C++ reader independently carry or verify that identity.
Because a target assigns exact mass to bins, it can only be used with a
reference and insert shape that can reach every positive bin.

## hg38 chr21 30x target-GC audit

The retained audit uses the UCSC hg38 primary chromosome
[`chr21.fa.gz`](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz).
Its official MD5 is `184df2bd9b812b6e6b6da16c6021369e`; the downloaded
bytes matched it and had SHA-256
`c979ca1e5065c2521a50773473e0d0cc018fd6f3e9bb3aa90493fe7b45d57d1b`.
The run fixes PE150, insert length 300, depth 30, seed 20260815, no de novo
mutation, and one C++ plus one Python worker:

```sh
bsreadsim run \
  -r /path/to/chr21.fa.gz \
  -o /path/to/chr21-30x-fastq \
  -d 30 --seed 20260815 --mutation-rate 0 \
  --read-length 150 --insert-size 300 \
  --coverage-profile data/experiments/wgbs-gc-target-mock.tsv \
  --error-rate 0 --workers 1 --core-workers 1 \
  --gzip-level 1 --prefix chr21_30x

python3 data/experiments/audit_wgbs_gc_target.py \
  --manifest /path/to/chr21-30x-fastq/chr21_30x.manifest.json \
  --output-tsv /path/to/chr21-30x-gc-bins.tsv \
  --output-json /path/to/chr21-30x-gc-summary.json
```

The audit does not calculate GC from bisulfite-converted FASTQ bases. It
reconstructs the completed manifest, verifies the current reference/profile
bytes, and replays the same run ID through protocol 2.0. The replayed stream
SHA-256 matched the FASTQ run exactly, so its fragment coordinates are the
coordinates that produced the FASTQ files.

For each bin, the retained `chr21-30x-gc-target-v1.tsv` records the requested
target, eligible starts, calibrated acceptance, and observed fragment
probability. `chr21-30x-gc-target-v1.json` records exact run identity,
accounting, distribution metrics, and replay timing. The files are generated
only after exact protocol replay and the full target-overlap checks pass.

The retained result emitted 4,670,998 fragments (9,341,996 reads). Target
versus observed overlap was `0.999009341`, total variation was `0.000990659`,
and the largest absolute per-bin difference was `0.000120237` (0.0120
percentage points). Expected versus observed proposal acceptance was
`0.552658070` versus `0.552628251`; the protocol stream SHA-256 matched the
FASTQ production manifest exactly.

On the recorded host, that gzip-level-1 production run took 116.27 s wall time
with 257,372 KiB peak RSS. R1 and R2 were 261,306,755 and 261,335,438 bytes.
These host measurements are descriptive, not a portable performance claim;
the paired uniform/profile benchmark uses a separately fixed execution lane.
