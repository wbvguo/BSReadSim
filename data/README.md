# Data workspace

All non-test experiments and performance benchmarks belong under this
directory. Correctness fixtures and regression checks remain under `tests/`;
they are part of the package gate, not research output.

- `experiments/` contains small, reviewable experiment drivers and configs.
- `benchmarks/` contains small, reviewable benchmark drivers and configs.
- raw references, downloaded inputs, trained model artifacts, run directories,
  FASTQ files, truth files, reports, and timing output stay untracked beneath
  those directories. Intentionally tiny, reviewed mock fixtures may be tracked
  beside their top-level experiment configs.

Every retained run should record the command, normalized configuration, input
and model SHA-256 digests, package/core versions, platform, seed, and output
manifest. Large inputs or results must use external artifact storage rather
than Git.
