# Output modes

`bsreadsim run` exposes exactly two modes:

- `--mode production` is the default. It emits FASTQ plus the reproducibility
  manifest and, unless truth BAM is enabled, requests no Truth-only protocol
  columns.
- `--mode debug` emits the same FASTQ and manifest plus canonical Truth JSONL
  1.2. It requests the sparse columns needed to reconstruct that artifact.

Truth is a run mode, not a user configuration field. `output.truth` is rejected
by the configuration schema. The selected internal policy is inserted during
normalization and is therefore recorded in the config SHA-256 and manifest.

`--truth-bam` is an independent output policy, not a third mode. It requests
the Full Truth coordinate projection, adds an unsorted SAM/BAM 1.6 artifact,
and leaves JSON Truth controlled by production/debug. Thus production with
truth BAM has Full Truth transport but no Truth JSONL. See
[truth-bam-v1.md](truth-bam-v1.md).

For the same biological inputs, models, seed, and software version, production
and debug must emit byte-identical decompressed FASTQ. Their protocol streams,
config identities, manifests, and artifact sets intentionally differ.

## FASTQ read identifiers

Read identifiers are assigned in the output formatter from fragment metadata
already emitted by HTSIM. They follow
`@<chr_id>:<c1>-<c4>:<read_num>/<pair_num>`; coordinates are one-based and
inclusive, `read_num` is the zero-based fragment ordinal, and `pair_num` is 1
or 2. The output filename prefix is not part of the identifier. See the
[read-name v1 contract](read-name-v1.md) for the exact coordinate and BAM
QNAME rules.

## FASTQ comment line

The third FASTQ line is exactly `+`. It carries no Truth data.

An earlier design placed a per-base `M/X/I/E` diagnostic string after `+`.
That output-symbol representation is not a lossless encoding of Truth JSONL:

- one symbol collapses methylation context/state, conversion outcome, variant
  class, insertion status, and sequencing-error precedence;
- it omits coordinates, event and site IDs, probabilities, source/allele
  identity, quality, and intermediate base states;
- constructing it would restore per-base annotation work to production and add
  roughly one output character per read base.

It is therefore neither compressed Truth nor a source from which Full Truth
can be restored. Any future diagnostic symbol output must be a separately
versioned debug feature and cannot change production FASTQ.

## Required gates

Any output-path change must prove all of the following:

1. production without truth BAM creates only R1, optional R2, and manifest roles;
2. debug additionally creates exactly one canonical Truth role;
3. truth BAM opt-in adds exactly one `truth_bam` role counted per mate;
4. protocol Truth-column policy matches the selected mode and BAM policy;
5. fixed-input, fixed-seed production/debug runs have byte-identical FASTQ;
6. core and Python fragment, mate, site, base, and per-contig counts reconcile;
7. output files are staged privately and the manifest is published last;
8. throughput is reported separately for production, debug, and truth BAM after the
   correctness gates pass.
