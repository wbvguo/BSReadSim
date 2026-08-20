# BSReadSim architecture

Status: current architecture for BSReadSim 0.3.

BSReadSim has one runtime path. The C++17 HTSIM process generates biological
fragment templates and writes the columnar binary protocol. BSReadSim's Python
layer validates that stream, applies bisulfite conversion and sequencing
effects, and transactionally publishes FASTQ plus a manifest. Full Truth is an
explicit extension of the same path: debug publishes JSONL, while the
independent truth-BAM option uses its coordinate projection for alignments.

```mermaid
flowchart LR
    C["CLI arguments"] --> P["Python: normalize, validate, resolve, hash"]
    P --> G["HTSIM C++: genome and fragment generation"]
    G --> W["Columnar protocol 2.0"]
    W --> S["BSReadSim Python: conversion, quality, errors"]
    S --> O["FASTQ + manifest"]
    S -. "debug only" .-> T["Truth JSONL"]
    S -. "--truth-bam" .-> B["HTSlib: unsorted truth BAM"]
```

Protocol 1 and the pre-0.2 native programs are available from the `v0.1.0`
tag; they are not compiled, imported, selected, or documented as current
runtime alternatives.

## Ownership

HTSIM owns:

- reference, VCF, CGmap/bedMethyl, ASM/ASM BED, technology BED, and coverage-profile parsing;
- haplotype construction and coordinate projection;
- methylation-site catalogs and probability assignment;
- WGBS, RRBS, and TBS fragment selection;
- deterministic fragment allocation and the protocol writer.

Python owns:

- public CLI and run-mode selection;
- configuration validation, defaults, path resolution, and immutable-input
  hashing (the RRBS scored-candidate exchange uses regenerated-row matching);
- core lifecycle and strict protocol decoding;
- fragment-level bisulfite conversion before mate derivation, quality
  generation, and sequencing errors;
- bounded process scheduling and ordered output publication;
- FASTQ, optional Truth JSONL/truth BAM, and the reproducibility manifest.

Scientific generation is not duplicated across languages. Python does not
reparse biological inputs, and C++ does not choose output paths, compression,
or Python worker scheduling.

## Source layout

The native tree is grouped by responsibility:

| Directory | Responsibility |
| --- | --- |
| `HTSIM/src/core` | argv validation, application assembly, generation loop |
| `HTSIM/src/reference` | FASTA loading and immutable contig storage |
| `HTSIM/src/variant` | variants, haplotypes, and coordinate projection |
| `HTSIM/src/methdb` | methylation sites, CGmap/bedMethyl and ASM/ASM BED inputs, and diploid MethDB catalogs |
| `HTSIM/src/model` | shared biological value types |
| `HTSIM/src/fragment` | fragment counts, allocation, insert lengths, and construction |
| `HTSIM/src/wgbs` | reference and haplotype-aware WGBS fragment selection |
| `HTSIM/src/rrbs` | restriction-site discovery, candidate catalog, and sampling |
| `HTSIM/src/tbs` | target BED projection, candidate catalog, and sampling |
| `HTSIM/src/protocol` | model-to-column adapter, ordered batch emission, and wire encoding |
| `HTSIM/src/shared` | RNG, distributions, SHA-256, and immutable text snapshots |

The Python package uses one module per runtime boundary:

| Module | Responsibility |
| --- | --- |
| `config` | schema validation and canonical configuration identity |
| `catalog` | launch the native RRBS candidate exporter; no motif parsing |
| `preparation` | seed materialization and immutable input snapshots |
| `core_argv` / `core_process` | C++ launch projection and process supervision |
| `protocol` / `protocol_adapter` | wire decoding and column-to-model projection |
| `postprocess` / `numpy_postprocess` | conversion, quality, and error stages |
| `process_pool` / `pipeline` | bounded scheduling and stage orchestration |
| `bam` | SAM 1.6 truth-alignment formatting from typed base projections |
| `output` / `manifest` | HTSlib-backed BAM streaming, transactional publication, and audit record |

Small compiled Python extensions accelerate protocol validation and common
post-processing loops. Their Python implementations remain the semantic
reference and are exercised by differential tests.

## Run sequence

1. The sole `bsreadsim run` command projects CLI arguments into an in-memory
   run document. The normalizer rejects unknown or invalid values, applies
   schema defaults, checks cross-field scientific rules, resolves paths, and
   freezes canonical JSON plus its SHA-256 identity for the manifest. No JSON
   run-configuration file is accepted. `--mode production` adds the internal
   no-Truth policy, while `--mode debug` adds Full Truth. Users cannot
   configure `output.truth` directly. `--truth-bam` independently requests
   Full Truth transport while leaving JSON publication controlled by the mode.
2. `prepare_run` materializes an omitted 64-bit seed once and hashes every
   reference, biological input, and model artifact. Declared model hashes are
   checked before any process starts or output directory is staged.
3. `core_argv` projects only C++-owned fields and verified paths/digests. The
   core independently validates argv syntax and the bytes it opens. Parsed
   argv and direct generator calls share one semantic validator for numeric
   ranges and technology combinations.
4. The core builds contig-level biological state, assigns stable fragment
   ordinals, and emits one header followed by ordered columnar batches and one
   trailer. Details are normative in [protocol.md](protocol.md).
5. Python validates preamble, frame length, CRC32C, column bounds, header
   identity, batch order, and trailer counts before consuming data.
6. With `workers=1`, batches use one reusable process-local slot. With more
   workers, the supervisor remains the sole stdout reader and dispatches
   authenticated batches through bounded shared-memory slots. Results are
   written strictly by fragment ordinal.
7. Output files are created in a private staging directory. Counts and hashes
   are computed there; the manifest is published last as the completion
   marker. Existing destinations are never overwritten.

Each validation above protects a different boundary. Configuration is not
fully revalidated in preparation or argv projection; the immutable canonical
identity prevents mutation between those stages. Input bytes, untrusted core
bytes, and final output accounting are still checked where they enter their
respective trust boundary.

## Output modes

Production is the default and, without truth BAM, requests no Truth-only
protocol columns. It creates R1, optional R2, and a manifest. Debug requests
sparse Full-Truth columns and additionally creates canonical Truth JSONL.
Both modes use the
same biological generation and RNG addresses, and must produce byte-identical
FASTQ for the same inputs, model, seed, and software version. See
[output-modes.md](output-modes.md).

Truth BAM is orthogonal to those modes. When selected, workers derive
reference-forward SAM records from the typed coordinate projection and stream
them through the same statically linked HTSlib in `htsim-core`; no temporary
SAM is published. BAM finalization is part of the same manifest-last output
transaction. Its normative contract is [truth-bam-v1.md](truth-bam-v1.md).

The third FASTQ line is exactly `+`. Historical output-symbol annotations were
lossy and would reintroduce per-base work, so they are not a compressed Truth
format and are not part of production output.

## RNG identity contract

The released RNG contract is `philox4x32-10+philox-domain-v2`. Philox remains
the stateless counter generator used for simulation draws; BLAKE2b is not part
of the RNG path. Domain keys are derived by one reserved Philox block:

```text
domain_key_entity = 0x4253522f4b455932  # numeric ASCII "BSR/KEY2"
domain_local      = (uint64(stage) << 32) | uint64(contig_index)
domain_block      = Philox4x32-10(master_seed,
                                  domain_key_entity,
                                  domain_local)
derived_key       = uint64(domain_block[0])
                    | (uint64(domain_block[1]) << 32)
```

`stage` is a frozen `uint32` enum:

| Value | Stage |
| ---: | --- |
| 0 | mutation |
| 1 | methylation-level |
| 2 | fragment |
| 3 | haplotype |
| 4 | site-state |
| 5 | library-orientation |
| 6 | conversion |
| 7 | quality |
| 8 | sequencing-error |

`contig_index` is the zero-based `uint32` ordinal in the verified reference
catalog and protocol header. Therefore changing reference contig order changes
the RNG stream even when names and sequences are otherwise identical. Names,
UTF-8 encoding, and string normalization never enter an RNG address.

HTSlib owns verified genomics parsing, decoded input access, and final BAM
serialization. It neither
defines random streams nor replaces the counter-based RNG: deterministic
worker/chunk invariance still requires explicit Philox addresses.

The ownership boundary is intentionally narrow. HTSlib 1.24 is fixed by the
top-level `HTSLIB/` gitlink and handles format detection, gzip/BGZF decoding,
VCF headers, records, alleles, and genotypes. BSReadSim still verifies raw
SHA-256 bytes through one stable descriptor and enforces its stricter FASTA,
single-sample VCF, BED6, coordinate, ordering, and variant-subset contracts.
Autotools runs only on a build-tree copy, so configuring never dirties the
submodule; the resulting core links HTSlib statically.

## Scientific invariants

- Reference intervals are zero-based, half-open `[begin, end)`; ordinal and
  offset arrays are monotone and bounds-checked.
- Complete diploid haplotypes and contig-level methylation catalogs are built
  before provider-specific fragment discovery.
- Methylation state and bisulfite conversion are physical fragment events.
  Conversion is attempted once per template offset before mate windows are
  derived, so overlapping mates observe the same converted molecule.
- Philox addresses stochastic decisions by numeric stage, contig index, fragment, mate,
  site, and cycle as applicable. Worker count, chunk size, and completion order
  cannot change generated reads.
- Truth omission, when neither JSON nor BAM requires the projection, removes
  transport, object construction, serialization, and I/O only. It cannot
  change fragment allocation, methylation probabilities, conversion, quality,
  or error draws.
- The manifest records normalized configuration, seed, input and output
  digests, component versions, stream digest, and reconciled core/Python
  counts.

Component-level algorithms remain versioned independently when their exact
definition affects reproducibility, for example
[`coverage-profile-target-v1`](coverage-profile-target-v1.md),
[`vcf-catalog-v1`](vcf-catalog-v1.md), and
[`sequencing-models-v1`](sequencing-models-v1.md). The `v1` in those names is
the component contract version, not the retired wire protocol.

## Failure and dependency policy

The core writes diagnostics only to stderr. Malformed CLI arguments, saved
config values, or core argv fail before protocol output; a malformed stream,
failed child, count mismatch, or output collision aborts the staging
transaction. Partial files and a complete manifest are never published
together.

Runtime dependencies are limited to NumPy and JSON Schema validation. JSON
Schema is an internal normalized contract, not a required user-authored file.
It keeps type/default/range rules declarative and shared by direct CLI and
saved-config runs. Optional scientific behavior must be implemented as a
reviewed, versioned component rather than a dynamically loaded plugin.
