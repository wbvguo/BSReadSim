# Generation-technology extension boundary

Status: current contract for adding a built-in sequencing technology.

WGBS, RRBS, and TBS are explicit C++ providers. A provider owns only its
candidate catalog, eligibility rules, allocation weights, and fragment sampling
policy. Shared C++ code owns reference/haplotype/MethDB construction, fragment
materialization, ordinals, protocol encoding, and accounting. Python receives
the same columnar stream for every provider and contains no second
technology-specific simulator.

## Requirements

A new technology must be delivered as a reviewed built-in provider with:

1. one mutually exclusive configuration section and normalized defaults;
2. preparation rules for each immutable input/model artifact and its SHA-256;
3. one deterministic Python-to-core argv projection;
4. a focused C++ catalog/planner/sampler component using checked `uint32_t` and
   `uint64_t` domains;
5. pre-preamble capability checks, stable RNG addresses, bounded rejection,
   and exact fragment/skipped counts;
6. component tests, argv tests, cross-language protocol tests, and fixed-seed
   end-to-end byte tests across worker and chunk settings; and
7. an explicit protocol-version decision if the existing columns cannot carry
   the required scientific values.

A provider cannot own FASTQ writing, quality/error models, methylation-state
draws, Python concurrency, or output publication. Dynamic generator and
methylation plugins are intentionally absent: adding scientific behavior
requires normal source review, versioning, and the complete validation gates.
