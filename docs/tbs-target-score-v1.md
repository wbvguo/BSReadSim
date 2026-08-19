# TBS target output-weight sampling v1

Status: normative optional C++ fragment-sampling contract.

The paper models target-specific capture efficiency by sampling TBS regions in
proportion to provided or empirically estimated output. In this contract the
BED6 score is the target-level `output_weight`: an expected relative share of
the observed target output, potentially estimated from aligned BAM read counts.
It is not named or interpreted as a pure molecular capture probability because
aligned depth can also contain library, PCR, bisulfite, and mapping effects.

The user-facing CLI selects this mode with `bsreadsim run --target-score`; no
JSON input or second coverage artifact is required. The normalized run identity
records `coverage.kind = "target-score"`. The mode is valid only with
`technology = "TBS"` because the verified BED snapshot owns both target
coordinates and output weights.

## Exact weight domain

Every BED score must be numerically integral and in `[0, UINT32_MAX]` when
output-weight mode is selected. Decimal spellings such as `1.0` and `1e3` are
accepted only when their parsed value is exactly an integer in that range.
Raw per-target BAM read counts satisfy this representation directly. Ratios
computed from those counts have the same sampling distribution when represented
by the original counts or another common exact integer scale. Multiplying every
positive weight by the same constant does not change the intended distribution.
A zero weight is valid metadata but gives that eligible target zero sampling
probability.

Target syntax and reference validation still follow
[tbs-catalog-v1.md](tbs-catalog-v1.md). All output weights are validated even
when their centered fragment is ineligible. The weight sum of eligible targets
on each contig must fit `uint32_t`. At least one eligible target in the complete
reference must have positive weight; otherwise planning fails before the
protocol header.

## Allocation and selection

For fixed-center sampling, a target contributes its output weight once when at
least one haplotype passes the coordinate and mate-ambiguity checks. For
positive center standard deviation, it contributes once when at least one
haplotype retains the target center and has a sequenceable start, matching the
displacement eligibility rule.

The output weight is aggregate target mass, not per-haplotype mass. It is never
multiplied by `popcount(eligible_haplotypes)`. After a target is selected, a
one-bit mask chooses its only available haplotype with probability one; mask
`both` chooses haplotype 0 or 1 with probability one half. No separate
haplotype weight is inferred from an unphased aggregate BAM.

The per-contig allocation weight is the exact sum of those eligible target
weights. C++ applies the existing Hamilton largest-remainder allocator to
these `uint32` contig weights, so the requested global read-pair count is
preserved exactly. Within a contig, eligible targets remain in canonical BED
order and form an inclusive `uint32` cumulative-weight array. For each output
candidate:

```text
key              = derive_key(seed, Stage::fragment, contig_index)
entity_ordinal   = the per-contig output candidate ordinal (uint64)
local_index      = 1
draw             = bounded_integer(key, entity_ordinal,
                                   local_index, total_weight)
selected target  = first cumulative weight strictly greater than draw
```

Haplotype selection uses `Stage::haplotype`, the contig-index-derived key, the
global fragment ordinal, and local index 0. With positive center standard
deviation, the target and haplotype are each selected once and all
normal-displacement retries remain attached to that pair. Retry addresses, the
attempt cap, and `skipped_fragment_count` remain exactly as defined in
[tbs-catalog-v1.md](tbs-catalog-v1.md). Chunk size therefore cannot change
target selection, haplotype selection, or fragment bytes.

## Width and bit boundary

BED coordinates, target ordinals, individual weights, per-contig cumulative
weights, and per-contig requested counts are checked `uint32_t`. Global
candidate/fragment ordinals, aggregate counts, and Philox counter halves are
`uint64_t`. Cross-contig allocation multiplies the checked `uint32` fragment
count by one checked `uint32` contig weight, and that product fits `uint64_t`.
No u96 coordinate, counter, or allocation value is used.

The component consumes typed `CaptureStrand` and `HaplotypeMask` values. A zero
eligibility bitset is represented by omission rather than by constructing an
invalid mask.
