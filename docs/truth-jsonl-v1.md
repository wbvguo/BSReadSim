# Truth JSONL 1.2

Status: normative optional BSReadSim 0.2 debug-output boundary.

When the effective `output.truth=full`, the Python writer emits one UTF-8 JSON
object per successfully processed fragment. The public command line selects
this policy only with `--mode debug`; its default production mode selects
`output.truth=none` and emits no Truth JSONL artifact. Emitted objects use sorted
keys, compact separators, and no non-finite JSON numbers. When gzip is selected,
the gzip header has zero modification time and no host filename. Fragment order
is the protocol fragment-ordinal order.

Version 1.1 added the required top-level `variant_events` array. Version 1.2
additionally preserves the complete protocol identity of every methylation site
alongside its sampled state. The manifest identifies the current format as
`bsreadsim-truth-jsonl` version `1.2`.

## Fragment fields

- `fragment_ordinal`, `contig`, and `haplotype` identify the C++ fragment.
- `fragment_conversion_mode` records the Python library-orientation decision.
- `site_states` records one sampled state per protocol methylation site. Every
  entry carries `site_index`, `template_offset`, `reference_pos`, `context`,
  `source`, `allele`, `probability`, and `methylated`; the typed context/source/
  allele enum names are preserved rather than inferred again in Python.
- `mates` records each final sequence and its complete per-base annotations.
- `variant_events` is the protocol event tuple copied through post-processing
  without filtering or reordering.

Each `variant_events` entry contains `event_id`, `kind`, `phased_haplotype`,
`reference_start`, `reference_end`, `ref_bases`, and `alt_bases`. `kind` is one
of `SNV`, `INSERTION`, or `DELETION`; bases are uppercase A/C/G/T/N strings.
Coordinates remain 0-based half-open. Insertions have an empty reference
interval and empty `ref_bases`; deletions have empty `alt_bases`. Event IDs are
fragment-local protocol IDs and may be referenced by mate annotations.

Every mate annotation records `read_offset`, `reference_pos`,
`variant_event_id`, `site_index`, the oriented/post-conversion/final base
encodings, methylation/conversion/error outcomes, and Phred score. The
`0xffffffff` variant ID remains the no-event sentinel. An inserted base uses
`reference_pos = -1`; a deleted base is absent and is represented only by its
fragment-level event.
