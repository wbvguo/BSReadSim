# Haplotype interval projection v1

Status: normative contract for applying the typed VCF catalog.

The projector is a pure C++ boundary between the verified per-contig VCF event
catalog and later fragment/methylation construction. It consumes original
reference bases, one zero-based selected haplotype, and a non-empty 0-based
half-open reference interval. It emits the haplotype template bases, one
original-reference position per base, one event id per base, and the typed
variant events needed by the fragment materializer. It performs no sampling
and no I/O.
The returned value also carries its `uint32` contig index and zero-based
haplotype, so later builders cannot accidentally pair otherwise valid arrays
with another contig or haplotype.
The validated catalog retains its contig index and reference length, so passing
it beside another materialized contig fails before projection. Event lookup is
`O(log V + K)` for `V` events on the contig and `K` events near the interval;
it does not rescan a whole human VCF for every fragment.

## Width and bit boundary

Reference coordinates and stable per-contig VCF event ordinals are `uint32_t`.
No u96 representation is used.

The catalog's `HaplotypeMask` uses bit value 1 for logical haplotype 0, value 2
for logical haplotype 1, and value 3 for both. The projector converts those
availability bits to the zero-based wire haplotype at this boundary.

## Event application

- An event is applied only when its two-bit mask contains the selected
  haplotype. Other events leave the reference unchanged and are not emitted.
- SNVs replace their one mapped reference base. That template base carries the
  stable per-contig event ordinal in `base_variant_indices`.
- Deletions remove their mapped reference bases. The deletion remains in
  `variants` even though no template base can point to it.
- Insertions are emitted immediately before `reference[anchor]`; their
  `reference_positions` are `-1`, and each inserted base points to the event.
- Event ids are the original zero-based ordinals in the validated per-contig
  catalog. They may therefore be sparse in one fragment, but remain stable
  across overlapping fragments and never use `0xffffffff`.
- A heterozygous event records the zero-based ALT haplotype in
  `phased_haplotype`. A homozygous ALT event uses `255` because phase is not
  applicable.

For an internal insertion anchor, the interval beginning at the anchor owns the
insertion; the interval ending there does not. An insertion at coordinate zero
belongs to an interval beginning at zero. A terminal insertion at contig length
belongs to an interval ending at contig length so projecting the complete
contig does not lose it.

An interval may not cut through an active SNV/deletion. A fully contained
deletion is legal, but projection fails if the selected events remove every
base because the protocol forbids an empty fragment template. Later sampling
must treat such intervals, and intervals whose projected template is shorter
than a requested read, as ineligible rather than truncating an event.
