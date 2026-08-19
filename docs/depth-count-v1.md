# WGBS depth-to-count contract v1

Status: normative C++ fragment-count boundary for WGBS.

The C++ generator, not Python, converts requested sequencing depth into the
number of fragments (read pairs in paired-end mode). Version 1 enables this
conversion for WGBS only. RRBS and TBS remain explicit-`read_pairs` modes until
their distinct effective-target-length denominators are frozen.

## Denominator

After complete pre-header eligibility planning, a contig contributes its full
verified reference length if it has at least one eligible fragment start on
either haplotype; a contig with zero eligible starts contributes zero. These
lengths are summed in FASTA order as checked `uint64`. Each contributing contig
length remains checked `uint32`.

For read length `R`, paired-end mode emits `2R` sequenced bases per fragment;
single-end mode emits `R`. The current paired-end generator requires equal
R1/R2 lengths.
The derived fragment count is:

```text
floor((binary64(effective_reference_bases) * binary64(depth))
      / binary64(R * mate_count))
```

The operation order above is normative and requires round-to-nearest floating
point. A zero result, non-finite intermediate, or result above `UINT32_MAX`
fails before the protocol preamble. The final count and every per-contig count
are `uint32`; the effective-length sum is `uint64`. No u96/u128 integer or
hidden Python count calculation is used.

This definition targets mean sequenced-base depth over reference contigs that
can actually emit fragments. Variant indels do not change the denominator:
FASTQ mate lengths remain fixed, while their typed reference-coordinate truth
is carried separately by the protocol.
