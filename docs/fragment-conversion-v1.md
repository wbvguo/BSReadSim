# Fragment-level bisulfite conversion v1

Status: normative Python post-processing boundary.

BSReadSim assigns methylation states to a physical DNA fragment, converts that
fragment with bisulfite, and only then derives reads from its declared mate
windows. R1 and R2 are not converted as independent molecules.

## Ordered transform

For one decoded protocol fragment, Python:

1. samples one binary state for every fragment-level methylation site;
2. selects the fragment's C-to-T or G-to-A conversion orientation;
3. visits the complete haplotype-derived template once and attempts conversion
   at every unmethylated target base;
4. slices the converted template into the declared mate windows and reverse
   complements a mate where required;
5. generates quality scores and then quality-conditioned sequencing errors.

A target C/G without a declared site is implicitly unmethylated, matching the
protocol boundary for context filtering.

## Randomness and overlapping mates

The counter address for conversion is:

```text
key             = derive_key(master_seed, CONVERSION, contig_index)
entity_ordinal  = fragment_ordinal
local_index     = template_offset
```

Each physical template base therefore receives at most one conversion draw.
If paired reads overlap, both mates project the same conversion event; a reverse
mate observes its reverse-complemented base and the opposite visible conversion
mode. Worker count, chunk size, mate order, and overlap length cannot duplicate
or move the draw.

Per-base truth records both the original oriented base and the oriented
post-conversion base. Conversion flags are physical fragment events projected
into each observing mate, while sequencing-error flags remain independent read
observations.
