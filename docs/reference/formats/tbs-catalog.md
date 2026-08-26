# TBS/WES/TS target BED

Use `bsreadsim run tbs`, `run wes`, or `run ts` with `--targets PATH` and a
plain or gzip-compressed BED6
file. Coordinates are zero-based and half-open.

## BED6 fields

Each data row contains exactly six tab-separated fields:

1. FASTA contig name;
2. interval start;
3. interval end;
4. non-empty target name;
5. finite non-negative score; and
6. strand `+`, `-`, or `.`.

Coordinates must satisfy `0 <= start < end <= contig length`. Empty lines,
`#` comments, and UCSC `track` or `browser` lines are ignored. Duplicate rows
with the same contig, interval, and strand are rejected. Contig names are
case-sensitive and must match the reference exactly.

The score is metadata during uniform sampling. Select `--sampling score` to use it
as a relative output weight; that mode has additional integer requirements
described in [TBS target-score sampling](../advanced/target-score-sampling.md).

## Fragment placement

Fragments are centered on each target, with optional displacement controlled
by `--fragment-center-stddev`. A target is ineligible when a complete requested
fragment cannot be placed on either haplotype or when an emitted mate exceeds
`--max-ambiguous-fraction`.

Uniform mode gives each eligible BED row equal target mass. Overlapping rows
remain separate targets even when they can produce the same fragment
coordinates. Strand controls capture orientation; `.` leaves orientation
unspecified.
