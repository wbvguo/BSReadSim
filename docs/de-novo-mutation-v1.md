# De novo mutation catalog v1

Status: normative deterministic random-variant boundary.

The C++ core owns random haplotype mutation. Python never scans the reference,
constructs haplotypes, or infers mutation flags. This component emits the same
typed `Event` representation used by the verified VCF path.

## Input and widths

The input is one verified, materialized FASTA contig and the normalized
`mutation.rate`, `indel_fraction`, `indel_extension_probability`, and
`homozygous_only` fields. All probabilities are finite binary64 values in
`[0,1]`. Only A/C/G/T positions are mutation candidates; N consumes no RNG
draw and never belongs to a generated event.

Contig coordinates, event intervals, event ordinals, and event counts are
checked `uint32`. The protocol's `0xffffffff` no-event sentinel is never a real
event ordinal. Master seeds, RNG entity ordinals, and aggregate counts are
`uint64`. No u96 coordinate or counter exists.

## Addressed sampling

The key is `derive_key(catalog_seed, Stage::mutation, contig_index)`. The entity is
the zero-based `uint32` reference position widened to `uint64`. Local indices
are frozen as follows:

| Local | Draw |
| ---: | --- |
| 0 | mutation Bernoulli |
| 1 | indel-versus-SNV Bernoulli |
| 2 | one of the other three SNV bases |
| 3 | deletion-versus-insertion Bernoulli |
| 4--6 | three possible indel-extension Bernoulli draws |
| 8--11 | inserted bases in emitted order |
| 16 | homozygous Bernoulli with the nearest binary64 value to one third |
| 17 | heterozygous haplotype Bernoulli |

Draws are stateless, so branch choice, contig chunking, and downstream worker
count cannot shift another position's random stream. If `homozygous_only` is
true the mask is always `both`; otherwise homozygous probability is one third,
and a heterozygous event selects haplotype 1 or 2 with equal probability. These
are the frozen `HaplotypeMask` values 3, 1, and 2, not protocol zero-based
haplotype numbers.

## Event construction

An accepted mutation becomes an SNV with probability `1-indel_fraction` or an
indel otherwise. An SNV chooses uniformly among the three bases unequal to REF.
An indel chooses deletion/insertion with equal probability. Its length starts
at one and is extended by independent Bernoulli draws, capped at the typed
short-indel limit of four bases.

A deletion begins at the candidate base, stops before an N or the contig end,
and advances the scan past its interval. An insertion is anchored immediately
after the candidate base. Because the
typed catalog forbids two events at one anchor, the immediately following
reference position is not a mutation candidate after an insertion. Terminal
insertions at `contig.length` are valid. Events are therefore canonical,
strictly ordered, non-overlapping, and accepted directly by `ContigVariants`.

## Released generator checkpoint

WGBS consumes the generated typed catalog with fixed or variable inserts and
uniform coverage, using the two-bit eligibility boundary in
[variable-haplotype-wgbs-v1.md](variable-haplotype-wgbs-v1.md). Fixed-insert target-GC is calibrated after haplotype construction; variable target-GC with variants remains gated.
RRBS and TBS remain gated until their fragment catalogs become haplotype-aware.
