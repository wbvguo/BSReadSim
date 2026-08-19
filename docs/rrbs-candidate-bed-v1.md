# RRBS candidate BED exchange v1

Status: normative optional C++/external-model exchange contract.

RRBS candidate discovery remains owned exclusively by C++. The exchange path
allows an external GAM or other predictor to attach one relative score to each
known fragment without reimplementing restriction-site or haplotype logic:

1. `bsreadsim catalog rrbs -r REF -o candidates.bed --cut-site 'C|CGG'`
   asks the native core to generate the canonical rows directly from CLI
   arguments;
2. an external program preserves every field except `score` and writes a
   scored BED; and
3. `bsreadsim run ... --technology RRBS --rrbs-candidates scored.bed
   --rrbs-score` consumes the scores. Omitting `--rrbs-score` validates the
   same rows but ignores scores completely.

No JSON configuration file is read or written by this workflow. Reference,
cut-site, read-layout, insert-range, ambiguity, VCF, mutation, and seed flags
that define the catalog must have the same effective values on the export and
run commands. Shared CLI defaults count as the same values.

BED is optional. Uniform RRBS without `rrbs.candidate_bed` retains the direct
in-memory path.

## No hash contract

The candidate BED contains no fingerprint, file digest, catalog digest, or
digest metadata. There is no `--rrbs-candidate-bed-sha256` core option, and
Python does not add this exchange file to the prepared-input hash list.

Instead, the sampling core regenerates the complete candidate catalog from the
current reference, variants, cut declarations, and fragment parameters. It
then requires an exact one-to-one match by candidate ID and all serialized
fixed fields. Missing, extra, duplicate, or changed candidate rows fail before
the protocol header. Input row order is irrelevant. This validation does not
turn the external BED into an alternative source of biological candidates;
the regenerated C++ candidate is always the object used to build a fragment.

Because unphased VCF records and de novo variants depend on the master seed, a
candidate BED used with either variant source requires the same explicit
`--seed` on both CLI commands. Reference-only geometry is seed independent.

## Ten-column format

The exporter writes a comment header followed by plain, tab-separated rows:

```text
#chrom start end candidate_id score strand haplotype_mask template_length gc_count restriction_site_count
```

The separators in the real file are tabs. Columns are:

1. `chrom`: exact reference contig name;
2. `start`: zero-based reference-envelope start;
3. `end`: zero-based half-open reference-envelope end; equality with `start`
   is permitted for a physical insertion-only envelope;
4. `candidate_id`: short deterministic row identity;
5. `score`: `.` or a finite non-negative decimal relative weight;
6. `strand`: `.`, because RRBS candidates have no capture strand;
7. `haplotype_mask`: `1`, `2`, or `3` using the shared two-bit mask;
8. `template_length`: physical haplotype fragment length;
9. `gc_count`: physical count of C and G bases in that template; and
10. `restriction_site_count`: endpoint plus internal recognized cut sites.

The native exporter initializes every score to `1`. `.` is accepted only when
coverage is uniform. Profile mode requires a numeric score for every row and
requires positive total mass somewhere in the reference. Scores are relative
and need not sum to one; multiplying all positive scores by one constant does
not change the intended distribution.

The v1 reader accepts plain text with LF or CRLF, blank lines, and comment lines
beginning with `#`. Data rows contain exactly ten fields and are limited to
1 MiB. Compressed BED is outside this first exchange contract.

## Candidate IDs

For an envelope represented by one row, the ID is:

```text
contig:start-end
```

If multiple physical candidates share that reference envelope, every member
receives a local base-36 suffix in canonical C++ catalog order:

```text
contig:start-end~0
contig:start-end~1
```

The suffix is only a compact disambiguator. It is not a global candidate
ordinal, is not used as an RNG address, and does not encode haplotype meaning.
`haplotype_mask` remains the authoritative biological field.

## Fair allocation and sampling

The fairness unit is one physical haplotype fragment copy. For candidate `i`,

```text
copy_mass_i = popcount(haplotype_mask_i)
```

so masks `1` and `2` contribute one copy and mask `3` contributes two. Uniform
RRBS uses `copy_mass_i`; profile RRBS uses
`score_i * copy_mass_i`. Per-contig masses drive the shared Hamilton
largest-remainder allocation, preserving the exact requested global fragment
count. A selected mask-3 row makes the later haplotype choice with probability
one half for each copy.

Within a profiled contig, candidates form a cumulative non-negative
floating-weight array in canonical catalog order. One addressed
`uniform01(seed, "fragment", contig, candidate_ordinal, 1)` draw selects the
first cumulative value strictly greater than the scaled draw. Candidate
ordinals continue across chunks, so profile and uniform selection remain chunk
independent and report zero skipped fragments.
