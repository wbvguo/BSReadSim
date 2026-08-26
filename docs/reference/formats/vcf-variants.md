# VCF variant sets

Use `--vcf PATH` to simulate predefined diploid variants. This automatically
disables de novo mutation generation.
Input may be plain-text or gzip-compressed VCF 4.2 or 4.3.

## Accepted subset

BSReadSim currently requires:

- exactly one named sample and ten tab-separated columns;
- rows ordered by FASTA contig order and nondecreasing one-based `POS`;
- one ALT allele and a diploid `GT` containing only alleles `0` and `1`;
- uppercase A, C, G, or T in the normalized event (a shared indel anchor may
  be `N`); and
- an SNV or a pure insertion/deletion of at most four bases after
  normalization.

Missing genotypes, multiallelic records, symbolic alleles, MNPs, complex
replacements, duplicate normalized events, and overlapping normalized events
are rejected. Two source rows may share `POS` when their normalized events are
distinct and non-overlapping—for example an SNV followed by a left-anchored
deletion. A `0/0` row is allowed but creates no variant. FILTER and INFO do not
silently remove records.

REF bases and coordinates are checked against the supplied FASTA. Contig names
must match exactly.

## Phasing and reproducibility

Phased `0|1` and `1|0` assignments are preserved; `1|1` applies ALT to both
haplotypes. Unphased heterozygous variants are assigned deterministically from
`--seed-phase`. Keep that seed fixed when regenerating a MethDB or an RRBS
candidate set from the same VCF.

VCF variants are supported by WGBS, RRBS, TBS, WGS, WES, and TS. Candidate fragments that
cannot be represented cleanly at a deletion or inserted boundary are omitted
rather than truncated. Annotated BAM records retain variant effects in their
CIGAR, sequence, and BSReadSim tags.

## Export generated variants

Save the exact de novo variant set as a deterministic, BGZF-compressed phased
VCF:

```bash
bsreadsim build variants \
  --reference reference.fa \
  --output variants.vcf.gz \
  --mutation-rate 0.001 \
  --seed-mut 7
```

The command accepts paths directly; BSReadSim computes input identities
internally. Use the same mutation seed and mutation parameters as a later run.
Alternatively, add the `--save-vcf` flag to a de novo run; BSReadSim writes
`OUTPUT/truth/PREFIX.variants.vcf.gz`. `build variants --vcf INPUT` can also
normalize and deterministically phase an existing VCF.
