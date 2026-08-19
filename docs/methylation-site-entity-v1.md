# Methylation site entity v1

Status: normative RNG identity contract for variant-aware methylation levels.

The `methylation-level` Philox counter already carries a 64-bit entity ordinal.
Variant-aware sites use that existing field; no u96 value and no additional
counter word are needed. The encoding is a typed C++ boundary.

```text
63                         60 59                                  0
+----------------------------+-------------------------------------+
| site entity kind (4 bits)  | kind-specific payload (60 bits)     |
+----------------------------+-------------------------------------+
```

The released kind tags are:

| Tag | Meaning | Payload |
| ---: | --- | --- |
| `0` | unchanged reference-baseline RNG identity | `uint32 reference_position` |
| `1` | variant-changed reference site, shared by both haplotypes | `uint32 reference_position` |
| `2` | variant-changed reference site on haplotype 0 | `uint32 reference_position` |
| `3` | variant-changed reference site on haplotype 1 | `uint32 reference_position` |
| `4` | inserted site shared by both haplotypes | `(uint32 event_ordinal << 2) | offset` |
| `5` | inserted site on haplotype 0 | `(uint32 event_ordinal << 2) | offset` |
| `6` | inserted site on haplotype 1 | `(uint32 event_ordinal << 2) | offset` |

Tags 7--15 and all unused payload bits are zero/reserved. Insertion offsets are
0--3 because the VCF v1 subset permits at most four inserted bases. Event
ordinal `0xffffffff` is forbidden because it is the protocol no-event sentinel.
The largest insertion payload needs only 34 bits. The largest reference payload
needs 32 bits.

Tag 0 deliberately encodes to the numeric reference position. Therefore an
unchanged site continues to use exactly the reference-only Beta RNG address and
probability. A base or context changed by a homozygous event uses tag 1; a
heterozygous ALT haplotype uses tag 2 or 3. Inserted sites use tag 4, 5, or 6.
Protocol allele ownership is assigned only after comparing both haplotypes. A
tag-0 site is `shared` when it is identical on both; if a heterozygous event
removes or changes the corresponding site on the other haplotype, the same
baseline RNG identity is labeled `reference_haplotype`. Tags 1 and 4 are
shared; tags 2, 3, 5, and 6 are `alternate_haplotype`.

The selected zero-based haplotype must be present in the design-deck
`HaplotypeMask`: numeric 1 contains haplotype 0, numeric 2 contains haplotype 1,
and numeric 3 contains both. Invalid masks, absent haplotypes, reserved tags,
reserved payload bits, and the no-event sentinel fail closed.

`beta_sampler::sample_beta_for_site` consumes this validated entity under
algorithm id `marsaglia-tsang-box-muller-beta-site-entity-v2`. Its Gamma,
Box-Muller, acceptance, and local-index algorithm is unchanged from the
reference-only v1 sampler. Because tag 0 is numerically the reference position,
unchanged sites produce the exact same binary32 result as v1. The new algorithm
identifier prevents callers from silently interpreting a tagged entity as an
ordinary reference coordinate.
