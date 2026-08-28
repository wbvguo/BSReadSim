# MethDB methylation profile

MethDB v2 is the immutable serialized form of a prepared BSReadSim
methylation world. It contains the reference baseline, the prepared variant
events, variant-induced reference and insertion overlays, ASM probabilities,
contexts, sources, and allele relationships. A later simulation therefore
needs only the same reference plus the MethDB; an external VCF is not part of
the reload contract.

`--save-methdb` writes the MethDB used by a run to
`OUTPUT/truth/PREFIX.methdb`. A later run loads that path with `--methdb`:

```bash
bsreadsim run wgbs -r reference.fa -o runs/create -n 100000 \
  --mutation-rate 0 --save-methdb \
  --seed-meth 7 --seed 42

bsreadsim run wgbs -r reference.fa -o runs/reuse -n 100000 \
  --methdb runs/create/truth/sim.methdb --seed 99
```

Loading a MethDB profile is mutually exclusive with VCF, de novo mutation,
CGmap, bedMethyl, ASM, ASM BED, and CGmap-pool inputs. The reference identity
must match. Build-time mutation, phasing, and methylation seeds are not reload
requirements because their resolved results are already stored in the file.

`--seed-meth` controls generated Beta probabilities. `--seed-mut` controls de
novo variants, and `--seed-phase` controls unphased VCF assignment. `--seed`
separately controls fragment selection, sampled states, conversion, quality,
and sequencing error.

MethDB is an opaque binary artifact. Do not edit or concatenate it; retain its
original bytes and the associated run manifest. To inspect it, export a
human-readable extended BED instead of modifying the MethDB file:

```bash
bsreadsim export methdb \
  -i runs/create/truth/sim.methdb \
  -o sim.methdb.bed.gz
```

The default is deterministic BGZF and the destination must end in `.bed.gz`.
For uncompressed text, use:

```bash
bsreadsim export methdb \
  -i runs/create/truth/sim.methdb \
  -o sim.methdb.bed \
  --no-compression
```

The exporter creates missing parent directories, refuses to overwrite an
existing file, validates the complete MethDB while decoding it, and streams
the result with bounded memory.

The file starts with the `methdb` magic and version byte `2`. Each contig has
an independently compressed baseline section. Diploid contigs add a logical
`VARIANT_LAYER`, represented by three independent subframes so a reader can
load only what it needs: prepared events, reference overlays, and insertion
overlays. ASM is a separate probability-only section. The footer stores a
verified canonical content root over metadata and the raw digest of every
section.

| Logical part | Physical section | Purpose |
| --- | --- | --- |
| Reference baseline | `BASELINE` | sorted reference sites and default probabilities |
| `VARIANT_LAYER` | `EVENTS` | normalized SNV/insertion/deletion events and haplotype masks |
| `VARIANT_LAYER` | `REF_OVERLAY` | only reference sites changed or removed by events, with causal event ordinals |
| `VARIANT_LAYER` | `INS_OVERLAY` | methylation sites created on inserted bases |
| ASM layer | `ASM` | linked allele-specific probability overrides |

Site coordinates/origin IDs are delta-encoded and row metadata is bit-packed.
Every section has its own uncompressed length and SHA-256 and is compressed as
an independent zlib frame. Consequently, opening a MethDB reads and validates
the directory but does not inflate every chromosome.

Methylation probabilities use a two-byte unsigned
normalized integer (uint16 UNORM): the writer stores
`round(probability * 65535)` and the reader divides that value by `65535` when
expanding it to float32. The resolution is therefore a constant `1/65535`
across `[0, 1]`, both endpoints are exact, and the maximum quantization error is
`1/131070`. The writer uses only a bounded compression buffer, rather than
retaining the serialized profile in memory.

At runtime, only the active contig is materialized. Reference baseline rows are
8 bytes each. Diploid lookup rows are packed into 64-bit words with the layout
`[key:32][probability_u16:16][metadata:16]` and split into shared/haplotype and
reference/insertion arrays. Query cursors advance monotonically through these
sorted arrays. Completed contigs are released, so MethDB RAM scales with the
largest active contig rather than the whole genome.

Without `--save-methdb`, these arrays are built in memory and discarded after
simulation. With `--save-methdb`, the same simulation process writes canonical
sections to a transactional sidecar while constructing the runtime catalog;
it does not run a second whole-genome MethDB build first.

This is the only supported MethDB representation; incompatible historical
files must be regenerated.

## Exported extended BED

The text export identifies itself with:

```text
#format	methdb-bed-v2
```

Lines beginning with `#` are metadata or insertion-origin records. Every other
line is one reference-backed methylation site and begins with the six standard
BED columns. Coordinates are zero-based, half-open, one-base intervals.

| Column | Name | Meaning |
| ---: | --- | --- |
| 1 | `chrom` | reference contig name |
| 2 | `chromStart` | zero-based site position |
| 3 | `chromEnd` | `chromStart + 1` |
| 4 | `name` | stable `methdb:SET:ORIGIN_ID` label |
| 5 | `score` | probability mapped to the BED integer range `0..1000` |
| 6 | `strand` | `+` for a cytosine context, `-` for a guanine context |
| 7 | `set` | `reference`, `shared`, `haplotype-1`, or `haplotype-2` |
| 8 | `origin_id` | MethDB origin identifier; equal to `chromStart` here |
| 9 | `origin_kind` | `reference` for a coordinate-backed row |
| 10 | `variant_event` | `.` for a coordinate-backed row |
| 11 | `insertion_offset` | `.` for a coordinate-backed row |
| 12 | `context` | `CG-C`, `CHG-C`, `CHH-C`, `CG-G`, `CHG-G`, or `CHH-G` |
| 13 | `source` | `beta`, `cgmap`, `pooled-cgmap`, or `asm` |
| 14 | `allele` | `shared`, `reference`, or `alternate` |
| 15 | `probability_u16` | exact stored uint16 UNORM value, `0..65535` |
| 16 | `probability` | decoded probability in `[0, 1]` |

Diploid exports include one `#variant` record per embedded prepared event:

```text
#variant	CONTIG	ORDINAL	START	END	KIND	REF	ALT	HAPLOTYPE_MASK	ID	SOURCE
```

MethDB can also contain methylation sites on inserted bases. The exporter does
not invent a one-base reference interval for those bases. It retains each such
site as a comment record:

```text
#insertion	CONTIG	SET	ORIGIN_ID	VARIANT_EVENT	INSERTION_OFFSET	CONTEXT	SOURCE	ALLELE	PROBABILITY_U16	PROBABILITY
```

This keeps the output acceptable to BED readers, which skip comment lines,
while preserving every MethDB site for human inspection. `#contig` records
retain contig order, length, reference SHA-256, and ploidy mode;
`#binding_sha256` identifies the reference binding, and
`#content_sha256` identifies the canonical section content independently of
compression bytes.
