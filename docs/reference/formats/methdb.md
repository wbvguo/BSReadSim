# MethDB methylation profile

MethDB is the immutable serialized form of a prepared BSReadSim methylation
profile, including site probabilities, contexts, sources, and allele
relationships. It allows multiple simulations to reuse exactly the same
prepared methylation profile.

`--save-methdb` writes the MethDB used by a run to
`OUTPUT/truth/PREFIX.methdb`. A later run loads that path with `--methdb`:

```bash
bsreadsim run wgbs -r reference.fa -o runs/create -n 100000 \
  --mutation-rate 0 --save-methdb \
  --seed-meth 7 --seed 42

bsreadsim run wgbs -r reference.fa -o runs/reuse -n 100000 \
  --methdb runs/create/truth/sim.methdb --mutation-rate 0 \
  --seed-meth 7 --seed 99
```

Loading a MethDB profile is mutually exclusive with CGmap, bedMethyl, ASM, ASM
BED, and CGmap-pool inputs. The reference identity and prepared variant-set
identity must match those used to create it.

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

The file starts with the `methdb` magic followed by a one-byte internal format
version. Site coordinates/origin IDs are delta-encoded, row metadata is
bit-packed, and the payload is streamed through zlib while it is written.
Methylation probabilities use a two-byte unsigned
normalized integer (uint16 UNORM): the writer stores
`round(probability * 65535)` and the reader divides that value by `65535` when
expanding it to float32. The resolution is therefore a constant `1/65535`
across `[0, 1]`, both endpoints are exact, and the maximum quantization error is
`1/131070`. The writer uses only a bounded compression buffer, rather than
retaining the serialized profile in memory.

This is an on-disk reduction. The current simulation loader still materializes
the validated site arrays in their compact in-memory structs, so peak RAM for
a whole-genome MethDB continues to scale with the number of methylation sites.
That memory is separate from transactional FASTQ/BAM output staging.

This is the only supported MethDB representation; incompatible historical
files must be regenerated.

## Exported extended BED

The text export identifies itself with:

```text
#format	methdb-bed
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

MethDB can also contain methylation sites on inserted bases. The compact
MethDB file stores their variant-event ordinal and insertion offset, but not the
event's standalone reference coordinate. The exporter therefore does not
invent a BED interval. It retains each such site as a comment record:

```text
#insertion	CONTIG	SET	ORIGIN_ID	VARIANT_EVENT	INSERTION_OFFSET	CONTEXT	SOURCE	ALLELE	PROBABILITY_U16	PROBABILITY
```

This keeps the output acceptable to BED readers, which skip comment lines,
while preserving every MethDB site for human inspection. `#contig` records
retain contig order, length, reference SHA-256, and ploidy mode;
`#binding_sha256` identifies the profile's input binding, and `#file_sha256`
identifies the source MethDB bytes.
