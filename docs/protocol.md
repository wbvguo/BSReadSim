# BSReadSim columnar core stream protocol

Status: normative protocol version `2.0`, the only wire format in BSReadSim
0.2. The key words MUST, MUST NOT, SHOULD, and MAY describe this contract.

## 1. Purpose and non-goals

The protocol transfers consecutive fragments in bounded columnar batches while
preserving these scientific rules:

- C++ constructs both complete haplotypes and the contig-level diploid MethDB
  before a technology provider fragments that contig.
- C++ emits the already selected haplotype template and its methylation
  probabilities; Python does not reconstruct either haplotype or the MethDB.
- Python samples one methylation state per unique fragment site and reuses that
  state in every mate that covers it.
- Python then performs bisulfite conversion, quality generation, sequencing
  error, and output in that order.
- WGBS, RRBS, and TBS use the same wire schema. Variants and indels add data to
  optional truth columns; they do not select a different protocol.

The goals are:

1. preserve the logical `Fragment`, `Mate`, `MethylationSite`, and truth
   contracts;
2. let Python expose authenticated read-only array views without first building
   a large graph of Python objects;
3. remove dense per-base reference-position, event-id, and mate-site-reference
   arrays from the common path;
4. keep the C++ packed representation private, so it can be optimized without
   changing Python or the scientific model; and
5. keep all coordinate and count widths within the checked 32/64-bit boundary;
   the protocol introduces no u96 value.

The protocol does **not** change the probability model, make empirical counts
part of the Bernoulli model, or expose a dynamic technology-plugin ABI.

Unless an example explicitly lists encoded bytes, its YAML/text/C++ notation is
a human-readable decoded view. It does not define a second JSON, YAML, or C++
wire representation.

## 2. Scientific ownership and stage order

The runtime boundary is:

```mermaid
flowchart LR
    R["Reference + variants"] --> H["C++: two complete packed haplotypes"]
    H --> D["C++: complete diploid MethDB"]
    D --> P["C++: WGBS / RRBS / TBS fragment planning"]
    P --> B["Protocol 2.0: selected fragment batches"]
    B --> S["Python: site-state model"]
    S --> C["Python: bisulfite conversion"]
    C --> Q["Python: quality"]
    Q --> E["Python: sequencing error"]
    E --> O["Python: FASTQ + optional Truth"]
```

The normative order is:

1. Validate the normalized config and immutable input identities.
2. Build both haplotypes for the current contig.
3. Discover the complete haplotype-aware methylation catalog and apply
   `Beta -> CGmap/bedMethyl -> ASM/ASM BED` precedence.
4. Let the selected technology provider sample a physical fragment.
5. Materialize that fragment from haplotype 0 or 1, including its unique
   methylation sites and optional sparse truth.
6. Emit consecutive fragments in a columnar batch.
7. In Python, sample site states before any mate-specific processing.
8. Project the shared states into oriented mates, perform conversion, generate
   quality, apply sequencing error, then serialize output.

This order is important. A context created by an insertion or by two bases
joined across a deletion is discovered on the complete haplotype, even when
the causal variant lies just outside the emitted fragment.

### Why two whole haplotypes are not sent to Python

C++ still constructs two logically complete haplotypes. A fragment carries:

- `haplotype = 0` or `1`;
- the exact selected `template_bases`;
- its original-reference envelope; and
- optional sparse variant projection.

Those values completely identify which haplotype supplied the fragment. Sending
both whole contig arrays for every Python worker would duplicate gigabytes of
state and would make Python repeat a C++ responsibility. It adds no information
needed by conversion, quality, error, FASTQ, or the site-state model.

Example: if haplotype 0 contains `ACGTACGT` and haplotype 1 contains
`ACATATCGT`, a row with `haplotype=1` and
`template_bases=ACATATCGT` is unambiguous. Full Truth additionally carries the
SNV and insertion events that explain the difference from the reference.

## 3. Width, coordinate, and base contracts

### 3.1 Integer widths

| Domain | Wire type | Rule |
|---|---:|---|
| base, enum, boolean | `u8` | Exact enumerated value |
| reserved/flag word | `u16` or `u32` | Unknown set bits are fatal |
| contig index | `u32` | Index into the header contig table |
| contig/haplotype coordinate | `u32` | 0-based; contig and physical haplotype length must fit |
| fragment ordinal | `u32` | Consecutive from zero; zero-extended to `u64` for the existing RNG address |
| per-batch counts and offsets | `u32` | Entire frame is at most 64 MiB |
| master seed | `u64` | Zero is valid |
| candidate ordinal, skipped count, totals | `u64` | Never narrowed during accounting |

The normalized requested fragment count is already checked `u32`, so a
fragment ordinal fits `u32`. The internal fragment-candidate ordinal remains
`u64` because rejected candidates may exceed emitted fragments. The Bernoulli
RNG address uses the ordinal zero-extended to `u64`.

Example: wire ordinal `4_000_000_000` is valid and becomes the RNG entity
`uint64(4_000_000_000)`. A contig length of `4_500_000_000` is rejected rather
than widened into a new coordinate type.

### 3.2 Coordinates

All intervals are 0-based and half-open.

- `reference_begin/reference_end` address the original FASTA.
- `haplotype_begin/haplotype_end` are private C++ physical coordinates and are
  not on the wire.
- `template_offset` addresses one selected fragment template.
- `read_offset` is derived after mate orientation.
- An inserted base has no original-reference position.

Example: `reference=[100,108)` contains eight reference bases. A one-base
insertion may produce a nine-base template without changing that reference
envelope.

### 3.3 Base encoding

Wire bases use one byte:

| Value | Base |
|---:|---|
| 0 | A |
| 1 | C |
| 2 | G |
| 3 | T |
| 4 | N |

C++ MAY store A/C/G/T internally as two bits plus an ambiguity bitmap. It MUST
expand the selected template to this typed byte encoding at the protocol
boundary.

Example: `ACGTN` is encoded as `[0, 1, 2, 3, 4]`.

An N MUST NOT be randomized independently for different fragments. The header
records one ambiguity policy:

- `PRESERVE_N=0`: the emitted template retains base code 4;
- `RESOLVE_ONCE=1`: C++ resolves the position once while constructing the
  complete haplotype, then every overlapping fragment reuses that base.

`RESOLVE_ONCE` may be enabled only with a separately frozen RNG contract for
that operation. The current implementation remains `PRESERVE_N`. When Full
Truth is enabled, `original_n_template_offsets` retains N provenance under
either policy.

## 4. Primitive binary encoding

- All integers and IEEE-754 binary32 values are little-endian.
- `u8`, `u16`, `u32`, `u64`, and binary32 have their usual fixed widths.
- `string` is `u32 byte_length` followed by valid UTF-8 bytes without NUL.
- `bytes32` is exactly 32 bytes, with no length prefix.
- Booleans are `u8` and MUST be 0 or 1.
- Probabilities MUST be finite binary32 values in `[0,1]`.
- Strings are limited to 1 MiB.
- A frame payload is limited to 64 MiB.
- Every payload ends with 0--3 zero bytes so its encoded length is divisible
  by four. Padding is included in the CRC and MUST be zero.

Batch columns do not carry individual array-length prefixes. Their lengths are
derived from the authenticated batch header and prefix offsets. This is what
allows direct array views.

Example: a `u32` value `300` is encoded as bytes `2c 01 00 00`. A binary32
probability `0.5` is encoded as `00 00 00 3f`.

## 5. Stream and frame envelope

### 5.1 Preamble

The stream begins with this 16-byte preamble:

| Field | Type | Value |
|---|---:|---|
| magic | 8 bytes | ASCII `BSRSTRM\0` |
| protocol_major | `u16` | `2` |
| protocol_minor | `u16` | `0` |
| preamble_flags | `u32` | `0` |

Example diagnostic view:

```text
magic="BSRSTRM\0", major=2, minor=0, flags=0
```

### 5.2 Frame envelope

Every frame uses this envelope:

| Field | Type | Meaning |
|---|---:|---|
| payload_length | `u32` | Payload bytes including final zero padding |
| frame_type | `u8` | `1` header, `2` fragment batch, `3` trailer, `255` error |
| frame_flags | `u8` | Type-specific flags |
| reserved | `u16` | Zero |
| sequence | `u64` | Consecutive from zero |
| payload | raw bytes | Type-specific payload |
| crc32c | `u32` | Castagnoli CRC over the complete envelope fields and payload |

The CRC input is the little-endian encoding of payload length, type, flags,
reserved, sequence, followed by payload. A normal stream is:

```text
preamble
header(sequence=0)
fragment_batch(sequence=1)
...
trailer(sequence=n)
EOF
```

Example envelope:

```text
payload_length=2048
frame_type=2
frame_flags=0x01
reserved=0
sequence=1
crc32c=CRC32C(encoded envelope fields || payload)
```

The trailer `stream_sha256` covers the exact preamble plus complete encoded
header and batch frames, including their CRCs, but excludes the trailer.

### 5.3 Version identity

Python accepts exactly major `2`, minor `0`, records that pair in the successful
manifest, and rejects every other preamble. The core has no protocol selector
or fallback and emits only this version. Protocol version is an execution
format identity, not a biological-model input, and does not appear in the run
configuration or any RNG address. Package and format versions are mapped in
[versioning.md](versioning.md).

## 6. Header frame

Header frame flags MUST be zero. Its fields appear in this order:

| Field | Type | Rule |
|---|---:|---|
| run_id | `string` | UUID text selected by Python |
| core_version | `string` | Released core semantic version |
| config_schema_version | `string` | Normalized config schema |
| rng_contract | `string` | Exact RNG contract identifier |
| master_seed | `u64` | Actual seed |
| normalized_config_sha256 | `bytes32` | Digest of canonical normalized config |
| technology | `u8` | `1` WGBS, `2` RRBS, `3` TBS |
| truth_columns | `u8` | `0` NONE, `1` FULL |
| mates_per_fragment | `u8` | `1` SE or `2` PE |
| base_encoding | `u8` | `1` ACGTN_U8 |
| ambiguity_policy | `u8` | `0` PRESERVE_N, `1` RESOLVE_ONCE |
| reserved | 3 bytes | All zero |
| read_length_r1 | `u32` | Positive |
| read_length_r2 | `u32` | Zero for SE; positive for PE |
| contig_count | `u32` | Positive |
| contigs | repeated entries | FASTA order |
| padding | 0--3 bytes | Zero |

Each contig entry is:

| Field | Type | Rule |
|---|---:|---|
| name | `string` | Non-empty and unique |
| length | `u32` | Positive |
| reference_sha256 | `bytes32` | Digest of uppercase reference sequence |

Technology is included for identity and diagnostics. Python MUST use the same
post-processing stages for all three values; it MUST NOT implement a second
RRBS/TBS fragment generator.

Example diagnostic form:

```yaml
run_id: "72d5658e-d074-4e7d-bd8c-1b1ce83cfe92"
core_version: "2.0.0"
config_schema_version: "1.0"
rng_contract: "philox4x32-10+philox-domain-v2"
master_seed: 42
normalized_config_sha256: "9f1c...e8a2"   # abbreviated only in this example
technology: WGBS
truth_columns: FULL
mates_per_fragment: 2
read_length_r1: 4
read_length_r2: 4
base_encoding: ACGTN_U8
ambiguity_policy: PRESERVE_N
contigs:
  - name: chrMini
    length: 1000
    reference_sha256: "4b72...0031"       # abbreviated only in this example
```

Real digest fields always contain all 32 bytes.

`truth_columns` is a transfer requirement derived from output policy. It MUST
be `FULL` for debug JSON Truth or truth BAM, and `NONE` only when neither
artifact needs the sparse projection. It is not a second protocol and cannot
be set directly in user configuration.

## 7. Fragment-batch frame

Frame type 2 contains one or more consecutive fragments. It uses one flag:

| Bit | Name | Meaning |
|---:|---|---|
| 0 | `TRUTH_COLUMNS_PRESENT` | The complete sparse truth section is appended |
| 1--7 | reserved | MUST be zero |

If the header says `truth_columns=FULL`, every batch MUST set bit 0. If it says
`NONE`, every batch MUST clear it. The flag means that truth **columns are
present**; it does not claim that a fragment contains a mutation.

### 7.1 Batch header

The payload begins with five `u32` values:

| Field | Meaning |
|---|---|
| first_fragment_ordinal | Ordinal of row zero |
| fragment_count (`F`) | Number of fragment rows |
| template_base_count (`B`) | Flat template-base elements |
| mate_count (`M`) | Flat mate rows |
| methylation_site_count (`S`) | Flat unique-site rows |

`F` MUST be positive. The batch represents ordinals
`first_fragment_ordinal ... first_fragment_ordinal + F - 1`, and that range
must follow the preceding batch without a gap.

Example:

```text
first_fragment_ordinal=42, F=2, B=14, M=4, S=7
```

This represents fragment ordinals 42 and 43. Ordinals are implicit rather than
repeated in a column.

### 7.2 Common column layout

After the batch header, columns occur in this exact order:

| Column | Element type | Elements |
|---|---:|---:|
| contig_indices | `u32` | `F` |
| reference_begins | `u32` | `F` |
| reference_ends | `u32` | `F` |
| template_offsets | `u32` | `F+1` |
| mate_offsets | `u32` | `F+1` |
| site_offsets | `u32` | `F+1` |
| mate_template_begins | `u32` | `M` |
| mate_template_ends | `u32` | `M` |
| site_template_offsets | `u32` | `S` |
| site_probabilities | binary32 | `S` |
| haplotypes | `u8` | `F` |
| capture_strands | `u8` | `F` |
| mate_indices | `u8` | `M` |
| mate_reverse_complements | `u8` | `M` |
| site_contexts | `u8` | `S` |
| site_sources | `u8` | `S` |
| site_alleles | `u8` | `S` |
| template_bases | `u8` | `B` |

The three prefix-offset columns begin at zero, are monotonically
non-decreasing, and end at `B`, `M`, and `S` respectively. All four-byte
columns precede the byte columns, so an aligned payload provides naturally
aligned read-only views.

The source and allele columns remain common because the built-in site-state and
post-processing stages consume one stable column layout in both output modes.

### 7.3 Fragment row

A logical fragment row is reconstructed from:

```text
ordinal          = first_fragment_ordinal + row
contig           = contig_indices[row]
reference        = [reference_begins[row], reference_ends[row])
haplotype        = haplotypes[row]
capture_strand   = capture_strands[row]
template         = template_bases[
                     template_offsets[row] : template_offsets[row + 1]]
mates            = mate columns[mate_offsets[row] : mate_offsets[row + 1]]
sites            = site columns[site_offsets[row] : site_offsets[row + 1]]
```

`capture_strand` uses `0=unknown`, `1=forward`, and `2=reverse`.

Two-row example used throughout this document:

| Ordinal | Contig | Haplotype | Reference | Template |
|---:|---:|---:|---|---|
| 42 | 0 | 0 | `[100,108)` | `ACGTCAGT` |
| 43 | 0 | 1 | `[200,205)` | `TTCCGA` |

The second template is six bases over a five-base reference envelope because
it contains one inserted G.

Its fragment columns are:

```text
contig_indices  = [0, 0]
reference_begins= [100, 200]
reference_ends  = [108, 205]
haplotypes      = [0, 1]
capture_strands = [0, 0]
```

### 7.4 Flat template storage

Templates are copied once into a flat byte column:

```text
template_offsets = [0, 8, 14]
template_bases   = [0,1,2,3,1,0,2,3, 3,3,1,1,2,0]
                    A C G T C A G T  T T C C G A
```

Prefix offsets describe storage adjacency, not genomic adjacency. If two
physical source slices overlap:

```text
fragment 10 source = haplotype[0:8]
fragment 11 source = haplotype[6:14]
```

their independent wire copies may use:

```text
template_offsets = [0, 8, 16]
```

The overlap at physical positions 6--7 is preserved in both templates. It does
not require overlapping ranges in the batch buffer.

### 7.5 Mate rows and derived site projection

Each fragment has one SE mate or two PE mates. Mate rows are grouped by
`mate_offsets` and ordered by `mate_index`.

For the two-row example:

```text
mate_offsets             = [0, 2, 4]
mate_indices             = [0, 1, 0, 1]
mate_reverse_complements = [0, 1, 0, 1]
mate_template_begins     = [0, 4, 0, 2]
mate_template_ends       = [4, 8, 4, 6]
```

Thus fragment 43 uses R1 `template[0:4]` and overlapping R2
`template[2:6]` in reverse-complement orientation.

The protocol does not transmit a dense mate-to-site reference array. Sites are
sorted by template offset, so Python derives membership:

```text
forward read_offset = site_template_offset - mate_template_begin
reverse read_offset = mate_template_end - 1 - site_template_offset
```

Example: fragment 43 has sites at template offsets `[2,3,4]`. R1 `[0,4)` sees
offsets `2,3` at read offsets `2,3`. Reverse R2 `[2,6)` sees all three at read
offsets `3,2,1`. Both mates refer to the same fragment-level site rows, so
Python samples each latent state once.

### 7.6 Methylation-site rows

Site rows are unique within a fragment, sorted by strictly increasing
`site_template_offset`, and get an implicit zero-based `site_index` within that
fragment's `site_offsets` range.

Context values are unchanged:

| Context | C-oriented | G-oriented |
|---|---:|---:|
| CG | 1 | 9 |
| CHG | 3 | 11 |
| CHH | 7 | 15 |

Source uses `1=CGmap-compatible profile`, `2=ASM`, `3=Beta`, `4=pooled
CGmap-compatible profile`. Thus bedMethyl preserves source value 1 or 4 and
ASM BED preserves source value 2 without changing the wire enum. Allele uses
`0=shared`, `1=reference haplotype`, `2=alternate haplotype`.

Example:

```text
site_offsets          = [0, 4, 7]
site_template_offsets = [1, 2, 4, 6,  2, 3, 4]
site_probabilities    = [.80,.75,.30,.35, .20,.90,.85]
site_contexts         = [1, 9, 3,11,  3, 1, 9]
site_sources          = [3, 3, 3, 3,  3, 3, 3]
site_alleles          = [0, 0, 0, 0,  0, 0, 2]
```

For fragment 42, `site_index=0` is template offset 1 with probability 0.80.
For fragment 43, site indices restart at zero; its inserted G is index 2,
context `CG_G=9`, source Beta, alternate-haplotype allele.

The complete contig MethDB owns these probabilities. Overlapping fragments
copy the same stored probability; neither C++ nor Python redraws a methylation
level per fragment.

## 8. Optional sparse Full-Truth section

When `TRUTH_COLUMNS_PRESENT` is set, zero padding aligns the current cursor to
four bytes and a truth section follows. It reconstructs the logical
`reference_positions`, `base_event_ids`, variant events, site reference
positions, and original-N provenance without storing three dense arrays in
every fragment.

### 8.1 Truth header

Five `u32` counts appear first:

| Field | Symbol |
|---|---|
| projection_run_count | `P` |
| variant_event_count | `V` |
| event_ref_base_count | `RB` |
| event_alt_base_count | `AB` |
| original_n_count | `N` |

Example for the two fragments above:

```text
P=3, V=1, RB=0, AB=1, N=0
```

### 8.2 Truth column layout

Columns then occur in this exact order:

| Column | Type | Elements |
|---|---:|---:|
| projection_offsets | `u32` | `F+1` |
| event_offsets | `u32` | `F+1` |
| original_n_offsets | `u32` | `F+1` |
| projection_template_begins | `u32` | `P` |
| projection_template_ends | `u32` | `P` |
| projection_reference_begins | `u32` | `P` |
| event_ids | `u32` | `V` |
| event_reference_begins | `u32` | `V` |
| event_reference_ends | `u32` | `V` |
| event_template_begins | `u32` | `V` |
| event_template_ends | `u32` | `V` |
| event_ref_offsets | `u32` | `V+1` |
| event_alt_offsets | `u32` | `V+1` |
| site_reference_positions | `u32` | `S` |
| original_n_template_offsets | `u32` | `N` |
| event_kinds | `u8` | `V` |
| event_phased_haplotypes | `u8` | `V` |
| event_ref_bases | `u8` | `RB` |
| event_alt_bases | `u8` | `AB` |
| padding | bytes | 0--3 zero bytes |

`site_reference_positions` uses `0xffffffff` for an inserted site. That
sentinel is valid because a contig is constrained to at most `UINT32_MAX`
bases, so valid positions end at `UINT32_MAX-1`. The wire view retains the
`u32` sentinel; the per-base Truth projection exposes it as signed `-1`.

### 8.3 Projection run

A projection run is the logical structure:

```text
ProjectionRun {
    template_begin: u32
    template_end: u32
    reference_begin: u32
}
```

For every `t` in `[template_begin,template_end)`:

```text
reference_position(t) = reference_begin + (t - template_begin)
```

Runs are sorted, non-empty, and do not overlap. Mapped template bases,
including an SNV base, occur in exactly one run. Inserted template bases occur
in no run and are instead covered by one insertion event.

Runs MUST be canonical and maximal within each fragment.  They are ordered by
strictly increasing `template_begin`; their reference intervals are also
strictly increasing and stay inside the fragment reference envelope.  Two
adjacent runs for which both template and reference coordinates continue by
one MUST be coalesced.  These requirements prevent two encoders from producing
different bytes for the same dense reference-position projection.

Example:

```text
projection_offsets          = [0, 1, 3]
projection_template_begins  = [0, 0, 5]
projection_template_ends    = [8, 4, 6]
projection_reference_begins = [100, 200, 204]
```

Fragment 42 is one linear run. Fragment 43 maps offsets 0--3 to references
200--203, skips inserted template offset 4, then maps offset 5 to reference
204.

This reconstructs dense per-base reference positions:

```text
fragment 43: [200, 201, 202, 203, -1, 204]
```

### 8.4 Variant event

A logical event row is:

```text
VariantEvent {
    event_id: u32
    kind: u8                    # 1 SNV, 2 insertion, 3 deletion
    phased_haplotype: u8        # 0, 1, or 255
    reference: [u32, u32)
    template: [u32, u32)
    ref_bases: byte slice
    alt_bases: byte slice
}
```

Event rows are grouped by `event_offsets`. Event IDs are stable per-contig typed
event ordinals, strictly increasing within each fragment, and unique within a
fragment.  A phased haplotype is either `255` or the fragment's selected
haplotype. `event_ref_offsets` and `event_alt_offsets` are ordinary prefix
arrays into the two flat base blobs.

Example insertion in fragment 43:

```text
event_offsets             = [0, 0, 1]
event_ids                 = [7]
event_kinds               = [2]          # insertion
event_phased_haplotypes   = [1]
event_reference_begins    = [204]
event_reference_ends      = [204]
event_template_begins     = [4]
event_template_ends       = [5]
event_ref_offsets         = [0, 0]
event_alt_offsets         = [0, 1]
event_ref_bases           = []
event_alt_bases           = [2]          # G
```

An SNV has equal non-zero reference and template spans. An insertion has an
empty reference interval and a non-empty template span. A deletion has a
non-empty reference interval and an empty template span at its physical
boundary.  Every reference span stays inside the fragment reference envelope;
every template boundary is in `[0, template_length]`.  For every event,
`len(ref_bases)` equals its reference span and `len(alt_bases)` equals its
template span.  Event bases use only A/C/G/T codes `0..3`.  Non-empty event
template spans MUST NOT overlap each other, so every reconstructed dense
`base_event_ids` element has at most one value.

The dense `base_event_ids` Truth field is reconstructed by assigning an event ID
to its non-empty template span. Every other template base gets
`0xffffffff`. A deletion remains represented even though it has no emitted
base.

### 8.5 Site reference positions

Mapped sites carry their original position; inserted sites carry the sentinel.

Example:

```text
site_reference_positions =
    [101,102,104,106, 202,203,0xffffffff]
```

The last site is the inserted G in fragment 43. Its event is found by locating
the insertion event whose template span contains that site's template offset.

### 8.6 Original-N provenance

`original_n_offsets` groups sorted template offsets that originated from a
reference/haplotype ambiguity.

Offsets are strictly increasing within each fragment, in range, and identify
exactly the template elements whose source base was ambiguous.  Under
`PRESERVE_N`, each identified template element MUST contain base code `4` and
every base-code-`4` element MUST be identified.  The first implementation
supports only `PRESERVE_N`; `RESOLVE_ONCE` remains gated on its separately
frozen RNG contract.

Example:

```text
original_n_offsets          = [0, 1]
original_n_template_offsets = [2]
template_bases              = [0, 1, 2, 3]   # A C G T after one-time resolution
```

This says template offset 2 was originally N even though `RESOLVE_ONCE` fixed
it to G. Under `PRESERVE_N`, the same template element would be base code 4.
Later fragment materialization and Python MUST NOT resolve it again.

### 8.7 Sparse-truth completeness

For every fragment, the union of:

1. projection-run template intervals; and
2. insertion-event template intervals

MUST cover every template offset exactly once. SNV event spans overlap mapped
projection runs because they add event provenance without changing coordinate
mapping. Deletion events have empty template spans. These rules let a strict
decoder reconstruct every per-base Truth annotation.

## 9. Trailer frame

Trailer frame flags MUST be zero. Fields appear in this order:

| Field | Type |
|---|---:|
| fragment_count | `u64` |
| fragment_batch_count | `u64` |
| mate_count | `u64` |
| template_base_count | `u64` |
| methylation_site_count | `u64` |
| skipped_fragment_count | `u64` |
| per_contig_fragment_count | `u32` |
| per_contig_fragment_counts | repeated `u64` |
| stream_sha256 | `bytes32` |
| padding | 0--3 zero bytes |

All values MUST match counts independently accumulated by Python.

Example:

```yaml
fragment_count: 2
fragment_batch_count: 1
mate_count: 4
template_base_count: 14
methylation_site_count: 7
skipped_fragment_count: 0
per_contig_fragment_counts: [2]
stream_sha256: "6a31...810c"   # abbreviated only in this example
```

## 10. Error frame

An error payload is:

| Field | Type |
|---|---:|
| error_code | `u32` |
| message | `string` |
| padding | 0--3 zero bytes |

It MAY occur only after a valid header, MUST be the last frame, and the core
MUST exit non-zero. Errors found during config/input/capability preflight MUST
instead occur before the preamble, producing zero stdout bytes.

Example:

```text
error_code=1204
message="fragment batch exceeds payload limit"
```

Python never uses the diagnostic string for control flow.


## 11. Python post-processing from a batch

For each fragment row, Python:

1. obtains its site range from `site_offsets`;
2. samples the built-in Bernoulli state once for those fragment-level sites;
3. validates exactly one state per site;
4. finds each mate's covered site interval using sorted
   `site_template_offsets`;
5. derives oriented read offsets using the formulas in section 7.5;
6. applies the shared states during bisulfite conversion;
7. generates quality;
8. samples the quality-conditioned sequencing-error category using the
   released cumulative-probability/search contract; and
9. formats FASTQ and, if requested, reconstructs full per-base truth from the
   sparse projection.

Example site-model input for fragment 42:

```text
ordinal: 42
contig: chrMini
template: ACGTCAGT
site_template_offsets: [1,2,4,6]
methylation_probabilities: [.80,.75,.30,.35]
contexts: [CG_C,CG_G,CHG_C,CHG_G]
inter-site template distances: [1,2,2]
```

Independent Bernoulli uses each probability directly. Alternative methylation
models are outside this protocol and are not dynamically loaded by the runtime.

## 12. Validation and fail-closed rules

Python MUST reject:

- a non-consecutive batch ordinal range;
- an offset vector that does not begin at zero, is not monotone, or does not
  end at its authenticated flat count;
- an invalid fragment coordinate, contig, haplotype, strand, base, mate, site,
  probability, context, source, or allele;
- mate rows not ordered as exactly one `0` for SE or `0,1` for PE;
- a mate slice outside its fragment template or with a length inconsistent with
  the header;
- non-increasing fragment-local site offsets or a site centered on a base
  incompatible with its C/G-oriented context;
- a truth flag inconsistent with the header policy;
- a malformed projection/event cover, inserted position without an insertion
  event, event bases inconsistent with its kind, or non-zero padding;
- any count disagreement among headers, prefix arrays, observed rows, and the
  trailer; or
- an error frame, truncated frame, CRC/digest mismatch, trailing bytes, or
  non-zero core exit status.

The decoder MUST validate the frame CRC before exposing array views. A batch
buffer MUST remain immutable and alive while any view or worker uses it.
Failure cancels workers, terminates and waits for the core, removes staged
outputs, and publishes no successful manifest.
