# C++ core launch contract

Status: normative argv boundary from the Python orchestrator to
`htsim-core`. Arguments are passed directly to `exec`; a shell is never
involved. Each option and value is a separate argv element.

The Python side MUST materialize schema defaults, choose the master seed, hash
inputs, and then pass the complete projection below.  The C++ side MUST reject
unknown, duplicate, missing, malformed, or cross-field-inconsistent options
before reading the reference or writing a protocol preamble.

## Execution-format envelope

Python passes `--truth-columns none|full` and
`--protocol-batch-fragments U32` on every launch. The batch limit is restricted
to `[1,64]`; a frame may contain fewer fragments when the payload-size bound is
reached. Both options control transport work and backpressure only. They cannot
alter the normalized scientific config, RNG addresses, fragment allocation,
haplotype choice, methylation catalogs, or probabilities.

A direct core invocation may omit them and receives `none` and `64`. There is
no protocol selector: the executable writes only protocol 2.0.

The exact helper envelope `htsim-core --sam-to-bam LEVEL` is reserved for the
Python output transaction. It reads SAM from stdin, writes BAM to stdout, and
accepts a BGZF compression level in `[0,9]`. HTSlib parses every header and
record and writes the final BGZF EOF block. This helper does not generate
fragments, accept biological options, or choose a destination path; malformed
SAM or an I/O/finalization error returns nonzero and aborts publication. See
[truth-bam-v1.md](truth-bam-v1.md).

`htsim-core rrbs-catalog [core contract options]` is the catalog-only envelope
used by the direct `bsreadsim catalog rrbs -r ... --cut-site ...` command. The
Python CLI constructs this projection in memory; it does not read or write a
JSON config. The native mode requires an RRBS uniform projection, forbids
`--rrbs-candidate-bed`, and writes the plain candidate BED rather than protocol
frames. All reference, variant, cut-site, and fragment-shape options otherwise
pass through the same validated core contract.

## Identity and inputs

Required options:

- `--run-id UUID`
- `--config-sha256 HEX64`
- `--seed U64_DECIMAL`
- `--reference PATH`
- `--reference-sha256 HEX64`
- `--technology WGBS|RRBS|TBS`

Optional inputs use paired path and digest options: `--vcf PATH
--vcf-sha256 HEX64`, `--cgmap PATH --cgmap-sha256 HEX64`, `--bed-methyl PATH
--bed-methyl-sha256 HEX64`, `--asm PATH --asm-sha256 HEX64`, and `--asm-bed
PATH --asm-bed-sha256 HEX64`. A path and its digest MUST appear together.
CGmap and bedMethyl are mutually exclusive; ASM and ASM BED are mutually
exclusive. Either ASM form requires `--vcf`.

## Fragment projection

Required options:

- `--paired-end true|false`
- `--read-length-1 U32`
- `--insert-min U32`
- `--insert-mean U32`
- `--insert-max U32`
- `--insert-stddev NONNEGATIVE_NUMBER`
- exactly one of `--depth POSITIVE_NUMBER` and `--read-pairs U32`
- `--max-ambiguous-fraction PROBABILITY`
- `--chunk-size U32`

`--core-workers U32` is optional and defaults to `1`; values are restricted to
`[1,64]`. It controls the bounded C++ protocol-frame or columnar-batch
preparation pool. Fragment ordinals, frame sequences, the stream SHA-256, and
stdout ownership remain single ordered contracts regardless of this value.

Paired-end runs also require `--read-length-2 U32`. The current generator
requires R1 and R2 to have equal lengths even though the wire
header carries them independently; single-end runs MUST omit
`--read-length-2`.
Every configured read length MUST be less than or equal to `--insert-min`, so
every sampled template can supply a complete mate without implicit adapter or
read-through sequence.

At the current released checkpoint, `--depth` is supported for WGBS only and
is converted to a fragment count entirely inside C++ according to
[depth-count-v1.md](depth-count-v1.md). RRBS and TBS require `--read-pairs`.

RRBS requires one or more repeated `--rrbs-cut-site MOTIF` options and accepts
an optional `--rrbs-candidate-bed PATH`. That path deliberately has no paired
digest option; C++ verifies it by regenerated-row matching according to
[rrbs-candidate-bed-v1.md](rrbs-candidate-bed-v1.md). RRBS forbids TBS options.
TBS requires `--tbs-bed PATH`, `--tbs-bed-sha256 HEX64`, and
`--tbs-center-stddev NONNEGATIVE_NUMBER`, and forbids RRBS cut sites.  WGBS
forbids both technology-specific projections.

At the current released capability checkpoint, TBS requires
`insert_min == insert_mean == insert_max` and `insert_stddev == 0`. WGBS may
use the addressed variable-insert contracts in
[variable-wgbs-v1.md](variable-wgbs-v1.md) and, for typed VCF/de novo catalogs,
[variable-haplotype-wgbs-v1.md](variable-haplotype-wgbs-v1.md). Typed VCF
catalogs retain optional CGmap and ASM overlays on that variable-span path.
Target-GC WGBS is a deliberately narrower fixed-insert, reference-only path.
TBS center dispersion may be any finite non-negative value and follows the
addressed normal/rejection contract in [tbs-catalog-v1.md](tbs-catalog-v1.md).

Coverage is projected with `--coverage uniform|profile|target-score`.
Artifact-backed WGBS profile coverage requires `--coverage-profile PATH
--coverage-profile-format FORMAT --coverage-profile-version VERSION
--coverage-profile-sha256 HEX64`, with format `tsv` and version
`wgbs-gc-target-v1`, following
[coverage-profile-target-v1.md](coverage-profile-target-v1.md).

RRBS profile coverage instead requires `--rrbs-candidate-bed PATH` and forbids
all four WGBS profile-artifact options. The BED score supplies each fragment's
relative probability as defined by
[rrbs-candidate-bed-v1.md](rrbs-candidate-bed-v1.md). RRBS uniform coverage may
omit the BED or supply it for exact round-trip validation while ignoring score.
`target-score` supports TBS only and consumes exact aggregate target
`output_weight` values from the verified BED6 input, following
[tbs-target-score-v1.md](tbs-target-score-v1.md). Invalid technology/mode
combinations fail before protocol output.

At the current VCF execution checkpoint, paired `--vcf` and `--vcf-sha256`
are accepted for WGBS, RRBS, and TBS. The core emits the strict catalog's SNVs,
insertions, and deletions. WGBS deletion-safe candidate planning follows
[variant-start-index-v1.md](variant-start-index-v1.md). RRBS restriction discovery
and TBS target projection run on the constructed haplotypes according to
[haplotype-fragmentation-v1.md](haplotype-fragmentation-v1.md).
The released VCF path also requires `--update-variant-boundaries true`; the
core never silently ignores a request for reference-only contexts.

At the current CGmap execution checkpoint, paired `--cgmap` and
`--cgmap-sha256` are accepted with any released technology/coverage
combination. The C++ core validates the complete input against the verified
reference before writing the protocol header, then overlays defined values on
reference-equivalent MethDB sites according to
[cgmap-profile-v1.md](cgmap-profile-v1.md).

Paired `--bed-methyl` and `--bed-methyl-sha256` select the explicit
BED9+2/BED9+9 parser in [bed-methyl-profile-v1.md](bed-methyl-profile-v1.md).
They have the same technology, overlay, pooling, and variant-equivalence gates
as CGmap; the parser is never inferred from a filename suffix.

At the current ASM execution checkpoint, paired `--asm` and `--asm-sha256`
require a paired VCF input and `--update-variant-boundaries true`. Before the
protocol preamble, every
ASM row must resolve to its exact typed, heterozygous VCF SNV and to one shared
reference-equivalent methylation site on both complete haplotypes. The core
then assigns `REF_METH` and `ALT_METH` by the frozen `HaplotypeMask` mapping and
emits typed source/allele metadata. The normative format and failure cases are
defined in [asm-profile-v1.md](asm-profile-v1.md); source precedence is
`ASM > CGmap/bedMethyl > Beta`.

Paired `--asm-bed` and `--asm-bed-sha256` select the BED6+6 serialization in
[asm-bed-profile-v1.md](asm-bed-profile-v1.md). It preserves the same exact VCF
SNV link and produces the same typed ASM overlay as an equivalent TSV profile.

## Mutation and methylation-level projection

Required options:

- `--mutation-rate PROBABILITY`
- `--indel-fraction PROBABILITY`
- `--indel-extension-probability PROBABILITY`
- `--homozygous-only true|false`
- `--collect-non-cpg true|false`
- `--cgmap-pool true|false`
- `--update-variant-boundaries true|false`
- `--beta-cg ALPHA,BETA`
- `--beta-chg ALPHA,BETA`
- `--beta-chh ALPHA,BETA`

Beta shape parameters MUST be finite and strictly positive.  Probabilities
MUST be finite values in `[0,1]`.  Paths are opaque non-empty argv strings;
existence and input digests were already checked by Python.  C++ MUST hash the
bytes it opens and compare them with the projected digests before accepting the
input, closing the preparation-to-launch replacement window.

`--cgmap-pool true` requires one paired `--cgmap` or `--bed-methyl` input and
enables the typed context-pool contract in
[cgmap-pool-v1.md](cgmap-pool-v1.md). It is valid for every released technology
and coverage combination. A profile with no defined probabilities fails before
protocol output. `false` retains the position-specific overlay in
[cgmap-profile-v1.md](cgmap-profile-v1.md) or
[bed-methyl-profile-v1.md](bed-methyl-profile-v1.md).

A positive `--mutation-rate` enables the deterministic typed catalog in
[de-novo-mutation-v1.md](de-novo-mutation-v1.md). The released execution path
allows every released technology, requires no VCF and
`--update-variant-boundaries true`, and retains each technology's coverage
gate. Generated SNVs, insertions, and deletions then use the same haplotype
projection, MethDB construction, provider-specific candidate catalog, and
protocol event representation as verified VCF events.

No Python-only publication or sequencing option (conversion, quality, error,
Python worker, output destination, or compression) may cross this boundary.
`--truth-columns` declares an
authenticated transport capability; it is not itself an output-publication
request.
