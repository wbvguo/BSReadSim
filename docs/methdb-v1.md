# Fixed MethDB snapshot v1

Status: normative catalog contract `bsreadsim-methdb-v1`.

MethDB is a fixed, normalized database of site probabilities and provenance.
It is generated before fragment sampling and is independent of the simulation
seed. `--methdb-seed` controls catalog-only randomness such as unphased VCF,
de novo variants, and generated Beta probabilities. `--seed` controls fragment
selection, sampled states, conversion, qualities, and sequencing errors.

`--save-methdb PATH` writes the exact catalog used by the run. `--methdb PATH`
loads that immutable snapshot and forbids CGmap, bedMethyl, ASM, ASM BED, and
CGmap-pool inputs. A recorded SHA-256 binds the bytes; reference identity and
catalog-defining variant identity are checked before protocol output.

The compact binary snapshot begins with magic `BSRMDB01`, binding digest,
contig count, then per-contig identity and rows. Reference-site rows store
position, binary32 probability, context, and source. Diploid rows additionally
store origin identity and allele. Integers are little-endian; rows are sorted;
truncation, trailing bytes, invalid enums, duplicates, and identity mismatch
fail closed.

MethDB stores probabilities, source, context, and allele. Per-fragment sampled
states and conversion outcomes belong to BAM `zt`/`zx`; variant definitions
belong to VCF or the snapshot binding. This separation prevents duplicated
details while preserving exact replay and position-based lookup.
