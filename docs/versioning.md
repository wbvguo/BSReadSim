# Versioning

BSReadSim package, wire protocol, configuration, RNG, and artifact formats have
independent version numbers. A change to one does not silently renumber the
others.

| Domain | Current version | Meaning |
| --- | --- | --- |
| Package | `0.3.0` | current unreleased software line |
| C++/Python wire | `2.0` | sole columnar stream accepted by this tree |
| Internal normalized run document | `1.0` | packaged validation schema; not a user input format |
| RNG | `philox4x32-10+philox-domain-v2` | deterministic address contract |
| Truth | `1.2` | canonical JSONL artifact format |
| Read names | `bsreadsim-read-name-v1` | FASTQ identifiers and truth-BAM QNAMEs |

## Release history

- `v0.1.0` preserves the old implementation and protocol 1 source.
- `v0.2.0` denotes the protocol-2-default release with FASTQ-only production,
  Full Truth under `--mode debug`, and an explicit protocol-1 compatibility
  path.
- The unreleased 0.3 line removes the protocol selector and legacy source,
  replaces string/BLAKE2b domain derivation with the frozen numeric-stage,
  numeric-contig Philox domain contract v2, and adds the versioned read-name
  and truth-BAM contracts. Fixed-seed and public-output bytes from 0.2 are
  intentionally not preserved.

The current 0.3 tree has no protocol selector, decoder fallback, or legacy
source copy. Exact 0.1 behavior requires checking out the `v0.1.0` tag. This is
intentional: the project has not been publicly released, so current code stays
small instead of carrying an unneeded compatibility branch.

## Change rules

- Compatible clarification of a format keeps its version; any changed bytes,
  interpretation, or required field increments that format's version.
- Unsupported protocol or artifact versions fail closed. There is no automatic
  downgrade.
- A successful manifest records every observed component version rather than
  inferring them from the package number.
- Performance changes must retain exact scientific counts and decompressed
  output bytes for a fixed seed before timing is accepted.
- Until the first public release, obsolete runtime paths are deleted and remain
  recoverable from Git history rather than receiving compatibility code.
