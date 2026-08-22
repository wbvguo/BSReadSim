# Versioning

BSReadSim package, wire protocol, configuration, RNG, manifest, read-name, and
BAM contracts have independent versions.

| Domain | Current version | Meaning |
| --- | --- | --- |
| Package | `0.4.0` | current unreleased software line |
| C++/Python wire | `2.0` | sole columnar stream accepted by this tree |
| Internal normalized run document | `1.1` | packaged validation schema; not a user input format |
| RNG | `philox4x32-10+philox-domain-v2` | deterministic address contract |
| Manifest | `1.1` | complete manifest-last audit record |
| Read names | `bsreadsim-read-name-v2` | variable-width hexadecimal fragment ordinal |
| annotated BAM | `bsreadsim-bam-v3` | fixed `zt`/`zr`, optional `zf`/`zx` |

## Release history

- `v0.1.0` preserves the old implementation and protocol 1 source.
- `v0.2.0` preserves the earlier protocol-2 output surface.
- The unreleased 0.3 line removed the protocol selector, legacy source,
  production/debug branch, and Details JSONL artifact. It uses the frozen numeric
  Philox domain contract v2 and annotated BAM v2 as the sole detailed details
  product.
- The unreleased 0.4 line adds fixed MethDB snapshots, target-region depth,
  haplotype-aware fixed-insert GC calibration, a state-model interface, and
  annotated BAM v3 fragment realization without changing wire protocol 2.0.

The current tree has no protocol fallback or compatibility serializer. Exact
historical behavior remains recoverable from Git history rather than runtime
branches.

## Change rules

- Compatible clarification keeps a contract version; changed bytes,
  interpretation, or required fields increment it.
- Unsupported protocol or artifact versions fail closed.
- A successful manifest records every observed component version.
- Performance changes must retain scientific counts and equivalent output
  semantics for a fixed seed before timing is accepted.
- Until the first public release, obsolete runtime paths are deleted rather
  than retained as compatibility code.
