# WGBS/WGS target GC distribution

Use `bsreadsim run wgbs --gc-profile PATH` or
`bsreadsim run wgs --gc-profile PATH` to request a fragment GC distribution.
The profile describes desired output, not a proposal distribution.

## Input

The file contains one finite probability in `[0,1]` per line. It must contain
at least two lines and sum to one within `1e-9`. Blank lines, comments, extra
fields, and surrounding whitespace are rejected.

The first line represents the lowest GC bin and the last line the highest.
Every eligible physical fragment is assigned to a bin using its actual length
and GC count. For GC count `g`, fragment length `F`, and `B` bins:

```text
bin = floor(g * (B - 1) / F + 0.5)
```

## Supported combinations

- A fixed insert (`--insert-mean N --insert-sd 0`) supports reference-only
  runs, VCF, fixed de novo mutation, and ASM.
- Reference-only variable inserts are supported using a mean-insert
  approximation; the manifest reports unreachable requested mass.
- Variable inserts combined with VCF or de novo variants are currently
  rejected because that calibration is not implemented.

A profile requesting positive mass for a bin with no eligible fragments is
rejected before output is published. Use an explicit `--seed` for comparable
runs; batching and worker settings do not change fixed-seed results.
