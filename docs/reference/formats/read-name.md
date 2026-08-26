# Read names

FASTQ identifiers and BAM QNAMEs share one fragment identity:

```text
<contig>:<start>-<end>:<ordinal-hex>
```

FASTQ appends `/1` or `/2`; BAM stores the fragment identity without that
suffix and uses SAM flags for mate identity.

- `contig` is the exact SAM-safe contig name.
- `start` and `end` are the one-based inclusive fragment envelope.
- `ordinal-hex` is the zero-based fragment ordinal in variable-width lowercase
  hexadecimal, with no fixed cardinality limit or `0x` prefix.

The identifier is deterministic for fixed simulation inputs and is independent
of worker count, batching, output format, and completion order.
