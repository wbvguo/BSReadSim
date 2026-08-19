# Protocol v2 frozen byte fixtures

Each `.hex` file is lowercase hexadecimal for one immutable binary fixture.
Whitespace is not part of the fixture.  The values were frozen after manually
checking the field offsets against `docs/protocol.md`; tests independently
recompute every frame CRC32C and the stream trailer SHA-256.

| Fixture | Bytes | SHA-256 | Checked content |
|---|---:|---|---|
| `header-none.hex` | 244 | `94dfd95954618c8b57b4b68e8ad63734af610514e9c169f16466293b8fbf28ea` | 16-byte v2 preamble and sequence-0, 208-byte Header payload |
| `batch-none.hex` | 156 | `57d85ec26bf7e7f92394a1b1847a6df99d127015c405c4f025d76bc9e6665ba4` | sequence-1, flags 0, `F=2 B=8 M=2 S=2` |
| `batch-full.hex` | 468 | `2618de26cdf85f69528e60ba100201989c963253bb7a301d57463b7d8691121c` | sequence-1, Truth flag, SNV/insertion/deletion/N columns |
| `trailer-none.hex` | 112 | `2a94417307f8acc1c4546e82cdba61afeb9ac5db552d68b89f01b8a2c910e32c` | sequence-2, one batch/two fragments, seven skipped |
| `error.hex` | 48 | `be10795217958261c7c4fd5c80a3f062fb210265c6bf3df9b05ed8962e550f5c` | sequence-1 terminal error 1204 |

The no-Truth stream is `header-none || batch-none || trailer-none`.  This
concatenation is also used to verify that the trailer digest covers the exact
preamble, Header frame, and Batch frame while excluding the Trailer frame.
