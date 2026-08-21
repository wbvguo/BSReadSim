# Protocol v2.1 frozen byte fixtures

Each `.hex` file is lowercase hexadecimal for one immutable binary fixture.
Whitespace is not part of the fixture.  The values were frozen after manually
checking the field offsets against `docs/protocol.md`; tests independently
recompute every frame CRC32C and the stream trailer SHA-256.

| Fixture | Bytes | SHA-256 | Checked content |
|---|---:|---|---|
| `header-none.hex` | 244 | `5bd464816d26450bbb801661077047092716bfafa1e81d06e924c957c726e61e` | 16-byte v2 preamble and sequence-0, 208-byte Header payload |
| `batch-none.hex` | 156 | `57d85ec26bf7e7f92394a1b1847a6df99d127015c405c4f025d76bc9e6665ba4` | sequence-1, flags 0, `F=2 B=8 M=2 S=2` |
| `batch-full.hex` | 496 | `7fca502060dae056c67a28d68225de2a23fb069b5212af74120137ea6eeca74c` | sequence-1, Details flag, variant IDs/sources and SNV/insertion/deletion/N columns |
| `trailer-none.hex` | 112 | `8065c6bdcc203aeb450b913d3dd1fbbeaa981d4722ea8f8b3cb4486deb5b337f` | sequence-2, one batch/two fragments, seven skipped |
| `error.hex` | 48 | `be10795217958261c7c4fd5c80a3f062fb210265c6bf3df9b05ed8962e550f5c` | sequence-1 terminal error 1204 |

The no-Details stream is `header-none || batch-none || trailer-none`.  This
concatenation is also used to verify that the trailer digest covers the exact
preamble, Header frame, and Batch frame while excluding the Trailer frame.
