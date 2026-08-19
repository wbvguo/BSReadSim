# Variant-aware fragment builder v1

Status: normative pure boundary from a typed haplotype projection to protocol
v1 `Fragment`.

The builder accepts one `ProjectedInterval` by value, one read-only diploid
methylation catalog, a global `uint64` fragment ordinal, capture strand, and
the read layout. It performs no sampling, file access, logging, or protocol
writing. A temporary projection transfers ownership of its arrays without
copying them.

The projection is self-identifying: its contig index, zero-based haplotype,
non-empty original-reference interval, template bases, reference positions,
per-base event ids, and typed events are validated together. The configured
insert length is the original-reference span. An indel may change the projected
template length, but that length must still fit `uint32` and must be at least
the read length. Partial variants, incomplete event bases, wrong phase, unknown
event ids, and malformed insertion/SNV/deletion shapes fail closed.

R1 is the first `read_length` projected bases. For paired-end data, R2 is the
last `read_length` projected bases and has `reverse_complement=true`. A mate's
reference envelope is the minimum and maximum mapped reference position in its
slice, not template-offset arithmetic. A mate containing only one insertion
uses that event's zero-width anchor. Methylation sites are attached from the
complete-haplotype diploid catalog, and every site in a mate slice receives a
site reference in final read orientation. Overlapping mates therefore share
one fragment-level latent methylation site, including inserted sites.

Payload preflight uses only `uint64` arithmetic. Its exact upper bound is the
protocol fixed fields plus projected base arrays, one possible methylation site
per projected base, oriented mate site references, and every typed event's REF
and ALT bytes. Contig-local coordinates and counts remain `uint32`; no u96
representation crosses this boundary.
