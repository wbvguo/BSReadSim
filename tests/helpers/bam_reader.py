"""Minimal standard-library BAM reader for behavior tests.

The reader decodes the public BAM contract only: header text, reference
dictionary, standard alignment fields, and typed auxiliary tags. It performs
no BSReadSim-specific interpretation, so tests built on it remain valid when
the simulator internals change.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
import struct

_CORE = struct.Struct("<iiIIiiii")
_CIGAR_OPERATIONS = "MIDNSHP=XB"
_SEQUENCE_CODES = "=ACMGRSVTWYHKDBN"
_SCALAR_TYPES = {
    "c": "<b", "C": "<B", "s": "<h", "S": "<H",
    "i": "<i", "I": "<I", "f": "<f",
}


@dataclass(frozen=True)
class BamRecord:
    query_name: str
    flag: int
    reference_id: int
    position: int  # zero-based leftmost reference position
    mapq: int
    cigar: tuple[tuple[int, str], ...]
    next_reference_id: int
    next_position: int
    template_length: int
    sequence: str  # reference-forward, as stored in BAM
    quality: tuple[int, ...]  # Phred values, reference-forward
    tags: dict

    @property
    def is_reverse(self) -> bool:
        return bool(self.flag & 0x10)

    @property
    def mate_number(self) -> int:
        if self.flag & 0x40:
            return 1
        if self.flag & 0x80:
            return 2
        return 1

    @property
    def reference_end(self) -> int:
        """Zero-based exclusive end of the aligned reference interval."""
        consumed = sum(length for length, op in self.cigar if op in "MDN=X")
        return self.position + consumed

    def aligned_pairs(self):
        """Yield (query_index, reference_position, operation) for M/I/D."""
        query = 0
        reference = self.position
        for length, op in self.cigar:
            if op in "M=X":
                for offset in range(length):
                    yield query + offset, reference + offset, "M"
                query += length
                reference += length
            elif op in "IS":
                for offset in range(length):
                    yield query + offset, None, op
                query += length
            elif op in "DN":
                for offset in range(length):
                    yield None, reference + offset, op
                reference += length


@dataclass(frozen=True)
class BamFile:
    header: str
    references: tuple[tuple[str, int], ...]
    records: tuple[BamRecord, ...]


def _parse_tags(data: bytes) -> dict:
    tags = {}
    cursor = 0
    while cursor < len(data):
        tag = data[cursor:cursor + 2].decode("ascii")
        kind = chr(data[cursor + 2])
        cursor += 3
        if kind == "A":
            tags[tag] = chr(data[cursor])
            cursor += 1
        elif kind in _SCALAR_TYPES:
            fmt = _SCALAR_TYPES[kind]
            tags[tag] = struct.unpack_from(fmt, data, cursor)[0]
            cursor += struct.calcsize(fmt)
        elif kind in "ZH":
            end = data.index(b"\x00", cursor)
            tags[tag] = data[cursor:end].decode("ascii")
            cursor = end + 1
        elif kind == "B":
            subtype = chr(data[cursor])
            count = struct.unpack_from("<i", data, cursor + 1)[0]
            cursor += 5
            fmt = "<{}{}".format(count, _SCALAR_TYPES[subtype][1])
            tags[tag] = struct.unpack_from(fmt, data, cursor)
            cursor += struct.calcsize(fmt)
        else:
            raise ValueError("unsupported BAM tag type {!r}".format(kind))
    return tags


def _parse_record(block: bytes) -> BamRecord:
    (
        reference_id, position, bin_mq_nl, flag_nc, sequence_length,
        next_reference_id, next_position, template_length,
    ) = _CORE.unpack_from(block)
    name_length = bin_mq_nl & 0xFF
    cigar_count = flag_nc & 0xFFFF
    cursor = _CORE.size
    query_name = block[cursor:cursor + name_length - 1].decode("ascii")
    cursor += name_length
    cigar = []
    for _ in range(cigar_count):
        encoded = struct.unpack_from("<I", block, cursor)[0]
        cigar.append((encoded >> 4, _CIGAR_OPERATIONS[encoded & 0xF]))
        cursor += 4
    packed = block[cursor:cursor + (sequence_length + 1) // 2]
    cursor += (sequence_length + 1) // 2
    sequence = "".join(
        _SEQUENCE_CODES[(packed[i // 2] >> (0 if i % 2 else 4)) & 0xF]
        for i in range(sequence_length)
    )
    quality = tuple(block[cursor:cursor + sequence_length])
    cursor += sequence_length
    return BamRecord(
        query_name=query_name,
        flag=flag_nc >> 16,
        reference_id=reference_id,
        position=position,
        mapq=(bin_mq_nl >> 8) & 0xFF,
        cigar=tuple(cigar),
        next_reference_id=next_reference_id,
        next_position=next_position,
        template_length=template_length,
        sequence=sequence,
        quality=quality,
        tags=_parse_tags(block[cursor:]),
    )


def read_bam(path: Path) -> BamFile:
    data = gzip.decompress(Path(path).read_bytes())
    if not data.startswith(b"BAM\x01"):
        raise ValueError("not a BAM file: {}".format(path))
    cursor = 4
    header_length = struct.unpack_from("<i", data, cursor)[0]
    cursor += 4
    header = data[cursor:cursor + header_length].decode("ascii")
    cursor += header_length
    reference_count = struct.unpack_from("<i", data, cursor)[0]
    cursor += 4
    references = []
    for _ in range(reference_count):
        name_length = struct.unpack_from("<i", data, cursor)[0]
        cursor += 4
        name = data[cursor:cursor + name_length - 1].decode("ascii")
        cursor += name_length
        length = struct.unpack_from("<i", data, cursor)[0]
        cursor += 4
        references.append((name, length))
    records = []
    while cursor < len(data):
        block_size = struct.unpack_from("<i", data, cursor)[0]
        cursor += 4
        records.append(_parse_record(data[cursor:cursor + block_size]))
        cursor += block_size
    return BamFile(header, tuple(references), tuple(records))
