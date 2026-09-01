"""The sole columnar stream codec exchanged with ``htsim-core``.

This module owns canonical bytes, strict semantic validation, and immutable
decoded batch views. Process supervision and output publication remain
elsewhere.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
import hashlib
import math
import re
import struct
import uuid
from collections.abc import Iterator, Sequence
from typing import BinaryIO

from ..rng import RNG_CONTRACT

from .._cext import (
    crc32c as _cext_crc32c,
    validate_protocol_batch_columns as _cext_validate_batch_columns,
)


MAGIC = b"BSRSTRM\x00"
PROTOCOL_MAJOR = 2
PROTOCOL_MINOR = 1
PROTOCOL_VERSION = "{}.{}".format(PROTOCOL_MAJOR, PROTOCOL_MINOR)
PREAMBLE_FLAGS = 0

MAX_STRING_BYTES = 1024 * 1024
MAX_FRAME_PAYLOAD = 64 * 1024 * 1024
NO_REFERENCE_POSITION = 0xFFFFFFFF


DETAILS_PRESENT = 0x01

_PREAMBLE = struct.Struct("<8sHHI")
_FRAME_ENVELOPE = struct.Struct("<IBBHQ")
_CRC = struct.Struct("<I")
_U8 = struct.Struct("<B")
_U32 = struct.Struct("<I")
_U64 = struct.Struct("<Q")
_F32 = struct.Struct("<f")
_SEMVER = re.compile(
    r"^(0|[1-9][0-9]*)\."
    r"(0|[1-9][0-9]*)\."
    r"(0|[1-9][0-9]*)"
    r"(?:-((?:0|[1-9][0-9]*|[0-9]*[A-Za-z-][0-9A-Za-z-]*)"
    r"(?:\.(?:0|[1-9][0-9]*|[0-9]*[A-Za-z-][0-9A-Za-z-]*))*))?"
    r"(?:\+([0-9A-Za-z-]+(?:\.[0-9A-Za-z-]+)*))?$"
)


def crc32c(data: bytes) -> int:
    """Return the Castagnoli CRC-32C of a bytes-like value."""
    if not isinstance(data, (bytes, bytearray, memoryview)):
        raise TypeError("crc32c data must be bytes-like")
    return _cext_crc32c(data)


class ProtocolError(ValueError):
    """The stream or a logical value violates the wire contract."""


class CoreReportedError(ProtocolError):
    """The core emitted a valid terminal error frame."""

    def __init__(self, error_code: int, message: str) -> None:
        self.error_code = error_code
        self.message = message
        super().__init__("core error {}: {}".format(error_code, message))


class FrameType(IntEnum):
    HEADER = 1
    FRAGMENT_BATCH = 2
    TRAILER = 3
    ERROR = 255


class Technology(IntEnum):
    WGBS = 1
    RRBS = 2
    TBS = 3
    WGS = 4
    WES = 5
    TS = 6


class BaseEncoding(IntEnum):
    ACGTN_U8 = 1


class AmbiguityPolicy(IntEnum):
    PRESERVE_N = 0


class CaptureStrand(IntEnum):
    UNKNOWN = 0
    FORWARD = 1
    REVERSE = 2


class VariantKind(IntEnum):
    SNV = 1
    INSERTION = 2
    DELETION = 3


class VariantSourceCode(IntEnum):
    VCF = 1
    DE_NOVO = 2
    ASM = 3


class MethylationContext(IntEnum):
    CG_C = 1
    CHG_C = 3
    CHH_C = 7
    CG_G = 9
    CHG_G = 11
    CHH_G = 15


class MethylationSource(IntEnum):
    CGMAP = 1
    ASM = 2
    BETA = 3
    POOLED_CGMAP = 4


class MethylationAllele(IntEnum):
    SHARED = 0
    REFERENCE_HAPLOTYPE = 1
    ALTERNATE_HAPLOTYPE = 2


@dataclass(frozen=True)
class Contig:
    name: str
    length: int
    reference_sha256: bytes


@dataclass(frozen=True)
class Header:
    run_id: str
    core_version: str
    rng_contract: str
    master_seed: int
    normalized_config_sha256: bytes
    technology: Technology
    has_details: bool
    mates_per_fragment: int
    base_encoding: BaseEncoding
    ambiguity_policy: AmbiguityPolicy
    read_length_r1: int
    read_length_r2: int
    contigs: tuple[Contig, ...]


@dataclass(frozen=True)
class FragmentDetails:
    projection_offsets: Sequence[int]
    variant_offsets: Sequence[int]
    original_n_offsets: Sequence[int]
    projection_template_starts: Sequence[int]
    projection_template_ends: Sequence[int]
    projection_reference_starts: Sequence[int]
    variant_indices: Sequence[int]
    variant_id_offsets: Sequence[int]
    variant_reference_starts: Sequence[int]
    variant_reference_ends: Sequence[int]
    variant_template_starts: Sequence[int]
    variant_template_ends: Sequence[int]
    variant_ref_offsets: Sequence[int]
    variant_alt_offsets: Sequence[int]
    site_reference_positions: Sequence[int]
    original_n_template_offsets: Sequence[int]
    variant_sources: Sequence[int]
    variant_kinds: Sequence[int]
    variant_phased_haplotypes: Sequence[int]
    variant_ids: Sequence[int]
    variant_ref_bases: Sequence[int]
    variant_alt_bases: Sequence[int]


@dataclass(frozen=True)
class FragmentBatch:
    first_fragment_ordinal: int
    contig_indices: Sequence[int]
    reference_starts: Sequence[int]
    reference_ends: Sequence[int]
    template_offsets: Sequence[int]
    mate_offsets: Sequence[int]
    site_offsets: Sequence[int]
    mate_template_starts: Sequence[int]
    mate_template_ends: Sequence[int]
    site_template_offsets: Sequence[int]
    site_probabilities: Sequence[float]
    haplotypes: Sequence[int]
    capture_strands: Sequence[int]
    mate_indices: Sequence[int]
    mate_reverse_complements: Sequence[int]
    site_contexts: Sequence[int]
    methylation_sources: Sequence[int]
    site_alleles: Sequence[int]
    template_bases: Sequence[int]
    details: FragmentDetails | None

    @property
    def fragment_count(self) -> int:
        return len(self.contig_indices)

    @property
    def template_base_count(self) -> int:
        return len(self.template_bases)

    @property
    def mate_count(self) -> int:
        return len(self.mate_indices)

    @property
    def methylation_site_count(self) -> int:
        return len(self.site_template_offsets)


@dataclass(frozen=True)
class DecodedBatchView(FragmentBatch):
    """A validated batch whose columns retain one immutable payload."""

    _payload: bytes = field(repr=False, compare=False)

    @property
    def raw_payload(self) -> memoryview:
        return memoryview(self._payload)


@dataclass(frozen=True)
class Trailer:
    fragment_count: int
    fragment_batch_count: int
    mate_count: int
    template_base_count: int
    methylation_site_count: int
    skipped_fragment_count: int
    per_contig_fragment_counts: tuple[int, ...]
    stream_sha256: bytes


@dataclass(frozen=True)
class ErrorFrame:
    error_code: int
    message: str


@dataclass(frozen=True)
class ProtocolStream:
    header: Header
    batches: tuple[DecodedBatchView, ...]
    trailer: Trailer


class _ArrayView(Sequence):
    """Read-only, zero-copy little-endian primitive array."""

    __slots__ = ("_view", "_item", "_count")

    def __init__(
        self,
        backing: memoryview,
        offset: int,
        count: int,
        item: struct.Struct,
    ) -> None:
        size = count * item.size
        self._view = backing[offset : offset + size].toreadonly()
        self._item = item
        self._count = count

    def __len__(self) -> int:
        return self._count

    def __getitem__(self, index):
        if isinstance(index, slice):
            return tuple(self[position] for position in range(*index.indices(self._count)))
        if isinstance(index, bool) or not isinstance(index, int):
            raise TypeError("array indices must be integers or slices")
        if index < 0:
            index += self._count
        if index < 0 or index >= self._count:
            raise IndexError("array index out of range")
        return self._item.unpack_from(self._view, index * self._item.size)[0]

    def __iter__(self):
        for index in range(self._count):
            yield self._item.unpack_from(self._view, index * self._item.size)[0]

    @property
    def raw(self) -> memoryview:
        return self._view


def _require_int(name: str, value: int, minimum: int, maximum: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ProtocolError("{} must be an integer".format(name))
    if value < minimum or value > maximum:
        raise ProtocolError(
            "{} must be in [{}, {}]".format(name, minimum, maximum)
        )


def _require_u8(name: str, value: int) -> None:
    _require_int(name, value, 0, 0xFF)


def _require_u32(name: str, value: int) -> None:
    _require_int(name, value, 0, 0xFFFFFFFF)


def _require_u64(name: str, value: int) -> None:
    _require_int(name, value, 0, 0xFFFFFFFFFFFFFFFF)


def _require_enum(name: str, value: int, enum_type) -> None:
    if isinstance(value, bool) or not isinstance(value, (int, IntEnum)):
        raise ProtocolError("{} must be a {}".format(name, enum_type.__name__))
    try:
        enum_type(int(value))
    except ValueError as error:
        raise ProtocolError(
            "{} is outside {}".format(name, enum_type.__name__)
        ) from error


def _require_digest(name: str, value: bytes) -> None:
    if not isinstance(value, bytes) or len(value) != 32:
        raise ProtocolError("{} must contain exactly 32 bytes".format(name))


def _validated_string(name: str, value: str, *, nonempty: bool = False) -> bytes:
    if not isinstance(value, str):
        raise ProtocolError("{} must be text".format(name))
    if nonempty and not value:
        raise ProtocolError("{} must not be empty".format(name))
    if "\x00" in value:
        raise ProtocolError("{} must not contain NUL".format(name))
    try:
        encoded = value.encode("utf-8")
    except UnicodeError as error:
        raise ProtocolError("{} must be valid UTF-8".format(name)) from error
    if len(encoded) > MAX_STRING_BYTES:
        raise ProtocolError("{} exceeds the 1 MiB string limit".format(name))
    return encoded


def _validate_header(header: Header) -> None:
    if not isinstance(header, Header):
        raise ProtocolError("header must be a protocol Header")
    _validated_string("header.run_id", header.run_id, nonempty=True)
    try:
        parsed_run_id = uuid.UUID(header.run_id)
    except (ValueError, AttributeError) as error:
        raise ProtocolError("header.run_id must be canonical UUID text") from error
    if str(parsed_run_id) != header.run_id:
        raise ProtocolError("header.run_id must be canonical lowercase UUID text")
    _validated_string("header.core_version", header.core_version, nonempty=True)
    if _SEMVER.fullmatch(header.core_version) is None:
        raise ProtocolError("header.core_version must be a semantic version")
    if header.rng_contract != RNG_CONTRACT:
        raise ProtocolError("unsupported RNG contract")
    _require_u64("header.master_seed", header.master_seed)
    _require_digest(
        "header.normalized_config_sha256", header.normalized_config_sha256
    )
    _require_enum("header.technology", header.technology, Technology)
    if not isinstance(header.has_details, bool):
        raise ProtocolError("header.has_details must be a boolean")
    if header.mates_per_fragment not in (1, 2):
        raise ProtocolError("header.mates_per_fragment must be 1 or 2")
    if header.base_encoding != BaseEncoding.ACGTN_U8:
        raise ProtocolError("unsupported base encoding")
    if header.ambiguity_policy != AmbiguityPolicy.PRESERVE_N:
        raise ProtocolError("unsupported ambiguity policy")
    _require_u32("header.read_length_r1", header.read_length_r1)
    _require_u32("header.read_length_r2", header.read_length_r2)
    if header.read_length_r1 == 0:
        raise ProtocolError("header.read_length_r1 must be positive")
    if header.mates_per_fragment == 1 and header.read_length_r2 != 0:
        raise ProtocolError("SE header.read_length_r2 must be zero")
    if header.mates_per_fragment == 2 and header.read_length_r2 == 0:
        raise ProtocolError("PE header.read_length_r2 must be positive")
    if not isinstance(header.contigs, tuple) or not header.contigs:
        raise ProtocolError("header.contigs must be a non-empty tuple")
    if len(header.contigs) > 0xFFFFFFFF:
        raise ProtocolError("header.contigs has too many entries")
    names = set()
    for index, contig in enumerate(header.contigs):
        if not isinstance(contig, Contig):
            raise ProtocolError(
                "header.contigs[{}] must be a Contig".format(index)
            )
        _validated_string("contig.name", contig.name, nonempty=True)
        if contig.name in names:
            raise ProtocolError("contig names must be unique")
        names.add(contig.name)
        _require_u32("contig.length", contig.length)
        if contig.length == 0:
            raise ProtocolError("contig.length must be positive")
        _require_digest("contig.reference_sha256", contig.reference_sha256)


class _Encoder:
    def __init__(self) -> None:
        self.data = bytearray()

    def u8(self, value: int) -> None:
        self.data.extend(_U8.pack(value))

    def u32(self, value: int) -> None:
        self.data.extend(_U32.pack(value))

    def u64(self, value: int) -> None:
        self.data.extend(_U64.pack(value))

    def raw(self, value: bytes) -> None:
        self.data.extend(value)

    def string(self, value: str) -> None:
        encoded = value.encode("utf-8")
        self.u32(len(encoded))
        self.raw(encoded)

    def array(self, values: Sequence, item: struct.Struct) -> None:
        for value in values:
            self.data.extend(item.pack(value))

    def align4(self) -> None:
        self.data.extend(b"\x00" * ((-len(self.data)) % 4))


class _Decoder:
    def __init__(self, payload: bytes) -> None:
        if not isinstance(payload, bytes):
            raise ProtocolError("protocol payload must be immutable bytes")
        self.payload = payload
        self.backing = memoryview(payload)
        self.cursor = 0

    def _take(self, size: int, name: str) -> memoryview:
        end = self.cursor + size
        if size < 0 or end > len(self.payload):
            raise ProtocolError("{} is truncated".format(name))
        value = self.backing[self.cursor:end]
        self.cursor = end
        return value

    def u8(self, name: str) -> int:
        return _U8.unpack(self._take(_U8.size, name))[0]

    def u32(self, name: str) -> int:
        return _U32.unpack(self._take(_U32.size, name))[0]

    def u64(self, name: str) -> int:
        return _U64.unpack(self._take(_U64.size, name))[0]

    def raw(self, size: int, name: str) -> bytes:
        return bytes(self._take(size, name))

    def string(self, name: str) -> str:
        length = self.u32(name + " length")
        if length > MAX_STRING_BYTES:
            raise ProtocolError("{} exceeds the 1 MiB string limit".format(name))
        raw = self.raw(length, name)
        if b"\x00" in raw:
            raise ProtocolError("{} must not contain NUL".format(name))
        try:
            return raw.decode("utf-8")
        except UnicodeError as error:
            raise ProtocolError("{} is not valid UTF-8".format(name)) from error

    def view(self, count: int, item: struct.Struct, name: str) -> _ArrayView:
        size = count * item.size
        start = self.cursor
        self._take(size, name)
        return _ArrayView(self.backing, start, count, item)

    def align4(self, name: str) -> None:
        size = (-self.cursor) % 4
        padding = self.raw(size, name)
        if any(padding):
            raise ProtocolError("{} must be zero".format(name))

    def finish_padding(self, name: str) -> None:
        remaining = len(self.payload) - self.cursor
        expected = (-self.cursor) % 4
        if remaining != expected or len(self.payload) % 4 != 0:
            raise ProtocolError("{} has non-canonical padding length".format(name))
        padding = self.raw(remaining, name + " padding")
        if any(padding):
            raise ProtocolError("{} padding must be zero".format(name))


def _require_length(name: str, values: Sequence, expected: int) -> None:
    try:
        observed = len(values)
    except TypeError as error:
        raise ProtocolError("{} must be a sequence".format(name)) from error
    if observed != expected:
        raise ProtocolError(
            "{} length {} does not equal {}".format(name, observed, expected)
        )


def _validate_u32_column(name: str, values: Sequence[int]) -> None:
    for index, value in enumerate(values):
        _require_u32("{}[{}]".format(name, index), value)


def _validate_u8_column(name: str, values: Sequence[int]) -> None:
    for index, value in enumerate(values):
        _require_u8("{}[{}]".format(name, index), value)


def _validate_prefix(
    name: str,
    values: Sequence[int],
    row_count: int,
    flat_count: int,
) -> None:
    _require_length(name, values, row_count + 1)
    previous = None
    for index, value in enumerate(values):
        _require_u32("{}[{}]".format(name, index), value)
        if previous is not None and value < previous:
            raise ProtocolError("{} must be monotone".format(name))
        previous = value
    if values[0] != 0 or values[-1] != flat_count:
        raise ProtocolError(
            "{} must begin at zero and end at {}".format(name, flat_count)
        )


def _validate_probability(name: str, value: float) -> None:
    if isinstance(value, bool):
        raise ProtocolError("{} must be a finite probability".format(name))
    try:
        converted = float(value)
    except (TypeError, ValueError, OverflowError) as error:
        raise ProtocolError("{} must be a finite probability".format(name)) from error
    if not math.isfinite(converted) or converted < 0.0 or converted > 1.0:
        raise ProtocolError("{} must be finite and in [0,1]".format(name))


def _validate_base_column(
    name: str, values: Sequence[int], *, allow_n: bool
) -> None:
    maximum = 4 if allow_n else 3
    for index, value in enumerate(values):
        _require_u8("{}[{}]".format(name, index), value)
        if value > maximum:
            raise ProtocolError("{} contains an invalid base code".format(name))


def _validate_details(
    batch: FragmentBatch,
    header: Header,
    details: FragmentDetails,
) -> None:
    if not isinstance(details, FragmentDetails):
        raise ProtocolError("batch.details must be FragmentDetails")
    fragment_count = batch.fragment_count
    site_count = batch.methylation_site_count
    projection_count = len(details.projection_template_starts)
    variant_count = len(details.variant_indices)
    variant_id_byte_count = len(details.variant_ids)
    ref_base_count = len(details.variant_ref_bases)
    alt_base_count = len(details.variant_alt_bases)
    original_n_count = len(details.original_n_template_offsets)
    for name, value in (
        ("projection_count", projection_count),
        ("variant_count", variant_count),
        ("variant_id_byte_count", variant_id_byte_count),
        ("ref_base_count", ref_base_count),
        ("alt_base_count", alt_base_count),
        ("original_n_count", original_n_count),
    ):
        _require_u32(name, value)

    _validate_prefix(
        "details.projection_offsets",
        details.projection_offsets,
        fragment_count,
        projection_count,
    )
    _validate_prefix(
        "details.variant_offsets", details.variant_offsets, fragment_count, variant_count
    )
    _validate_prefix(
        "details.original_n_offsets",
        details.original_n_offsets,
        fragment_count,
        original_n_count,
    )
    _validate_prefix(
        "details.variant_id_offsets",
        details.variant_id_offsets,
        variant_count,
        variant_id_byte_count,
    )
    _validate_prefix(
        "details.variant_ref_offsets",
        details.variant_ref_offsets,
        variant_count,
        ref_base_count,
    )
    _validate_prefix(
        "details.variant_alt_offsets",
        details.variant_alt_offsets,
        variant_count,
        alt_base_count,
    )
    for name, values, expected in (
        ("projection_template_ends", details.projection_template_ends, projection_count),
        (
            "projection_reference_starts",
            details.projection_reference_starts,
            projection_count,
        ),
        ("variant_reference_starts", details.variant_reference_starts, variant_count),
        ("variant_reference_ends", details.variant_reference_ends, variant_count),
        ("variant_template_starts", details.variant_template_starts, variant_count),
        ("variant_template_ends", details.variant_template_ends, variant_count),
        ("variant_sources", details.variant_sources, variant_count),
        ("variant_kinds", details.variant_kinds, variant_count),
        (
            "variant_phased_haplotypes",
            details.variant_phased_haplotypes,
            variant_count,
        ),
        ("site_reference_positions", details.site_reference_positions, site_count),
    ):
        _require_length("details." + name, values, expected)
    for name, values in (
        ("projection_template_starts", details.projection_template_starts),
        ("projection_template_ends", details.projection_template_ends),
        ("projection_reference_starts", details.projection_reference_starts),
        ("variant_indices", details.variant_indices),
        ("variant_id_offsets", details.variant_id_offsets),
        ("variant_reference_starts", details.variant_reference_starts),
        ("variant_reference_ends", details.variant_reference_ends),
        ("variant_template_starts", details.variant_template_starts),
        ("variant_template_ends", details.variant_template_ends),
        ("site_reference_positions", details.site_reference_positions),
        ("original_n_template_offsets", details.original_n_template_offsets),
    ):
        _validate_u32_column("details." + name, values)
    _validate_u8_column("details.variant_kinds", details.variant_kinds)
    _validate_u8_column("details.variant_sources", details.variant_sources)
    _validate_u8_column(
        "details.variant_phased_haplotypes", details.variant_phased_haplotypes
    )
    _validate_base_column(
        "details.variant_ref_bases", details.variant_ref_bases, allow_n=False
    )
    _validate_base_column(
        "details.variant_alt_bases", details.variant_alt_bases, allow_n=False
    )

    for row in range(fragment_count):
        template_flat_begin = batch.template_offsets[row]
        template_length = batch.template_offsets[row + 1] - template_flat_begin
        reference_begin = batch.reference_starts[row]
        reference_end = batch.reference_ends[row]
        haplotype = batch.haplotypes[row]
        projection_cover = bytearray(template_length)
        insertion_cover = bytearray(template_length)
        event_cover = bytearray(template_length)
        mapped_positions = [NO_REFERENCE_POSITION] * template_length

        projection_begin = details.projection_offsets[row]
        projection_end = details.projection_offsets[row + 1]
        previous_template_end = None
        previous_reference_begin = None
        previous_reference_end = None
        for index in range(projection_begin, projection_end):
            template_begin = details.projection_template_starts[index]
            template_end = details.projection_template_ends[index]
            mapped_begin = details.projection_reference_starts[index]
            if template_begin >= template_end or template_end > template_length:
                raise ProtocolError("projection run has an invalid template interval")
            mapped_end = mapped_begin + (template_end - template_begin)
            if mapped_end > 0xFFFFFFFF:
                raise ProtocolError("projection run reference interval overflows u32")
            if mapped_begin < reference_begin or mapped_end > reference_end:
                raise ProtocolError("projection run exceeds its reference envelope")
            if previous_template_end is not None:
                if template_begin < previous_template_end:
                    raise ProtocolError("projection runs overlap or are unordered")
                if mapped_begin < previous_reference_end:
                    raise ProtocolError("projection references are not increasing")
                if (
                    template_begin == previous_template_end
                    and mapped_begin == previous_reference_end
                ):
                    raise ProtocolError("projection runs are not maximal")
                if mapped_begin <= previous_reference_begin:
                    raise ProtocolError("projection references are not strictly ordered")
            for local_offset in range(template_begin, template_end):
                if projection_cover[local_offset]:
                    raise ProtocolError("projection runs cover a template base twice")
                projection_cover[local_offset] = 1
                mapped_positions[local_offset] = mapped_begin + (
                    local_offset - template_begin
                )
            previous_template_end = template_end
            previous_reference_begin = mapped_begin
            previous_reference_end = mapped_end

        event_begin = details.variant_offsets[row]
        event_end = details.variant_offsets[row + 1]
        previous_variant_index = None
        for index in range(event_begin, event_end):
            variant_index = details.variant_indices[index]
            if variant_index == NO_REFERENCE_POSITION:
                raise ProtocolError("event id uses the reserved sentinel")
            if previous_variant_index is not None and variant_index <= previous_variant_index:
                raise ProtocolError("event ids must be strictly increasing")
            previous_variant_index = variant_index
            try:
                VariantSourceCode(details.variant_sources[index])
            except ValueError as error:
                raise ProtocolError("variant source is invalid") from error
            id_start = details.variant_id_offsets[index]
            id_end = details.variant_id_offsets[index + 1]
            if id_start == id_end:
                raise ProtocolError("variant ID must not be empty")
            try:
                bytes(details.variant_ids[id_start:id_end]).decode("utf-8")
            except UnicodeDecodeError as error:
                raise ProtocolError("variant ID must be valid UTF-8") from error
            kind_value = details.variant_kinds[index]
            _require_enum("event.kind", kind_value, VariantKind)
            kind = VariantKind(kind_value)
            phased = details.variant_phased_haplotypes[index]
            if phased not in (255, haplotype):
                raise ProtocolError("event phased haplotype disagrees with fragment")
            event_reference_begin = details.variant_reference_starts[index]
            event_reference_end = details.variant_reference_ends[index]
            event_template_begin = details.variant_template_starts[index]
            event_template_end = details.variant_template_ends[index]
            if (
                event_reference_begin > event_reference_end
                or event_reference_begin < reference_begin
                or event_reference_end > reference_end
            ):
                raise ProtocolError("event reference interval exceeds its fragment")
            if (
                event_template_begin > event_template_end
                or event_template_end > template_length
            ):
                raise ProtocolError("event template interval exceeds its fragment")
            ref_begin = details.variant_ref_offsets[index]
            ref_end = details.variant_ref_offsets[index + 1]
            alt_begin = details.variant_alt_offsets[index]
            alt_end = details.variant_alt_offsets[index + 1]
            reference_span = event_reference_end - event_reference_begin
            template_span = event_template_end - event_template_begin
            if ref_end - ref_begin != reference_span:
                raise ProtocolError("event REF bases disagree with its reference span")
            if alt_end - alt_begin != template_span:
                raise ProtocolError("event ALT bases disagree with its template span")
            if kind is VariantKind.SNV:
                if reference_span == 0 or reference_span != template_span:
                    raise ProtocolError("SNV event has an invalid shape")
            elif kind is VariantKind.INSERTION:
                if reference_span != 0 or template_span == 0:
                    raise ProtocolError("insertion event has an invalid shape")
            elif reference_span == 0 or template_span != 0:
                raise ProtocolError("deletion event has an invalid shape")

            if kind is VariantKind.INSERTION:
                previous_mapped = next(
                    (
                        mapped_positions[offset]
                        for offset in range(event_template_begin - 1, -1, -1)
                        if mapped_positions[offset] != NO_REFERENCE_POSITION
                    ),
                    None,
                )
                next_mapped = next(
                    (
                        mapped_positions[offset]
                        for offset in range(event_template_end, template_length)
                        if mapped_positions[offset] != NO_REFERENCE_POSITION
                    ),
                    None,
                )
                if (
                    (previous_mapped is None and event_reference_begin != reference_begin)
                    or (
                        previous_mapped is not None
                        and previous_mapped >= event_reference_begin
                    )
                    or (next_mapped is None and event_reference_begin != reference_end)
                    or (
                        next_mapped is not None
                        and next_mapped < event_reference_begin
                    )
                ):
                    raise ProtocolError("insertion anchor disagrees with projection")
            elif kind is VariantKind.DELETION:
                previous_mapped = next(
                    (
                        mapped_positions[offset]
                        for offset in range(event_template_begin - 1, -1, -1)
                        if mapped_positions[offset] != NO_REFERENCE_POSITION
                    ),
                    None,
                )
                next_mapped = next(
                    (
                        mapped_positions[offset]
                        for offset in range(event_template_begin, template_length)
                        if mapped_positions[offset] != NO_REFERENCE_POSITION
                    ),
                    None,
                )
                if (
                    (previous_mapped is not None and previous_mapped >= event_reference_begin)
                    or (next_mapped is not None and next_mapped < event_reference_end)
                    or any(
                        event_reference_begin <= position < event_reference_end
                        for position in mapped_positions
                        if position != NO_REFERENCE_POSITION
                    )
                ):
                    raise ProtocolError("deletion boundary disagrees with projection")

            for relative in range(template_span):
                local_offset = event_template_begin + relative
                if event_cover[local_offset]:
                    raise ProtocolError("event template spans overlap")
                event_cover[local_offset] = 1
                if (
                    batch.template_bases[template_flat_begin + local_offset]
                    != details.variant_alt_bases[alt_begin + relative]
                ):
                    raise ProtocolError("event ALT bases disagree with the template")
                if kind is VariantKind.SNV:
                    if mapped_positions[local_offset] != event_reference_begin + relative:
                        raise ProtocolError("SNV projection disagrees with its event")
                elif kind is VariantKind.INSERTION:
                    if projection_cover[local_offset]:
                        raise ProtocolError("insertion overlaps a projection run")
                    insertion_cover[local_offset] = 1

        for local_offset in range(template_length):
            if projection_cover[local_offset] + insertion_cover[local_offset] != 1:
                raise ProtocolError(
                    "projection and insertion runs do not cover the template exactly"
                )

        site_begin = batch.site_offsets[row]
        site_end = batch.site_offsets[row + 1]
        for index in range(site_begin, site_end):
            local_offset = batch.site_template_offsets[index]
            if details.site_reference_positions[index] != mapped_positions[local_offset]:
                raise ProtocolError("site reference position disagrees with projection")

        original_n_begin = details.original_n_offsets[row]
        original_n_end = details.original_n_offsets[row + 1]
        observed_n = bytearray(template_length)
        previous_n = None
        for index in range(original_n_begin, original_n_end):
            local_offset = details.original_n_template_offsets[index]
            if local_offset >= template_length:
                raise ProtocolError("original-N offset exceeds its template")
            if previous_n is not None and local_offset <= previous_n:
                raise ProtocolError("original-N offsets must be strictly increasing")
            previous_n = local_offset
            observed_n[local_offset] = 1
            if batch.template_bases[template_flat_begin + local_offset] != 4:
                raise ProtocolError("PRESERVE_N provenance does not point to N")
        for local_offset in range(template_length):
            is_n = batch.template_bases[template_flat_begin + local_offset] == 4
            if bool(observed_n[local_offset]) != is_n:
                raise ProtocolError("PRESERVE_N provenance is incomplete")


def _validate_batch_python(
    batch: FragmentBatch,
    header: Header,
    *,
    expected_first_ordinal: int | None = None,
) -> None:
    if not isinstance(batch, FragmentBatch):
        raise ProtocolError("batch must be a FragmentBatch")
    _validate_header(header)
    fragment_count = batch.fragment_count
    template_count = batch.template_base_count
    mate_count = batch.mate_count
    site_count = batch.methylation_site_count
    for name, value in (
        ("batch.fragment_count", fragment_count),
        ("batch.template_base_count", template_count),
        ("batch.mate_count", mate_count),
        ("batch.methylation_site_count", site_count),
    ):
        _require_u32(name, value)
    if fragment_count == 0:
        raise ProtocolError("fragment batch must contain at least one row")
    _require_u32("batch.first_fragment_ordinal", batch.first_fragment_ordinal)
    if batch.first_fragment_ordinal + fragment_count > 0xFFFFFFFF:
        raise ProtocolError("batch fragment ordinal range exceeds u32")
    if (
        expected_first_ordinal is not None
        and batch.first_fragment_ordinal != expected_first_ordinal
    ):
        raise ProtocolError("fragment batch ordinal range is not consecutive")

    for name, values, expected in (
        ("reference_starts", batch.reference_starts, fragment_count),
        ("reference_ends", batch.reference_ends, fragment_count),
        ("haplotypes", batch.haplotypes, fragment_count),
        ("capture_strands", batch.capture_strands, fragment_count),
        ("mate_template_starts", batch.mate_template_starts, mate_count),
        ("mate_template_ends", batch.mate_template_ends, mate_count),
        ("mate_reverse_complements", batch.mate_reverse_complements, mate_count),
        ("site_probabilities", batch.site_probabilities, site_count),
        ("site_contexts", batch.site_contexts, site_count),
        ("methylation_sources", batch.methylation_sources, site_count),
        ("site_alleles", batch.site_alleles, site_count),
    ):
        _require_length("batch." + name, values, expected)
    _validate_prefix(
        "batch.template_offsets",
        batch.template_offsets,
        fragment_count,
        template_count,
    )
    _validate_prefix(
        "batch.mate_offsets", batch.mate_offsets, fragment_count, mate_count
    )
    _validate_prefix(
        "batch.site_offsets", batch.site_offsets, fragment_count, site_count
    )
    for name, values in (
        ("contig_indices", batch.contig_indices),
        ("reference_starts", batch.reference_starts),
        ("reference_ends", batch.reference_ends),
        ("mate_template_starts", batch.mate_template_starts),
        ("mate_template_ends", batch.mate_template_ends),
        ("site_template_offsets", batch.site_template_offsets),
    ):
        _validate_u32_column("batch." + name, values)
    for name, values in (
        ("haplotypes", batch.haplotypes),
        ("capture_strands", batch.capture_strands),
        ("mate_indices", batch.mate_indices),
        ("mate_reverse_complements", batch.mate_reverse_complements),
        ("site_contexts", batch.site_contexts),
        ("methylation_sources", batch.methylation_sources),
        ("site_alleles", batch.site_alleles),
    ):
        _validate_u8_column("batch." + name, values)
    _validate_base_column("batch.template_bases", batch.template_bases, allow_n=True)

    cytosine_contexts = {
        MethylationContext.CG_C,
        MethylationContext.CHG_C,
        MethylationContext.CHH_C,
    }
    guanine_contexts = {
        MethylationContext.CG_G,
        MethylationContext.CHG_G,
        MethylationContext.CHH_G,
    }
    for row in range(fragment_count):
        contig_index = batch.contig_indices[row]
        if contig_index >= len(header.contigs):
            raise ProtocolError("fragment contig index is outside the header")
        reference_begin = batch.reference_starts[row]
        reference_end = batch.reference_ends[row]
        if (
            reference_begin >= reference_end
            or reference_end > header.contigs[contig_index].length
        ):
            raise ProtocolError("fragment reference envelope is invalid")
        if batch.haplotypes[row] not in (0, 1):
            raise ProtocolError("fragment haplotype must be zero or one")
        _require_enum(
            "fragment.capture_strand", batch.capture_strands[row], CaptureStrand
        )
        template_begin = batch.template_offsets[row]
        template_end = batch.template_offsets[row + 1]
        template_length = template_end - template_begin
        if template_length == 0:
            raise ProtocolError("fragment template must not be empty")

        mate_begin = batch.mate_offsets[row]
        mate_end = batch.mate_offsets[row + 1]
        if mate_end - mate_begin != header.mates_per_fragment:
            raise ProtocolError("fragment mate count disagrees with header")
        for local_mate, index in enumerate(range(mate_begin, mate_end)):
            if batch.mate_indices[index] != local_mate:
                raise ProtocolError("mate rows must be ordered as 0 or 0,1")
            if batch.mate_reverse_complements[index] not in (0, 1):
                raise ProtocolError("mate reverse-complement flag must be boolean")
            begin = batch.mate_template_starts[index]
            end = batch.mate_template_ends[index]
            expected_length = (
                header.read_length_r1 if local_mate == 0 else header.read_length_r2
            )
            if begin >= end or end > template_length or end - begin != expected_length:
                raise ProtocolError("mate template slice disagrees with header")

        site_begin = batch.site_offsets[row]
        site_end = batch.site_offsets[row + 1]
        previous_site = None
        for index in range(site_begin, site_end):
            local_offset = batch.site_template_offsets[index]
            if local_offset >= template_length:
                raise ProtocolError("site offset exceeds its template")
            if previous_site is not None and local_offset <= previous_site:
                raise ProtocolError("site offsets must be strictly increasing")
            previous_site = local_offset
            _validate_probability(
                "site.probability", batch.site_probabilities[index]
            )
            _require_enum(
                "site.context", batch.site_contexts[index], MethylationContext
            )
            _require_enum("site.source", batch.methylation_sources[index], MethylationSource)
            _require_enum("site.allele", batch.site_alleles[index], MethylationAllele)
            context = MethylationContext(batch.site_contexts[index])
            base = batch.template_bases[template_begin + local_offset]
            if (context in cytosine_contexts and base != 1) or (
                context in guanine_contexts and base != 2
            ):
                raise ProtocolError("site context is incompatible with template base")

    if header.has_details:
        if batch.details is None:
            raise ProtocolError("header requires Full-Details batch columns")
        _validate_details(batch, header, batch.details)
    elif batch.details is not None:
        raise ProtocolError("header forbids Full-Details batch columns")


_cext_header_cache = None


def _validated_cext_contig_lengths(header: Header) -> tuple[int, ...]:
    global _cext_header_cache
    if _cext_header_cache is not None and _cext_header_cache[0] is header:
        return _cext_header_cache[1]
    values = tuple(contig.length for contig in header.contigs)
    _cext_header_cache = (header, values)
    return values


def _validate_batch(
    batch: FragmentBatch,
    header: Header,
    *,
    expected_first_ordinal: int | None = None,
) -> None:
    """Validate a batch without weakening the reference protocol contract.

    Regular Python-owned batches use the reference validator.  Immutable
    decoded views may use the C-extension validator, which independently checks the
    same column lengths, prefixes, enums, ranges, projections, variant events,
    methylation mappings, and preserved-N provenance directly on their
    little-endian buffers.
    """

    if not isinstance(batch, DecodedBatchView):
        _validate_batch_python(
            batch,
            header,
            expected_first_ordinal=expected_first_ordinal,
        )
        return

    if not isinstance(batch, FragmentBatch):
        raise ProtocolError("batch must be a FragmentBatch")
    _validate_header(header)
    if header.has_details:
        if batch.details is None:
            raise ProtocolError("header requires Full-Details batch columns")
    elif batch.details is not None:
        raise ProtocolError("header forbids Full-Details batch columns")

    common_columns = (
        batch.contig_indices.raw,
        batch.reference_starts.raw,
        batch.reference_ends.raw,
        batch.template_offsets.raw,
        batch.mate_offsets.raw,
        batch.site_offsets.raw,
        batch.mate_template_starts.raw,
        batch.mate_template_ends.raw,
        batch.site_template_offsets.raw,
        batch.site_probabilities.raw,
        batch.haplotypes.raw,
        batch.capture_strands.raw,
        batch.mate_indices.raw,
        batch.mate_reverse_complements.raw,
        batch.site_contexts.raw,
        batch.methylation_sources.raw,
        batch.site_alleles.raw,
        batch.template_bases.raw,
    )
    has_details = None
    if batch.details is not None:
        details = batch.details
        has_details = (
            details.projection_offsets.raw,
            details.variant_offsets.raw,
            details.original_n_offsets.raw,
            details.projection_template_starts.raw,
            details.projection_template_ends.raw,
            details.projection_reference_starts.raw,
            details.variant_indices.raw,
            details.variant_id_offsets.raw,
            details.variant_reference_starts.raw,
            details.variant_reference_ends.raw,
            details.variant_template_starts.raw,
            details.variant_template_ends.raw,
            details.variant_ref_offsets.raw,
            details.variant_alt_offsets.raw,
            details.site_reference_positions.raw,
            details.original_n_template_offsets.raw,
            details.variant_sources.raw,
            details.variant_kinds.raw,
            details.variant_phased_haplotypes.raw,
            details.variant_ids.raw,
            details.variant_ref_bases.raw,
            details.variant_alt_bases.raw,
        )
    try:
        _cext_validate_batch_columns(
            common_columns,
            has_details,
            batch.first_fragment_ordinal,
            _validated_cext_contig_lengths(header),
            header.mates_per_fragment,
            header.read_length_r1,
            header.read_length_r2,
            expected_first_ordinal,
            int(header.has_details),
        )
    except (BufferError, OverflowError, TypeError, ValueError) as error:
        raise ProtocolError(str(error)) from error




def _decode_header_payload(payload: bytes) -> Header:
    decoder = _Decoder(payload)
    run_id = decoder.string("header.run_id")
    core_version = decoder.string("header.core_version")
    rng_contract = decoder.string("header.rng_contract")
    master_seed = decoder.u64("header.master_seed")
    normalized_config_sha256 = decoder.raw(32, "header.normalized_config_sha256")
    technology = decoder.u8("header.technology")
    has_details = decoder.u8("header.has_details")
    mates_per_fragment = decoder.u8("header.mates_per_fragment")
    base_encoding = decoder.u8("header.base_encoding")
    ambiguity_policy = decoder.u8("header.ambiguity_policy")
    if any(decoder.raw(3, "header.reserved")):
        raise ProtocolError("header reserved bytes must be zero")
    read_length_r1 = decoder.u32("header.read_length_r1")
    read_length_r2 = decoder.u32("header.read_length_r2")
    contig_count = decoder.u32("header.contig_count")
    contigs = []
    for index in range(contig_count):
        contigs.append(
            Contig(
                name=decoder.string("contig[{}].name".format(index)),
                length=decoder.u32("contig[{}].length".format(index)),
                reference_sha256=decoder.raw(
                    32, "contig[{}].reference_sha256".format(index)
                ),
            )
        )
    decoder.finish_padding("header payload")
    header = Header(
        run_id=run_id,
        core_version=core_version,
        rng_contract=rng_contract,
        master_seed=master_seed,
        normalized_config_sha256=normalized_config_sha256,
        technology=(
            Technology(technology)
            if technology in tuple(item.value for item in Technology)
            else technology
        ),
        has_details=bool(has_details) if has_details in (0, 1) else has_details,
        mates_per_fragment=mates_per_fragment,
        base_encoding=(
            BaseEncoding(base_encoding) if base_encoding == 1 else base_encoding
        ),
        ambiguity_policy=(
            AmbiguityPolicy(ambiguity_policy)
            if ambiguity_policy == 0
            else ambiguity_policy
        ),
        read_length_r1=read_length_r1,
        read_length_r2=read_length_r2,
        contigs=tuple(contigs),
    )
    _validate_header(header)
    return header




def _decode_batch_payload(
    payload: bytes,
    frame_flags: int,
    header: Header,
    *,
    expected_first_ordinal: int | None = None,
) -> DecodedBatchView:
    if frame_flags & ~DETAILS_PRESENT:
        raise ProtocolError("fragment-batch frame contains unknown flags")
    frame_has_details = bool(frame_flags & DETAILS_PRESENT)
    if frame_has_details != header.has_details:
        raise ProtocolError("fragment-batch details flag disagrees with header")
    decoder = _Decoder(payload)
    first_fragment_ordinal = decoder.u32("batch.first_fragment_ordinal")
    fragment_count = decoder.u32("batch.fragment_count")
    template_count = decoder.u32("batch.template_base_count")
    mate_count = decoder.u32("batch.mate_count")
    site_count = decoder.u32("batch.methylation_site_count")
    contig_indices = decoder.view(fragment_count, _U32, "batch.contig_indices")
    reference_starts = decoder.view(fragment_count, _U32, "batch.reference_starts")
    reference_ends = decoder.view(fragment_count, _U32, "batch.reference_ends")
    template_offsets = decoder.view(fragment_count + 1, _U32, "batch.template_offsets")
    mate_offsets = decoder.view(fragment_count + 1, _U32, "batch.mate_offsets")
    site_offsets = decoder.view(fragment_count + 1, _U32, "batch.site_offsets")
    mate_template_starts = decoder.view(
        mate_count, _U32, "batch.mate_template_starts"
    )
    mate_template_ends = decoder.view(
        mate_count, _U32, "batch.mate_template_ends"
    )
    site_template_offsets = decoder.view(
        site_count, _U32, "batch.site_template_offsets"
    )
    site_probabilities = decoder.view(site_count, _F32, "batch.site_probabilities")
    haplotypes = decoder.view(fragment_count, _U8, "batch.haplotypes")
    capture_strands = decoder.view(fragment_count, _U8, "batch.capture_strands")
    mate_indices = decoder.view(mate_count, _U8, "batch.mate_indices")
    mate_reverse_complements = decoder.view(
        mate_count, _U8, "batch.mate_reverse_complements"
    )
    site_contexts = decoder.view(site_count, _U8, "batch.site_contexts")
    methylation_sources = decoder.view(site_count, _U8, "batch.methylation_sources")
    site_alleles = decoder.view(site_count, _U8, "batch.site_alleles")
    template_bases = decoder.view(template_count, _U8, "batch.template_bases")

    details = None
    if frame_has_details:
        decoder.align4("batch common-column alignment")
        projection_count = decoder.u32("details.projection_run_count")
        variant_count = decoder.u32("details.variant_count")
        variant_id_byte_count = decoder.u32("details.variant_id_byte_count")
        ref_base_count = decoder.u32("details.variant_ref_base_count")
        alt_base_count = decoder.u32("details.variant_alt_base_count")
        original_n_count = decoder.u32("details.original_n_count")
        details = FragmentDetails(
            projection_offsets=decoder.view(
                fragment_count + 1, _U32, "details.projection_offsets"
            ),
            variant_offsets=decoder.view(
                fragment_count + 1, _U32, "details.variant_offsets"
            ),
            original_n_offsets=decoder.view(
                fragment_count + 1, _U32, "details.original_n_offsets"
            ),
            projection_template_starts=decoder.view(
                projection_count, _U32, "details.projection_template_starts"
            ),
            projection_template_ends=decoder.view(
                projection_count, _U32, "details.projection_template_ends"
            ),
            projection_reference_starts=decoder.view(
                projection_count, _U32, "details.projection_reference_starts"
            ),
            variant_indices=decoder.view(variant_count, _U32, "details.variant_indices"),
            variant_id_offsets=decoder.view(
                variant_count + 1, _U32, "details.variant_id_offsets"
            ),
            variant_reference_starts=decoder.view(
                variant_count, _U32, "details.variant_reference_starts"
            ),
            variant_reference_ends=decoder.view(
                variant_count, _U32, "details.variant_reference_ends"
            ),
            variant_template_starts=decoder.view(
                variant_count, _U32, "details.variant_template_starts"
            ),
            variant_template_ends=decoder.view(
                variant_count, _U32, "details.variant_template_ends"
            ),
            variant_ref_offsets=decoder.view(
                variant_count + 1, _U32, "details.variant_ref_offsets"
            ),
            variant_alt_offsets=decoder.view(
                variant_count + 1, _U32, "details.variant_alt_offsets"
            ),
            site_reference_positions=decoder.view(
                site_count, _U32, "details.site_reference_positions"
            ),
            original_n_template_offsets=decoder.view(
                original_n_count, _U32, "details.original_n_template_offsets"
            ),
            variant_sources=decoder.view(
                variant_count, _U8, "details.variant_sources"
            ),
            variant_kinds=decoder.view(variant_count, _U8, "details.variant_kinds"),
            variant_phased_haplotypes=decoder.view(
                variant_count, _U8, "details.variant_phased_haplotypes"
            ),
            variant_ids=decoder.view(
                variant_id_byte_count, _U8, "details.variant_ids"
            ),
            variant_ref_bases=decoder.view(
                ref_base_count, _U8, "details.variant_ref_bases"
            ),
            variant_alt_bases=decoder.view(
                alt_base_count, _U8, "details.variant_alt_bases"
            ),
        )
    decoder.finish_padding("fragment-batch payload")
    batch = DecodedBatchView(
        first_fragment_ordinal=first_fragment_ordinal,
        contig_indices=contig_indices,
        reference_starts=reference_starts,
        reference_ends=reference_ends,
        template_offsets=template_offsets,
        mate_offsets=mate_offsets,
        site_offsets=site_offsets,
        mate_template_starts=mate_template_starts,
        mate_template_ends=mate_template_ends,
        site_template_offsets=site_template_offsets,
        site_probabilities=site_probabilities,
        haplotypes=haplotypes,
        capture_strands=capture_strands,
        mate_indices=mate_indices,
        mate_reverse_complements=mate_reverse_complements,
        site_contexts=site_contexts,
        methylation_sources=methylation_sources,
        site_alleles=site_alleles,
        template_bases=template_bases,
        details=details,
        _payload=payload,
    )
    _validate_batch(
        batch, header, expected_first_ordinal=expected_first_ordinal
    )
    return batch


def decode_batch_payload(
    payload: bytes,
    header: Header,
    *,
    expected_first_ordinal: int | None = None,
) -> DecodedBatchView:
    """Decode one already authenticated batch payload.

    This entry point is for a supervisor that has already verified the frame
    envelope and CRC32C before copying the immutable payload across a process
    boundary.  The payload is still decoded and semantically validated in the
    worker.  Details-column presence is derived from the accepted header, so a
    caller cannot weaken that contract independently.
    """

    if not isinstance(payload, bytes):
        raise ProtocolError("protocol batch payload must be immutable bytes")
    frame_flags = (
        DETAILS_PRESENT
        if header.has_details
        else 0
    )
    return _decode_batch_payload(
        payload,
        frame_flags,
        header,
        expected_first_ordinal=expected_first_ordinal,
    )


def _validate_trailer(trailer: Trailer, contig_count: int) -> None:
    if not isinstance(trailer, Trailer):
        raise ProtocolError("trailer must be a protocol Trailer")
    for name, value in (
        ("fragment_count", trailer.fragment_count),
        ("fragment_batch_count", trailer.fragment_batch_count),
        ("mate_count", trailer.mate_count),
        ("template_base_count", trailer.template_base_count),
        ("methylation_site_count", trailer.methylation_site_count),
        ("skipped_fragment_count", trailer.skipped_fragment_count),
    ):
        _require_u64("trailer." + name, value)
    if not isinstance(trailer.per_contig_fragment_counts, tuple):
        raise ProtocolError("trailer per-contig counts must be a tuple")
    if len(trailer.per_contig_fragment_counts) != contig_count:
        raise ProtocolError("trailer per-contig count cardinality is invalid")
    for index, value in enumerate(trailer.per_contig_fragment_counts):
        _require_u64("trailer.per_contig[{}]".format(index), value)
    if sum(trailer.per_contig_fragment_counts) != trailer.fragment_count:
        raise ProtocolError("trailer per-contig counts do not sum to fragments")
    _require_digest("trailer.stream_sha256", trailer.stream_sha256)




def _decode_trailer_payload(payload: bytes, contig_count: int) -> Trailer:
    decoder = _Decoder(payload)
    fragment_count = decoder.u64("trailer.fragment_count")
    fragment_batch_count = decoder.u64("trailer.fragment_batch_count")
    mate_count = decoder.u64("trailer.mate_count")
    template_base_count = decoder.u64("trailer.template_base_count")
    methylation_site_count = decoder.u64("trailer.methylation_site_count")
    skipped_fragment_count = decoder.u64("trailer.skipped_fragment_count")
    per_contig_count = decoder.u32("trailer.per_contig_fragment_count")
    per_contig = tuple(
        decoder.u64("trailer.per_contig[{}]".format(index))
        for index in range(per_contig_count)
    )
    stream_sha256 = decoder.raw(32, "trailer.stream_sha256")
    decoder.finish_padding("trailer payload")
    trailer = Trailer(
        fragment_count=fragment_count,
        fragment_batch_count=fragment_batch_count,
        mate_count=mate_count,
        template_base_count=template_base_count,
        methylation_site_count=methylation_site_count,
        skipped_fragment_count=skipped_fragment_count,
        per_contig_fragment_counts=per_contig,
        stream_sha256=stream_sha256,
    )
    _validate_trailer(trailer, contig_count)
    return trailer


def _encode_error_payload(error: ErrorFrame) -> bytes:
    if not isinstance(error, ErrorFrame):
        raise ProtocolError("error must be an ErrorFrame")
    _require_u32("error.error_code", error.error_code)
    _validated_string("error.message", error.message)
    encoder = _Encoder()
    encoder.u32(error.error_code)
    encoder.string(error.message)
    encoder.align4()
    return bytes(encoder.data)


def _decode_error_payload(payload: bytes) -> ErrorFrame:
    decoder = _Decoder(payload)
    error = ErrorFrame(
        error_code=decoder.u32("error.error_code"),
        message=decoder.string("error.message"),
    )
    decoder.finish_padding("error payload")
    _encode_error_payload(error)
    return error


@dataclass(frozen=True)
class _DecodedFrame:
    frame_type: FrameType
    frame_flags: int
    sequence: int
    payload: bytes
    encoded: bytes


def _read_exact(source: BinaryIO, size: int, name: str) -> bytes:
    output = bytearray()
    while len(output) < size:
        chunk = source.read(size - len(output))
        if chunk is None:
            raise ProtocolError("{} read made no progress".format(name))
        if not isinstance(chunk, bytes):
            raise ProtocolError("{} source returned non-bytes data".format(name))
        if not chunk:
            raise ProtocolError("{} is truncated".format(name))
        output.extend(chunk)
    return bytes(output)


def _read_frame(source: BinaryIO, expected_sequence: int) -> _DecodedFrame:
    envelope = _read_exact(source, _FRAME_ENVELOPE.size, "frame envelope")
    payload_length, raw_type, flags, reserved, sequence = _FRAME_ENVELOPE.unpack(
        envelope
    )
    if payload_length > MAX_FRAME_PAYLOAD:
        raise ProtocolError("frame exceeds the 64 MiB payload limit")
    if reserved != 0:
        raise ProtocolError("frame reserved field must be zero")
    if sequence != expected_sequence:
        raise ProtocolError("frame sequences are not consecutive")
    try:
        frame_type = FrameType(raw_type)
    except ValueError as error:
        raise ProtocolError("frame type is unknown") from error
    payload = _read_exact(source, payload_length, "frame payload")
    encoded_crc = _read_exact(source, _CRC.size, "frame CRC32C")
    observed_crc = _CRC.unpack(encoded_crc)[0]
    checksummed = envelope + payload
    if observed_crc != crc32c(checksummed):
        raise ProtocolError("frame CRC32C mismatch")
    return _DecodedFrame(
        frame_type=frame_type,
        frame_flags=flags,
        sequence=sequence,
        payload=payload,
        encoded=checksummed + encoded_crc,
    )


class ProtocolReader:
    """Strict streaming reader that exposes a batch only after frame CRC."""

    def __init__(
        self,
        source: BinaryIO,
        *,
        expected_header: Header | None = None,
        expected_skipped_fragment_count: int | None = None,
    ) -> None:
        if not hasattr(source, "read"):
            raise ValueError("source must be a binary reader")
        if expected_header is not None:
            _validate_header(expected_header)
        if expected_skipped_fragment_count is not None:
            _require_u64(
                "expected_skipped_fragment_count",
                expected_skipped_fragment_count,
            )
        self._source = source
        self._expected_skipped = expected_skipped_fragment_count
        self._sequence = 0
        self._next_ordinal = 0
        self._digest = hashlib.sha256()
        self._fragment_count = 0
        self._batch_count = 0
        self._mate_count = 0
        self._template_count = 0
        self._site_count = 0
        self._iterated = False
        self._terminal = False
        self.trailer = None  # type: Trailer | None

        preamble = _read_exact(source, _PREAMBLE.size, "protocol preamble")
        magic, major, minor, flags = _PREAMBLE.unpack(preamble)
        if magic != MAGIC:
            raise ProtocolError("protocol magic is invalid")
        if (major, minor) != (PROTOCOL_MAJOR, PROTOCOL_MINOR):
            raise ProtocolError("unsupported protocol version {}.{}".format(major, minor))
        if flags != PREAMBLE_FLAGS:
            raise ProtocolError("protocol preamble flags are unsupported")
        self._digest.update(preamble)

        frame = _read_frame(source, self._sequence)
        if frame.frame_type is not FrameType.HEADER or frame.frame_flags != 0:
            raise ProtocolError("protocol stream must begin with an unflagged header")
        self._sequence += 1
        self._digest.update(frame.encoded)
        self.header = _decode_header_payload(frame.payload)
        if expected_header is not None and self.header != expected_header:
            raise ProtocolError("core header disagrees with the Python projection")
        self._per_contig = [0] * len(self.header.contigs)

    def __iter__(self) -> Iterator[DecodedBatchView]:
        if self._iterated:
            raise ProtocolError("protocol reader can be iterated only once")
        self._iterated = True
        while True:
            frame = _read_frame(self._source, self._sequence)
            self._sequence += 1
            if frame.frame_type is FrameType.FRAGMENT_BATCH:
                self._digest.update(frame.encoded)
                batch = _decode_batch_payload(
                    frame.payload,
                    frame.frame_flags,
                    self.header,
                    expected_first_ordinal=self._next_ordinal,
                )
                self._accept_batch(batch)
                yield batch
                continue
            if frame.frame_type is FrameType.ERROR:
                if frame.frame_flags != 0:
                    raise ProtocolError("error frame flags must be zero")
                error = _decode_error_payload(frame.payload)
                self._require_eof("error frame")
                self._terminal = True
                raise CoreReportedError(error.error_code, error.message)
            if frame.frame_type is not FrameType.TRAILER or frame.frame_flags != 0:
                raise ProtocolError("unexpected frame in protocol stream")
            trailer = _decode_trailer_payload(
                frame.payload, len(self.header.contigs)
            )
            self._validate_trailer_counts(trailer)
            self._require_eof("trailer")
            self.trailer = trailer
            self._terminal = True
            return

    def _accept_batch(self, batch: DecodedBatchView) -> None:
        self._next_ordinal += batch.fragment_count
        self._fragment_count += batch.fragment_count
        self._batch_count += 1
        self._mate_count += batch.mate_count
        self._template_count += batch.template_base_count
        self._site_count += batch.methylation_site_count
        for contig_index in batch.contig_indices:
            self._per_contig[contig_index] += 1

    def _validate_trailer_counts(self, trailer: Trailer) -> None:
        observed = (
            self._fragment_count,
            self._batch_count,
            self._mate_count,
            self._template_count,
            self._site_count,
        )
        claimed = (
            trailer.fragment_count,
            trailer.fragment_batch_count,
            trailer.mate_count,
            trailer.template_base_count,
            trailer.methylation_site_count,
        )
        if claimed != observed:
            raise ProtocolError("trailer aggregate counts disagree with batches")
        if trailer.per_contig_fragment_counts != tuple(self._per_contig):
            raise ProtocolError("trailer per-contig counts disagree with batches")
        if trailer.stream_sha256 != self._digest.digest():
            raise ProtocolError("trailer stream SHA-256 mismatch")
        if (
            self._expected_skipped is not None
            and trailer.skipped_fragment_count != self._expected_skipped
        ):
            raise ProtocolError("trailer skipped-fragment count is unexpected")

    def _require_eof(self, name: str) -> None:
        trailing = self._source.read(1)
        if trailing is None:
            raise ProtocolError("{} EOF check made no progress".format(name))
        if not isinstance(trailing, bytes):
            raise ProtocolError("{} EOF check returned non-bytes data".format(name))
        if trailing:
            raise ProtocolError("bytes follow the terminal {}".format(name))

    def read_all(self) -> ProtocolStream:
        batches = tuple(self)
        if self.trailer is None:
            raise ProtocolError("protocol stream ended without a trailer")
        return ProtocolStream(self.header, batches, self.trailer)

    def validate_core_exit_status(self, exit_status: int) -> None:
        if not self._terminal or self.trailer is None:
            raise ProtocolError("protocol trailer must be accepted first")
        if isinstance(exit_status, bool) or not isinstance(exit_status, int):
            raise ProtocolError("core exit status must be an integer")
        if exit_status != 0:
            raise ProtocolError(
                "core exited with non-zero status {} after its trailer".format(
                    exit_status
                )
            )


__all__ = [
    "AmbiguityPolicy",
    "BaseEncoding",
    "CaptureStrand",
    "Contig",
    "CoreReportedError",
    "DecodedBatchView",
    "ErrorFrame",
    "FragmentBatch",
    "Header",
    "MAGIC",
    "MAX_FRAME_PAYLOAD",
    "MAX_STRING_BYTES",
    "MethylationAllele",
    "MethylationContext",
    "MethylationSource",
    "NO_REFERENCE_POSITION",
    "PREAMBLE_FLAGS",
    "PROTOCOL_MAJOR",
    "PROTOCOL_MINOR",
    "PROTOCOL_VERSION",
    "ProtocolError",
    "ProtocolReader",
    "ProtocolStream",
    "RNG_CONTRACT",
    "DETAILS_PRESENT",
    "Technology",
    "Trailer",
    "FragmentDetails",
    "VariantKind",
    "decode_batch_payload",
]
