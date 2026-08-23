"""Test-only reference writer and in-memory protocol conveniences."""

from __future__ import annotations

from collections.abc import Iterable
from typing import BinaryIO
import hashlib
import io

from bsreadsim.htsim import protocol as _protocol
from bsreadsim.htsim.protocol import (
    ErrorFrame,
    FragmentBatch,
    FrameType,
    Header,
    ProtocolError,
    ProtocolReader,
    ProtocolStream,
    Trailer,
)


def _encode_frame(
    frame_type: FrameType,
    frame_flags: int,
    sequence: int,
    payload: bytes,
) -> bytes:
    _protocol._require_enum("frame_type", frame_type, FrameType)
    _protocol._require_u8("frame_flags", frame_flags)
    _protocol._require_u64("frame sequence", sequence)
    if not isinstance(payload, bytes):
        raise ProtocolError("frame payload must be bytes")
    if len(payload) > _protocol.MAX_FRAME_PAYLOAD:
        raise ProtocolError("frame exceeds the 64 MiB payload limit")
    envelope = _protocol._FRAME_ENVELOPE.pack(
        len(payload), int(frame_type), frame_flags, 0, sequence
    )
    checksummed = envelope + payload
    return checksummed + _protocol._CRC.pack(_protocol.crc32c(checksummed))


def _write_all(sink: BinaryIO, data: bytes) -> None:
    view = memoryview(data)
    cursor = 0
    while cursor < len(view):
        written = sink.write(view[cursor:])
        if written is None:
            raise OSError("protocol sink write made no progress")
        if isinstance(written, bool) or not isinstance(written, int) or written <= 0:
            raise OSError("protocol sink write made no progress")
        if written > len(view) - cursor:
            raise OSError("protocol sink reported an oversized write")
        cursor += written


class ProtocolWriter:
    """Canonical reference writer retained only for protocol tests."""

    def __init__(self, sink: BinaryIO) -> None:
        if not hasattr(sink, "write"):
            raise ValueError("sink must be a binary writer")
        self._sink = sink
        self._header = None  # type: Header | None
        self._sequence = 0
        self._next_ordinal = 0
        self._digest = hashlib.sha256()
        self._fragment_count = 0
        self._batch_count = 0
        self._mate_count = 0
        self._template_count = 0
        self._site_count = 0
        self._per_contig = None  # type: list | None
        self._terminal = False
        self._failed = False

    @property
    def header(self) -> Header | None:
        return self._header

    def _require_open(self) -> None:
        if self._failed:
            raise ProtocolError("protocol writer is poisoned")
        if self._terminal:
            raise ProtocolError("protocol writer is already terminal")

    def write_header(self, header: Header) -> None:
        self._require_open()
        if self._header is not None:
            self._failed = True
            raise ProtocolError("protocol header was written twice")
        try:
            payload = _encode_header_payload(header)
            preamble = _protocol._PREAMBLE.pack(
                _protocol.MAGIC,
                _protocol.PROTOCOL_MAJOR,
                _protocol.PROTOCOL_MINOR,
                _protocol.PREAMBLE_FLAGS,
            )
            frame = _encode_frame(FrameType.HEADER, 0, 0, payload)
            _write_all(self._sink, preamble)
            _write_all(self._sink, frame)
        except Exception:
            self._failed = True
            raise
        self._header = header
        self._sequence = 1
        self._per_contig = [0] * len(header.contigs)
        self._digest.update(preamble)
        self._digest.update(frame)

    def write_batch(self, batch: FragmentBatch) -> None:
        self._require_open()
        if self._header is None or self._per_contig is None:
            self._failed = True
            raise ProtocolError("protocol header must precede batches")
        try:
            flags, payload = _encode_batch_payload(
                batch,
                self._header,
                expected_first_ordinal=self._next_ordinal,
            )
            frame = _encode_frame(
                FrameType.FRAGMENT_BATCH, flags, self._sequence, payload
            )
            _write_all(self._sink, frame)
        except Exception:
            self._failed = True
            raise
        self._sequence += 1
        self._next_ordinal += batch.fragment_count
        self._fragment_count += batch.fragment_count
        self._batch_count += 1
        self._mate_count += batch.mate_count
        self._template_count += batch.template_base_count
        self._site_count += batch.methylation_site_count
        for contig_index in batch.contig_indices:
            self._per_contig[contig_index] += 1
        self._digest.update(frame)

    def finish(self, *, skipped_fragment_count: int = 0) -> Trailer:
        self._require_open()
        if self._header is None or self._per_contig is None:
            self._failed = True
            raise ProtocolError("protocol header must precede the trailer")
        try:
            _protocol._require_u64("skipped_fragment_count", skipped_fragment_count)
            trailer = Trailer(
                fragment_count=self._fragment_count,
                fragment_batch_count=self._batch_count,
                mate_count=self._mate_count,
                template_base_count=self._template_count,
                methylation_site_count=self._site_count,
                skipped_fragment_count=skipped_fragment_count,
                per_contig_fragment_counts=tuple(self._per_contig),
                stream_sha256=self._digest.digest(),
            )
            payload = _encode_trailer_payload(
                trailer, len(self._header.contigs)
            )
            frame = _encode_frame(FrameType.TRAILER, 0, self._sequence, payload)
            _write_all(self._sink, frame)
        except Exception:
            self._failed = True
            raise
        self._sequence += 1
        self._terminal = True
        return trailer

    def write_error(self, error: ErrorFrame) -> None:
        self._require_open()
        if self._header is None:
            self._failed = True
            raise ProtocolError("protocol header must precede an error frame")
        try:
            payload = _protocol._encode_error_payload(error)
            frame = _encode_frame(FrameType.ERROR, 0, self._sequence, payload)
            _write_all(self._sink, frame)
        except Exception:
            self._failed = True
            raise
        self._sequence += 1
        self._terminal = True


def encode_stream(
    header: Header,
    batches: Iterable[FragmentBatch],
    *,
    skipped_fragment_count: int = 0,
) -> bytes:
    output = io.BytesIO()
    writer = ProtocolWriter(output)
    writer.write_header(header)
    for batch in batches:
        writer.write_batch(batch)
    writer.finish(skipped_fragment_count=skipped_fragment_count)
    return output.getvalue()


def read_stream(
    source: BinaryIO | bytes | bytearray | memoryview,
    *,
    expected_header: Header | None = None,
    expected_skipped_fragment_count: int | None = None,
    core_exit_status: int = 0,
) -> ProtocolStream:
    binary_source = (
        io.BytesIO(bytes(source))
        if isinstance(source, (bytes, bytearray, memoryview))
        else source
    )
    reader = ProtocolReader(
        binary_source,
        expected_header=expected_header,
        expected_skipped_fragment_count=expected_skipped_fragment_count,
    )
    stream = reader.read_all()
    reader.validate_core_exit_status(core_exit_status)
    return stream


__all__ = ["ProtocolWriter", "encode_stream", "read_stream"]


# Test-only protocol reference encoders.
from bsreadsim.htsim import protocol as _protocol_reference
FragmentBatch = _protocol_reference.FragmentBatch
Header = _protocol_reference.Header
MAX_FRAME_PAYLOAD = _protocol_reference.MAX_FRAME_PAYLOAD
ProtocolError = _protocol_reference.ProtocolError
DETAILS_PRESENT = _protocol_reference.DETAILS_PRESENT
Trailer = _protocol_reference.Trailer
_Encoder = _protocol_reference._Encoder
_F32 = _protocol_reference._F32
_U32 = _protocol_reference._U32
_U64 = _protocol_reference._U64
_U8 = _protocol_reference._U8
_validate_batch = _protocol_reference._validate_batch
_validate_header = _protocol_reference._validate_header
_validate_trailer = _protocol_reference._validate_trailer

def _encode_header_payload(header: Header) -> bytes:
    _validate_header(header)
    encoder = _Encoder()
    encoder.string(header.run_id)
    encoder.string(header.core_version)
    encoder.string(header.config_schema_version)
    encoder.string(header.rng_contract)
    encoder.u64(header.master_seed)
    encoder.raw(header.normalized_config_sha256)
    encoder.u8(int(header.technology))
    encoder.u8(int(header.has_details))
    encoder.u8(header.mates_per_fragment)
    encoder.u8(int(header.base_encoding))
    encoder.u8(int(header.ambiguity_policy))
    encoder.raw(b"\x00\x00\x00")
    encoder.u32(header.read_length_r1)
    encoder.u32(header.read_length_r2)
    encoder.u32(len(header.contigs))
    for contig in header.contigs:
        encoder.string(contig.name)
        encoder.u32(contig.length)
        encoder.raw(contig.reference_sha256)
    encoder.align4()
    return bytes(encoder.data)

def _encode_batch_payload(
    batch: FragmentBatch,
    header: Header,
    *,
    expected_first_ordinal: int | None = None,
) -> tuple[int, bytes]:
    _validate_batch(
        batch, header, expected_first_ordinal=expected_first_ordinal
    )
    encoder = _Encoder()
    for value in (
        batch.first_fragment_ordinal,
        batch.fragment_count,
        batch.template_base_count,
        batch.mate_count,
        batch.methylation_site_count,
    ):
        encoder.u32(value)
    for values in (
        batch.contig_indices,
        batch.reference_starts,
        batch.reference_ends,
        batch.template_offsets,
        batch.mate_offsets,
        batch.site_offsets,
        batch.mate_template_starts,
        batch.mate_template_ends,
        batch.site_template_offsets,
    ):
        encoder.array(values, _U32)
    encoder.array(batch.site_probabilities, _F32)
    for values in (
        batch.haplotypes,
        batch.capture_strands,
        batch.mate_indices,
        batch.mate_reverse_complements,
        batch.site_contexts,
        batch.methylation_sources,
        batch.site_alleles,
        batch.template_bases,
    ):
        encoder.array(values, _U8)

    flags = 0
    if batch.details is not None:
        flags = DETAILS_PRESENT
        details = batch.details
        encoder.align4()
        for value in (
            len(details.projection_template_starts),
            len(details.variant_indices),
            len(details.variant_ids),
            len(details.variant_ref_bases),
            len(details.variant_alt_bases),
            len(details.original_n_template_offsets),
        ):
            encoder.u32(value)
        for values in (
            details.projection_offsets,
            details.variant_offsets,
            details.original_n_offsets,
            details.projection_template_starts,
            details.projection_template_ends,
            details.projection_reference_starts,
            details.variant_indices,
            details.variant_id_offsets,
            details.variant_reference_starts,
            details.variant_reference_ends,
            details.variant_template_starts,
            details.variant_template_ends,
            details.variant_ref_offsets,
            details.variant_alt_offsets,
            details.site_reference_positions,
            details.original_n_template_offsets,
        ):
            encoder.array(values, _U32)
        for values in (
            details.variant_sources,
            details.variant_kinds,
            details.variant_phased_haplotypes,
            details.variant_ids,
            details.variant_ref_bases,
            details.variant_alt_bases,
        ):
            encoder.array(values, _U8)
    encoder.align4()
    payload = bytes(encoder.data)
    if len(payload) > MAX_FRAME_PAYLOAD:
        raise ProtocolError("fragment batch exceeds the 64 MiB payload limit")
    return flags, payload

def _encode_trailer_payload(trailer: Trailer, contig_count: int) -> bytes:
    _validate_trailer(trailer, contig_count)
    encoder = _Encoder()
    for value in (
        trailer.fragment_count,
        trailer.fragment_batch_count,
        trailer.mate_count,
        trailer.template_base_count,
        trailer.methylation_site_count,
        trailer.skipped_fragment_count,
    ):
        encoder.u64(value)
    encoder.u32(len(trailer.per_contig_fragment_counts))
    encoder.array(trailer.per_contig_fragment_counts, _U64)
    encoder.raw(trailer.stream_sha256)
    encoder.align4()
    return bytes(encoder.data)
