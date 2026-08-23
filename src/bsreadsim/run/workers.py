"""Buffer-backed batch primitives for fragment processing.

The processing core accepts a caller-owned writable memoryview. Multi-process
adapters keep payloads and formatted output in bounded parent-owned shared
memory while only small descriptors cross the control queue. The inline
adapter uses the same core with one reusable local allocation.
"""

from __future__ import annotations
from contextlib import suppress

from dataclasses import dataclass
from multiprocessing import shared_memory

from ..output.fastq import format_fragment_records_trusted
from ..process.batch import FragmentSummary
from ..process.config import ProcessConfig
from ..htsim.protocol import Header, decode_batch_payload
from ..process.fragment import (
    decode_common_numpy_batch,
    decode_fragments,
    encode_fastq_batch,
    generate_columnar_reads,
    process_fragment_batch,
    supports_common_processing,
)
from ..output.bam import (
    format_sam_batch,
    format_sam_columns,
)


_MINIMUM_SLOT_BYTES = 1024 * 1024
_MAXIMUM_INITIAL_OUTPUT_BYTES = 32 * 1024 * 1024
_ALLOCATION_QUANTUM = 1024 * 1024


@dataclass(frozen=True)
class PayloadSlice:
    expected_ordinal: int
    offset: int
    length: int


@dataclass(frozen=True)
class BufferRegion:
    offset: int
    length: int


@dataclass(frozen=True)
class WorkerBatchResult:
    first_ordinal: int
    record_lengths: tuple[tuple[int, int], ...]
    alignment_record_lengths: tuple[int, ...]
    read1: BufferRegion | None
    read2: BufferRegion | None
    alignment: BufferRegion | None
    summary: FragmentSummary
    required_capacity: int = 0

    @property
    def needs_resize(self) -> bool:
        return self.required_capacity > 0


_worker_header = None  # type: Header | None
_worker_config = None  # type: ProcessConfig | None
_worker_paired_end = None  # type: bool | None
_worker_include_details = None  # type: bool | None
_worker_include_alignment = None  # type: bool | None
_worker_include_fragment_summary = None  # type: bool | None
_worker_include_fragment_realization = None  # type: bool | None


def initialize_process_worker(
    header: Header,
    config: ProcessConfig,
    paired_end: bool,
    include_details: bool = True,
    include_alignment: bool = False,
    include_fragment_summary: bool = False,
    include_fragment_realization: bool = False,
) -> None:
    """Initialize immutable per-run worker state exactly once per process."""

    global _worker_header
    global _worker_config
    global _worker_paired_end
    global _worker_include_details
    global _worker_include_alignment
    global _worker_include_fragment_summary
    global _worker_include_fragment_realization
    if not isinstance(include_details, bool):
        raise RuntimeError("process worker include_details must be a boolean")
    if not isinstance(include_alignment, bool):
        raise RuntimeError("process worker include_alignment must be a boolean")
    if not isinstance(include_fragment_summary, bool):
        raise RuntimeError("process worker fragment-summary policy must be boolean")
    if not isinstance(include_fragment_realization, bool):
        raise RuntimeError("process worker fragment-realization policy must be boolean")
    if include_fragment_realization and not include_fragment_summary:
        raise RuntimeError("fragment realization requires fragment summary")
    if include_fragment_summary and not include_alignment:
        raise RuntimeError("fragment summary requires BAM output")
    if include_alignment and not include_details:
        raise RuntimeError("formatted details outputs require Full Details columns")
    if not isinstance(header, Header):
        raise RuntimeError("process worker header must use the protocol contract")
    _worker_header = header
    _worker_config = config
    _worker_paired_end = paired_end
    _worker_include_details = include_details
    _worker_include_alignment = include_alignment
    _worker_include_fragment_summary = include_fragment_summary
    _worker_include_fragment_realization = include_fragment_realization


def process_shared_batch(
    shared_name: str,
    payload_slices: tuple[PayloadSlice, ...],
    output_offset: int,
) -> WorkerBatchResult:
    """Decode, process, and format one shared-memory payload batch."""

    state = (
        _worker_header,
        _worker_config,
        _worker_paired_end,
        _worker_include_details,
        _worker_include_alignment,
        _worker_include_fragment_summary,
        _worker_include_fragment_realization,
    )
    if (
        state[0] is None
        or state[1] is None
        or state[2] is None
        or state[3] is None
        or state[4] is None
        or state[5] is None
        or state[6] is None
    ):
        raise RuntimeError("fragment process worker was not initialized")
    memory = shared_memory.SharedMemory(name=shared_name, create=False)
    buffer = memory.buf
    try:
        return process_payload_batch(
            buffer,
            payload_slices,
            output_offset,
            header=state[0],
            config=state[1],
            paired_end=state[2],
            include_details=state[3],
            include_alignment=state[4],
            include_fragment_summary=state[5],
            include_fragment_realization=state[6],
        )
    finally:
        buffer.release()
        memory.close()


def process_payload_batch(
    buffer: memoryview,
    payload_slices: tuple[PayloadSlice, ...],
    output_offset: int,
    *,
    header: Header,
    config: ProcessConfig,
    paired_end: bool,
    include_details: bool,
    include_alignment: bool = False,
    include_fragment_summary: bool = False,
    include_fragment_realization: bool = False,
) -> WorkerBatchResult:
    """Revalidate, process, and format one authenticated batch payload."""

    if not isinstance(buffer, memoryview) or buffer.readonly:
        raise RuntimeError("payload batch buffer must be a writable memoryview")
    if len(payload_slices) != 1:
        raise RuntimeError("one frame payload must occupy one process batch")
    if not isinstance(header, Header):
        raise RuntimeError("process batch requires a protocol header")
    if header.has_details is not include_details:
        raise RuntimeError("process batch details policy disagrees with its header")
    if (
        not isinstance(include_alignment, bool)
        or not isinstance(include_fragment_summary, bool)
        or not isinstance(include_fragment_realization, bool)
    ):
        raise RuntimeError("process batch formatted details policies must be booleans")
    if include_alignment and not include_details:
        raise RuntimeError("process batch details output requires Full Details columns")
    if include_fragment_summary and not include_alignment:
        raise RuntimeError("fragment summary requires BAM output")
    if include_fragment_realization and not include_fragment_summary:
        raise RuntimeError("fragment realization requires fragment summary")
    descriptor = payload_slices[0]
    if (
        descriptor.offset < 0
        or descriptor.length <= 0
        or descriptor.offset + descriptor.length > output_offset
    ):
        raise RuntimeError("process batch payload is outside its input region")
    payload_view = buffer[
        descriptor.offset : descriptor.offset + descriptor.length
    ]
    try:
        # Shared memory is mutable and owned by the parent.  Freeze one bounded
        # frame payload before constructing retained immutable decoded views.
        payload = bytes(payload_view)
    finally:
        payload_view.release()
    batch = decode_batch_payload(
        payload,
        header,
        expected_first_ordinal=descriptor.expected_ordinal,
    )

    common_path = supports_common_processing(config)
    processed_common = None
    if common_path:
        common_batch = decode_common_numpy_batch(batch, header)
        if include_alignment:
            formatted_common = None
            processed_common = generate_columnar_reads(
                common_batch,
                config,
                include_fragment_summary=include_fragment_summary,
                include_fragment_realization=include_fragment_realization,
            )
        else:
            formatted_common = encode_fastq_batch(
                common_batch,
                config,
                paired_end,
            )
        processed_fragments = ()
    else:
        formatted_common = None
        fragments = decode_fragments(batch, header)
        contig_names = tuple(
            header.contigs[fragment.contig_index].name
            for fragment in fragments
        )
        processed_fragments = process_fragment_batch(
            fragments,
            contig_names,
            config,
            compact_base_states=include_details,
            include_details=include_details,
            include_fragment_realization=include_fragment_realization,
        )

    observed_count = (
        processed_common.fragment_count
        if processed_common is not None
        else (
            len(formatted_common.record_lengths)
            if formatted_common is not None
            else len(processed_fragments)
        )
    )
    if observed_count != batch.fragment_count:
        raise RuntimeError("process changed the batch fragment count")
    read1_parts = []  # type: list[bytes]
    read2_parts = []  # type: list[bytes]
    alignment_parts = []  # type: list[bytes]
    record_lengths = []  # type: list[tuple[int, int]]
    alignment_record_lengths = []  # type: list[int]
    formatted = []
    if processed_common is not None:
        record_lengths.extend(((0, 0),) * processed_common.fragment_count)
        alignment_payload, alignment_lengths = format_sam_columns(
            batch,
            header,
            processed_common,
            paired_end=paired_end,
            include_fragment_summary=include_fragment_summary,
            include_fragment_realization=include_fragment_realization,
        )
        alignment_parts.append(alignment_payload)
        alignment_record_lengths.extend(alignment_lengths)
    elif formatted_common is not None:
        read1_parts.append(formatted_common.read1)
        if formatted_common.read2 is not None:
            read2_parts.append(formatted_common.read2)
        record_lengths.extend(formatted_common.record_lengths)
        if include_alignment:
            alignment_payload, alignment_lengths = format_sam_columns(
                batch,
                header,
                formatted_common,
                paired_end=paired_end,
                include_fragment_summary=include_fragment_summary,
                include_fragment_realization=include_fragment_realization,
            )
            alignment_parts.append(alignment_payload)
            alignment_record_lengths.extend(alignment_lengths)
    else:
        if include_alignment:
            alignment_payload, alignment_lengths = format_sam_batch(
                processed_fragments,
                paired_end=paired_end,
                read_group_id=header.run_id,
                contig_lengths=tuple(
                    header.contigs[contig_index].length
                    for contig_index in batch.contig_indices
                ),
                include_fragment_summary=include_fragment_summary,
                include_fragment_realization=include_fragment_realization,
            )
            alignment_parts.append(alignment_payload)
            alignment_record_lengths.extend(alignment_lengths)
            record_lengths.extend(((0, 0),) * len(processed_fragments))
        else:
            for processed in processed_fragments:
                formatted.append(
                    format_fragment_records_trusted(
                        processed,
                        paired_end=paired_end,
                    )
                )
            for read1, read2 in formatted:
                if paired_end != (read2 is not None):
                    raise RuntimeError("process result changed mate cardinality")
                read1_parts.append(read1)
                if read2 is not None:
                    read2_parts.append(read2)
                record_lengths.append(
                    (
                        len(read1),
                        0 if read2 is None else len(read2),
                    )
                )

    read1_size = sum(len(value) for value in read1_parts)
    read2_size = sum(len(value) for value in read2_parts)
    alignment_size = sum(len(value) for value in alignment_parts)
    read1_offset = output_offset
    read2_offset = read1_offset + read1_size
    alignment_offset = read2_offset + read2_size
    required_capacity = alignment_offset + alignment_size
    per_contig = [0 for _ in header.contigs]
    for contig_index in batch.contig_indices:
        per_contig[contig_index] += 1
    summary = FragmentSummary(
        fragment_count=batch.fragment_count,
        mate_count=batch.mate_count,
        template_base_count=batch.template_base_count,
        methylation_site_count=batch.methylation_site_count,
        per_contig_fragment_counts=tuple(per_contig),
    )
    if required_capacity > len(buffer):
        return WorkerBatchResult(
            first_ordinal=batch.first_fragment_ordinal,
            record_lengths=tuple(record_lengths),
            alignment_record_lengths=tuple(alignment_record_lengths),
            read1=None,
            read2=None,
            alignment=None,
            summary=summary,
            required_capacity=required_capacity,
        )

    _copy_parts(buffer, read1_offset, read1_parts)
    _copy_parts(buffer, read2_offset, read2_parts)
    if alignment_parts:
        _copy_parts(buffer, alignment_offset, alignment_parts)
    return WorkerBatchResult(
        first_ordinal=batch.first_fragment_ordinal,
        record_lengths=tuple(record_lengths),
        alignment_record_lengths=tuple(alignment_record_lengths),
        read1=(BufferRegion(read1_offset, read1_size) if read1_parts else None),
        read2=(BufferRegion(read2_offset, read2_size) if read2_parts else None),
        alignment=(
            BufferRegion(alignment_offset, alignment_size)
            if include_alignment
            else None
        ),
        summary=summary,
    )


def _copy_parts(buffer: memoryview, offset: int, parts: list[bytes]) -> None:
    cursor = offset
    for value in parts:
        end = cursor + len(value)
        buffer[cursor:end] = value
        cursor = end


def _slot_requirements(
    payloads: tuple[tuple[int, bytes], ...],
    minimum_capacity: int,
) -> tuple[int, int]:
    input_size = sum(len(payload) for _, payload in payloads)
    output_offset = _align(input_size, 64)
    estimated_output = min(
        _MAXIMUM_INITIAL_OUTPUT_BYTES,
        max(
            _MINIMUM_SLOT_BYTES,
            input_size * 16 + len(payloads) * 4096,
        ),
    )
    required = max(minimum_capacity, output_offset + estimated_output)
    return required, output_offset


def _copy_payloads(
    buffer: memoryview,
    payloads: tuple[tuple[int, bytes], ...],
) -> tuple[PayloadSlice, ...]:
    descriptors = []
    cursor = 0
    for expected_ordinal, payload in payloads:
        end = cursor + len(payload)
        buffer[cursor:end] = payload
        descriptors.append(
            PayloadSlice(
                expected_ordinal=expected_ordinal,
                offset=cursor,
                length=len(payload),
            )
        )
        cursor = end
    return tuple(descriptors)


class SharedBatchSlot:
    """One reusable parent-owned input/output shared-memory allocation."""

    def __init__(self, slot_id: int) -> None:
        self.slot_id = slot_id
        self._memory = None  # type: shared_memory.SharedMemory | None
        self._capacity = 0

    @property
    def name(self) -> str:
        if self._memory is None:
            raise RuntimeError("shared batch slot is not allocated")
        return self._memory.name

    @property
    def capacity(self) -> int:
        return self._capacity

    def prepare(
        self,
        payloads: tuple[tuple[int, bytes], ...],
        *,
        minimum_capacity: int = 0,
    ) -> tuple[tuple[PayloadSlice, ...], int]:
        if not payloads:
            raise RuntimeError("shared batch slot requires at least one payload")
        required, output_offset = _slot_requirements(payloads, minimum_capacity)
        self._ensure_capacity(required)
        assert self._memory is not None
        descriptors = _copy_payloads(self._memory.buf, payloads)
        return descriptors, output_offset

    def view(self, region: BufferRegion) -> memoryview:
        if self._memory is None:
            raise RuntimeError("shared batch slot is not allocated")
        if (
            region.offset < 0
            or region.length < 0
            or region.offset + region.length > self._capacity
        ):
            raise RuntimeError("shared batch output region is outside the slot")
        return self._memory.buf[region.offset : region.offset + region.length]

    def buffer_view(self) -> memoryview:
        """Expose the complete allocation to an inline benchmark executor."""
        if self._memory is None:
            raise RuntimeError("shared batch slot is not allocated")
        return self._memory.buf[:]

    def close(self) -> None:
        memory = self._memory
        self._memory = None
        self._capacity = 0
        if memory is None:
            return
        memory.close()
        with suppress(FileNotFoundError):
            memory.unlink()

    def _ensure_capacity(self, required: int) -> None:
        if required <= self._capacity:
            return
        if required <= 0:
            raise RuntimeError("shared batch capacity must be positive")
        self.close()
        capacity = _align(max(required, _MINIMUM_SLOT_BYTES), _ALLOCATION_QUANTUM)
        self._memory = shared_memory.SharedMemory(create=True, size=capacity)
        self._capacity = capacity


class LocalBatchSlot:
    """One reusable in-process input/output allocation for inline batches."""

    def __init__(self, slot_id: int) -> None:
        self.slot_id = slot_id
        self._buffer = bytearray()

    @property
    def capacity(self) -> int:
        return len(self._buffer)

    def prepare(
        self,
        payloads: tuple[tuple[int, bytes], ...],
        *,
        minimum_capacity: int = 0,
    ) -> tuple[tuple[PayloadSlice, ...], int]:
        if not payloads:
            raise RuntimeError("local batch slot requires at least one payload")
        required, output_offset = _slot_requirements(payloads, minimum_capacity)
        self._ensure_capacity(required)
        buffer = memoryview(self._buffer)
        try:
            descriptors = _copy_payloads(buffer, payloads)
        finally:
            buffer.release()
        return descriptors, output_offset

    def buffer_view(self) -> memoryview:
        if not self._buffer:
            raise RuntimeError("local batch slot is not allocated")
        return memoryview(self._buffer)

    def view(self, region: BufferRegion) -> memoryview:
        if (
            region.offset < 0
            or region.length < 0
            or region.offset + region.length > self.capacity
        ):
            raise RuntimeError("local batch output region is outside the slot")
        return memoryview(self._buffer)[
            region.offset : region.offset + region.length
        ]

    def close(self) -> None:
        self._buffer = bytearray()

    def _ensure_capacity(self, required: int) -> None:
        if required <= self.capacity:
            return
        if required <= 0:
            raise RuntimeError("local batch capacity must be positive")
        capacity = _align(max(required, _MINIMUM_SLOT_BYTES), _ALLOCATION_QUANTUM)
        self._buffer = bytearray(capacity)


def _align(value: int, quantum: int) -> int:
    return ((value + quantum - 1) // quantum) * quantum


__all__ = [
    "BufferRegion",
    "LocalBatchSlot",
    "PayloadSlice",
    "SharedBatchSlot",
    "WorkerBatchResult",
    "initialize_process_worker",
    "process_payload_batch",
    "process_shared_batch",
]
