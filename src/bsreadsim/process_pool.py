"""Buffer-backed batch primitives for fragment post-processing.

The processing core accepts a caller-owned writable memoryview. Multi-process
adapters keep payloads and formatted output in bounded parent-owned shared
memory while only small descriptors cross the control queue. The inline
adapter uses the same core with one reusable local allocation.
"""

from __future__ import annotations

from dataclasses import dataclass
from multiprocessing import shared_memory
from typing import List, Optional, Tuple

from .bam import format_truth_sam_fragment
from .output import (
    NATIVE_TRUTH_JSON_AVAILABLE,
    _format_fragment_records_trusted,
)
from .postprocess import (
    PostprocessConfig,
    process_fragment_batch,
)
from .model import FragmentSummary
from .protocol import Header, TruthMode, decode_batch_payload
from .protocol_adapter import _decode_common_numpy_batch, decode_fragments


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
    record_lengths: Tuple[Tuple[int, int, int], ...]
    alignment_record_lengths: Tuple[int, ...]
    read1: Optional[BufferRegion]
    read2: Optional[BufferRegion]
    truth: Optional[BufferRegion]
    alignment: Optional[BufferRegion]
    summary: FragmentSummary
    required_capacity: int = 0

    @property
    def needs_resize(self) -> bool:
        return self.required_capacity > 0


_worker_header = None  # type: Optional[Header]
_worker_config = None  # type: Optional[PostprocessConfig]
_worker_paired_end = None  # type: Optional[bool]
_worker_include_truth = None  # type: Optional[bool]
_worker_emit_truth = None  # type: Optional[bool]
_worker_include_alignment = None  # type: Optional[bool]


def initialize_process_worker(
    header: Header,
    config: PostprocessConfig,
    paired_end: bool,
    include_truth: bool = True,
    emit_truth: Optional[bool] = None,
    include_alignment: bool = False,
) -> None:
    """Initialize immutable per-run worker state exactly once per process."""

    global _worker_header
    global _worker_config
    global _worker_paired_end
    global _worker_include_truth
    global _worker_emit_truth
    global _worker_include_alignment
    if not isinstance(include_truth, bool):
        raise RuntimeError("process worker include_truth must be a boolean")
    if emit_truth is None:
        emit_truth = include_truth
    if not isinstance(emit_truth, bool):
        raise RuntimeError("process worker emit_truth must be a boolean")
    if not isinstance(include_alignment, bool):
        raise RuntimeError("process worker include_alignment must be a boolean")
    if (emit_truth or include_alignment) and not include_truth:
        raise RuntimeError("formatted truth outputs require Full Truth columns")
    if not isinstance(header, Header):
        raise RuntimeError("process worker header must use the protocol contract")
    _worker_header = header
    _worker_config = config
    _worker_paired_end = paired_end
    _worker_include_truth = include_truth
    _worker_emit_truth = emit_truth
    _worker_include_alignment = include_alignment


def process_shared_batch(
    shared_name: str,
    payload_slices: Tuple[PayloadSlice, ...],
    output_offset: int,
) -> WorkerBatchResult:
    """Decode, process, and format one shared-memory payload batch."""

    state = (
        _worker_header,
        _worker_config,
        _worker_paired_end,
        _worker_include_truth,
        _worker_emit_truth,
        _worker_include_alignment,
    )
    if (
        state[0] is None
        or state[1] is None
        or state[2] is None
        or state[3] is None
        or state[4] is None
        or state[5] is None
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
            include_truth=state[3],
            emit_truth=state[4],
            include_alignment=state[5],
        )
    finally:
        buffer.release()
        memory.close()


def process_payload_batch(
    buffer: memoryview,
    payload_slices: Tuple[PayloadSlice, ...],
    output_offset: int,
    *,
    header: Header,
    config: PostprocessConfig,
    paired_end: bool,
    include_truth: bool,
    emit_truth: Optional[bool] = None,
    include_alignment: bool = False,
) -> WorkerBatchResult:
    """Revalidate, process, and format one authenticated batch payload."""

    if not isinstance(buffer, memoryview) or buffer.readonly:
        raise RuntimeError("payload batch buffer must be a writable memoryview")
    if len(payload_slices) != 1:
        raise RuntimeError("one frame payload must occupy one process batch")
    if not isinstance(header, Header):
        raise RuntimeError("process batch requires a protocol header")
    if (header.truth_columns is TruthMode.FULL) != include_truth:
        raise RuntimeError("process batch truth policy disagrees with its header")
    if emit_truth is None:
        emit_truth = include_truth
    if not isinstance(emit_truth, bool) or not isinstance(include_alignment, bool):
        raise RuntimeError("process batch formatted truth policies must be booleans")
    if (emit_truth or include_alignment) and not include_truth:
        raise RuntimeError("process batch truth output requires Full Truth columns")
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

    from .numpy_postprocess import (
        format_common_fastq_batch,
        supports_common_postprocess,
    )

    common_native = supports_common_postprocess(
        config,
        include_truth=include_truth,
    )
    if common_native:
        formatted_common = format_common_fastq_batch(
            _decode_common_numpy_batch(batch, header),
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
            compact_annotations=(include_truth and NATIVE_TRUTH_JSON_AVAILABLE),
            include_truth=include_truth,
        )

    observed_count = (
        len(formatted_common.record_lengths)
        if formatted_common is not None
        else len(processed_fragments)
    )
    if observed_count != batch.fragment_count:
        raise RuntimeError("post-process changed the batch fragment count")
    read1_parts = []  # type: List[bytes]
    read2_parts = []  # type: List[bytes]
    truth_parts = []  # type: List[bytes]
    alignment_parts = []  # type: List[bytes]
    record_lengths = []  # type: List[Tuple[int, int, int]]
    alignment_record_lengths = []  # type: List[int]
    formatted = []
    if formatted_common is not None:
        read1_parts.append(formatted_common.read1)
        if formatted_common.read2 is not None:
            read2_parts.append(formatted_common.read2)
        record_lengths.extend(formatted_common.record_lengths)
    else:
        for index, processed in enumerate(processed_fragments):
            formatted.append(
                _format_fragment_records_trusted(
                    processed,
                    paired_end=paired_end,
                    include_truth=emit_truth,
                )
            )
            if include_alignment:
                contig_index = batch.contig_indices[index]
                alignment_records = format_truth_sam_fragment(
                    processed,
                    paired_end=paired_end,
                    read_group_id=header.run_id,
                    contig_length=header.contigs[contig_index].length,
                )
                alignment_parts.extend(alignment_records)
                alignment_record_lengths.extend(
                    len(record) for record in alignment_records
                )
        for read1, read2, truth in formatted:
            if paired_end != (read2 is not None):
                raise RuntimeError("process result changed mate cardinality")
            if emit_truth != (truth is not None):
                raise RuntimeError("process result changed Truth presence")
            read1_parts.append(read1)
            if read2 is not None:
                read2_parts.append(read2)
            if truth is not None:
                truth_parts.append(truth)
            record_lengths.append(
                (
                    len(read1),
                    0 if read2 is None else len(read2),
                    0 if truth is None else len(truth),
                )
            )

    read1_size = sum(len(value) for value in read1_parts)
    read2_size = sum(len(value) for value in read2_parts)
    truth_size = sum(len(value) for value in truth_parts)
    alignment_size = sum(len(value) for value in alignment_parts)
    read1_offset = output_offset
    read2_offset = read1_offset + read1_size
    truth_offset = read2_offset + read2_size
    alignment_offset = truth_offset + truth_size
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
            truth=None,
            alignment=None,
            summary=summary,
            required_capacity=required_capacity,
        )

    _copy_parts(buffer, read1_offset, read1_parts)
    _copy_parts(buffer, read2_offset, read2_parts)
    if truth_parts:
        _copy_parts(buffer, truth_offset, truth_parts)
    if alignment_parts:
        _copy_parts(buffer, alignment_offset, alignment_parts)
    return WorkerBatchResult(
        first_ordinal=batch.first_fragment_ordinal,
        record_lengths=tuple(record_lengths),
        alignment_record_lengths=tuple(alignment_record_lengths),
        read1=BufferRegion(read1_offset, read1_size),
        read2=(BufferRegion(read2_offset, read2_size) if paired_end else None),
        truth=(BufferRegion(truth_offset, truth_size) if emit_truth else None),
        alignment=(
            BufferRegion(alignment_offset, alignment_size)
            if include_alignment
            else None
        ),
        summary=summary,
    )


def _copy_parts(buffer: memoryview, offset: int, parts: List[bytes]) -> None:
    cursor = offset
    for value in parts:
        end = cursor + len(value)
        buffer[cursor:end] = value
        cursor = end


class SharedBatchSlot:
    """One reusable parent-owned input/output shared-memory allocation."""

    def __init__(self, slot_id: int) -> None:
        self.slot_id = slot_id
        self._memory = None  # type: Optional[shared_memory.SharedMemory]
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
        payloads: Tuple[Tuple[int, bytes], ...],
        *,
        minimum_capacity: int = 0,
    ) -> Tuple[Tuple[PayloadSlice, ...], int]:
        if not payloads:
            raise RuntimeError("shared batch slot requires at least one payload")
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
        self._ensure_capacity(required)
        assert self._memory is not None
        descriptors = []
        cursor = 0
        buffer = self._memory.buf
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
        return tuple(descriptors), output_offset

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
        try:
            memory.unlink()
        except FileNotFoundError:
            pass

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
        payloads: Tuple[Tuple[int, bytes], ...],
        *,
        minimum_capacity: int = 0,
    ) -> Tuple[Tuple[PayloadSlice, ...], int]:
        if not payloads:
            raise RuntimeError("local batch slot requires at least one payload")
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
        self._ensure_capacity(required)
        descriptors = []
        cursor = 0
        buffer = memoryview(self._buffer)
        try:
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
        finally:
            buffer.release()
        return tuple(descriptors), output_offset

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
