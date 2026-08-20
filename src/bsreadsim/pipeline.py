"""Bounded orchestration for the C++ generator and Python post-processor.

The module is the sole runtime owner of the component boundary:

* C++ owns reference parsing, genomic catalogs, fragment allocation, and the
  versioned protocol stream.
* Python validates that stream identity, applies pure fragment-level stages,
  and publishes FASTQ/truth/manifest artifacts transactionally.

For a single worker, every policy uses the raw-payload batch core inline from
one reusable local slot. Eligible FASTQ-only blocks take its native column and
NumPy branch; other policies take the byte-equivalent typed-object fallback
inside the same core. With multiple workers, the supervisor remains the sole
stdout reader and dispatches CRC-verified payload batches through bounded
shared-memory slots. Worker processes decode, post-process, and format those
batches; only small descriptors cross the multiprocessing control queue. A
single writer publishes completed slots strictly by fragment ordinal, so
worker count and completion order cannot affect output bytes.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass
import multiprocessing
import os
from pathlib import Path
import pickle
import queue
import threading
from typing import (
    Deque,
    Dict,
    Mapping,
    Optional,
    Tuple,
    Union,
)
import uuid

from . import __version__
from .bam import build_truth_sam_header
from .config import (
    LoadedRunConfig,
    WGBS_GC_PROFILE_FORMAT,
    WGBS_GC_PROFILE_VERSION,
    normalize_run_config,
)
from .core_argv import build_core_argv
from .core_process import CoreProcess
from .manifest import (
    CompleteManifest,
    build_complete_manifest,
    validate_header_projection,
)
from .output import (
    NATIVE_TRUTH_JSON_AVAILABLE,
    OutputConfig,
    OutputSummary,
    OutputTransaction,
    TruthBamConfig,
)
from .process_pool import (
    LocalBatchSlot,
    PayloadSlice,
    SharedBatchSlot,
    WorkerBatchResult,
    initialize_process_worker,
    process_payload_batch,
    process_shared_batch,
)
from .postprocess import (
    PostprocessConfig,
    PostprocessError,
    UniformError,
    UniformQuality,
)
from .preparation import (
    EntropySource,
    PreparationError,
    PreparedRun,
    prepare_run,
    snapshot_prepared_file,
)
from .model import FragmentSummary
from .protocol import (
    Header,
    Trailer,
)
from .runtime import (
    CoreExecutableError,
    resolve_core_executable as _resolve_core_executable,
)
from .sequencing_models import (
    MAX_MODEL_BYTES,
    QUALITY_CONFUSION_FORMAT,
    QUALITY_CONFUSION_VERSION,
    QUALITY_MARKOV_FORMAT,
    QUALITY_MARKOV_VERSION,
    QualityConfusionModel,
    QualityMarkovModel,
    SequencingModelError,
    parse_quality_confusion,
    parse_quality_markov,
)


PathLike = Union[str, os.PathLike]

_PROCESS_BATCH_MAX_FRAGMENTS = 64
_FULL_TRUTH_BATCH_FRAGMENTS = 8
_PROCESS_RESULT_POLL_SECONDS = 0.05

class PipelineError(RuntimeError):
    """The requested run is outside the implemented component boundary."""


@dataclass(frozen=True)
class RunResult:
    """Published result of one completely verified simulation run."""

    run_id: str
    manifest_path: Path
    manifest: CompleteManifest
    outputs: OutputSummary


@dataclass(frozen=True)
class _QueuedFormattedBatch:
    slot: SharedBatchSlot
    result: WorkerBatchResult
    release_queue: "queue.Queue[Tuple[int, int]]"


@dataclass
class _PendingProcessBatch:
    slot: Union[LocalBatchSlot, SharedBatchSlot]
    payloads: Tuple[Tuple[int, bytes], ...]
    descriptors: Tuple[PayloadSlice, ...]
    output_offset: int
    result: object
    fragment_count: int


class _SharedSlotWindow:
    """Bound running, completed, and writer-buffered process work together."""

    def __init__(self, slots: Tuple[SharedBatchSlot, ...], maximum: int) -> None:
        self.slots = {slot.slot_id: slot for slot in slots}
        self.free = deque(slots)  # type: Deque[SharedBatchSlot]
        self.busy = {}  # type: Dict[int, int]
        self.maximum = maximum
        self.in_flight = 0
        self.releases = queue.Queue(  # type: queue.Queue[Tuple[int, int]]
            maxsize=len(slots)
        )

    def can_submit(self, fragment_count: int) -> bool:
        return bool(self.free) and self.in_flight + fragment_count <= self.maximum

    def acquire(self, fragment_count: int) -> SharedBatchSlot:
        if not self.can_submit(fragment_count):
            raise PipelineError("shared process window capacity was exceeded")
        slot = self.free.popleft()
        if slot.slot_id in self.busy:
            raise PipelineError("shared process slot was acquired twice")
        self.busy[slot.slot_id] = fragment_count
        self.in_flight += fragment_count
        return slot

    def release(self, slot_id: int, fragment_count: int) -> None:
        expected = self.busy.pop(slot_id, None)
        if expected is None:
            raise PipelineError("shared process slot was released twice")
        if expected != fragment_count:
            raise PipelineError("shared process slot release count changed")
        slot = self.slots.get(slot_id)
        if slot is None:
            raise PipelineError("unknown shared process slot was released")
        self.in_flight -= fragment_count
        if self.in_flight < 0:
            raise PipelineError("shared process window count became negative")
        self.free.append(slot)

    def drain_releases(self) -> None:
        while True:
            try:
                slot_id, fragment_count = self.releases.get_nowait()
            except queue.Empty:
                return
            self.release(slot_id, fragment_count)

    def wait_for_release(self, writer: "_AsyncFragmentWriter") -> None:
        while True:
            writer._raise_if_failed()
            try:
                slot_id, fragment_count = self.releases.get(
                    timeout=_PROCESS_RESULT_POLL_SECONDS
                )
            except queue.Empty:
                continue
            self.release(slot_id, fragment_count)
            self.drain_releases()
            return

    def require_fully_released(self) -> None:
        self.drain_releases()
        if self.busy or self.in_flight != 0 or len(self.free) != len(self.slots):
            raise PipelineError("shared process slots were not fully released")


class _AsyncFragmentWriter:
    """Single-owner output writer behind a small bounded handoff queue."""

    _PUT_TIMEOUT_SECONDS = 0.05

    def __init__(self, output: OutputTransaction, capacity: int) -> None:
        self._output = output
        self._queue = queue.Queue(maxsize=capacity)
        self._sentinel = object()
        self._cancelled = threading.Event()
        self._error = None  # type: Optional[BaseException]
        self._closed = False
        self._thread = threading.Thread(
            target=self._run,
            name="bsreadsim-output-writer",
            daemon=True,
        )

    def __enter__(self) -> "_AsyncFragmentWriter":
        self._thread.start()
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        self.close(cancel=exc_type is not None)
        return False

    def write_formatted_batch(
        self,
        slot: SharedBatchSlot,
        result: WorkerBatchResult,
        release_queue: "queue.Queue[Tuple[int, int]]",
    ) -> None:
        if self._closed:
            raise PipelineError("output writer is already closed")
        item = _QueuedFormattedBatch(slot, result, release_queue)
        while True:
            self._raise_if_failed()
            try:
                self._queue.put(item, timeout=self._PUT_TIMEOUT_SECONDS)
                return
            except queue.Full:
                continue

    def close(self, *, cancel: bool = False) -> None:
        if self._closed:
            return
        if cancel:
            self._cancelled.set()
        if self._thread.is_alive():
            while self._thread.is_alive():
                try:
                    self._queue.put(
                        self._sentinel,
                        timeout=self._PUT_TIMEOUT_SECONDS,
                    )
                    break
                except queue.Full:
                    continue
        self._thread.join()
        self._closed = True
        if not cancel:
            self._raise_if_failed()

    def _run(self) -> None:
        try:
            while True:
                item = self._queue.get()
                if item is self._sentinel:
                    return
                if isinstance(item, _QueuedFormattedBatch):
                    try:
                        if not self._cancelled.is_set():
                            _write_shared_formatted_batch(
                                self._output,
                                item.slot,
                                item.result,
                            )
                    finally:
                        item.release_queue.put(
                            (item.slot.slot_id, len(item.result.record_lengths))
                        )
                else:
                    raise PipelineError("output writer received an invalid item")
        except BaseException as error:
            self._error = error

    def _raise_if_failed(self) -> None:
        if self._error is not None:
            raise PipelineError("output writer failed") from self._error


def _write_shared_formatted_batch(
    output: OutputTransaction,
    slot: Union[LocalBatchSlot, SharedBatchSlot],
    result: WorkerBatchResult,
) -> None:
    truth_enabled = output.config.truth == "full"
    alignment_enabled = output.config.truth_bam is not None
    if result.needs_resize or result.read1 is None:
        raise PipelineError("process worker returned an incomplete output batch")
    if truth_enabled != (result.truth is not None):
        raise PipelineError("process worker truth output disagrees with policy")
    if alignment_enabled != (result.alignment is not None):
        raise PipelineError(
            "process worker truth alignment output disagrees with policy"
        )
    views = []
    try:
        read1 = slot.view(result.read1)
        views.append(read1)
        truth = None
        if result.truth is not None:
            truth = slot.view(result.truth)
            views.append(truth)
        read2 = None
        if result.read2 is not None:
            read2 = slot.view(result.read2)
            views.append(read2)
        alignment = None
        if result.alignment is not None:
            alignment = slot.view(result.alignment)
            views.append(alignment)
        output.write_formatted_batch(
            result.first_ordinal,
            result.record_lengths,
            read1,
            read2,
            truth,
            alignment_record_lengths=result.alignment_record_lengths,
            alignment=alignment,
        )
    finally:
        for view in views:
            view.release()


def resolve_core_executable(value: Optional[PathLike] = None) -> Path:
    try:
        return _resolve_core_executable(value)
    except CoreExecutableError as error:
        raise PipelineError(str(error)) from error


def run_document(
    document: Mapping[str, object],
    *,
    base_directory: PathLike = ".",
    core_executable: Optional[PathLike] = None,
    run_id: Optional[str] = None,
    entropy: Optional[EntropySource] = None,
    mode: str = "production",
) -> RunResult:
    """Normalize and execute one in-memory command-line configuration."""
    loaded = normalize_run_config(document, base_directory, mode=mode)
    prepared = (
        prepare_run(loaded)
        if entropy is None
        else prepare_run(loaded, entropy=entropy)
    )
    return run_prepared(
        prepared,
        core_executable=core_executable,
        run_id=run_id,
    )


def run_prepared(
    prepared: PreparedRun,
    *,
    core_executable: Optional[PathLike] = None,
    run_id: Optional[str] = None,
) -> RunResult:
    """Run one prepared configuration through the production data path."""
    if not isinstance(prepared, PreparedRun):
        raise PipelineError("prepared must be a PreparedRun")
    config = prepared.config
    if not isinstance(config, LoadedRunConfig) or config.master_seed is None:
        raise PipelineError("prepared config must contain a materialized seed")

    normalized = config.normalized
    _require_released_capabilities(normalized)
    postprocess_config = _build_postprocess_config(prepared)
    executable = resolve_core_executable(core_executable)
    effective_run_id = str(uuid.uuid4()) if run_id is None else run_id
    fragments = _mapping(normalized, "fragments")
    execution = _mapping(normalized, "execution")
    output = _mapping(normalized, "output")
    emit_truth = output["truth"] == "full"
    include_alignment = bool(output["truth_bam"])
    include_truth = emit_truth or include_alignment
    truth_columns = "full" if include_truth else "none"
    argv = build_core_argv(
        prepared,
        effective_run_id,
        executable,
        truth_columns=truth_columns,
        protocol_batch_fragments=_protocol_batch_fragment_limit(
            truth_columns=truth_columns,
            max_in_flight=execution["max_in_flight_fragments"],
        ),
    )
    if execution["workers"] > 1:
        _require_picklable_process_state(postprocess_config)

    core = CoreProcess(
        argv,
        read_length=fragments["read_length_1"],
        paired_end=fragments["paired_end"],
        expected_skipped_fragment_count=_expected_skipped_fragment_count(
            normalized
        ),
    )
    with core:
        _validate_core_header(
            prepared,
            core.header,
            expected_run_id=effective_run_id,
        )
        truth_bam_config = None
        if include_alignment:
            truth_bam_config = TruthBamConfig(
                writer_argv=(
                    str(executable),
                    "--sam-to-bam",
                    str(output["gzip_level"]),
                ),
                sam_header=build_truth_sam_header(
                    core.header,
                    sample_name=output["prefix"],
                    program_version=__version__,
                ),
                references=tuple(
                    (contig.name, contig.length) for contig in core.header.contigs
                ),
                read_group_id=core.header.run_id,
            )
        output_config = OutputConfig(
            directory=Path(output["directory"]),
            prefix=output["prefix"],
            paired_end=fragments["paired_end"],
            compression=output["compression"],
            gzip_level=output["gzip_level"],
            truth=output["truth"],
            truth_bam=truth_bam_config,
        )
        with OutputTransaction(output_config) as transaction:
            if execution["workers"] == 1:
                _consume_batches_inline(
                    core,
                    core.header,
                    postprocess_config,
                    transaction,
                    paired_end=fragments["paired_end"],
                    include_truth=include_truth,
                    emit_truth=emit_truth,
                    include_alignment=include_alignment,
                    max_in_flight=execution["max_in_flight_fragments"],
                )
            else:
                _consume_batches(
                    core,
                    core.header,
                    postprocess_config,
                    transaction,
                    workers=execution["workers"],
                    max_in_flight=execution["max_in_flight_fragments"],
                    paired_end=fragments["paired_end"],
                    include_truth=include_truth,
                    emit_truth=emit_truth,
                    include_alignment=include_alignment,
                )
            output_summary = transaction.finalize()
            manifest = build_complete_manifest(
                prepared,
                core.header,
                core.trailer,
                output_summary,
            )
            transaction.commit(manifest.canonical_json)

    manifest_path = Path(output["directory"]) / "{}.manifest.json".format(
        output["prefix"]
    )
    return RunResult(
        run_id=effective_run_id,
        manifest_path=manifest_path,
        manifest=manifest,
        outputs=output_summary,
    )


def _consume_batches_inline(
    batches: CoreProcess,
    header: Header,
    config: PostprocessConfig,
    output: OutputTransaction,
    *,
    paired_end: bool,
    include_truth: bool,
    emit_truth: bool,
    include_alignment: bool,
    max_in_flight: int,
) -> None:
    """Process authenticated batches inline through one bounded slot."""

    slot = LocalBatchSlot(0)
    summary = _FragmentSummaryAccumulator(len(header.contigs))
    try:
        for decoded in batches.iter_batches():
            if decoded.fragment_count > max_in_flight:
                raise PipelineError(
                    "protocol batch exceeds the configured in-flight bound"
                )
            payloads = ((
                decoded.first_fragment_ordinal,
                bytes(decoded.raw_payload),
            ),)
            descriptors, output_offset = slot.prepare(payloads)
            while True:
                buffer = slot.buffer_view()
                try:
                    result = process_payload_batch(
                        buffer,
                        descriptors,
                        output_offset,
                        header=header,
                        config=config,
                        paired_end=paired_end,
                        include_truth=include_truth,
                        emit_truth=emit_truth,
                        include_alignment=include_alignment,
                    )
                finally:
                    buffer.release()
                pending = _PendingProcessBatch(
                    slot=slot,
                    payloads=payloads,
                    descriptors=descriptors,
                    output_offset=output_offset,
                    result=result,
                    fragment_count=decoded.fragment_count,
                )
                _validate_process_batch_result(pending, result)
                if not result.needs_resize:
                    break
                descriptors, output_offset = slot.prepare(
                    payloads,
                    minimum_capacity=result.required_capacity,
                )
            summary.add(result.summary)
            _write_shared_formatted_batch(output, slot, result)
        _validate_processed_summary(batches.trailer, summary.build())
    finally:
        slot.close()


def _consume_batches(
    batches: CoreProcess,
    header: Header,
    config: PostprocessConfig,
    output: OutputTransaction,
    *,
    workers: int,
    max_in_flight: int,
    paired_end: bool,
    include_truth: bool,
    emit_truth: bool,
    include_alignment: bool,
) -> None:
    """Dispatch whole frame payloads through ordered shared-memory slots."""

    worker_count = min(workers, max_in_flight)
    slot_count = min(max_in_flight, max(1, worker_count * 2))
    context = _safe_process_context()
    slots = tuple(SharedBatchSlot(index) for index in range(slot_count))
    window = _SharedSlotWindow(slots, max_in_flight)
    pending = deque()  # type: Deque[_PendingProcessBatch]
    summary = _FragmentSummaryAccumulator(len(header.contigs))
    pool = None
    successful = False
    try:
        pool = context.Pool(
            processes=worker_count,
            initializer=initialize_process_worker,
            initargs=(
                header,
                config,
                paired_end,
                include_truth,
                emit_truth,
                include_alignment,
            ),
        )
        with _AsyncFragmentWriter(output, min(slot_count, 4)) as writer:
            for decoded in batches.iter_batches():
                fragment_count = decoded.fragment_count
                if fragment_count > max_in_flight:
                    raise PipelineError(
                        "protocol batch exceeds the configured in-flight bound"
                    )
                while not window.can_submit(fragment_count):
                    if pending:
                        _collect_process_batch(
                            pending,
                            pool,
                            window,
                            writer,
                            summary,
                        )
                    else:
                        if not isinstance(writer, _AsyncFragmentWriter):
                            raise PipelineError(
                                "process window stalled without pending work"
                            )
                        window.wait_for_release(writer)

                payloads = ((
                    decoded.first_fragment_ordinal,
                    bytes(decoded.raw_payload),
                ),)
                slot = window.acquire(fragment_count)
                try:
                    descriptors, output_offset = slot.prepare(payloads)
                    pending.append(
                        _submit_process_batch(
                            pool,
                            slot,
                            payloads,
                            descriptors,
                            output_offset,
                            fragment_count=fragment_count,
                        )
                    )
                except Exception:
                    window.release(slot.slot_id, fragment_count)
                    raise
                window.drain_releases()

            while pending:
                _collect_process_batch(
                    pending,
                    pool,
                    window,
                    writer,
                    summary,
                )

        window.require_fully_released()
        _validate_processed_summary(batches.trailer, summary.build())
        successful = True
    finally:
        if pool is not None:
            if successful:
                pool.close()
            else:
                pool.terminate()
            pool.join()
        for slot in slots:
            slot.close()


def _require_picklable_process_state(config: PostprocessConfig) -> None:
    try:
        pickle.dumps(config, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as error:
        raise PipelineError(
            "multi-process post-processing state is not serializable"
        ) from error


def _protocol_batch_fragment_limit(
    *,
    truth_columns: str,
    max_in_flight: int,
) -> int:
    """Bound Full Truth object lifetimes without shrinking common batches."""

    maximum = min(_PROCESS_BATCH_MAX_FRAGMENTS, max_in_flight)
    if truth_columns == "full":
        return min(_FULL_TRUTH_BATCH_FRAGMENTS, maximum)
    return maximum


def _safe_process_context():
    # The core supervisor and ordered-output path already own live threads.
    # ``spawn`` neither forks those threads nor needs the forkserver's private
    # Unix-domain listener, which is unavailable in some sandboxed/HPC jobs.
    return multiprocessing.get_context("spawn")


def _submit_process_batch(
    pool,
    slot: SharedBatchSlot,
    payloads: Tuple[Tuple[int, bytes], ...],
    descriptors: Tuple[PayloadSlice, ...],
    output_offset: int,
    *,
    fragment_count: Optional[int] = None,
) -> _PendingProcessBatch:
    result = pool.apply_async(
        process_shared_batch,
        (slot.name, descriptors, output_offset),
    )
    expected_count = len(payloads) if fragment_count is None else fragment_count
    return _PendingProcessBatch(
        slot=slot,
        payloads=payloads,
        descriptors=descriptors,
        output_offset=output_offset,
        result=result,
        fragment_count=expected_count,
    )


def _collect_process_batch(
    pending: Deque[_PendingProcessBatch],
    pool,
    window: _SharedSlotWindow,
    writer: _AsyncFragmentWriter,
    summary: "_FragmentSummaryAccumulator",
) -> None:
    if not pending:
        raise PipelineError("cannot collect an empty process queue")
    batch = pending[0]
    while True:
        if isinstance(writer, _AsyncFragmentWriter):
            writer._raise_if_failed()
        window.drain_releases()
        try:
            result = batch.result.get(timeout=_PROCESS_RESULT_POLL_SECONDS)
        except multiprocessing.TimeoutError:
            continue
        if not isinstance(result, WorkerBatchResult):
            raise PipelineError("process worker returned an unsupported result")
        _validate_process_batch_result(batch, result)
        if not result.needs_resize:
            break
        descriptors, output_offset = batch.slot.prepare(
            batch.payloads,
            minimum_capacity=result.required_capacity,
        )
        batch.descriptors = descriptors
        batch.output_offset = output_offset
        batch.result = pool.apply_async(
            process_shared_batch,
            (batch.slot.name, descriptors, output_offset),
        )

    pending.popleft()
    summary.add(result.summary)
    writer.write_formatted_batch(
        batch.slot,
        result,
        window.releases,
    )
    window.drain_releases()


def _validate_process_batch_result(
    batch: _PendingProcessBatch,
    result: WorkerBatchResult,
) -> None:
    expected_count = batch.fragment_count
    expected_first = batch.payloads[0][0]
    if (
        isinstance(result.first_ordinal, bool)
        or not isinstance(result.first_ordinal, int)
        or result.first_ordinal != expected_first
    ):
        raise PipelineError("process worker changed a batch ordinal")
    if len(result.record_lengths) != expected_count:
        raise PipelineError("process worker changed a batch fragment count")
    if any(
        not isinstance(lengths, tuple)
        or len(lengths) != 3
        or any(
            isinstance(value, bool)
            or not isinstance(value, int)
            or value < 0
            for value in lengths
        )
        for lengths in result.record_lengths
    ):
        raise PipelineError("process worker returned invalid record lengths")
    if (
        isinstance(result.summary.fragment_count, bool)
        or not isinstance(result.summary.fragment_count, int)
        or result.summary.fragment_count != expected_count
    ):
        raise PipelineError("process worker summary changed a batch fragment count")
    if any(
        isinstance(value, bool) or not isinstance(value, int) or value < 0
        for value in (
            result.summary.mate_count,
            result.summary.template_base_count,
            result.summary.methylation_site_count,
        )
    ):
        raise PipelineError("process worker returned an invalid batch summary")
    if not isinstance(result.alignment_record_lengths, tuple) or any(
        isinstance(length, bool) or not isinstance(length, int) or length <= 0
        for length in result.alignment_record_lengths
    ):
        raise PipelineError(
            "process worker returned invalid truth alignment lengths"
        )
    if len(result.alignment_record_lengths) not in (
        0,
        result.summary.mate_count,
    ):
        raise PipelineError("process worker changed the truth alignment mate count")
    if result.needs_resize:
        if result.required_capacity <= batch.slot.capacity:
            raise PipelineError("process worker requested a non-growing shared slot")
        if (
            result.read1 is not None
            or result.read2 is not None
            or result.truth is not None
            or result.alignment is not None
        ):
            raise PipelineError("resize result unexpectedly exposed output regions")
        return
    if result.read1 is None:
        raise PipelineError("process worker omitted an output region")
    if result.required_capacity != 0:
        raise PipelineError("process worker returned conflicting resize state")
    truth_lengths = tuple(lengths[2] for lengths in result.record_lengths)
    if (result.truth is None) != all(length == 0 for length in truth_lengths):
        raise PipelineError("process worker returned conflicting truth regions")
    if (result.alignment is None) != (
        len(result.alignment_record_lengths) == 0
    ):
        raise PipelineError(
            "process worker returned conflicting truth alignment regions"
        )
    if result.alignment is not None:
        if sum(result.alignment_record_lengths) != result.alignment.length:
            raise PipelineError(
                "process worker truth alignment sizes disagree with records"
            )
    regions = [result.read1]
    if result.truth is not None:
        regions.append(result.truth)
    if result.read2 is not None:
        regions.append(result.read2)
    if result.alignment is not None:
        regions.append(result.alignment)
    ordered = sorted(regions, key=lambda region: region.offset)
    previous_end = batch.output_offset
    for region in ordered:
        if region.offset < previous_end or region.length < 0:
            raise PipelineError("process worker output regions overlap")
        if region.offset + region.length > batch.slot.capacity:
            raise PipelineError("process worker output region exceeds its shared slot")
        previous_end = region.offset + region.length


class _FragmentSummaryAccumulator:
    def __init__(self, contig_count: int) -> None:
        self.fragment_count = 0
        self.mate_count = 0
        self.template_base_count = 0
        self.methylation_site_count = 0
        self.per_contig_fragment_counts = [0 for _ in range(contig_count)]

    def add(self, summary: FragmentSummary) -> None:
        if len(summary.per_contig_fragment_counts) != len(
            self.per_contig_fragment_counts
        ):
            raise PipelineError("process worker summary changed contig cardinality")
        if any(
            isinstance(value, bool) or not isinstance(value, int) or value < 0
            for value in summary.per_contig_fragment_counts
        ):
            raise PipelineError("process worker returned invalid per-contig counts")
        if sum(summary.per_contig_fragment_counts) != summary.fragment_count:
            raise PipelineError("process worker per-contig counts changed batch size")
        self.fragment_count += summary.fragment_count
        self.mate_count += summary.mate_count
        self.template_base_count += summary.template_base_count
        self.methylation_site_count += summary.methylation_site_count
        for index, count in enumerate(summary.per_contig_fragment_counts):
            self.per_contig_fragment_counts[index] += count

    def build(self) -> FragmentSummary:
        return FragmentSummary(
            fragment_count=self.fragment_count,
            mate_count=self.mate_count,
            template_base_count=self.template_base_count,
            methylation_site_count=self.methylation_site_count,
            per_contig_fragment_counts=tuple(self.per_contig_fragment_counts),
        )


def _validate_processed_summary(
    trailer: Trailer,
    summary: FragmentSummary,
) -> None:
    """Cross-check independently processed aggregates before publication."""

    if not isinstance(trailer, Trailer):
        raise PipelineError("protocol processing did not retain a trailer")
    observed = (
        summary.fragment_count,
        summary.mate_count,
        summary.template_base_count,
        summary.methylation_site_count,
        summary.per_contig_fragment_counts,
    )
    declared = (
        trailer.fragment_count,
        trailer.mate_count,
        trailer.template_base_count,
        trailer.methylation_site_count,
        trailer.per_contig_fragment_counts,
    )
    if observed != declared:
        raise PipelineError(
            "protocol trailer disagrees with worker-observed aggregates"
        )


def _validate_core_header(
    prepared: PreparedRun,
    header: Header,
    *,
    expected_run_id: str,
) -> None:
    if not isinstance(header, Header):
        raise PipelineError("core did not emit a protocol header")
    if header.run_id != expected_run_id:
        raise PipelineError("core run_id disagrees with the requested run")
    if header.core_version != __version__:
        raise PipelineError(
            "core version {} is incompatible with Python {}".format(
                header.core_version, __version__
            )
        )
    validate_header_projection(prepared, header)


def _expected_skipped_fragment_count(
    config: Mapping[str, object]
) -> Optional[int]:
    coverage = _mapping(config, "coverage")
    if coverage["kind"] == "profile" and config["technology"] == "WGBS":
        return None
    fragments = _mapping(config, "fragments")
    if config["technology"] == "WGBS" and _uses_variable_insert(fragments):
        return None
    if (
        config["technology"] == "TBS"
        and _mapping(config, "tbs")["fragment_center_stddev"] > 0
    ):
        return None
    return 0


def _require_released_capabilities(config: Mapping[str, object]) -> None:
    """Mirror the released capability gate before starting the core."""
    technology = config["technology"]
    if technology not in ("WGBS", "RRBS", "TBS"):
        raise PipelineError("technology is outside the current pipeline contract")
    inputs = _mapping(config, "inputs")
    supported_inputs = {"vcf", "cgmap", "bed_methyl", "asm", "asm_bed"}
    unsupported_inputs = sorted(set(inputs) - supported_inputs)
    if unsupported_inputs:
        raise PipelineError(
            "the current pipeline does not support input(s): {}".format(
                ", ".join(unsupported_inputs)
            )
        )
    if ("asm" in inputs or "asm_bed" in inputs) and "vcf" not in inputs:
        raise PipelineError("ASM v1 generation requires a VCF input")
    mutation = _mapping(config, "mutation")
    coverage = _mapping(config, "coverage")
    if mutation["rate"] != 0 and "vcf" in inputs:
        raise PipelineError(
            "VCF and de novo mutation generation are mutually exclusive"
        )
    if coverage["kind"] == "profile":
        if "artifact" in coverage:
            if technology != "WGBS":
                raise PipelineError(
                    "artifact-backed profile coverage currently supports WGBS only"
                )
            artifact = _mapping(coverage, "artifact")
            if (
                artifact["format"] != WGBS_GC_PROFILE_FORMAT
                or artifact["version"] != WGBS_GC_PROFILE_VERSION
            ):
                raise PipelineError(
                    "unsupported WGBS coverage profile format or version"
                )
            if "vcf" in inputs or mutation["rate"] != 0:
                raise PipelineError(
                    "target GC profile requires reference-only WGBS"
                )
        elif technology == "RRBS":
            rrbs = _mapping(config, "rrbs")
            if "candidate_bed" not in rrbs:
                raise PipelineError(
                    "RRBS profile coverage requires a candidate BED"
                )
        else:
            raise PipelineError(
                "profile coverage without an artifact currently supports RRBS only"
            )
    elif coverage["kind"] == "target-score":
        if technology != "TBS":
            raise PipelineError("target-score coverage currently supports TBS only")
    elif coverage["kind"] != "uniform":
        raise PipelineError("coverage kind is outside the released contract")

    fragments = _mapping(config, "fragments")
    has_depth = "depth" in fragments
    has_read_pairs = "read_pairs" in fragments
    if has_depth == has_read_pairs:
        raise PipelineError(
            "exactly one of fragments.depth and fragments.read_pairs is required"
        )
    if has_depth and technology != "WGBS":
        raise PipelineError("depth-to-count v1 currently supports WGBS only")
    variable_insert = _uses_variable_insert(fragments)
    if technology == "TBS" and variable_insert:
        raise PipelineError("the TBS baseline requires one fixed insert length")
    if technology == "TBS":
        tbs = _mapping(config, "tbs")
        if tbs["fragment_center_stddev"] < 0:
            raise PipelineError(
                "TBS fragment_center_stddev must be non-negative"
            )

    methylation = _mapping(config, "methylation")
    if methylation["cgmap_pool"] and not (
        "cgmap" in inputs or "bed_methyl" in inputs
    ):
        raise PipelineError(
            "cgmap_pool=true requires a CGmap or bedMethyl input"
        )
    if (
        "vcf" in inputs or mutation["rate"] != 0
    ) and not methylation["update_variant_boundaries"]:
        raise PipelineError(
            "variant generation requires update_variant_boundaries=true"
        )

    sequencing = _mapping(config, "sequencing")
    quality = _mapping(sequencing, "quality")
    if quality["kind"] == "markov":
        _require_model_contract(
            quality,
            "quality",
            QUALITY_MARKOV_FORMAT,
            QUALITY_MARKOV_VERSION,
        )
    elif quality["kind"] != "uniform":
        raise PipelineError("quality kind is outside the released contract")
    error = _mapping(sequencing, "error")
    if error["kind"] == "quality-confusion":
        _require_model_contract(
            error,
            "error",
            QUALITY_CONFUSION_FORMAT,
            QUALITY_CONFUSION_VERSION,
        )
    elif error["kind"] != "uniform":
        raise PipelineError("error kind is outside the released contract")

def _mapping(parent: Mapping[str, object], field: str) -> Mapping[str, object]:
    value = parent[field]
    if not isinstance(value, Mapping):
        raise PipelineError("{} must be an object".format(field))
    return value


def _require_model_contract(
    declaration: Mapping[str, object],
    label: str,
    expected_format: str,
    expected_version: str,
) -> None:
    artifact = _mapping(declaration, "artifact")
    if (
        artifact["format"] != expected_format
        or artifact["version"] != expected_version
    ):
        raise PipelineError(
            "unsupported {} model format or version".format(label)
        )


def _build_postprocess_config(prepared: PreparedRun) -> PostprocessConfig:
    normalized = prepared.config.normalized
    master_seed = prepared.config.master_seed
    if master_seed is None:
        raise PipelineError("prepared config must contain a materialized seed")
    sequencing = _mapping(normalized, "sequencing")
    quality_declaration = _mapping(sequencing, "quality")
    error_declaration = _mapping(sequencing, "error")

    try:
        if quality_declaration["kind"] == "uniform":
            quality = UniformQuality(quality_declaration["phred"])
        else:
            quality = parse_quality_markov(
                _snapshot_model_artifact(
                    prepared,
                    "model.quality",
                    _mapping(quality_declaration, "artifact"),
                )
            )

        if error_declaration["kind"] == "uniform":
            error = UniformError(error_declaration["rate"])
        else:
            error = parse_quality_confusion(
                _snapshot_model_artifact(
                    prepared,
                    "model.error",
                    _mapping(error_declaration, "artifact"),
                )
            )
        _require_model_compatibility(quality, error)
        return PostprocessConfig(
            master_seed=master_seed,
            directional=sequencing["directional"],
            conversion_rate=sequencing["conversion_rate"],
            quality=quality,
            error=error,
        )
    except (PreparationError, SequencingModelError, PostprocessError) as error_value:
        raise PipelineError(
            "cannot load sequencing model: {}".format(error_value)
        ) from error_value


def _snapshot_model_artifact(
    prepared: PreparedRun,
    role: str,
    artifact: Mapping[str, object],
) -> bytes:
    try:
        identity = prepared.file_for_role(role)
    except KeyError as error:
        raise PipelineError("prepared run is missing {}".format(role)) from error
    expected_path = Path(str(artifact["path"]))
    expected_sha256 = str(artifact["sha256"])
    if (
        identity.role != role
        or identity.path != expected_path
        or identity.sha256 != expected_sha256
        or identity.declared_sha256 != expected_sha256
    ):
        raise PipelineError(
            "prepared {} identity disagrees with normalized config".format(role)
        )
    return snapshot_prepared_file(identity, maximum_size=MAX_MODEL_BYTES)


def _require_model_compatibility(quality: object, error: object) -> None:
    if not isinstance(error, QualityConfusionModel):
        return
    available = set(error.quality_scores)
    if isinstance(quality, UniformQuality):
        required = {quality.phred}
    elif isinstance(quality, QualityMarkovModel):
        required = set(quality.quality_scores)
    else:
        raise PipelineError("quality policy is outside the released contract")
    missing = sorted(required - available)
    if missing:
        raise PipelineError(
            "error model is missing quality score(s): {}".format(
                ", ".join(str(value) for value in missing)
            )
        )


def _uses_variable_insert(fragments: Mapping[str, object]) -> bool:
    return (
        fragments["insert_min"] != fragments["insert_mean"]
        or fragments["insert_min"] != fragments["insert_max"]
        or fragments["insert_stddev"] != 0
    )


__all__ = [
    "PipelineError",
    "RunResult",
    "resolve_core_executable",
    "run_document",
    "run_prepared",
]
