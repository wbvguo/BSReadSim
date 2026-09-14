"""One simulation run: htsim fragments -> reads (in worker processes) -> files and manifest.

htsim owns the reference, variants, methylation profile, and fragment
sampling, and streams fragment batches. Python realizes methylation states,
bisulfite conversion, qualities, and errors batch by batch (in worker
processes when threads allow), writes the outputs in fragment order, and
publishes them together with the manifest.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator, Sequence
from concurrent.futures import ProcessPoolExecutor
from concurrent.futures.process import BrokenProcessPool
from contextlib import contextmanager
from dataclasses import dataclass
import multiprocessing
import os
from pathlib import Path
import queue
import shutil
import threading
import uuid

from .errors import BSReadSimError
from .batch import FragmentBatch
from .htsim import HtsimCore, HtsimProcess
from .input import describe_inputs
from .output import FastqFormat, RunFiles, SamFormat, build_manifest
from .reads import ReadSimulator
from .settings import Settings


@dataclass(frozen=True)
class RunResult:
    manifest_path: Path
    manifest: dict


@dataclass(frozen=True)
class ThreadPlan:
    """How the requested threads are split between htsim and Python."""

    workers: int  # Python read-processing processes
    core: int  # htsim threads
    bam: int  # BAM compression threads


def thread_plan(threads: int, output_format: str) -> ThreadPlan:
    """Split ``--threads`` between read processing, htsim, and BAM compression.

    Python read processing costs the most. htsim needs about a quarter of the
    threads, or half for plain FASTQ, whose reads are cheap to format; BAM
    output gives just under half to BGZF compression, whose threads are busy
    about a third of the time (at the default --gzip-level 4). The main
    process, which relays every batch, and htsim's SAM parsing come on top
    (2M WGBS reads as BAM: 4.9 busy CPUs at 4 threads, 7.6 at 8, 12.4 at 16).
    """
    if output_format == "fastq":
        core = max(1, (threads + 1) // 2)
    else:
        core = max(1, threads // 4)
    bam = min(64, max(0, threads // 2 - 1)) if output_format == "bam" else 0  # htsim allows 64
    return ThreadPlan(workers=max(1, threads - core - bam // 3), core=core, bam=bam)


def run(settings: Settings, argv: Sequence[str],
        core: str | os.PathLike | None = None) -> RunResult:
    """Simulate the reads of ``settings``; ``argv`` is the command as received."""
    core = HtsimCore.find(core)
    simulator = ReadSimulator.from_settings(settings)
    inputs = describe_inputs(settings)
    plan = thread_plan(settings.threads, settings.format)
    run_id = str(uuid.uuid4())
    prefix = settings.prefix

    with RunFiles(settings.output, prefix, paired_end=settings.paired_end,
                  format=settings.format) as files:
        methdb_output = None
        if settings.save_methdb:
            methdb_output = files.truth_path(prefix + ".methdb")
            if settings.methdb is not None:
                shutil.copyfile(settings.methdb, methdb_output)
                methdb_output = None
            files.add_truth("truth.methdb", prefix + ".methdb")
        if settings.save_vcf:
            vcf = files.truth_path(prefix + ".variants.vcf.gz")
            vcf.unlink(missing_ok=True)
            core.build(settings, "variants", vcf)
            files.add_truth("truth.vcf", vcf.name)

        with core.stream(settings, details=settings.alignments, threads=plan.core,
                         site_kmers=simulator.methylation.site_kmers,
                         methdb_output=methdb_output) as htsim:
            header = htsim.header
            if settings.alignments:
                output = SamFormat.of_run(header, settings, run_id)
                files.start_alignments(output.header(),
                                       core.sam_to_bam(settings.gzip_level, plan.bam))
            else:
                output = FastqFormat(
                    contig_names=tuple(contig.name for contig in header.contigs),
                    gzip_level=settings.gzip_level if settings.format == "fastq.gz" else None)
            _write_batches(htsim, BatchJob(simulator, output), files, plan.workers)
        outputs = files.finish()
        manifest = build_manifest(settings, argv, header, htsim.summary, inputs, outputs,
                                  run_id=run_id, methylation_model=simulator.methylation.name)
        files.publish(manifest)
    return RunResult(files.manifest_path, manifest)


# ---------------------------------------------------------------- batches

@dataclass(frozen=True)
class FormattedBatch:
    fragment_count: int
    template_base_count: int
    site_count: int
    records: dict[str, bytes]  # by output file role


@dataclass(frozen=True)
class BatchJob:
    """What becomes of every batch of a run (in a worker): its reads, formatted."""

    reads: ReadSimulator
    output: FastqFormat | SamFormat

    def __call__(self, payload: bytes) -> FormattedBatch:
        batch = FragmentBatch.decode(payload)
        return FormattedBatch(
            fragment_count=batch.size, template_base_count=len(batch.bases),
            site_count=len(batch.site_positions),
            records=self.output.records(batch, self.reads.simulate(batch)))


def _write_batches(htsim: HtsimProcess, job: BatchJob, files: RunFiles, workers: int) -> None:
    """Process every streamed batch and write it in fragment order."""
    totals = [0, 0, 0]
    for batch in _formatted_batches(htsim.batches(), job, workers):
        files.write(batch.fragment_count, batch.records)
        totals[0] += batch.fragment_count
        totals[1] += batch.template_base_count
        totals[2] += batch.site_count
    summary = htsim.summary
    if totals != [summary.fragment_count, summary.template_base_count,
                  summary.methylation_site_count]:
        raise BSReadSimError("processed fragments disagree with the htsim summary")


def _formatted_batches(payloads: Iterable[bytes], job: Callable[[bytes], FormattedBatch],
                       workers: int) -> Iterator[FormattedBatch]:
    """Formatted batches in input order, using up to ``workers`` processes.

    The processes compute. In this process, a thread reads the stream and submits each
    batch while the caller writes the results; it runs at most two batches a worker ahead.
    """
    if workers <= 1:
        for payload in payloads:
            yield job(payload)
        return
    submitted: queue.Queue = queue.Queue()  # futures in batch order, then None or an error
    room = threading.Semaphore(2 * workers)
    stop = threading.Event()

    def submit_all(pool: ProcessPoolExecutor) -> None:
        try:
            for payload in payloads:
                room.acquire()
                if stop.is_set():
                    return
                submitted.put(pool.submit(_format_in_worker, payload))
            submitted.put(None)
        except BaseException as error:  # raised by the caller's loop below
            submitted.put(error)

    # spawn: this process already runs threads that must not be forked
    with _single_threaded_blas(), ProcessPoolExecutor(
            workers, mp_context=multiprocessing.get_context("spawn"),
            initializer=_start_worker, initargs=(job,)) as pool:
        reader = threading.Thread(target=submit_all, args=(pool,), name="bsreadsim-reader",
                                  daemon=True)
        reader.start()
        try:
            while (item := submitted.get()) is not None:
                if isinstance(item, BaseException):
                    raise item
                try:
                    batch = item.result()
                except BrokenProcessPool as error:
                    raise BSReadSimError("a worker process stopped unexpectedly") from error
                room.release()
                yield batch
        finally:
            stop.set()
            room.release(2 * workers)  # a reader waiting for room sees stop
            reader.join()
            pool.shutdown(cancel_futures=True)


_job: Callable[[bytes], FormattedBatch] | None = None  # in a worker process: its job


def _start_worker(job: Callable[[bytes], FormattedBatch]) -> None:
    global _job  # a worker gets the job once: it may hold a large model
    _job = job


def _format_in_worker(payload: bytes) -> FormattedBatch:
    return _job(payload)


@contextmanager
def _single_threaded_blas() -> Iterator[None]:
    """Workers do no linear algebra: keep OpenBLAS from starting a thread per CPU in each."""
    name = "OPENBLAS_NUM_THREADS"
    if name in os.environ:
        yield
        return
    os.environ[name] = "1"
    try:
        yield
    finally:
        del os.environ[name]
