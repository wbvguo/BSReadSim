"""Tests for the bounded local/shared batch processing boundary."""

from dataclasses import replace
import multiprocessing
from pathlib import Path
import sys
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim.postprocess import UniformPostprocessConfig  # noqa: E402
from bsreadsim.process_pool import (  # noqa: E402
    LocalBatchSlot,
    SharedBatchSlot,
    initialize_process_worker,
    process_payload_batch,
    process_shared_batch,
)
from bsreadsim.protocol import (  # noqa: E402
    ProtocolError,
    TruthMode,
    encode_stream,
    read_stream,
)
from tests.test_protocol import (  # noqa: E402
    make_full_truth_batch,
    make_header,
)


def decoded_payload(*, truth: bool):
    header = make_header(truth=truth)
    source = make_full_truth_batch()
    batch = source if truth else replace(source, truth=None)
    decoded = read_stream(encode_stream(header, (batch,))).batches[0]
    return header, decoded


def run_inline(
    slot,
    header,
    decoded,
    *,
    include_truth: bool,
    emit_truth=None,
    include_alignment: bool = False,
):
    descriptors, output_offset = slot.prepare(
        ((decoded.first_fragment_ordinal, bytes(decoded.raw_payload)),)
    )
    buffer = slot.buffer_view()
    try:
        result = process_payload_batch(
            buffer,
            descriptors,
            output_offset,
            header=header,
            config=UniformPostprocessConfig(
                master_seed=7,
                conversion_rate=1.0,
                quality_phred=30,
                error_rate=0.0,
            ),
            paired_end=False,
            include_truth=include_truth,
            emit_truth=emit_truth,
            include_alignment=include_alignment,
        )
    finally:
        buffer.release()
    return descriptors, output_offset, result


def region_bytes(slot, region):
    if region is None:
        return None
    view = slot.view(region)
    try:
        return bytes(view)
    finally:
        view.release()


class ProcessPoolTests(unittest.TestCase):
    def test_local_and_shared_slots_use_the_same_processing_core(self) -> None:
        header, decoded = decoded_payload(truth=True)
        local = LocalBatchSlot(0)
        shared = SharedBatchSlot(1)
        try:
            _, _, local_result = run_inline(
                local, header, decoded, include_truth=True
            )
            _, _, shared_result = run_inline(
                shared, header, decoded, include_truth=True
            )
            self.assertEqual(local_result, shared_result)
            for name in ("read1", "read2", "truth"):
                self.assertEqual(
                    region_bytes(local, getattr(local_result, name)),
                    region_bytes(shared, getattr(shared_result, name)),
                )
        finally:
            local.close()
            shared.close()

    def test_fastq_only_common_columns_match_full_truth_fastq(self) -> None:
        full_header, full_batch = decoded_payload(truth=True)
        no_truth_header, no_truth_batch = decoded_payload(truth=False)
        full_slot = LocalBatchSlot(0)
        fastq_slot = LocalBatchSlot(1)
        try:
            _, _, full = run_inline(
                full_slot, full_header, full_batch, include_truth=True
            )
            _, _, fastq = run_inline(
                fastq_slot, no_truth_header, no_truth_batch, include_truth=False
            )
            self.assertEqual(
                region_bytes(full_slot, full.read1),
                region_bytes(fastq_slot, fastq.read1),
            )
            self.assertIsNotNone(full.truth)
            self.assertIsNone(fastq.truth)
            self.assertTrue(all(lengths[2] == 0 for lengths in fastq.record_lengths))
            self.assertEqual(full.summary, fastq.summary)
        finally:
            full_slot.close()
            fastq_slot.close()

    def test_bam_only_worker_retains_projection_without_json_truth(self) -> None:
        header, decoded = decoded_payload(truth=True)
        slot = LocalBatchSlot(0)
        try:
            _, _, result = run_inline(
                slot,
                header,
                decoded,
                include_truth=True,
                emit_truth=False,
                include_alignment=True,
            )
            alignment = region_bytes(slot, result.alignment)
            self.assertIsNone(result.truth)
            self.assertIsNotNone(alignment)
            self.assertEqual(
                len(result.alignment_record_lengths),
                result.summary.mate_count,
            )
            self.assertEqual(
                sum(result.alignment_record_lengths),
                len(alignment),
            )
            self.assertTrue(all(lengths[2] == 0 for lengths in result.record_lengths))
            self.assertEqual(len(alignment.splitlines()), result.summary.mate_count)
        finally:
            slot.close()

    def test_spawned_worker_reads_one_authenticated_shared_batch(self) -> None:
        header, decoded = decoded_payload(truth=False)
        slot = SharedBatchSlot(0)
        pool = None
        completed = False
        config = UniformPostprocessConfig(master_seed=7)
        try:
            descriptors, output_offset = slot.prepare(
                ((decoded.first_fragment_ordinal, bytes(decoded.raw_payload)),)
            )
            pool = multiprocessing.get_context("spawn").Pool(
                processes=1,
                initializer=initialize_process_worker,
                initargs=(header, config, False, False),
            )
            result = pool.apply_async(
                process_shared_batch,
                (slot.name, descriptors, output_offset),
            ).get(timeout=10)
            pool.close()
            pool.join()
            completed = True

            self.assertFalse(result.needs_resize)
            self.assertEqual(result.summary.fragment_count, decoded.fragment_count)
            self.assertIsNotNone(result.read1)
            self.assertIsNone(result.read2)
            self.assertIsNone(result.truth)
        finally:
            if pool is not None and not completed:
                pool.terminate()
                pool.join()
            slot.close()

    def test_worker_revalidates_payload_before_postprocessing(self) -> None:
        header, decoded = decoded_payload(truth=False)
        slot = LocalBatchSlot(0)
        try:
            descriptors, output_offset = slot.prepare(
                ((decoded.first_fragment_ordinal, bytes(decoded.raw_payload)),)
            )
            invalid = (replace(descriptors[0], length=1),)
            buffer = slot.buffer_view()
            try:
                with self.assertRaisesRegex(ProtocolError, "truncated"):
                    process_payload_batch(
                        buffer,
                        invalid,
                        output_offset,
                        header=header,
                        config=UniformPostprocessConfig(master_seed=7),
                        paired_end=False,
                        include_truth=False,
                    )
            finally:
                buffer.release()
        finally:
            slot.close()

    def test_truth_policy_must_match_the_protocol_header(self) -> None:
        header, decoded = decoded_payload(truth=False)
        slot = LocalBatchSlot(0)
        try:
            with self.assertRaisesRegex(Exception, "truth|Truth"):
                run_inline(slot, header, decoded, include_truth=True)
        finally:
            slot.close()


if __name__ == "__main__":
    unittest.main()
