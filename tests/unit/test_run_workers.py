"""Tests for the bounded local/shared batch processing boundary."""

from dataclasses import replace
import gzip
import multiprocessing
import unittest


from tests.helpers.process_support import UniformProcessConfig
from bsreadsim.output.bam import format_sam_batch
from bsreadsim.process.fragment import (
    decode_fragments,
    process_fragment_batch,
)
from bsreadsim.run.workers import (
    LocalBatchSlot,
    SharedBatchSlot,
    initialize_process_worker,
    process_payload_batch,
    process_shared_batch,
)
from bsreadsim.htsim.protocol import (
    ProtocolError,
)
from tests.helpers.stream_support import encode_stream, read_stream
from tests.unit.test_protocol import (
    make_full_annotation_batch,
    make_header,
)


def decoded_payload(*, details: bool):
    header = make_header(details=details)
    source = make_full_annotation_batch()
    batch = source if details else replace(source, details=None)
    decoded = read_stream(encode_stream(header, (batch,))).batches[0]
    return header, decoded


def run_inline(
    slot,
    header,
    decoded,
    *,
    include_details: bool,
    include_alignment: bool = False,
    include_fragment_summary: bool = False,
    fastq_gzip_level: int | None = None,
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
            config=UniformProcessConfig(
                master_seed=7,
                conversion_rate=1.0,
                quality_phred=30,
                error_rate=0.0,
            ),
            paired_end=False,
            include_details=include_details,
            include_alignment=include_alignment,
            include_fragment_summary=include_fragment_summary,
            fastq_gzip_level=fastq_gzip_level,
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
    def test_fastq_batch_can_be_encoded_as_a_deterministic_gzip_member(self) -> None:
        header, decoded = decoded_payload(details=False)
        raw_slot = LocalBatchSlot(0)
        compressed_slot = LocalBatchSlot(1)
        try:
            _, _, raw = run_inline(
                raw_slot, header, decoded, include_details=False
            )
            _, _, compressed = run_inline(
                compressed_slot,
                header,
                decoded,
                include_details=False,
                fastq_gzip_level=6,
            )
            self.assertFalse(raw.precompressed_fastq)
            self.assertTrue(compressed.precompressed_fastq)
            self.assertEqual(compressed.record_lengths, raw.record_lengths)
            self.assertEqual(
                gzip.decompress(region_bytes(compressed_slot, compressed.read1)),
                region_bytes(raw_slot, raw.read1),
            )
        finally:
            raw_slot.close()
            compressed_slot.close()

    def test_local_and_shared_slots_use_the_same_processing_core(self) -> None:
        header, decoded = decoded_payload(details=True)
        local = LocalBatchSlot(0)
        shared = SharedBatchSlot(1)
        try:
            _, _, local_result = run_inline(
                local, header, decoded, include_details=True
            )
            _, _, shared_result = run_inline(
                shared, header, decoded, include_details=True
            )
            self.assertEqual(local_result, shared_result)
            for name in ("read1", "read2", "alignment"):
                self.assertEqual(
                    region_bytes(local, getattr(local_result, name)),
                    region_bytes(shared, getattr(shared_result, name)),
                )
        finally:
            local.close()
            shared.close()

    def test_fastq_only_common_columns_match_full_annotation_fastq(self) -> None:
        full_header, full_batch = decoded_payload(details=True)
        without_annotations_header, without_annotations_batch = decoded_payload(details=False)
        full_slot = LocalBatchSlot(0)
        fastq_slot = LocalBatchSlot(1)
        try:
            _, _, full = run_inline(
                full_slot, full_header, full_batch, include_details=True
            )
            _, _, fastq = run_inline(
                fastq_slot, without_annotations_header, without_annotations_batch, include_details=False
            )
            self.assertEqual(
                region_bytes(full_slot, full.read1),
                region_bytes(fastq_slot, fastq.read1),
            )
            self.assertTrue(all(len(lengths) == 2 for lengths in fastq.record_lengths))
            self.assertEqual(full.summary, fastq.summary)
        finally:
            full_slot.close()
            fastq_slot.close()

    def test_bam_only_worker_retains_full_annotation_projection(self) -> None:
        header, decoded = decoded_payload(details=True)
        slot = LocalBatchSlot(0)
        try:
            _, _, result = run_inline(
                slot,
                header,
                decoded,
                include_details=True,
                include_alignment=True,
            )
            alignment = region_bytes(slot, result.alignment)
            self.assertIsNotNone(alignment)
            self.assertEqual(
                len(result.alignment_record_lengths),
                result.summary.mate_count,
            )
            self.assertEqual(
                sum(result.alignment_record_lengths),
                len(alignment),
            )
            self.assertIsNone(result.read1)
            self.assertIsNone(result.read2)
            self.assertTrue(
                all(lengths == (0, 0) for lengths in result.record_lengths)
            )
            self.assertEqual(len(alignment.splitlines()), result.summary.mate_count)
        finally:
            slot.close()

    def test_fragment_summary_columnar_output_matches_typed_reference(self) -> None:
        header, decoded = decoded_payload(details=True)
        slot = LocalBatchSlot(0)
        baseline_slot = LocalBatchSlot(1)
        config = UniformProcessConfig(
            master_seed=7,
            conversion_rate=1.0,
            quality_phred=30,
            error_rate=0.0,
        )
        try:
            _, _, result = run_inline(
                slot,
                header,
                decoded,
                include_details=True,
                include_alignment=True,
                include_fragment_summary=True,
            )
            observed = region_bytes(slot, result.alignment)
            _, _, baseline_result = run_inline(
                baseline_slot,
                header,
                decoded,
                include_details=True,
                include_alignment=True,
            )
            baseline = region_bytes(baseline_slot, baseline_result.alignment)
            fragments = decode_fragments(decoded, header)
            processed = process_fragment_batch(
                fragments,
                tuple(
                    header.contigs[fragment.contig_index].name
                    for fragment in fragments
                ),
                config,
                compact_base_states=True,
                include_details=True,
            )
            expected, expected_lengths = format_sam_batch(
                processed,
                paired_end=False,
                read_group_id=header.run_id,
                contig_lengths=tuple(
                    header.contigs[fragment.contig_index].length
                    for fragment in fragments
                ),
                include_fragment_summary=True,
            )
            self.assertEqual(result.summary.mate_count, len(expected_lengths))
            observed_lines = observed.splitlines()
            expected_lines = expected.splitlines()
            observed_zf = tuple(
                next(field for field in line.split(b"\t") if field.startswith(b"zf:B:S,"))
                for line in observed_lines
            )
            expected_zf = tuple(
                next(field for field in line.split(b"\t") if field.startswith(b"zf:B:S,"))
                for line in expected_lines
            )
            self.assertEqual(observed_zf, expected_zf)
            without_zf = b"".join(
                b"\t".join(
                    field
                    for field in line.split(b"\t")
                    if not field.startswith(b"zf:B:S,")
                )
                + b"\n"
                for line in observed_lines
            )
            self.assertEqual(without_zf, baseline)
        finally:
            slot.close()
            baseline_slot.close()

    def test_spawned_worker_reads_one_authenticated_shared_batch(self) -> None:
        header, decoded = decoded_payload(details=False)
        slot = SharedBatchSlot(0)
        pool = None
        completed = False
        config = UniformProcessConfig(master_seed=7)
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
        finally:
            if pool is not None and not completed:
                pool.terminate()
                pool.join()
            slot.close()

    def test_worker_revalidates_payload_before_processing(self) -> None:
        header, decoded = decoded_payload(details=False)
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
                        config=UniformProcessConfig(master_seed=7),
                        paired_end=False,
                        include_details=False,
                    )
            finally:
                buffer.release()
        finally:
            slot.close()

    def test_annotation_policy_must_match_the_protocol_header(self) -> None:
        header, decoded = decoded_payload(details=False)
        slot = LocalBatchSlot(0)
        try:
            with self.assertRaisesRegex(Exception, "details|Details"):
                run_inline(slot, header, decoded, include_details=True)
        finally:
            slot.close()


if __name__ == "__main__":
    unittest.main()
