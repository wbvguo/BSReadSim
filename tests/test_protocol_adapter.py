"""Scientific and FASTQ equivalence tests for protocol adapters."""

from dataclasses import replace
from pathlib import Path
import struct
import sys
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim.numpy_postprocess import (  # noqa: E402
    format_common_fastq_batch,
    process_common_fastq_batch,
    supports_common_postprocess,
)
from bsreadsim.output import (  # noqa: E402
    _fastq_record_fields,
    _format_fragment_records_trusted,
)
from bsreadsim.postprocess import (  # noqa: E402
    PostprocessError,
    UniformPostprocessConfig,
    process_fragment,
)
from bsreadsim.protocol import (  # noqa: E402
    TruthMode,
    ProtocolError,
    decode_batch_payload,
    encode_stream,
    read_stream,
)
from bsreadsim.protocol_adapter import (  # noqa: E402
    _decode_common_numpy_batch,
    decode_fragments,
)
from tests.test_protocol import (  # noqa: E402
    make_full_truth_batch,
    make_header,
)


class ProtocolTypedAdapterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.header = make_header(truth=True)
        self.batch = read_stream(
            encode_stream(self.header, (make_full_truth_batch(),))
        ).batches[0]

    def test_full_truth_reconstructs_the_scientific_model(self) -> None:
        fragments = decode_fragments(self.batch, self.header)
        self.assertEqual(
            tuple(fragments[1].reference_positions),
            (20, 21, -1, 22, 23),
        )
        self.assertEqual(
            tuple(fragments[1].base_event_ids),
            (0xFFFFFFFF, 0xFFFFFFFF, 2, 0xFFFFFFFF, 0xFFFFFFFF),
        )
        self.assertEqual(fragments[1].methylation_sites[1].reference_pos, -1)
        self.assertEqual(
            tuple(
                (item.read_offset, item.site_index)
                for item in fragments[1].mates[0].site_refs
            ),
            ((1, 0), (2, 1)),
        )
        self.assertEqual(fragments[2].variant_events[0].alt_bases, b"")

    def test_fastq_only_columns_build_a_compact_processing_model(self) -> None:
        no_truth_header = replace(self.header, truth_columns=TruthMode.NONE)
        no_truth_batch = replace(make_full_truth_batch(), truth=None)
        decoded = read_stream(
            encode_stream(no_truth_header, (no_truth_batch,))
        ).batches[0]
        fragments = decode_fragments(decoded, no_truth_header)
        self.assertEqual(len(fragments), decoded.fragment_count)
        self.assertEqual(
            fragments[1].reference_positions,
            (-1,) * len(fragments[1].template_bases),
        )
        self.assertEqual(fragments[1].variant_events, ())
        self.assertTrue(
            all(site.reference_pos == -1 for site in fragments[1].methylation_sites)
        )

    def test_authenticated_payload_worker_entry_revalidates_semantics(self) -> None:
        decoded = decode_batch_payload(
            bytes(self.batch.raw_payload),
            self.header,
            expected_first_ordinal=0,
        )
        self.assertEqual(tuple(decoded.template_bases), tuple(self.batch.template_bases))
        corrupted = bytearray(self.batch.raw_payload)
        struct.pack_into("<I", corrupted, 44, 0)
        with self.assertRaisesRegex(ProtocolError, "reference envelope"):
            decode_batch_payload(bytes(corrupted), self.header)


class ProtocolCommonColumnTests(unittest.TestCase):
    def test_direct_common_path_matches_typed_fastq_bytes(self) -> None:
        full_header = make_header(truth=True)
        source_batch = make_full_truth_batch()
        full_batch = read_stream(
            encode_stream(full_header, (source_batch,))
        ).batches[0]
        fragments = decode_fragments(full_batch, full_header)

        no_truth_header = replace(full_header, truth_columns=TruthMode.NONE)
        no_truth_batch = read_stream(
            encode_stream(
                no_truth_header,
                (replace(source_batch, truth=None),),
            )
        ).batches[0]
        config = UniformPostprocessConfig(
            master_seed=0x123456789ABCDEF0,
            directional=False,
            conversion_rate=0.73,
            quality_phred=31,
            error_rate=0.19,
        )
        self.assertTrue(
            supports_common_postprocess(config, include_truth=False)
        )
        direct = process_common_fastq_batch(
            _decode_common_numpy_batch(no_truth_batch, no_truth_header),
            config,
        )
        formatted_direct = format_common_fastq_batch(
            _decode_common_numpy_batch(no_truth_batch, no_truth_header),
            config,
            False,
        )
        typed = tuple(
            process_fragment(
                fragment,
                full_header.contigs[fragment.contig_index].name,
                config,
                include_truth=False,
            )
            for fragment in fragments
        )
        direct_fastq = tuple(
            (
                _fastq_record_fields(
                    fragment.contig_name,
                    fragment.reference_start,
                    fragment.reference_end,
                    fragment.fragment_ordinal,
                    fragment.mates[0].mate_index,
                    fragment.mates[0].sequence,
                    fragment.mates[0].quality,
                ).encode("utf-8"),
                None,
            )
            for fragment in direct
        )
        typed_fastq = tuple(
            _format_fragment_records_trusted(
                fragment, paired_end=False, include_truth=False
            )[:2]
            for fragment in typed
        )
        self.assertEqual(direct_fastq, typed_fastq)
        self.assertEqual(
            formatted_direct.read1,
            b"".join(read1 for read1, _ in direct_fastq),
        )
        self.assertIsNone(formatted_direct.read2)
        self.assertEqual(
            formatted_direct.record_lengths,
            tuple((len(read1), 0, 0) for read1, _ in direct_fastq),
        )
        self.assertEqual(
            direct[1].mates[0].sequence,
            typed[1].mates[0].sequence,
        )

    def test_common_adapter_rejects_full_truth(self) -> None:
        header = make_header(truth=True)
        batch = read_stream(
            encode_stream(header, (make_full_truth_batch(),))
        ).batches[0]
        with self.assertRaisesRegex(PostprocessError, "TruthMode.NONE"):
            _decode_common_numpy_batch(batch, header)


if __name__ == "__main__":
    unittest.main()
