"""Scientific and FASTQ equivalence tests for protocol adapters."""

from dataclasses import replace
from pathlib import Path
import struct
import sys
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from process_support import UniformProcessConfig, process_fragment  # noqa: E402

from bsreadsim.process.fragment import (  # noqa: E402
    encode_fastq_batch,
    supports_common_processing,
)
from bsreadsim.output.fastq import (  # noqa: E402
    format_fragment_records_trusted,
)
from bsreadsim.process import (  # noqa: E402
    ProcessError,
)
from bsreadsim.native.protocol import (  # noqa: E402
    ProtocolError,
    decode_batch_payload,
)
from stream_support import encode_stream, read_stream  # noqa: E402
from bsreadsim.process.fragment import (  # noqa: E402
    decode_common_numpy_batch,
    decode_fragments,
)
from tests.test_protocol import (  # noqa: E402
    make_full_annotation_batch,
    make_header,
)


class ProtocolTypedAdapterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.header = make_header(details=True)
        self.batch = read_stream(
            encode_stream(self.header, (make_full_annotation_batch(),))
        ).batches[0]

    def test_full_annotation_reconstructs_the_scientific_model(self) -> None:
        fragments = decode_fragments(self.batch, self.header)
        self.assertEqual(
            tuple(fragments[1].reference_positions),
            (20, 21, -1, 22, 23),
        )
        self.assertEqual(
            tuple(fragments[1].base_variant_indices),
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
        self.assertEqual(fragments[2].variants[0].alt_bases, b"")

    def test_fastq_only_columns_build_a_compact_processing_model(self) -> None:
        without_annotations_header = replace(self.header, has_details=False)
        without_annotations_batch = replace(make_full_annotation_batch(), details=None)
        decoded = read_stream(
            encode_stream(without_annotations_header, (without_annotations_batch,))
        ).batches[0]
        fragments = decode_fragments(decoded, without_annotations_header)
        self.assertEqual(len(fragments), decoded.fragment_count)
        self.assertEqual(
            fragments[1].reference_positions,
            (-1,) * len(fragments[1].template_bases),
        )
        self.assertEqual(fragments[1].variants, ())
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
    def test_common_formatter_matches_typed_fastq_bytes(self) -> None:
        full_header = make_header(details=True)
        source_batch = make_full_annotation_batch()
        full_batch = read_stream(
            encode_stream(full_header, (source_batch,))
        ).batches[0]
        fragments = decode_fragments(full_batch, full_header)

        without_annotations_header = replace(full_header, has_details=False)
        without_annotations_batch = read_stream(
            encode_stream(
                without_annotations_header,
                (replace(source_batch, details=None),),
            )
        ).batches[0]
        config = UniformProcessConfig(
            master_seed=0x123456789ABCDEF0,
            directional=False,
            conversion_rate=0.73,
            quality_phred=31,
            error_rate=0.19,
        )
        self.assertTrue(
            supports_common_processing(config)
        )
        formatted_direct = encode_fastq_batch(
            decode_common_numpy_batch(without_annotations_batch, without_annotations_header),
            config,
            False,
        )
        typed = tuple(
            process_fragment(
                fragment,
                full_header.contigs[fragment.contig_index].name,
                config,
                include_details=False,
            )
            for fragment in fragments
        )
        typed_fastq = tuple(
            format_fragment_records_trusted(fragment, paired_end=False)[:2]
            for fragment in typed
        )
        self.assertEqual(
            formatted_direct.read1,
            b"".join(read1 for read1, _ in typed_fastq),
        )
        self.assertIsNone(formatted_direct.read2)
        self.assertEqual(
            formatted_direct.record_lengths,
            tuple((len(read1), 0) for read1, _ in typed_fastq),
        )

    def test_common_adapter_rejects_full_annotations(self) -> None:
        header = make_header(details=True)
        batch = read_stream(
            encode_stream(header, (make_full_annotation_batch(),))
        ).batches[0]
        decoded = decode_common_numpy_batch(batch, header)
        self.assertEqual(decoded.fragment_count, batch.fragment_count)
        self.assertEqual(decoded.mate_count, batch.mate_count)


if __name__ == "__main__":
    unittest.main()
