"""Contract tests for the columnar protocol stream."""

from dataclasses import replace
import hashlib
import io
from pathlib import Path
import struct
import sys
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim.protocol import (  # noqa: E402
    AmbiguityPolicy,
    BaseEncoding,
    CONFIG_SCHEMA_VERSION,
    CaptureStrand,
    Contig,
    CoreReportedError,
    ErrorFrame,
    FragmentBatch,
    Header,
    MethylationAllele,
    MethylationContext,
    MethylationSource,
    NO_REFERENCE_POSITION,
    ProtocolError,
    ProtocolReader,
    ProtocolWriter,
    RNG_CONTRACT,
    Technology,
    TruthColumns,
    TruthMode,
    VariantKind,
    encode_stream,
    read_stream,
)
from bsreadsim.protocol import crc32c  # noqa: E402


PREAMBLE_SIZE = struct.calcsize("<8sHHI")
FRAME_ENVELOPE = struct.Struct("<IBBHQ")
CRC = struct.Struct("<I")
FIXTURE_ROOT = REPOSITORY_ROOT / "tests/fixtures/protocol-v2"


def frozen_fixture(name: str) -> bytes:
    return bytes.fromhex((FIXTURE_ROOT / (name + ".hex")).read_text(encoding="ascii"))


def make_header(*, truth: bool) -> Header:
    return Header(
        run_id="00000000-0000-4000-8000-000000000002",
        core_version="2.0.0-alpha.1",
        config_schema_version=CONFIG_SCHEMA_VERSION,
        rng_contract=RNG_CONTRACT,
        master_seed=0x0123456789ABCDEF,
        normalized_config_sha256=hashlib.sha256(b"protocol-v2-config").digest(),
        technology=Technology.WGBS,
        truth_columns=TruthMode.FULL if truth else TruthMode.NONE,
        mates_per_fragment=1,
        base_encoding=BaseEncoding.ACGTN_U8,
        ambiguity_policy=AmbiguityPolicy.PRESERVE_N,
        read_length_r1=3,
        read_length_r2=0,
        contigs=(
            Contig(
                name="chrMini",
                length=100,
                reference_sha256=hashlib.sha256(b"ACGT" * 25).digest(),
            ),
        ),
    )


def make_no_truth_batch() -> FragmentBatch:
    return FragmentBatch(
        first_fragment_ordinal=0,
        contig_indices=(0, 0),
        reference_begins=(10, 20),
        reference_ends=(14, 24),
        template_offsets=(0, 4, 8),
        mate_offsets=(0, 1, 2),
        site_offsets=(0, 1, 2),
        mate_template_begins=(0, 1),
        mate_template_ends=(3, 4),
        site_template_offsets=(1, 1),
        site_probabilities=(0.25, 0.75),
        haplotypes=(0, 1),
        capture_strands=(CaptureStrand.UNKNOWN, CaptureStrand.FORWARD),
        mate_indices=(0, 0),
        mate_reverse_complements=(0, 0),
        site_contexts=(MethylationContext.CG_C, MethylationContext.CG_G),
        site_sources=(MethylationSource.BETA, MethylationSource.CGMAP),
        site_alleles=(
            MethylationAllele.SHARED,
            MethylationAllele.ALTERNATE_HAPLOTYPE,
        ),
        template_bases=(0, 1, 2, 3, 3, 2, 1, 0),
        truth=None,
    )


def make_full_truth_batch() -> FragmentBatch:
    truth = TruthColumns(
        projection_offsets=(0, 1, 3, 5),
        event_offsets=(0, 1, 2, 3),
        original_n_offsets=(0, 0, 0, 1),
        projection_template_begins=(0, 0, 3, 0, 2),
        projection_template_ends=(4, 2, 5, 2, 3),
        projection_reference_begins=(10, 20, 22, 30, 33),
        event_ids=(1, 2, 3),
        event_reference_begins=(11, 22, 32),
        event_reference_ends=(12, 22, 33),
        event_template_begins=(1, 2, 2),
        event_template_ends=(2, 3, 2),
        event_ref_offsets=(0, 1, 1, 2),
        event_alt_offsets=(0, 1, 2, 2),
        site_reference_positions=(12, 21, NO_REFERENCE_POSITION, 31),
        original_n_template_offsets=(0,),
        event_kinds=(VariantKind.SNV, VariantKind.INSERTION, VariantKind.DELETION),
        event_phased_haplotypes=(1, 1, 1),
        event_ref_bases=(1, 2),
        event_alt_bases=(3, 2),
    )
    return FragmentBatch(
        first_fragment_ordinal=0,
        contig_indices=(0, 0, 0),
        reference_begins=(10, 20, 30),
        reference_ends=(14, 24, 34),
        template_offsets=(0, 4, 9, 12),
        mate_offsets=(0, 1, 2, 3),
        site_offsets=(0, 1, 3, 4),
        mate_template_begins=(0, 0, 0),
        mate_template_ends=(3, 3, 3),
        site_template_offsets=(2, 1, 2, 1),
        site_probabilities=(0.5, 0.25, 0.75, 1.0),
        haplotypes=(1, 1, 1),
        capture_strands=(0, 0, 0),
        mate_indices=(0, 0, 0),
        mate_reverse_complements=(0, 0, 0),
        site_contexts=(
            MethylationContext.CG_G,
            MethylationContext.CG_C,
            MethylationContext.CG_G,
            MethylationContext.CG_C,
        ),
        site_sources=(MethylationSource.BETA,) * 4,
        site_alleles=(
            MethylationAllele.ALTERNATE_HAPLOTYPE,
            MethylationAllele.SHARED,
            MethylationAllele.ALTERNATE_HAPLOTYPE,
            MethylationAllele.SHARED,
        ),
        template_bases=(0, 3, 2, 3, 0, 1, 2, 3, 0, 4, 1, 3),
        truth=truth,
    )


def make_unequal_pe_batch() -> FragmentBatch:
    return FragmentBatch(
        first_fragment_ordinal=0,
        contig_indices=(0,),
        reference_begins=(40,),
        reference_ends=(45,),
        template_offsets=(0, 5),
        mate_offsets=(0, 2),
        site_offsets=(0, 0),
        mate_template_begins=(0, 2),
        mate_template_ends=(2, 5),
        site_template_offsets=(),
        site_probabilities=(),
        haplotypes=(0,),
        capture_strands=(0,),
        mate_indices=(0, 1),
        mate_reverse_complements=(0, 1),
        site_contexts=(),
        site_sources=(),
        site_alleles=(),
        template_bases=(0, 1, 2, 3, 0),
        truth=None,
    )


def frame_bounds(stream: bytes):
    bounds = []
    cursor = PREAMBLE_SIZE
    while cursor < len(stream):
        if len(stream) - cursor < FRAME_ENVELOPE.size:
            break
        payload_length = struct.unpack_from("<I", stream, cursor)[0]
        payload_start = cursor + FRAME_ENVELOPE.size
        payload_end = payload_start + payload_length
        end = payload_end + CRC.size
        if end > len(stream):
            break
        bounds.append((cursor, payload_start, payload_end, end))
        cursor = end
    return bounds


def rewrite_frame(stream: bytes, frame_index: int, mutate) -> bytes:
    changed = bytearray(stream)
    start, payload_start, payload_end, _ = frame_bounds(stream)[frame_index]
    mutate(changed, start, payload_start, payload_end)
    CRC.pack_into(changed, payload_end, crc32c(bytes(changed[start:payload_end])))
    return bytes(changed)


class ProtocolRoundTripTests(unittest.TestCase):
    def test_crc_primitive_remains_castagnoli(self) -> None:
        self.assertEqual(crc32c(b"123456789"), 0xE3069283)

    def test_no_truth_batch_round_trips_as_immutable_views(self) -> None:
        header = make_header(truth=False)
        encoded = encode_stream(
            header, (make_no_truth_batch(),), skipped_fragment_count=7
        )
        decoded = read_stream(
            encoded,
            expected_header=header,
            expected_skipped_fragment_count=7,
        )

        self.assertEqual(decoded.header, header)
        self.assertEqual(len(decoded.batches), 1)
        batch = decoded.batches[0]
        self.assertEqual(tuple(batch.template_offsets), (0, 4, 8))
        self.assertEqual(tuple(batch.template_bases), (0, 1, 2, 3, 3, 2, 1, 0))
        self.assertIsNone(batch.truth)
        self.assertTrue(batch.raw_payload.readonly)
        self.assertEqual(decoded.trailer.fragment_count, 2)
        self.assertEqual(decoded.trailer.fragment_batch_count, 1)
        self.assertEqual(decoded.trailer.per_contig_fragment_counts, (2,))

    def test_full_truth_batch_covers_snv_insertion_deletion_and_n(self) -> None:
        header = make_header(truth=True)
        decoded = read_stream(encode_stream(header, (make_full_truth_batch(),)))
        batch = decoded.batches[0]
        self.assertIsNotNone(batch.truth)
        truth = batch.truth
        self.assertEqual(tuple(truth.event_kinds), (1, 2, 3))
        self.assertEqual(
            tuple(truth.site_reference_positions),
            (12, 21, NO_REFERENCE_POSITION, 31),
        )
        self.assertEqual(tuple(truth.original_n_template_offsets), (0,))
        self.assertEqual(tuple(batch.site_probabilities), (0.5, 0.25, 0.75, 1.0))
        self.assertTrue(truth.event_ids.raw.readonly)

    def test_header_only_stream_has_zero_batch_trailer(self) -> None:
        decoded = read_stream(encode_stream(make_header(truth=False), ()))
        self.assertEqual(decoded.batches, ())
        self.assertEqual(decoded.trailer.fragment_count, 0)
        self.assertEqual(decoded.trailer.fragment_batch_count, 0)

    def test_stream_digest_excludes_trailer(self) -> None:
        encoded = encode_stream(make_header(truth=False), (make_no_truth_batch(),))
        trailer_start = frame_bounds(encoded)[-1][0]
        decoded = read_stream(encoded)
        self.assertEqual(
            decoded.trailer.stream_sha256,
            hashlib.sha256(encoded[:trailer_start]).digest(),
        )

    def test_header_carries_distinct_paired_end_read_lengths(self) -> None:
        header = replace(
            make_header(truth=False),
            mates_per_fragment=2,
            read_length_r1=2,
            read_length_r2=3,
        )
        decoded = read_stream(encode_stream(header, (make_unequal_pe_batch(),)))
        self.assertEqual(tuple(decoded.batches[0].mate_indices), (0, 1))
        self.assertEqual(decoded.trailer.mate_count, 2)

    def test_error_frame_is_terminal_and_rejects_the_run(self) -> None:
        output = io.BytesIO()
        writer = ProtocolWriter(output)
        writer.write_header(make_header(truth=False))
        writer.write_error(ErrorFrame(1204, "batch exceeds limit"))
        with self.assertRaisesRegex(CoreReportedError, "1204") as raised:
            read_stream(output.getvalue())
        self.assertEqual(raised.exception.error_code, 1204)

    def test_empty_error_diagnostic_is_valid_but_still_terminal(self) -> None:
        output = io.BytesIO()
        writer = ProtocolWriter(output)
        writer.write_header(make_header(truth=False))
        writer.write_error(ErrorFrame(1, ""))
        with self.assertRaises(CoreReportedError) as raised:
            read_stream(output.getvalue())
        self.assertEqual(raised.exception.message, "")


class ProtocolGoldenFixtureTests(unittest.TestCase):
    def test_frozen_binary_lengths_and_sha256_are_manifested(self) -> None:
        expected = {
            "header-none": (
                244,
                "94dfd95954618c8b57b4b68e8ad63734af610514e9c169f16466293b8fbf28ea",
            ),
            "batch-none": (
                156,
                "57d85ec26bf7e7f92394a1b1847a6df99d127015c405c4f025d76bc9e6665ba4",
            ),
            "batch-full": (
                468,
                "2618de26cdf85f69528e60ba100201989c963253bb7a301d57463b7d8691121c",
            ),
            "trailer-none": (
                112,
                "2a94417307f8acc1c4546e82cdba61afeb9ac5db552d68b89f01b8a2c910e32c",
            ),
            "error": (
                48,
                "be10795217958261c7c4fd5c80a3f062fb210265c6bf3df9b05ed8962e550f5c",
            ),
        }
        for name, (length, digest) in expected.items():
            with self.subTest(fixture=name):
                fixture = frozen_fixture(name)
                self.assertEqual(len(fixture), length)
                self.assertEqual(hashlib.sha256(fixture).hexdigest(), digest)

    def test_reference_encoder_matches_every_frozen_frame(self) -> None:
        no_truth_stream = encode_stream(
            make_header(truth=False),
            (make_no_truth_batch(),),
            skipped_fragment_count=7,
        )
        no_truth_bounds = frame_bounds(no_truth_stream)
        self.assertEqual(
            no_truth_stream[: no_truth_bounds[0][3]],
            frozen_fixture("header-none"),
        )
        self.assertEqual(
            no_truth_stream[no_truth_bounds[1][0] : no_truth_bounds[1][3]],
            frozen_fixture("batch-none"),
        )
        self.assertEqual(
            no_truth_stream[no_truth_bounds[2][0] : no_truth_bounds[2][3]],
            frozen_fixture("trailer-none"),
        )

        full_stream = encode_stream(
            make_header(truth=True), (make_full_truth_batch(),)
        )
        full_bounds = frame_bounds(full_stream)
        self.assertEqual(
            full_stream[full_bounds[1][0] : full_bounds[1][3]],
            frozen_fixture("batch-full"),
        )

        error_output = io.BytesIO()
        error_writer = ProtocolWriter(error_output)
        error_writer.write_header(make_header(truth=False))
        error_writer.write_error(ErrorFrame(1204, "batch exceeds limit"))
        error_stream = error_output.getvalue()
        error_bounds = frame_bounds(error_stream)
        self.assertEqual(
            error_stream[error_bounds[1][0] : error_bounds[1][3]],
            frozen_fixture("error"),
        )

    def test_frozen_envelopes_have_checked_lengths_types_sequences_and_crc(self) -> None:
        expected = {
            "batch-none": (136, 2, 0, 1, 156),
            "batch-full": (448, 2, 1, 1, 468),
            "trailer-none": (92, 3, 0, 2, 112),
            "error": (28, 255, 0, 1, 48),
        }
        for name, fields in expected.items():
            with self.subTest(fixture=name):
                fixture = frozen_fixture(name)
                payload_length, frame_type, flags, reserved, sequence = (
                    FRAME_ENVELOPE.unpack_from(fixture)
                )
                self.assertEqual(
                    (payload_length, frame_type, flags, sequence, len(fixture)),
                    fields,
                )
                self.assertEqual(reserved, 0)
                payload_end = FRAME_ENVELOPE.size + payload_length
                self.assertEqual(
                    CRC.unpack_from(fixture, payload_end)[0],
                    crc32c(fixture[:payload_end]),
                )

    def test_frozen_no_truth_stream_decodes_and_authenticates(self) -> None:
        encoded = b"".join(
            frozen_fixture(name)
            for name in ("header-none", "batch-none", "trailer-none")
        )
        decoded = read_stream(
            encoded,
            expected_header=make_header(truth=False),
            expected_skipped_fragment_count=7,
        )
        self.assertEqual(decoded.trailer.fragment_count, 2)
        self.assertEqual(tuple(decoded.batches[0].template_offsets), (0, 4, 8))


class ProtocolSemanticTests(unittest.TestCase):
    def assert_batch_rejected(
        self, batch: FragmentBatch, pattern: str, *, truth: bool
    ) -> None:
        with self.assertRaisesRegex(ProtocolError, pattern):
            encode_stream(make_header(truth=truth), (batch,))

    def test_common_prefix_and_context_are_strict(self) -> None:
        batch = make_no_truth_batch()
        cases = (
            (replace(batch, template_offsets=(1, 4, 8)), "begin at zero"),
            (replace(batch, mate_offsets=(0, 2, 2)), "mate count"),
            (replace(batch, site_template_offsets=(1, 9)), "site offset"),
            (
                replace(
                    batch,
                    site_contexts=(MethylationContext.CG_G, MethylationContext.CG_G),
                ),
                "incompatible",
            ),
            (replace(batch, site_probabilities=(float("nan"), 0.75)), "finite"),
        )
        for invalid, pattern in cases:
            with self.subTest(pattern=pattern):
                self.assert_batch_rejected(invalid, pattern, truth=False)

    def test_projection_runs_must_be_maximal(self) -> None:
        batch = make_full_truth_batch()
        truth = batch.truth
        split = replace(
            truth,
            projection_offsets=(0, 2, 4, 6),
            projection_template_begins=(0, 2, 0, 3, 0, 2),
            projection_template_ends=(2, 4, 2, 5, 2, 3),
            projection_reference_begins=(10, 12, 20, 22, 30, 33),
        )
        self.assert_batch_rejected(
            replace(batch, truth=split), "not maximal", truth=True
        )

    def test_sparse_truth_cover_and_event_shapes_fail_closed(self) -> None:
        batch = make_full_truth_batch()
        truth = batch.truth
        cases = (
            (
                replace(truth, projection_template_ends=(4, 3, 5, 2, 3)),
                "projection|insertion|cover",
            ),
            (
                replace(truth, event_alt_offsets=(0, 1, 1, 2)),
                "ALT bases disagree",
            ),
            (
                replace(
                    truth,
                    event_offsets=(0, 2, 2, 3),
                    event_ids=(2, 1, 3),
                ),
                "strictly increasing",
            ),
            (
                replace(truth, site_reference_positions=(12, 21, 22, 31)),
                "site reference position",
            ),
            (
                replace(
                    truth,
                    event_reference_begins=(11, 21, 32),
                    event_reference_ends=(12, 21, 33),
                ),
                "insertion anchor",
            ),
            (
                replace(
                    truth,
                    event_reference_begins=(11, 22, 31),
                    event_reference_ends=(12, 22, 32),
                ),
                "deletion boundary",
            ),
            (
                replace(truth, original_n_offsets=(0, 0, 0, 0)),
                "original_n_offsets|PRESERVE_N provenance",
            ),
        )
        for invalid_truth, pattern in cases:
            with self.subTest(pattern=pattern):
                self.assert_batch_rejected(
                    replace(batch, truth=invalid_truth), pattern, truth=True
                )

    def test_truth_policy_is_stream_wide(self) -> None:
        self.assert_batch_rejected(
            replace(make_no_truth_batch(), truth=make_full_truth_batch().truth),
            "forbids",
            truth=False,
        )
        self.assert_batch_rejected(
            replace(make_full_truth_batch(), truth=None), "requires", truth=True
        )

    def test_writer_is_poisoned_after_a_component_failure(self) -> None:
        output = io.BytesIO()
        writer = ProtocolWriter(output)
        writer.write_header(make_header(truth=False))
        invalid = replace(make_no_truth_batch(), first_fragment_ordinal=1)
        with self.assertRaisesRegex(ProtocolError, "ordinal"):
            writer.write_batch(invalid)
        with self.assertRaisesRegex(ProtocolError, "poisoned"):
            writer.finish()


class ProtocolCorruptionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.encoded = encode_stream(
            make_header(truth=False), (make_no_truth_batch(),)
        )

    def assert_rejected(self, stream: bytes, pattern: str) -> None:
        with self.assertRaisesRegex(ProtocolError, pattern):
            read_stream(stream)

    def test_crc_precedes_batch_view_exposure(self) -> None:
        corrupted = bytearray(self.encoded)
        _, payload_start, _, _ = frame_bounds(self.encoded)[1]
        corrupted[payload_start + 20] ^= 1
        self.assert_rejected(bytes(corrupted), "CRC32C")

    def test_unknown_batch_flag_is_rejected_with_valid_crc(self) -> None:
        corrupted = rewrite_frame(
            self.encoded,
            1,
            lambda data, start, _payload_start, _payload_end: data.__setitem__(
                start + 5, 0x80
            ),
        )
        self.assert_rejected(corrupted, "unknown flags|truth flag")

    def test_valid_crc_cannot_hide_an_invalid_prefix(self) -> None:
        def mutate(data, _start, payload_start, _payload_end):
            # Five u32 counts + three F-sized u32 columns precede template_offsets.
            struct.pack_into("<I", data, payload_start + 44, 1)

        self.assert_rejected(
            rewrite_frame(self.encoded, 1, mutate), "begin at zero"
        )

    def test_nonzero_padding_is_rejected_with_valid_crc(self) -> None:
        def mutate(data, _start, _payload_start, payload_end):
            data[payload_end - 1] = 1

        self.assert_rejected(
            rewrite_frame(self.encoded, 1, mutate), "padding must be zero"
        )

    def test_missing_padding_is_rejected_even_when_length_and_crc_agree(self) -> None:
        start, _payload_start, _payload_end, end = frame_bounds(self.encoded)[1]
        frame = bytearray(self.encoded[start:end])
        payload_length = struct.unpack_from("<I", frame)[0]
        shortened_payload = bytes(
            frame[FRAME_ENVELOPE.size : FRAME_ENVELOPE.size + payload_length - 1]
        )
        struct.pack_into("<I", frame, 0, payload_length - 1)
        shortened_envelope = bytes(frame[:FRAME_ENVELOPE.size])
        shortened_frame = (
            shortened_envelope
            + shortened_payload
            + CRC.pack(crc32c(shortened_envelope + shortened_payload))
        )
        changed = self.encoded[:start] + shortened_frame + self.encoded[end:]
        self.assert_rejected(changed, "non-canonical padding length")

    def test_trailer_count_and_digest_are_independently_checked(self) -> None:
        def count(data, _start, payload_start, _payload_end):
            struct.pack_into("<Q", data, payload_start, 3)

        self.assert_rejected(
            rewrite_frame(self.encoded, 2, count), "aggregate counts|sum"
        )

        def digest(data, _start, _payload_start, payload_end):
            data[payload_end - 1] ^= 0x80

        self.assert_rejected(
            rewrite_frame(self.encoded, 2, digest), "SHA-256"
        )

    def test_truncation_and_trailing_bytes_never_succeed(self) -> None:
        for removed in (1, 4, 20, len(self.encoded) // 2):
            with self.subTest(removed=removed):
                self.assert_rejected(self.encoded[:-removed], "truncated")
        self.assert_rejected(self.encoded + b"x", "follow the terminal trailer")

    def test_nonzero_core_status_rejects_complete_stream(self) -> None:
        self.assert_rejected_core_status(9)

    def assert_rejected_core_status(self, status: int) -> None:
        with self.assertRaisesRegex(ProtocolError, "non-zero status"):
            read_stream(self.encoded, core_exit_status=status)

    def test_header_projection_mismatch_is_fatal(self) -> None:
        expected = replace(make_header(truth=False), master_seed=7)
        self.assert_rejected_expected_header(expected)

    def assert_rejected_expected_header(self, header: Header) -> None:
        with self.assertRaisesRegex(ProtocolError, "Python projection"):
            read_stream(self.encoded, expected_header=header)


if __name__ == "__main__":
    unittest.main()
