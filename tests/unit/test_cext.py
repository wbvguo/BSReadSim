"""Equivalence checks for the private C extension."""

from __future__ import annotations

import math
import struct
import unittest


from bsreadsim import _cext
from bsreadsim.rng import _bernoulli_unchecked
from bsreadsim.rng import _philox4x32_10_unchecked
from bsreadsim.rng import (
    RNGStage,
    _u64_unchecked,
    derive_key,
)


def _crc32c_reference(value: bytes) -> int:
    """Return an independent bitwise CRC32C oracle for C-extension tests."""

    crc = 0xFFFFFFFF
    for byte in value:
        crc ^= byte
        for _ in range(8):
            crc = (crc >> 1) ^ (0x82F63B78 if crc & 1 else 0)
    return (~crc) & 0xFFFFFFFF


class CExtensionEquivalenceTests(unittest.TestCase):
    def test_batched_site_bernoulli_matches_reference(self) -> None:
        from tests.unit.test_process_stages import make_fragment  # pylint: disable=import-outside-toplevel

        fragment = make_fragment(paired_end=True, ordinal=73)
        key = derive_key(991, RNGStage.SITE_STATE, fragment.contig_index)
        expected = tuple(
            _bernoulli_unchecked(
                key,
                fragment.fragment_ordinal,
                site.site_index,
                site.methylation_probability,
            )
            for site in fragment.methylation_sites
        )
        self.assertEqual(
            _cext.sample_bernoulli_sites(
                fragment.methylation_sites,
                key,
                fragment.fragment_ordinal,
            ),
            expected,
        )

    def test_uniform_error_mate_matches_reference(self) -> None:
        alternatives = (
            (1, 2, 3),
            (0, 2, 3),
            (0, 1, 3),
            (0, 1, 2),
        )
        bases = bytes((0, 1, 2, 3, 4, 3, 2, 1, 0)) * 19
        key = derive_key(123456, RNGStage.SEQUENCING_ERROR, 6)
        entity = 0xFEDCBA9876543210
        for mate_index in (0, 1):
            for probability in (0.0, 0.005, 0.5, 1.0):
                expected_bases = bytearray(bases)
                expected_flags = bytearray(len(bases))
                for offset, base in enumerate(bases):
                    if base == 4:
                        continue
                    local_index = (mate_index << 32) | offset
                    if _bernoulli_unchecked(
                        key,
                        entity,
                        local_index,
                        probability,
                    ):
                        choice = (
                            _u64_unchecked(key, entity, local_index, 1) * 3
                        ) >> 64
                        expected_bases[offset] = alternatives[base][choice]
                        expected_flags[offset] = 1
                self.assertEqual(
                    _cext.apply_uniform_errors(
                        bases,
                        key,
                        entity,
                        mate_index,
                        probability,
                    ),
                    (bytes(expected_bases), bytes(expected_flags)),
                )

    def test_uniform_error_cext_boundary_rejects_invalid_values(self) -> None:
        with self.assertRaisesRegex(ValueError, "base code"):
            _cext.apply_uniform_errors(b"\x05", 0, 0, 0, 0.1)
        with self.assertRaisesRegex(ValueError, "probability"):
            _cext.apply_uniform_errors(b"\x00", 0, 0, 0, math.nan)
        with self.assertRaisesRegex(OverflowError, "mate_index"):
            _cext.apply_uniform_errors(b"\x00", 0, 0, 1 << 32, 0.1)

    def test_philox_u64_and_bernoulli_match_reference_words(self) -> None:
        cases = (
            (0, 0, 0),
            (1, 2, 3),
            ((1 << 64) - 1, (1 << 64) - 1, (1 << 64) - 1),
            (0x0123456789ABCDEF, 17, (1 << 32) | 149),
        )
        for key, entity, local in cases:
            block = _philox4x32_10_unchecked(key, entity, local)
            for pair in (0, 1):
                with self.subTest(key=key, entity=entity, local=local, pair=pair):
                    expected = block[pair * 2] | (block[pair * 2 + 1] << 32)
                    self.assertEqual(_cext.u64(key, entity, local, pair), expected)
                    uniform = math.ldexp(float(expected >> 11), -53)
                    self.assertEqual(
                        _cext.bernoulli(key, entity, local, 0.375, pair),
                        uniform < 0.375,
                    )

        with self.assertRaisesRegex(ValueError, "pair"):
            _cext.u64(0, 0, 0, 2)

    def test_packed_philox_pairs_match_scalar_reference(self) -> None:
        key = 0x0123456789ABCDEF  # deterministic RNG key
        entities = (0, 1, 17, 1 << 32, (1 << 64) - 1)
        local_indices = (0, 3, (1 << 32) | 149, 7, (1 << 64) - 1)
        packed_entities = struct.pack("={}Q".format(len(entities)), *entities)
        packed_locals = struct.pack(
            "={}Q".format(len(local_indices)), *local_indices
        )
        observed = _cext.philox_pairs(
            key,
            packed_entities,
            packed_locals,
        )
        expected = tuple(
            struct.pack(
                "={}Q".format(len(entities)),
                *(
                    _cext.u64(key, entity, local, pair)
                    for entity, local in zip(
                        entities, local_indices, strict=True
                    )
                ),
            )
            for pair in (0, 1)
        )
        self.assertEqual(observed, expected)
        self.assertEqual(_cext.philox_pairs(key, b"", b""), (b"", b""))
        with self.assertRaisesRegex(ValueError, "same length"):
            _cext.philox_pairs(key, packed_entities, packed_locals[:-8])
        with self.assertRaisesRegex(ValueError, "packed uint64"):
            _cext.philox_pairs(key, b"x", b"y")

    def test_fastq_batch_formatter_matches_reference_bytes(self) -> None:
        def record(contig, start, end, ordinal, mate_number, sequence, quality):
            left = start + 1
            right = end if end > start else left
            return "@{}:{}-{}:{:x}/{}\n{}\n+\n{}\n".format(
                contig,
                left,
                right,
                ordinal,
                mate_number,
                sequence.decode("ascii"),
                chr(quality) * len(sequence),
            ).encode("utf-8")

        contigs = ("chr1", "chrTwo", "chr雪")
        starts = (0, 10, (1 << 64) - 2)
        ends = (5, 10, (1 << 64) - 1)
        ordinals = (0, 9, (1 << 64) - 1)
        read_length = 5
        for paired_end in (False, True):
            mates_per_fragment = 2 if paired_end else 1
            sequences = tuple(
                "".join(
                    "ACGTN"[(row + cycle) % 5]
                    for cycle in range(read_length)
                ).encode("ascii")
                for row in range(len(ordinals) * mates_per_fragment)
            )
            mate_starts = tuple(
                start
                for start in starts
                for _ in range(mates_per_fragment)
            )
            mate_ends = tuple(
                end
                for end in ends
                for _ in range(mates_per_fragment)
            )
            observed = _cext.format_fastq_batch(
                contigs,
                struct.pack("<{}Q".format(len(ordinals)), *ordinals),
                struct.pack(
                    "<{}Q".format(len(ordinals) + 1),
                    *(index * mates_per_fragment
                      for index in range(len(ordinals) + 1)),
                ),
                bytes(
                    index % mates_per_fragment
                    for index in range(len(sequences))
                ),
                struct.pack("<{}Q".format(len(mate_starts)), *mate_starts),
                struct.pack("<{}Q".format(len(mate_ends)), *mate_ends),
                b"".join(sequences),
                read_length,
                73,
                int(paired_end),
            )
            expected_read1 = b"".join(
                record(
                    contigs[index],
                    starts[index],
                    ends[index],
                    ordinal,
                    1,
                    sequences[index * mates_per_fragment],
                    73,
                )
                for index, ordinal in enumerate(ordinals)
            )
            expected_read2 = (
                b"".join(
                    record(
                        contigs[index],
                        starts[index],
                        ends[index],
                        ordinal,
                        2,
                        sequences[index * mates_per_fragment + 1],
                        73,
                    )
                    for index, ordinal in enumerate(ordinals)
                )
                if paired_end
                else None
            )
            expected_lengths = tuple(
                (
                    len(record(
                        contigs[index],
                        starts[index],
                        ends[index],
                        ordinal,
                        1,
                        b"A" * read_length,
                        73,
                    )),
                    (
                        len(record(
                            contigs[index],
                            starts[index],
                            ends[index],
                            ordinal,
                            2,
                            b"A" * read_length,
                            73,
                        ))
                        if paired_end
                        else 0
                    ),
                )
                for index, ordinal in enumerate(ordinals)
            )
            with self.subTest(paired_end=paired_end):
                self.assertEqual(
                    observed,
                    (expected_read1, expected_read2, expected_lengths),
                )

    def test_fastq_batch_formatter_rejects_malformed_columns(self) -> None:
        good = (
            ("chr1",),
            struct.pack("<Q", 0),
            struct.pack("<2Q", 0, 1),
            b"\x00",
            struct.pack("<Q", 10),
            struct.pack("<Q", 15),
            b"ACGTN",
            5,
            73,
            0,
        )
        mutations = (
            (("bad names",) + good[1:], "contig"),
            (("chr\u0080",) + good[1:], "contig"),
            (good[:1] + (b"x",) + good[2:], "ordinals"),
            (good[:2] + (struct.pack("<2Q", 0, 2),) + good[3:], "offsets"),
            (good[:3] + (b"\x01",) + good[4:], "mate rows"),
            (
                good[:4]
                + (struct.pack("<Q", 16), struct.pack("<Q", 15))
                + good[6:],
                "envelope",
            ),
            (good[:6] + (b"ACGTX",) + good[7:], "non-ACGTN"),
            (good[:7] + (0,) + good[8:], "read length"),
            (good[:8] + (32,) + good[9:], "quality"),
            (good[:9] + (2,), "paired-end"),
        )
        for arguments, message in mutations:
            with self.subTest(message=message):
                with self.assertRaisesRegex((ValueError, TypeError), message):
                    _cext.format_fastq_batch(*arguments)

    def test_crc32c_matches_python_reference_across_release_threshold(self) -> None:
        for length in (0, 1, 7, 8, 31, 4095, 4096, 4097, 65536):
            value = bytes((index * 37 + 11) & 0xFF for index in range(length))
            with self.subTest(length=length):
                expected = _crc32c_reference(value)
                self.assertEqual(_cext.crc32c(value), expected)
                self.assertEqual(_cext.crc32c(bytearray(value)), expected)
                framed = b"x" + value + b"y"
                self.assertEqual(_cext.crc32c(memoryview(framed)[1:-1]), expected)
