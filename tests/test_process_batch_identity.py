"""Contract tests for FASTQ identifiers and BAM query names."""

from pathlib import Path
import sys
import unittest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim.process.batch import (  # noqa: E402
    READ_NAME_CONTRACT,
    ReadNameError,
    format_fragment_identifier,
    fragment_identifier_coordinates,
)


class ReadNameTests(unittest.TestCase):
    def test_contract_uses_variable_width_lowercase_hexadecimal(self) -> None:
        self.assertEqual(READ_NAME_CONTRACT, "bsreadsim-read-name-v2")
        self.assertEqual(
            format_fragment_identifier("chr1", 0, 1, 0x2A),
            "chr1:1-1:2a",
        )

    def test_regular_fragment_uses_one_based_inclusive_coordinates(self) -> None:
        self.assertEqual(
            format_fragment_identifier("chr1", 100, 108, 7),
            "chr1:101-108:7",
        )
        self.assertEqual(
            format_fragment_identifier(
                "chr1", 100, 108, 7, pair_number=2
            ),
            "chr1:101-108:7/2",
        )

    def test_zero_width_insertion_uses_its_one_based_anchor(self) -> None:
        self.assertEqual(fragment_identifier_coordinates(20, 20), (21, 21))
        self.assertEqual(
            format_fragment_identifier("chrI", 20, 20, 0, pair_number=1),
            "chrI:21-21:0/1",
        )

    def test_invalid_fields_fail_closed(self) -> None:
        failures = (
            (("", 0, 1, 0), "contig"),
            (("chr 1", 0, 1, 0), "contig"),
            (("chr\u0080", 0, 1, 0), "contig"),
            (("chr1", 2, 1, 0), "reversed"),
            (("chr1", (1 << 64) - 1, (1 << 64) - 1, 0), "converted"),
            (("chr1", 0, 1, -1), "ordinal"),
        )
        for arguments, message in failures:
            with self.subTest(arguments=arguments):
                with self.assertRaisesRegex(ReadNameError, message):
                    format_fragment_identifier(*arguments)

        with self.assertRaisesRegex(ReadNameError, "pair"):
            format_fragment_identifier("chr1", 0, 1, 0, pair_number=3)
        with self.assertRaisesRegex(ReadNameError, "pair"):
            format_fragment_identifier("chr1", 0, 1, 0, pair_number=True)


if __name__ == "__main__":
    unittest.main()
