from dataclasses import replace
import unittest

import bsreadsim.output.bam as bam
from tests.test_output_session import processed


class BamBatchTests(unittest.TestCase):
    def expected(self, fragments, *, paired_end):
        records = []
        for fragment in fragments:
            records.extend(
                bam.format_sam_fragment(
                    fragment,
                    paired_end=paired_end,
                    read_group_id="run",
                    contig_length=1000,
                )
            )
        return b"".join(records), tuple(len(record) for record in records)

    def cases(self, *, paired_end):
        regular = processed(paired_end=paired_end)
        compact = processed(
            paired_end=paired_end,
            compact_base_states=True,
        )
        return (
            regular,
            replace(compact, fragment_ordinal=compact.fragment_ordinal + 1),
        )


    def test_batch_matches_fragment_reference_byte_for_byte(self):
        for paired_end in (False, True):
            fragments = self.cases(paired_end=paired_end)
            with self.subTest(paired_end=paired_end):
                self.assertEqual(
                    bam.format_sam_batch(
                        fragments,
                        paired_end=paired_end,
                        read_group_id="run",
                        contig_lengths=(1000, 1000),
                    ),
                    self.expected(fragments, paired_end=paired_end),
                )


if __name__ == "__main__":
    unittest.main()
