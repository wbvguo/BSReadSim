"""Contract tests for SAM formatting behind annotated BAM output."""

from dataclasses import replace
import unittest


from tests.helpers.process_support import UniformProcessConfig, process_fragment

from bsreadsim.output.bam import (
    BamError,
    build_sam_header,
    format_sam_fragment,
)
from bsreadsim.process import (
    BaseState,
    CaptureStrand,
    ConversionMode,
    ProcessedFragment,
    ProcessedMate,
)
from bsreadsim.htsim.protocol import Contig
from tests.unit.test_process_stages import make_fragment
from tests.unit.test_protocol import make_header


RUN_ID = "00000000-0000-4000-8000-000000000002"


def annotation(offset: int, position: int) -> BaseState:
    return BaseState(
        read_offset=offset,
        reference_pos=position,
        variant_index=0xFFFFFFFF,
        site_index=None,
        methylated=None,
        oriented_base=0,
        post_conversion_base=0,
        final_base=0,
        conversion_attempted=False,
        conversion_succeeded=False,
        sequencing_error=False,
        quality_phred=30,
    )


def fragment_for(
    positions,
    *,
    sequence="ACGTN",
    quality="ABCDE",
    reverse=False,
    reference_start=10,
    reference_end=16,
) -> ProcessedFragment:
    mate = ProcessedMate(
        mate_index=0,
        reverse_complement=reverse,
        conversion_mode=ConversionMode.C2T,
        reference_start=reference_start,
        reference_end=reference_end,
        sequence=sequence,
        quality=quality,
        base_states=tuple(
            annotation(offset, position)
            for offset, position in enumerate(positions)
        ),
    )
    return ProcessedFragment(
        fragment_ordinal=7,
        contig_name="chrMini",
        reference_start=reference_start,
        reference_end=reference_end,
        haplotype=1,
        fragment_conversion_mode=ConversionMode.C2T,
        variants=(),
        site_states=(),
        mates=(mate,),
    )


def fields(record: bytes):
    return record.decode("ascii").rstrip("\n").split("\t")


def tag_fields(record_fields):
    return {field[:2]: field for field in record_fields[11:]}


class BamFormattingTests(unittest.TestCase):
    def test_header_declares_unsorted_sam_16_dictionary_and_provenance(self) -> None:
        value = build_sam_header(
            make_header(details=True),
            sample_name="sample",
            program_version="1.2.3",
        ).decode("ascii")

        self.assertTrue(value.startswith("@HD\tVN:1.6\tSO:unsorted\n"))
        self.assertIn("@SQ\tSN:chrMini\tLN:100\n", value)
        self.assertIn("@RG\tID:{}\tSM:sample\n".format(RUN_ID), value)
        self.assertIn("@PG\tID:bsreadsim\tPN:bsreadsim\tVN:1.2.3\n", value)
        self.assertIn("MAPQ 60 denotes simulated origin", value)
        self.assertIn("BSREADSIM_ZT=state64", value)
        self.assertIn(
            "BSREADSIM_ZS=informative-strand-conversion-v1;REQUIRED=1",
            value,
        )
        self.assertIn("BSREADSIM_ZR=u16x12;REQUIRED=1", value)
        self.assertIn("BSREADSIM_ZF=u16x12;ENABLED=0", value)

    def test_indel_projection_forms_query_complete_cigar(self) -> None:
        record = format_sam_fragment(
            fragment_for((10, 11, -1, 12, 15)),
            paired_end=False,
            read_group_id=RUN_ID,
            contig_length=100,
        )[0]
        value = fields(record)

        self.assertEqual(value[:11], [
            "chrMini:11-16:7", "0", "chrMini", "11", "60",
            "2M1I1M2D1M", "*", "0", "0", "ACGTN", "ABCDE",
        ])
        tags = tag_fields(value)
        self.assertEqual(tags["RG"], "RG:Z:" + RUN_ID)
        self.assertEqual(tags["AS"], "AS:i:5")
        self.assertEqual(tags["zs"], "zs:Z:W_C2T")
        self.assertTrue(tags["zt"].startswith("zt:Z:"))
        self.assertTrue(tags["zr"].startswith("zr:B:S,"))
        self.assertNotIn("PG", tags)

    def test_reverse_record_uses_reference_forward_sequence_and_quality(self) -> None:
        record = format_sam_fragment(
            fragment_for(
                (15, 12, -1, 11, 10),
                reverse=True,
            ),
            paired_end=False,
            read_group_id=RUN_ID,
            contig_length=100,
        )[0]
        value = fields(record)

        self.assertEqual(value[1], "16")
        self.assertEqual(value[3], "11")
        self.assertEqual(value[5], "2M1I1M2D1M")
        self.assertEqual(value[9], "NACGT")
        self.assertEqual(value[10], "EDCBA")

    def test_pure_insertion_is_mapped_at_its_reference_anchor(self) -> None:
        record = format_sam_fragment(
            fragment_for(
                (-1, -1, -1),
                sequence="ACG",
                quality="ABC",
                reference_start=20,
                reference_end=20,
            ),
            paired_end=False,
            read_group_id=RUN_ID,
            contig_length=100,
        )[0]
        value = fields(record)

        self.assertEqual(value[3], "21")
        self.assertEqual(value[5], "3I")

    def test_paired_records_have_standard_flags_mate_fields_and_tlen(self) -> None:
        processed = process_fragment(
            make_fragment(paired_end=True),
            "chrMini",
            UniformProcessConfig(
                master_seed=7,
                conversion_rate=1,
                error_rate=0,
                quality_phred=30,
            ),
        )
        first, second = tuple(
            fields(record)
            for record in format_sam_fragment(
                processed,
                paired_end=True,
                read_group_id="run",
                contig_length=1000,
            )
        )

        self.assertEqual(first[:9], [
            "chrMini:101-108:0", "99", "chrMini", "101", "60", "5M", "=", "103", "7",
        ])
        self.assertEqual(second[:9], [
            "chrMini:101-108:0", "147", "chrMini", "103", "60", "5M", "=", "101", "-7",
        ])
        self.assertEqual(tag_fields(first)["MC"], "MC:Z:5M")
        self.assertEqual(tag_fields(second)["MC"], "MC:Z:5M")
        self.assertEqual(tag_fields(first)["zs"], "zs:Z:W_C2T")
        self.assertEqual(tag_fields(second)["zs"], "zs:Z:W_G2A")
        self.assertEqual(second[9], "TGTCG")

    def test_crick_fragment_has_readable_strand_conversion_tags(self) -> None:
        source = make_fragment(
            paired_end=True,
            capture_strand=CaptureStrand.REVERSE,
        )
        left, right = source.mates
        original_bottom = replace(
            source,
            mates=(
                replace(right, mate_index=0),
                replace(left, mate_index=1),
            ),
        )
        processed = process_fragment(
            original_bottom,
            "chrMini",
            UniformProcessConfig(
                master_seed=7,
                conversion_rate=1,
                error_rate=0,
                quality_phred=30,
            ),
        )
        first, second = tuple(
            tag_fields(fields(record))
            for record in format_sam_fragment(
                processed,
                paired_end=True,
                read_group_id="run",
                contig_length=1000,
            )
        )

        self.assertEqual(first["zs"], "zs:Z:C_C2T")
        self.assertEqual(second["zs"], "zs:Z:C_G2A")

    def test_rich_tags_are_required_and_fragment_summary_is_optional(self) -> None:
        fragment = fragment_for((10, 11, -1, 12, 15))
        regular = tag_fields(fields(format_sam_fragment(
            fragment,
            paired_end=False,
            read_group_id=RUN_ID,
            contig_length=100,
        )[0]))
        summarized = tag_fields(fields(format_sam_fragment(
            fragment,
            paired_end=False,
            read_group_id=RUN_ID,
            contig_length=100,
            include_fragment_summary=True,
        )[0]))

        self.assertEqual(len(regular["zt"].removeprefix("zt:Z:")), 5)
        self.assertEqual(len(regular["zr"].split(",")), 13)
        self.assertEqual(regular["zs"], "zs:Z:W_C2T")
        self.assertNotIn("zf", regular)
        self.assertEqual(len(summarized["zf"].split(",")), 13)

    def test_reference_dictionary_rejects_names_other_tools_cannot_parse(self) -> None:
        header = replace(
            make_header(details=True),
            contigs=(
                Contig(
                    name="bad,name",
                    length=100,
                    reference_sha256=b"0" * 32,
                ),
            ),
        )
        with self.assertRaisesRegex(BamError, "reference name"):
            build_sam_header(
                header,
                sample_name="sample",
                program_version="1.2.3",
            )

    def test_query_name_enforces_sam_character_and_length_limits(self) -> None:
        maximum = replace(fragment_for((10, 11, -1, 12, 15)), contig_name="a" * 246)
        record = format_sam_fragment(
            maximum,
            paired_end=False,
            read_group_id=RUN_ID,
            contig_length=100,
        )[0]
        self.assertEqual(len(fields(record)[0]), 254)

        invalid_names = ("a" * 247, "@chrMini")
        for contig_name in invalid_names:
            with self.subTest(contig_name=contig_name[:12]):
                fragment = replace(maximum, contig_name=contig_name)
                with self.assertRaisesRegex(BamError, "query name"):
                    format_sam_fragment(
                        fragment,
                        paired_end=False,
                        read_group_id=RUN_ID,
                        contig_length=100,
                    )


if __name__ == "__main__":
    unittest.main()
