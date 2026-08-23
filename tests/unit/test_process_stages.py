"""Contract tests for dependency-free fragment processing."""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
import json
import unittest
from unittest.mock import patch


from tests.helpers.process_support import UniformProcessConfig, process_fragment
import tests.helpers.process_support as process_support_module

from bsreadsim.process import (
    PROCESS_STAGE_ORDER,
    ConversionMode,
    ProcessConfig,
    ProcessError,
    UniformError,
    UniformQuality,
    process_fragment_batch,
)
import bsreadsim.process.bisulfite as bisulfite_module
from bsreadsim.process.batch import (
    CaptureStrand,
    Fragment,
    Mate,
    MethylationAllele,
    MethylationContext,
    MethylationSite,
    MethylationSource,
    NO_VARIANT_INDEX,
    SiteReference,
    Variant,
    VariantKind,
    VariantSource,
)
from bsreadsim.rng import RNGStage, derive_key
from bsreadsim.process.sequencing import (
    parse_quality_confusion,
    parse_quality_markov,
)


def advanced_quality_model():
    mate = {
        "initial_counts": [[1, 3]] * 5,
        "transition_counts": [[1, 3], [3, 1]],
    }
    return parse_quality_markov(
        json.dumps(
            {
                "schema": "quality-markov-v1",
                "quality_scores": [10, 30],
                "mates": [mate, mate],
            },
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    )


def advanced_error_model():
    identity = [
        [1, 0, 0, 0],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
    ]
    rotate = [
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1],
        [1, 0, 0, 0],
    ]
    mate = {"base_transition_counts": [identity, rotate]}
    return parse_quality_confusion(
        json.dumps(
            {
                "schema": "quality-confusion-v1",
                "quality_scores": [10, 30],
                "mates": [mate, mate],
            },
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    )


def make_fragment(
    *,
    paired_end: bool,
    ordinal: int = 0,
    capture_strand: CaptureStrand = CaptureStrand.FORWARD,
    probabilities=(1.0, 1.0),
) -> Fragment:
    sites = (
        MethylationSite(
            site_index=0,
            template_offset=1,
            reference_pos=101,
            context=MethylationContext.CG_C,
            methylation_source=MethylationSource.CGMAP,
            allele=MethylationAllele.SHARED,
            methylation_probability=probabilities[0],
        ),
        MethylationSite(
            site_index=1,
            template_offset=5,
            reference_pos=105,
            context=MethylationContext.CHG_C,
            methylation_source=MethylationSource.BETA,
            allele=MethylationAllele.SHARED,
            methylation_probability=probabilities[1],
        ),
    )
    mates = (
        Mate(
            mate_index=0,
            reverse_complement=False,
            template_start=0,
            template_end=5,
            reference_start=100,
            reference_end=105,
            site_refs=(SiteReference(read_offset=1, site_index=0),),
        ),
    )
    if paired_end:
        mates += (
            Mate(
                mate_index=1,
                reverse_complement=True,
                template_start=2,
                template_end=7,
                reference_start=102,
                reference_end=107,
                site_refs=(SiteReference(read_offset=1, site_index=1),),
            ),
        )
    return Fragment(
        fragment_ordinal=ordinal,
        contig_index=0,
        haplotype=1,
        capture_strand=capture_strand,
        reference_start=100,
        reference_end=108,
        template_bases=bytes((0, 1, 1, 2, 3, 1, 2, 0)),
        reference_positions=tuple(range(100, 108)),
        base_variant_indices=(NO_VARIANT_INDEX,) * 8,
        variants=(),
        methylation_sites=sites,
        mates=mates,
    )


def make_explicit_overlap_fragment(probability=0.0) -> Fragment:
    site = MethylationSite(
        site_index=0,
        template_offset=0,
        reference_pos=300,
        context=MethylationContext.CG_C,
        methylation_source=MethylationSource.BETA,
        allele=MethylationAllele.SHARED,
        methylation_probability=probability,
    )
    return Fragment(
        fragment_ordinal=0,
        contig_index=0,
        haplotype=0,
        capture_strand=CaptureStrand.FORWARD,
        reference_start=300,
        reference_end=301,
        template_bases=bytes((1,)),
        reference_positions=(300,),
        base_variant_indices=(NO_VARIANT_INDEX,),
        variants=(),
        methylation_sites=(site,),
        mates=(
            Mate(
                mate_index=0,
                reverse_complement=False,
                template_start=0,
                template_end=1,
                reference_start=300,
                reference_end=301,
                site_refs=(SiteReference(read_offset=0, site_index=0),),
            ),
            Mate(
                mate_index=1,
                reverse_complement=True,
                template_start=0,
                template_end=1,
                reference_start=300,
                reference_end=301,
                site_refs=(SiteReference(read_offset=0, site_index=0),),
            ),
        ),
    )


def make_site_free_fragment(
    *,
    template_bases=(0, 1, 2, 3, 4),
    capture_strand: CaptureStrand = CaptureStrand.FORWARD,
    reverse_complement: bool = False,
    paired_end: bool = False,
    ordinal: int = 0,
) -> Fragment:
    length = len(template_bases)
    mates = (
        Mate(
            mate_index=0,
            reverse_complement=reverse_complement,
            template_start=0,
            template_end=length,
            reference_start=200,
            reference_end=200 + length,
            site_refs=(),
        ),
    )
    if paired_end:
        mates += (
            Mate(
                mate_index=1,
                reverse_complement=True,
                template_start=0,
                template_end=length,
                reference_start=200,
                reference_end=200 + length,
                site_refs=(),
            ),
        )
    return Fragment(
        fragment_ordinal=ordinal,
        contig_index=0,
        haplotype=0,
        capture_strand=capture_strand,
        reference_start=200,
        reference_end=200 + length,
        template_bases=bytes(template_bases),
        reference_positions=tuple(range(200, 200 + length)),
        base_variant_indices=(NO_VARIANT_INDEX,) * length,
        variants=(),
        methylation_sites=(),
        mates=mates,
    )


class UniformProcessTests(unittest.TestCase):
    def test_stage_order_is_explicit(self) -> None:
        self.assertEqual(
            PROCESS_STAGE_ORDER,
            (
                "methylation",
                "fragment-orientation",
                "fragment-conversion",
                "mate-derivation",
                "quality",
                "sequencing-error",
                "format-ready",
            ),
        )

    def test_fragment_variants_cross_process_unchanged(self) -> None:
        events = (
            Variant(
                index=7,
                id="varsim_7",
                source=VariantSource.DE_NOVO,
                kind=VariantKind.DELETION,
                phased_haplotype=1,
                reference_start=103,
                reference_end=105,
                ref_bases=bytes((3, 1)),
                alt_bases=b"",
            ),
        )
        fragment = replace(
            make_fragment(paired_end=False), variants=events
        )
        result = process_fragment(
            fragment,
            "chr1",
            UniformProcessConfig(
                master_seed=7,
                conversion_rate=1.0,
                error_rate=0.0,
            ),
        )

        self.assertEqual(result.variants, events)

    def test_fastq_only_path_omits_base_states_without_changing_reads(self) -> None:
        fragment = make_fragment(paired_end=True)
        config = UniformProcessConfig(
            master_seed=7,
            conversion_rate=1.0,
            error_rate=0.0,
            quality_phred=30,
        )
        full = process_fragment(fragment, "chr1", config)
        fastq_only = process_fragment(
            fragment,
            "chr1",
            config,
            include_details=False,
        )

        self.assertEqual(
            tuple((mate.sequence, mate.quality) for mate in fastq_only.mates),
            tuple((mate.sequence, mate.quality) for mate in full.mates),
        )
        self.assertTrue(all(mate.base_states == () for mate in fastq_only.mates))
        self.assertEqual(fastq_only.site_states, ())

    def test_single_end_probability_boundaries(self) -> None:
        result = process_fragment(
            make_fragment(paired_end=False, probabilities=(0.0, 1.0)),
            "chr1",
            UniformProcessConfig(
                master_seed=7,
                conversion_rate=1.0,
                error_rate=0.0,
                quality_phred=30,
            ),
        )

        self.assertEqual(tuple(state.methylated for state in result.site_states), (False, True))
        self.assertEqual(
            (
                result.site_states[0].template_offset,
                result.site_states[0].reference_pos,
                result.site_states[0].context,
                result.site_states[0].methylation_source,
                result.site_states[0].allele,
            ),
            (
                1,
                101,
                MethylationContext.CG_C,
                MethylationSource.CGMAP,
                MethylationAllele.SHARED,
            ),
        )
        self.assertEqual(len(result.mates), 1)
        self.assertEqual(result.mates[0].sequence, "ATTGT")
        self.assertEqual(result.mates[0].quality, "?" * 5)
        implicit = result.mates[0].base_states[2]
        self.assertIsNone(implicit.site_index)
        self.assertIs(implicit.methylated, False)
        self.assertTrue(implicit.conversion_attempted)
        self.assertTrue(implicit.conversion_succeeded)

    def test_site_free_c_g_and_reverse_complement_targets_convert(self) -> None:
        config = UniformProcessConfig(
            master_seed=11,
            conversion_rate=1.0,
            error_rate=0.0,
        )
        cases = (
            (
                make_site_free_fragment(),
                "ATGTN",
                ConversionMode.C2T,
                (1,),
                (2,),
            ),
            (
                make_site_free_fragment(capture_strand=CaptureStrand.REVERSE),
                "ACATN",
                ConversionMode.G2A,
                (2,),
                (1,),
            ),
            (
                make_site_free_fragment(
                    template_bases=(0, 1, 1, 2, 3),
                    reverse_complement=True,
                ),
                "ACAAT",
                ConversionMode.G2A,
                (2, 3),
                (1,),
            ),
        )

        for fragment, sequence, mode, target_offsets, non_target_offsets in cases:
            with self.subTest(sequence=sequence, mode=mode):
                result = process_fragment(fragment, "chrImplicit", config)
                mate = result.mates[0]
                self.assertEqual(mate.sequence, sequence)
                self.assertEqual(mate.conversion_mode, mode)
                self.assertEqual(result.site_states, ())
                for offset in target_offsets:
                    annotation = mate.base_states[offset]
                    self.assertIsNone(annotation.site_index)
                    self.assertIs(annotation.methylated, False)
                    self.assertTrue(annotation.conversion_attempted)
                    self.assertTrue(annotation.conversion_succeeded)
                for offset in non_target_offsets:
                    annotation = mate.base_states[offset]
                    self.assertIsNone(annotation.site_index)
                    self.assertIsNone(annotation.methylated)
                    self.assertFalse(annotation.conversion_attempted)

    def test_explicit_methylated_site_is_protected_beside_implicit_target(self) -> None:
        result = process_fragment(
            make_fragment(paired_end=False, probabilities=(1.0, 1.0)),
            "chr1",
            UniformProcessConfig(
                master_seed=7,
                conversion_rate=1.0,
                error_rate=0.0,
            ),
        )

        self.assertEqual(result.mates[0].sequence, "ACTGT")
        explicit = result.mates[0].base_states[1]
        implicit = result.mates[0].base_states[2]
        self.assertEqual((implicit.site_index, implicit.methylated), (None, False))
        self.assertTrue(implicit.conversion_succeeded)
        self.assertEqual((explicit.site_index, explicit.methylated), (0, True))
        self.assertFalse(explicit.conversion_attempted)
        self.assertFalse(explicit.conversion_succeeded)

    def test_overlapping_mates_share_one_latent_site_state(self) -> None:
        fragment = make_explicit_overlap_fragment()
        config = UniformProcessConfig(
            master_seed=7,
            conversion_rate=1.0,
            error_rate=0.0,
        )
        calls = []
        real_bernoulli = process_support_module.bernoulli

        def traced_bernoulli(key, entity_ordinal, local_index, probability):
            calls.append((key, entity_ordinal, local_index, probability))
            return real_bernoulli(key, entity_ordinal, local_index, probability)

        # The test-only scalar adapter exposes the exact RNG addresses. C
        # extension equivalence has a separate differential test in test_cext.py.
        with patch.object(
            process_support_module,
            "bernoulli",
            side_effect=traced_bernoulli,
        ), patch.object(
            bisulfite_module,
            "bernoulli",
            side_effect=traced_bernoulli,
        ):
            result = process_fragment(fragment, "chr1", config)

        self.assertEqual(tuple(state.methylated for state in result.site_states), (False,))
        self.assertEqual(tuple(mate.sequence for mate in result.mates), ("T", "A"))
        mate_1_site = result.mates[0].base_states[0]
        mate_2_site = result.mates[1].base_states[0]
        self.assertEqual(mate_1_site.site_index, 0)
        self.assertEqual(mate_2_site.site_index, 0)
        self.assertIs(mate_1_site.methylated, mate_2_site.methylated)
        self.assertTrue(mate_1_site.conversion_succeeded)
        self.assertTrue(mate_2_site.conversion_succeeded)
        self.assertEqual(result.mates[1].conversion_mode, ConversionMode.G2A)
        site_key = derive_key(config.master_seed, RNGStage.SITE_STATE, 0)
        conversion_key = derive_key(config.master_seed, RNGStage.CONVERSION, 0)
        relevant_calls = [
            call for call in calls if call[0] in (site_key, conversion_key)
        ]
        self.assertEqual(
            relevant_calls,
            [
                (site_key, 0, 0, 0.0),
                (conversion_key, 0, 0, 1.0),
            ],
        )

    def test_site_free_overlap_converts_each_fragment_base_once(self) -> None:
        fragment = make_site_free_fragment(
            template_bases=(1, 2, 1, 2),
            paired_end=True,
            ordinal=17,
        )
        config = UniformProcessConfig(
            master_seed=19,
            conversion_rate=1.0,
            error_rate=0.0,
        )
        calls = []
        real_bernoulli = process_support_module.bernoulli

        def traced_bernoulli(key, entity_ordinal, local_index, probability):
            calls.append((key, entity_ordinal, local_index, probability))
            return real_bernoulli(key, entity_ordinal, local_index, probability)

        with patch.object(bisulfite_module, "bernoulli", side_effect=traced_bernoulli):
            result = process_fragment(fragment, "chrOverlap", config)

        self.assertEqual(tuple(mate.sequence for mate in result.mates), ("TGTG", "CACA"))
        for mate in result.mates:
            implicit_targets = [
                annotation
                for annotation in mate.base_states
                if annotation.conversion_attempted
            ]
            self.assertEqual(len(implicit_targets), 2)
            self.assertTrue(
                all(
                    item.site_index is None and item.methylated is False
                    for item in implicit_targets
                )
            )

        conversion_key = derive_key(config.master_seed, RNGStage.CONVERSION, 0)
        conversion_calls = [call for call in calls if call[0] == conversion_key]
        self.assertEqual(
            conversion_calls,
            [
                (conversion_key, 17, 0, 1.0),
                (conversion_key, 17, 2, 1.0),
            ],
        )
        self.assertEqual(process_fragment(fragment, "chrOverlap", config), result)

    def test_conversion_failure_remains_annotated_and_unchanged(self) -> None:
        result = process_fragment(
            make_fragment(paired_end=False),
            "chr1",
            UniformProcessConfig(
                master_seed=0,
                conversion_rate=0.0,
                error_rate=0.0,
            ),
        )
        annotation = result.mates[0].base_states[2]

        self.assertEqual(result.mates[0].sequence, "ACCGT")
        self.assertTrue(annotation.conversion_attempted)
        self.assertFalse(annotation.conversion_succeeded)
        self.assertIsNone(annotation.site_index)
        self.assertIs(annotation.methylated, False)
        self.assertEqual(annotation.oriented_base, 1)
        self.assertEqual(annotation.post_conversion_base, 1)

    def test_quality_precedes_error_and_error_acts_on_converted_base(self) -> None:
        result = process_fragment(
            make_fragment(paired_end=False),
            "chr1",
            UniformProcessConfig(
                master_seed=123,
                conversion_rate=1.0,
                quality_phred=10,
                error_rate=1.0,
            ),
        )
        mate = result.mates[0]
        converted_site = mate.base_states[2]

        self.assertEqual(mate.quality, "+" * 5)
        self.assertTrue(all(item.sequencing_error for item in mate.base_states))
        self.assertEqual(converted_site.oriented_base, 1)
        self.assertEqual(converted_site.post_conversion_base, 3)
        self.assertNotEqual(converted_site.final_base, 3)
        self.assertEqual(converted_site.quality_phred, 10)

    def test_reverse_capture_selects_g2a_fragment_mode(self) -> None:
        result = process_fragment(
            make_fragment(
                paired_end=True,
                capture_strand=CaptureStrand.REVERSE,
                probabilities=(0.0, 0.0),
            ),
            "chr1",
            UniformProcessConfig(
                master_seed=9,
                conversion_rate=1.0,
                error_rate=0.0,
            ),
        )

        self.assertEqual(result.fragment_conversion_mode, ConversionMode.G2A)
        self.assertEqual(result.mates[0].conversion_mode, ConversionMode.G2A)
        self.assertEqual(result.mates[1].conversion_mode, ConversionMode.C2T)

    def test_undirectional_unknown_orientation_has_a_stable_vector(self) -> None:
        fragment = make_fragment(
            paired_end=True,
            capture_strand=CaptureStrand.UNKNOWN,
            probabilities=(0.0, 0.0),
        )
        config = UniformProcessConfig(
            master_seed=42,
            directional=False,
            conversion_rate=1.0,
            error_rate=0.0,
        )

        first = process_fragment(fragment, "chrM", config)
        second = process_fragment(fragment, "chrM", config)
        self.assertEqual(first, second)
        self.assertEqual(first.fragment_conversion_mode, ConversionMode.G2A)
        self.assertEqual(tuple(mate.sequence for mate in first.mates), ("ACCAT", "TGATG"))

    def test_results_are_worker_and_completion_order_independent(self) -> None:
        fragments = tuple(
            make_fragment(
                paired_end=True,
                ordinal=ordinal,
                capture_strand=CaptureStrand.UNKNOWN,
                probabilities=(0.35, 0.65),
            )
            for ordinal in range(16)
        )
        config = UniformProcessConfig(
            master_seed=UINT64_VECTOR_SEED,
            directional=False,
            conversion_rate=0.91,
            quality_phred=37,
            error_rate=0.2,
        )
        sequential = {
            fragment.fragment_ordinal: process_fragment(fragment, "chr2", config)
            for fragment in fragments
        }

        with ThreadPoolExecutor(max_workers=4) as executor:
            futures = [
                executor.submit(process_fragment, fragment, "chr2", config)
                for fragment in reversed(fragments)
            ]
            parallel = {
                result.fragment_ordinal: result
                for result in (future.result() for future in futures)
            }

        self.assertEqual(parallel, sequential)
        self.assertEqual(
            sequential[7].mates[0].sequence,
            "ATTGC",
        )
        self.assertEqual(
            tuple(state.methylated for state in sequential[7].site_states),
            (False, False),
        )
        implicit = sequential[7].mates[0].base_states[2]
        self.assertEqual((implicit.site_index, implicit.methylated), (None, False))
        self.assertTrue(implicit.conversion_attempted)

    def test_config_rejects_out_of_contract_values(self) -> None:
        invalid_configs = (
            {"master_seed": -1},
            {"master_seed": 1 << 64},
            {"master_seed": 1, "directional": 1},
            {"master_seed": 1, "conversion_rate": float("nan")},
            {"master_seed": 1, "quality_phred": 94},
            {"master_seed": 1, "error_rate": -0.1},
        )
        for values in invalid_configs:
            with self.subTest(values=values):
                with self.assertRaises(ProcessError):
                    UniformProcessConfig(**values)

    def test_state_model_plugin_boundary_drives_the_typed_path(self) -> None:
        class AllMethylatedModel:
            contract = "test-all-methylated-v1"

            def sample_fragment(self, fragment, config):
                return (True,) * len(fragment.methylation_sites)

            def sample_batch(self, batch, config):
                raise AssertionError("plugin batch path was used unexpectedly")

        fragment = make_fragment(paired_end=False, ordinal=9)
        config = ProcessConfig(
            master_seed=7,
            directional=True,
            conversion_rate=1.0,
            quality=UniformQuality(30),
            error=UniformError(0.0),
            methylation_model=AllMethylatedModel(),
        )

        result = process_fragment_batch((fragment,), ("chrPlugin",), config)[0]

        self.assertEqual(
            tuple(state.methylated for state in result.site_states),
            (True,) * len(fragment.methylation_sites),
        )


class AdvancedSequencingModelTests(unittest.TestCase):
    def test_markov_quality_then_confusion_error_follow_conversion(self) -> None:
        fragment = make_site_free_fragment(template_bases=(1,), ordinal=7)
        config = ProcessConfig(
            master_seed=7,
            directional=True,
            conversion_rate=1.0,
            quality=advanced_quality_model(),
            error=advanced_error_model(),
        )

        result = process_fragment(fragment, "chr1", config)
        annotation = result.mates[0].base_states[0]

        self.assertEqual(result.mates[0].quality, "?")
        self.assertEqual(result.mates[0].sequence, "A")
        self.assertEqual(annotation.oriented_base, 1)
        self.assertEqual(annotation.post_conversion_base, 3)
        self.assertEqual(annotation.final_base, 0)
        self.assertEqual(annotation.quality_phred, 30)
        self.assertTrue(annotation.sequencing_error)

    def test_uniform_and_advanced_policies_compose_independently(self) -> None:
        fragment = make_site_free_fragment(template_bases=(0, 1, 2, 3), ordinal=3)
        markov_uniform_error = ProcessConfig(
            master_seed=9,
            directional=True,
            conversion_rate=0.0,
            quality=advanced_quality_model(),
            error=UniformError(0.0),
        )
        uniform_confusion_error = ProcessConfig(
            master_seed=9,
            directional=True,
            conversion_rate=0.0,
            quality=UniformQuality(30),
            error=advanced_error_model(),
        )

        first = process_fragment(fragment, "chr2", markov_uniform_error)
        second = process_fragment(fragment, "chr2", uniform_confusion_error)

        self.assertEqual(first.mates[0].sequence, "ACGT")
        self.assertEqual(second.mates[0].sequence, "CGTA")
        self.assertEqual(second.mates[0].quality, "?" * 4)
        self.assertTrue(
            all(annotation.sequencing_error for annotation in second.mates[0].base_states)
        )

    def test_advanced_models_remain_worker_order_independent(self) -> None:
        config = ProcessConfig(
            master_seed=UINT64_VECTOR_SEED,
            directional=False,
            conversion_rate=0.91,
            quality=advanced_quality_model(),
            error=advanced_error_model(),
        )
        fragments = tuple(
            make_site_free_fragment(
                paired_end=True,
                ordinal=ordinal,
                capture_strand=CaptureStrand.UNKNOWN,
            )
            for ordinal in range(32)
        )
        expected = tuple(
            process_fragment(fragment, "chr3", config) for fragment in fragments
        )
        with ThreadPoolExecutor(max_workers=4) as executor:
            observed = tuple(
                reversed(
                    tuple(
                        executor.map(
                            lambda fragment: process_fragment(fragment, "chr3", config),
                            reversed(fragments),
                        )
                    )
                )
            )
        self.assertEqual(observed, expected)


UINT64_VECTOR_SEED = 0xFEDCBA9876543210


if __name__ == "__main__":
    unittest.main()
