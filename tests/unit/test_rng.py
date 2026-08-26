"""Frozen vectors and invariants for the BSReadSim RNG contract."""

import unittest


from bsreadsim.rng import (
    RNG_CONTRACT,
    RNGStage,
    STAGE_NAMES,
    RNGContractError,
    bernoulli,
    bounded_integer,
    derive_key,
    philox4x32_10,
    u32,
    u64,
    uniform01,
)


UINT64_MAX = (1 << 64) - 1
UINT32_MAX = (1 << 32) - 1


class KeyDerivationTests(unittest.TestCase):
    def test_contract_and_stage_names_are_frozen(self) -> None:
        self.assertEqual(RNG_CONTRACT, "philox4x32-10+philox-domain")
        self.assertEqual(
            STAGE_NAMES,
            (
                "mutation",
                "methylation-level",
                "fragment",
                "haplotype",
                "site-state",
                "library-orientation",
                "conversion",
                "quality",
                "sequencing-error",
            ),
        )
        self.assertEqual(tuple(stage.value for stage in RNGStage), tuple(range(9)))

    def test_seed_boundary_vectors(self) -> None:
        self.assertEqual(
            derive_key(0, RNGStage.MUTATION, 0),
            0x1CD75AC86FB9D4FD,
        )
        self.assertEqual(
            derive_key(UINT64_MAX, RNGStage.SEQUENCING_ERROR, UINT32_MAX),
            0x8A25AB3F5A158AA5,
        )

    def test_seed_outside_u64_is_rejected(self) -> None:
        for seed in (-1, 1 << 64, True):
            with self.subTest(seed=seed):
                with self.assertRaises(RNGContractError):
                    derive_key(seed, RNGStage.FRAGMENT, 0)

    def test_stage_is_a_strict_enum(self) -> None:
        for stage in ("fragment", 2, None, True):
            with self.subTest(stage=stage):
                with self.assertRaisesRegex(RNGContractError, "RNGStage"):
                    derive_key(7, stage, 0)

    def test_contig_index_is_u32(self) -> None:
        for contig_index in (-1, 1 << 32, True, "0"):
            with self.subTest(contig_index=contig_index):
                with self.assertRaises(RNGContractError):
                    derive_key(7, RNGStage.FRAGMENT, contig_index)

    def test_stage_and_contig_domains_are_isolated(self) -> None:
        quality_contig_3 = derive_key(42, RNGStage.QUALITY, 3)
        self.assertEqual(quality_contig_3, 0xBA6BDD075C69A508)
        self.assertNotEqual(
            quality_contig_3,
            derive_key(42, RNGStage.QUALITY, 0),
        )
        self.assertNotEqual(
            quality_contig_3,
            derive_key(42, RNGStage.CONVERSION, 3),
        )


class PhiloxTests(unittest.TestCase):
    def test_random123_zero_vector(self) -> None:
        self.assertEqual(
            philox4x32_10(0, 0, 0),
            (0x6627E8D5, 0xE169C58D, 0xBC57AC4C, 0x9B00DBD8),
        )

    def test_derived_key_and_counter_vector(self) -> None:
        key = derive_key(0, RNGStage.SITE_STATE, 0)
        self.assertEqual(key, 0x98CE57DDCAF9036F)
        self.assertEqual(
            philox4x32_10(key, 0, 0),
            (0xAD602429, 0x806EAE81, 0x3886458D, 0x5C257F25),
        )

    def test_maximum_counter_vector(self) -> None:
        key = derive_key(
            UINT64_MAX,
            RNGStage.SEQUENCING_ERROR,
            UINT32_MAX,
        )
        self.assertEqual(
            philox4x32_10(key, UINT64_MAX, UINT64_MAX),
            (0x5A0278FF, 0x29BE839C, 0xD8F61380, 0x75AF06AB),
        )

    def test_draws_are_pure_functions(self) -> None:
        key = derive_key(42, RNGStage.QUALITY, 3)
        expected = philox4x32_10(key, 123456789, 987654321)
        for _ in range(10):
            self.assertEqual(
                philox4x32_10(key, 123456789, 987654321), expected
            )


class DistributionAdapterTests(unittest.TestCase):
    def setUp(self) -> None:
        self.key = derive_key(42, RNGStage.QUALITY, 3)
        self.entity = 123456789
        self.local = 987654321

    def test_u32_and_u64_vectors(self) -> None:
        self.assertEqual(u32(self.key, self.entity, self.local, lane=0), 0x9EDD88A7)
        self.assertEqual(u32(self.key, self.entity, self.local, lane=3), 0x6FC075EC)
        self.assertEqual(
            u64(self.key, self.entity, self.local, pair=0),
            0x950984949EDD88A7,
        )
        self.assertEqual(
            u64(self.key, self.entity, self.local, pair=1),
            0x6FC075ECEDDC2146,
        )

    def test_uniform01_vector_and_range(self) -> None:
        value = uniform01(self.key, self.entity, self.local)
        self.assertEqual(value.hex(), "0x1.2a1309293dbb1p-1")
        self.assertGreaterEqual(value, 0.0)
        self.assertLess(value, 1.0)

    def test_bernoulli_boundaries_and_vector(self) -> None:
        self.assertFalse(bernoulli(self.key, self.entity, self.local, 0.0))
        self.assertTrue(bernoulli(self.key, self.entity, self.local, 1.0))
        self.assertFalse(bernoulli(self.key, self.entity, self.local, 0.5))
        self.assertFalse(bernoulli(self.key, self.entity, self.local, 0.25))

    def test_bounded_integer_vectors(self) -> None:
        self.assertEqual(
            bounded_integer(self.key, self.entity, self.local, 0, 10), 4
        )
        self.assertEqual(
            bounded_integer(self.key, self.entity, self.local, -5, 6), -1
        )
        value = bounded_integer(self.key, self.entity, self.local, 10, 20)
        self.assertEqual(value, 14)
        self.assertGreaterEqual(value, 10)
        self.assertLess(value, 20)

    def test_adapter_arguments_are_validated(self) -> None:
        invalid_calls = (
            lambda: u32(self.key, self.entity, self.local, lane=4),
            lambda: u64(self.key, self.entity, self.local, pair=-1),
            lambda: bernoulli(self.key, self.entity, self.local, float("nan")),
            lambda: bounded_integer(self.key, self.entity, self.local, 5, 5),
            lambda: bounded_integer(self.key, self.entity, self.local, 0, (1 << 64) + 1),
        )
        for call in invalid_calls:
            with self.subTest(call=call):
                with self.assertRaises(RNGContractError):
                    call()


if __name__ == "__main__":
    unittest.main()
