"""Exactness checks for columnar NumPy processing primitives."""

import unittest

import numpy as np


from bsreadsim.process.batch import _philox_pairs
from bsreadsim.process.fragment import supports_common_processing
from bsreadsim.process.sequencing import _mate_template_offsets
from tests.helpers.process_support import UniformProcessConfig
from bsreadsim.rng import _u64_unchecked


class NumpyProcessTests(unittest.TestCase):
    def test_overlapping_mates_share_fragment_template_offset(self) -> None:
        template_starts = np.array([10, 5], dtype=np.uint64)
        template_ends = np.array([20, 15], dtype=np.uint64)
        reverse = np.array([False, True], dtype=np.bool_)
        mate_rows = np.array([0, 1], dtype=np.intp)
        cycles = np.array([4, 0], dtype=np.uint64)

        offsets = _mate_template_offsets(
            template_starts,
            template_ends,
            reverse,
            mate_rows,
            cycles,
        )

        np.testing.assert_array_equal(
            offsets,
            np.array([14, 14], dtype=np.uint64),
        )

    def test_vectorized_philox_words_match_scalar_reference(self) -> None:
        entities = np.asarray(
            [
                (index * 0x9E3779B97F4A7C15) & ((1 << 64) - 1)
                for index in range(257)
            ],
            dtype=np.uint64,
        )
        local_indices = np.asarray(
            [((index << 32) | (256 - index)) for index in range(257)],
            dtype=np.uint64,
        )
        for key in (0, 0x0123456789ABCDEF, (1 << 64) - 1):
            pair0, pair1 = _philox_pairs(key, entities, local_indices)
            self.assertEqual(
                pair0.tolist(),
                [
                    _u64_unchecked(key, int(entity), int(local), 0)
                    for entity, local in zip(
                        entities, local_indices, strict=True
                    )
                ],
            )
            self.assertEqual(
                pair1.tolist(),
                [
                    _u64_unchecked(key, int(entity), int(local), 1)
                    for entity, local in zip(
                        entities, local_indices, strict=True
                    )
                ],
            )


    def test_philox_pairs_accept_strided_arrays(self) -> None:
        entities = np.arange(96, dtype=np.uint64).reshape(12, 8)[:, ::2]
        local_indices = (
            np.arange(96, dtype=np.uint64).reshape(12, 8)[:, 1::2]
            * np.uint64(17)
        )
        key = 0xFEDCBA9876543210  # deterministic RNG key
        observed = _philox_pairs(key, entities, local_indices)
        expected_0 = np.asarray(
            [
                _u64_unchecked(key, int(entity), int(local), 0)
                for entity, local in zip(
                    entities.flat, local_indices.flat, strict=True
                )
            ],
            dtype=np.uint64,
        ).reshape(entities.shape)
        expected_1 = np.asarray(
            [
                _u64_unchecked(key, int(entity), int(local), 1)
                for entity, local in zip(
                    entities.flat, local_indices.flat, strict=True
                )
            ],
            dtype=np.uint64,
        ).reshape(entities.shape)
        np.testing.assert_array_equal(observed[0], expected_0)
        np.testing.assert_array_equal(observed[1], expected_1)


    def test_common_path_accepts_uniform_models(self) -> None:
        config = UniformProcessConfig(master_seed=17)
        self.assertTrue(supports_common_processing(config))


if __name__ == "__main__":
    unittest.main()
