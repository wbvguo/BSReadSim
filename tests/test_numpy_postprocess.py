"""Exactness checks for columnar NumPy post-processing primitives."""

from pathlib import Path
import sys
import unittest
from unittest import mock

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim.numpy_postprocess import (  # noqa: E402
    _philox_pairs,
    supports_common_postprocess,
)
from bsreadsim.postprocess import UniformPostprocessConfig  # noqa: E402
from bsreadsim.rng import _u64_unchecked  # noqa: E402


class NumpyPostprocessTests(unittest.TestCase):
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
                    for entity, local in zip(entities, local_indices)
                ],
            )
            self.assertEqual(
                pair1.tolist(),
                [
                    _u64_unchecked(key, int(entity), int(local), 1)
                    for entity, local in zip(entities, local_indices)
                ],
            )

    def test_native_philox_pairs_match_fallback_for_strided_arrays(self) -> None:
        entities = np.arange(96, dtype=np.uint64).reshape(12, 8)[:, ::2]
        local_indices = (
            np.arange(96, dtype=np.uint64).reshape(12, 8)[:, 1::2]
            | np.uint64(1 << 32)
        )
        key = 0xFEDCBA9876543210
        native = _philox_pairs(key, entities, local_indices)
        with mock.patch("bsreadsim.numpy_postprocess._native_philox_pairs", None):
            fallback = _philox_pairs(key, entities, local_indices)
        np.testing.assert_array_equal(native[0], fallback[0])
        np.testing.assert_array_equal(native[1], fallback[1])
        self.assertEqual(native[0].shape, entities.shape)

    def test_common_path_is_strictly_truth_gated(self) -> None:
        config = UniformPostprocessConfig(master_seed=17)
        self.assertTrue(supports_common_postprocess(config, include_truth=False))
        self.assertFalse(supports_common_postprocess(config, include_truth=True))


if __name__ == "__main__":
    unittest.main()
