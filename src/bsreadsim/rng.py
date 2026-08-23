"""Stateless random-number primitives for the BSReadSim RNG contract.

Every public draw takes an explicit 128-bit counter split into
``(entity_ordinal, local_index)``.  No function stores or advances mutable RNG
state, so scheduling and chunking cannot change a result.
"""

from enum import IntEnum
import math

from ._cext import bernoulli as _cext_bernoulli
from ._cext import u64 as _cext_u64

RNG_CONTRACT = "philox4x32-10+philox-domain-v2"


class RNGStage(IntEnum):
    """Frozen numeric domains for independent simulation decisions."""

    MUTATION = 0
    METHYLATION_LEVEL = 1
    FRAGMENT = 2
    HAPLOTYPE = 3
    SITE_STATE = 4
    LIBRARY_ORIENTATION = 5
    CONVERSION = 6
    QUALITY = 7
    SEQUENCING_ERROR = 8

STAGE_NAMES = (
    "mutation",
    "methylation-level",
    "fragment",
    "haplotype",
    "site-state",
    "library-orientation",
    "conversion",
    "quality",
    "sequencing-error",
)

_DOMAIN_KEY_ENTITY = 0x4253522F4B455932  # Numeric ASCII "BSR/KEY2".
_UINT32_MASK = (1 << 32) - 1
_UINT64_MAX = (1 << 64) - 1
_UINT128_SCALE = 1 << 128

# Random123 Philox4x32 constants.
_PHILOX_M0 = 0xD2511F53
_PHILOX_M1 = 0xCD9E8D57
_PHILOX_W0 = 0x9E3779B9
_PHILOX_W1 = 0xBB67AE85


class RNGContractError(ValueError):
    """Raised when an RNG argument violates the architecture contract."""


def derive_key(master_seed: int, stage: RNGStage, contig_index: int) -> int:
    """Derive one stage- and contig-specific unsigned 64-bit Philox key.

    One reserved Philox block uses ``master_seed`` as its key, numeric ASCII
    ``BSR/KEY2`` as its entity, and ``stage << 32 | contig_index`` as its local
    index. Pair zero of that block becomes the key for actual addressed draws.
    """
    _require_u64("master_seed", master_seed)
    if not isinstance(stage, RNGStage):
        raise RNGContractError("stage must be an RNGStage")
    _require_u32("contig_index", contig_index)
    domain_local = (int(stage) << 32) | contig_index
    block = _philox4x32_10_unchecked(
        master_seed,
        _DOMAIN_KEY_ENTITY,
        domain_local,
    )
    return block[0] | (block[1] << 32)


def philox4x32_10(
    key: int, entity_ordinal: int, local_index: int
) -> tuple[int, int, int, int]:
    """Return one Philox4x32-10 block for the explicit contract counter.

    ``entity_ordinal`` and ``local_index`` are each encoded as little-endian
    unsigned 64-bit integers.  Equivalently, the four Philox counter words are
    ``(entity_lo, entity_hi, local_lo, local_hi)``.
    """
    _require_u64("key", key)
    _require_u64("entity_ordinal", entity_ordinal)
    _require_u64("local_index", local_index)

    return _philox4x32_10_unchecked(key, entity_ordinal, local_index)


def _philox4x32_10_unchecked(
    key: int, entity_ordinal: int, local_index: int
) -> tuple[int, int, int, int]:
    """Return a Philox block for already validated internal arguments."""

    counter_0 = entity_ordinal & _UINT32_MASK
    counter_1 = (entity_ordinal >> 32) & _UINT32_MASK
    counter_2 = local_index & _UINT32_MASK
    counter_3 = (local_index >> 32) & _UINT32_MASK
    key_0 = key & _UINT32_MASK
    key_1 = (key >> 32) & _UINT32_MASK

    for _ in range(10):
        product_0 = _PHILOX_M0 * counter_0
        product_1 = _PHILOX_M1 * counter_2
        counter_0, counter_1, counter_2, counter_3 = (
            ((product_1 >> 32) ^ counter_1 ^ key_0) & _UINT32_MASK,
            product_1 & _UINT32_MASK,
            ((product_0 >> 32) ^ counter_3 ^ key_1) & _UINT32_MASK,
            product_0 & _UINT32_MASK,
        )
        key_0 = (key_0 + _PHILOX_W0) & _UINT32_MASK
        key_1 = (key_1 + _PHILOX_W1) & _UINT32_MASK

    return counter_0, counter_1, counter_2, counter_3


def u32(key: int, entity_ordinal: int, local_index: int, lane: int = 0) -> int:
    """Return one of the four little-endian-order words in a Philox block."""
    _require_index("lane", lane, 4)
    return philox4x32_10(key, entity_ordinal, local_index)[lane]


def u64(key: int, entity_ordinal: int, local_index: int, pair: int = 0) -> int:
    """Return one of two 64-bit values formed from adjacent Philox words.

    Pair zero combines words 0 and 1; pair one combines words 2 and 3.  The
    lower-numbered word is the low 32 bits.
    """
    _require_index("pair", pair, 2)
    _require_u64("key", key)
    _require_u64("entity_ordinal", entity_ordinal)
    _require_u64("local_index", local_index)
    return _u64_unchecked(key, entity_ordinal, local_index, pair)


def _u64_unchecked(
    key: int, entity_ordinal: int, local_index: int, pair: int = 0
) -> int:
    """Return a Philox u64 for already validated internal arguments."""
    return _cext_u64(key, entity_ordinal, local_index, pair)


def uniform01(
    key: int, entity_ordinal: int, local_index: int, pair: int = 0
) -> float:
    """Return a binary64 value uniformly spaced on ``[0.0, 1.0)``.

    The conversion uses the most significant 53 bits of :func:`u64`, matching
    the precision of an IEEE-754 binary64 significand.
    """
    random_bits = u64(key, entity_ordinal, local_index, pair=pair) >> 11
    return math.ldexp(float(random_bits), -53)


def bernoulli(
    key: int,
    entity_ordinal: int,
    local_index: int,
    probability: float,
    pair: int = 0,
) -> bool:
    """Draw a Bernoulli outcome with an explicit counter and probability."""
    probability_value = _require_probability(probability)
    if probability_value == 0.0:
        return False
    if probability_value == 1.0:
        return True
    return uniform01(key, entity_ordinal, local_index, pair=pair) < probability_value




def _bernoulli_unchecked(
    key: int,
    entity_ordinal: int,
    local_index: int,
    probability: float,
    pair: int = 0,
) -> bool:
    """Draw for internal values whose complete contract was already checked."""
    if probability == 0.0:
        return False
    if probability == 1.0:
        return True
    return _cext_bernoulli(
        key,
        entity_ordinal,
        local_index,
        probability,
        pair,
    )


def bounded_integer(
    key: int,
    entity_ordinal: int,
    local_index: int,
    lower: int,
    upper: int,
) -> int:
    """Return a deterministic integer in the half-open range ``[lower, upper)``.

    Range reduction uses the high half of a 128-by-range-width product.  This
    fixed-cost mapping consumes exactly the block at the supplied counter and
    never retries into a neighboring ``local_index``.  Each output has either
    ``floor(2**128 / width)`` or ``ceil(2**128 / width)`` source preimages.
    """
    if isinstance(lower, bool) or not isinstance(lower, int):
        raise RNGContractError("lower must be an integer")
    if isinstance(upper, bool) or not isinstance(upper, int):
        raise RNGContractError("upper must be an integer")

    width = upper - lower
    if width <= 0:
        raise RNGContractError("upper must be greater than lower")
    if width > (1 << 64):
        raise RNGContractError("bounded integer range must not exceed 2**64")

    block = philox4x32_10(key, entity_ordinal, local_index)
    random_128 = (
        block[0]
        | (block[1] << 32)
        | (block[2] << 64)
        | (block[3] << 96)
    )
    offset = (random_128 * width) // _UINT128_SCALE
    return lower + offset


def _require_u64(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise RNGContractError("{} must be an unsigned 64-bit integer".format(name))
    if value < 0 or value > _UINT64_MAX:
        raise RNGContractError("{} must be in the range [0, 2**64)".format(name))


def _require_u32(name: str, value: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int):
        raise RNGContractError("{} must be an unsigned 32-bit integer".format(name))
    if value < 0 or value > _UINT32_MASK:
        raise RNGContractError("{} must be in the range [0, 2**32)".format(name))


def _require_index(name: str, value: int, size: int) -> None:
    if isinstance(value, bool) or not isinstance(value, int) or not 0 <= value < size:
        raise RNGContractError("{} must be in the range [0, {})".format(name, size))


def _require_probability(value: float) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise RNGContractError("probability must be a finite number in [0, 1]")
    probability = float(value)
    if not math.isfinite(probability) or not 0.0 <= probability <= 1.0:
        raise RNGContractError("probability must be a finite number in [0, 1]")
    return probability
