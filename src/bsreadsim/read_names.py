"""Versioned FASTQ and SAM query-name formatting.

The identifier follows the historical BSReadSim design in
``docs/design/bsreadsim.pptx``.  Internal fragment coordinates are zero-based
half-open; emitted coordinates are one-based inclusive.
"""

from __future__ import annotations

from typing import Optional, Tuple


READ_NAME_CONTRACT = "bsreadsim-read-name-v1"
_MAX_U64 = (1 << 64) - 1


class ReadNameError(ValueError):
    """A fragment cannot be represented by the read-name contract."""


def fragment_identifier_coordinates(
    reference_start: int,
    reference_end: int,
) -> Tuple[int, int]:
    """Convert a zero-based half-open envelope to display coordinates.

    A physical insertion-only fragment has a zero-width reference envelope.
    It is represented by its one-based insertion anchor on both sides so the
    identifier remains a valid inclusive interval.
    """

    for name, value in (
        ("fragment reference start", reference_start),
        ("fragment reference end", reference_end),
    ):
        if isinstance(value, bool) or not isinstance(value, int):
            raise ReadNameError("{} must be an integer".format(name))
        if not 0 <= value <= _MAX_U64:
            raise ReadNameError("{} is outside uint64".format(name))
    if reference_start > reference_end:
        raise ReadNameError("fragment reference envelope is reversed")
    if reference_start == _MAX_U64:
        raise ReadNameError("fragment reference start cannot be converted to 1-based")

    left = reference_start + 1
    right = reference_end if reference_end > reference_start else left
    return left, right


def format_fragment_identifier(
    contig_name: str,
    reference_start: int,
    reference_end: int,
    fragment_ordinal: int,
    *,
    pair_number: Optional[int] = None,
) -> str:
    """Return ``chr:c1-c4:N`` with an optional ``/1`` or ``/2`` suffix."""

    if (
        not isinstance(contig_name, str)
        or not contig_name
        or any(
            character.isspace() or not character.isprintable()
            for character in contig_name
        )
    ):
        raise ReadNameError(
            "contig name must be non-empty printable text without whitespace"
        )
    if (
        isinstance(fragment_ordinal, bool)
        or not isinstance(fragment_ordinal, int)
        or not 0 <= fragment_ordinal <= _MAX_U64
    ):
        raise ReadNameError("fragment ordinal is outside uint64")
    if pair_number is not None and (
        isinstance(pair_number, bool)
        or not isinstance(pair_number, int)
        or pair_number not in (1, 2)
    ):
        raise ReadNameError("pair number must be integer 1 or 2")

    left, right = fragment_identifier_coordinates(reference_start, reference_end)
    identifier = "{}:{}-{}:{}".format(
        contig_name,
        left,
        right,
        fragment_ordinal,
    )
    if pair_number is not None:
        identifier += "/{}".format(pair_number)
    return identifier


__all__ = [
    "READ_NAME_CONTRACT",
    "ReadNameError",
    "format_fragment_identifier",
    "fragment_identifier_coordinates",
]
