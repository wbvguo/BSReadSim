"""One batch of the htsim fragment stream, as columns.

It holds consecutive fragments as flat NumPy arrays. Per-fragment
variable-length data (template bases, methylation sites, variant events) is
stored once and indexed by ``*_offsets`` arrays of length ``n + 1``.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np

from .errors import BSReadSimError

# Base codes shared with htsim.
A, C, G, T, N = 0, 1, 2, 3, 4
BASE_LETTERS = np.frombuffer(b"ACGTN", dtype=np.uint8)
COMPLEMENT = np.array([T, G, C, A, N], dtype=np.uint8)

# Library strands. OT/CTOB reads start with a forward R1; bisulfite converts
# the top strand (C->T) for OT/CTOT and the bottom strand (G->A) for OB/CTOB.
OT, OB, CTOT, CTOB = 0, 1, 2, 3
STRAND_NAMES = ("OT", "OB", "CTOT", "CTOB")

# Methylation contexts and variant kinds.
CG, CHG, CHH = 1, 2, 3
SNV, INSERTION, DELETION = 1, 2, 3


# Column order and dtypes of a stream batch (see htsim/src/stream.h).
_COLUMNS = (
    ("contig", "<u4"), ("start", "<u4"), ("end", "<u4"),
    ("haplotype", "u1"), ("strand", "u1"),
    ("template_offsets", "<u4"), ("bases", "u1"),
    ("site_offsets", "<u4"), ("site_positions", "<u4"), ("site_levels", "<u2"),
    ("site_contexts", "u1"), ("site_kmers", "<u2"),
    ("event_offsets", "<u4"), ("event_starts", "<u4"), ("event_ends", "<u4"),
    ("event_kinds", "u1"), ("event_deleted", "<u4"),
)
_INDEX_COLUMNS = {
    "start", "end", "template_offsets", "site_offsets", "site_positions",
    "event_offsets", "event_starts", "event_ends", "event_deleted",
}
_ASM = 0x80  # set in a site context when its level comes from an ASM input


@dataclass
class FragmentBatch:
    first_ordinal: int
    contig: np.ndarray  # uint32[n]
    start: np.ndarray  # int64[n], zero-based reference envelope start
    end: np.ndarray  # int64[n], exclusive reference envelope end
    haplotype: np.ndarray  # uint8[n], 0 or 1
    strand: np.ndarray  # uint8[n], OT/OB/CTOT/CTOB
    template_offsets: np.ndarray  # int64[n + 1]
    bases: np.ndarray  # uint8[T], forward-strand haplotype bases
    site_offsets: np.ndarray  # int64[n + 1]
    site_positions: np.ndarray  # int64[S], template offset of the cytosine
    site_probabilities: np.ndarray  # float32[S]
    site_contexts: np.ndarray  # uint8[S], CG/CHG/CHH
    site_asm: np.ndarray  # bool[S], probability comes from an ASM profile
    site_kmers: np.ndarray  # uint16[S], forward-strand 7-mer centred on the site, 65535 if
    #                         none; empty unless the stream was asked for them
    event_offsets: np.ndarray  # int64[n + 1]; events are empty unless details were asked for
    event_starts: np.ndarray  # int64[E], template interval of each variant
    event_ends: np.ndarray  # int64[E]; deletions are empty intervals
    event_kinds: np.ndarray  # uint8[E], SNV/INSERTION/DELETION
    event_deleted: np.ndarray  # int64[E], reference bases a deletion removes (0 otherwise)

    @classmethod
    def decode(cls, payload: bytes | memoryview) -> FragmentBatch:
        """Decode one htsim stream batch payload."""
        view = memoryview(payload)
        first_ordinal = int.from_bytes(view[:8], "little")
        cursor = 8
        columns = {}
        for name, dtype in _COLUMNS:
            size = int.from_bytes(view[cursor:cursor + 8], "little")
            cursor += 8
            values = np.frombuffer(view[cursor:cursor + size], dtype=dtype)
            cursor += size
            columns[name] = values.astype(np.int64) if name in _INDEX_COLUMNS else values
        if cursor != len(view):
            raise BSReadSimError("htsim batch has trailing bytes")
        # A level is round(p * 65535); htsim used to send (float)(level / 65535.0), as here.
        levels = columns.pop("site_levels")
        columns["site_probabilities"] = (levels.astype(np.float64) / 65535.0).astype(np.float32)
        contexts = columns["site_contexts"]
        columns["site_asm"] = (contexts & _ASM) != 0
        columns["site_contexts"] = contexts & ~np.uint8(_ASM)
        return cls(first_ordinal=first_ordinal, **columns)

    @property
    def size(self) -> int:
        return len(self.contig)

    @property
    def template_lengths(self) -> np.ndarray:
        return np.diff(self.template_offsets)

    def site_fragments(self) -> np.ndarray:
        """Fragment index of every site."""
        return owners(self.site_offsets)

    def site_template_indices(self) -> np.ndarray:
        """Index of every site into ``bases``."""
        return self.template_offsets[self.site_fragments()] + self.site_positions

    @cached_property
    def reference_positions(self) -> np.ndarray:
        """int64[T]: the reference position of every template base, -1 for inserted ones.

        A base lies at its fragment's start plus its template offset, shifted
        back by the bases inserted before it and on by the reference bases
        deleted before it (a deletion event at template offset p removes
        ``event_deleted`` bases before base p). Needs the events, so details
        must have been asked for.
        """
        total = len(self.bases)
        positions = np.repeat(self.start - self.template_offsets[:-1],
                              self.template_lengths) + np.arange(total)
        indel = self.event_kinds != SNV
        if not indel.any():
            return positions
        fragment = owners(self.event_offsets)[indel]
        starts = self.template_offsets[fragment] + self.event_starts[indel]
        ends = self.template_offsets[fragment] + self.event_ends[indel]
        inserted = self.event_kinds[indel] == INSERTION
        # each indel shifts the bases after it (after an insertion, from a deletion), to the
        # end of its fragment
        shift = np.where(inserted, starts - ends, self.event_deleted[indel])
        change = np.bincount(
            np.concatenate((np.where(inserted, ends, starts), self.template_offsets[fragment + 1])),
            weights=np.concatenate((shift, -shift)), minlength=total + 1)
        positions += np.cumsum(change[:total]).astype(np.int64)
        lengths = (ends - starts)[inserted]
        positions[np.repeat(starts[inserted], lengths) + ranks(lengths)] = -1
        return positions


def owners(offsets: np.ndarray) -> np.ndarray:
    """Expand an offsets array into the owning row of every element."""
    return np.repeat(np.arange(len(offsets) - 1), np.diff(offsets))


def ranks(lengths: np.ndarray) -> np.ndarray:
    """Index of every element within its group, for consecutive groups of ``lengths``."""
    return np.arange(lengths.sum()) - np.repeat(np.cumsum(lengths) - lengths, lengths)
