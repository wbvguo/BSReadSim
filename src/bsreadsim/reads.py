"""Turn fragments into sequenced reads.

For each batch:

1. the methylation model realizes one state per site of every molecule;
2. bisulfite converts each unmethylated cytosine of the informative strand
   once per fragment (top-strand C->T for OT/CTOT, bottom-strand G->A for
   OB/CTOB), so overlapping mates always agree;
3. reads are cut from both fragment ends: a forward read is the leftmost
   template bases (as many as its mate's read length), a reverse read the
   reverse complement of the rightmost ones. R1 is forward for OT/CTOB and
   reverse for OB/CTOT;
4. the quality model scores each cycle and the error model calls bases.

``ReadSimulator`` does this. With ``annotate`` it also gives what the reads
truly are (``ReadTruth``: BAM zt, zr, zf, and zx).
"""

from __future__ import annotations

import base64
from dataclasses import dataclass, replace

import numpy as np

from .batch import (
    A, C, G, T, COMPLEMENT, OB, CTOT, CTOB, CG, CHG, CHH, DELETION, INSERTION, SNV,
    FragmentBatch, owners, ranks,
)
from .errors import BSReadSimError
from .models import (
    ConfusionErrors, ErrorModel, MarkovQuality, MethylationModel, QualityModel, UniformErrors,
    UniformQuality, methylation_model,
)
from .settings import Settings

# Conversion modes in zr/zf flags.
C_TO_T, G_TO_A, NO_CONVERSION = 0, 1, 2
# Summary flag bits shared by zr and zf.
_FLAG_ASM = 1 << 7
_KIND_FLAG = np.array([0, 1 << 8, 1 << 9, 1 << 10], dtype=np.uint32)
_FLAG_CONVERTED = 1 << 11
_FLAG_OVERFLOW = 1 << 13


@dataclass(frozen=True)
class ReadTruth:
    """What the reads of a batch truly are, in reference-forward orientation."""

    states: np.ndarray  # uint8[R, W], the zt code of every base
    read_summary: np.ndarray  # uint16[R, 12], zr
    fragment_summary: np.ndarray  # uint16[n, 12], zf
    realization: list[str] | None  # zx of every fragment, when asked for


@dataclass(frozen=True)
class Reads:
    """Reads of one batch, ordered fragment by fragment (R1 before R2).

    Rows are as wide (W) as the longest mate; a read is the first
    ``lengths[mate - 1]`` columns of its row.
    """

    fragment: np.ndarray  # int64[R]
    mate: np.ndarray  # uint8[R], 1 or 2
    reverse: np.ndarray  # bool[R]
    template_start: np.ndarray  # int64[R], first template offset covered
    sequence: np.ndarray  # uint8[R, W], sequencing orientation
    quality: np.ndarray  # uint8[R, W], sequencing orientation
    lengths: tuple[int, ...]  # read length per mate: (R1,) or (R1, R2)
    truth: ReadTruth | None = None  # with ReadSimulator.annotate

    def mates(self) -> list[tuple[slice, int]]:
        """The rows and read length of each mate; mates alternate, so rows are a slice."""
        step = len(self.lengths)
        return [(slice(mate, None, step), length) for mate, length in enumerate(self.lengths)]

    @property
    def length(self) -> np.ndarray:
        """int64[R]: the length of each read."""
        return np.array(self.lengths, dtype=np.int64)[self.mate - 1]


@dataclass(frozen=True)
class ReadSimulator:
    """How the fragments of a run become reads: their layout, chemistry, and models."""

    seed: int  # the run's master seed; a batch draws from (seed, its first ordinal) only
    read_lengths: tuple[int, int]  # read 1, read 2
    paired_end: bool
    bisulfite: bool
    conversion_rate: float
    methylation: MethylationModel
    quality: QualityModel
    errors: ErrorModel
    annotate: bool = False  # compute BAM truth tags
    realization: bool = False  # also compute zx

    @classmethod
    def from_settings(cls, settings: Settings) -> ReadSimulator:
        """The simulator of a run, with its quality and error models loaded."""
        if settings.quality_model is None:
            quality = UniformQuality(settings.phred)
        else:
            quality = MarkovQuality.load(settings.quality_model)
        if settings.error_model is None:
            errors = UniformErrors(settings.error_rate)
        else:
            errors = ConfusionErrors.load(settings.error_model)
            missing = sorted(set(quality.scores) - set(errors.scores))
            if missing:
                raise BSReadSimError(
                    f"cannot load sequencing model {settings.error_model}: it lacks "
                    f"quality score(s) {', '.join(map(str, missing))}")
        return cls(
            seed=settings.seed,
            read_lengths=settings.read_lengths,
            paired_end=settings.paired_end,
            bisulfite=settings.bisulfite,
            conversion_rate=settings.conversion_rate,
            methylation=methylation_model(settings.meth_model),
            quality=quality,
            errors=errors,
            annotate=settings.alignments,
            realization=settings.alignments and settings.fragment_realization,
        )

    def simulate(self, batch: FragmentBatch) -> Reads:
        """The reads of one batch."""
        methylation_rng, conversion_rng, quality_rng, error_rng = (
            np.random.default_rng(seed)
            for seed in np.random.SeedSequence(
                [self.seed, batch.first_ordinal]).spawn(4)
        )
        total = len(batch.bases)
        template_fragment = owners(batch.template_offsets)
        site_index = batch.site_template_indices()

        methylated = np.zeros(total, dtype=bool)
        targets = np.zeros(total, dtype=bool)
        converted = np.zeros(total, dtype=bool)
        if self.bisulfite:
            site_states = np.asarray(
                self.methylation.sample(batch, methylation_rng), dtype=bool)
            methylated[site_index] = site_states
            bottom = np.isin(batch.strand, (OB, CTOB))[template_fragment]
            targets = np.where(bottom, batch.bases == G, batch.bases == C)
            attempts = np.flatnonzero(targets & ~methylated)
            converted[attempts] = (
                conversion_rng.random(len(attempts)) < self.conversion_rate)
        template = np.where(
            converted, np.where(batch.bases == C, T, A), batch.bases).astype(np.uint8)

        # Read layout: R1 at the left end for OT/CTOB, at the right end otherwise. A row
        # spans the longest mate's bases from its read's outer end.
        mates_per_fragment = 2 if self.paired_end else 1
        lengths = self.read_lengths[:mates_per_fragment]
        width = max(lengths)
        fragment = np.repeat(np.arange(batch.size), mates_per_fragment)
        mate = np.tile(np.arange(1, mates_per_fragment + 1), batch.size).astype(np.uint8)
        r1_reverse = np.isin(batch.strand, (OB, CTOT))
        reverse = r1_reverse[fragment] ^ (mate == 2)
        template_length = batch.template_lengths[fragment]
        template_start = np.where(reverse, template_length - np.array(lengths)[mate - 1], 0)
        row_start = np.where(reverse, template_length - width, 0)
        index = (batch.template_offsets[fragment] + row_start)[:, None] + np.arange(width)

        true_bases = _orient(template[index], reverse)
        quality = self.quality.sample(mate, width, quality_rng)
        sequence = self.errors.apply(true_bases, quality, mate, error_rng)
        reads = Reads(fragment, mate, reverse, template_start, sequence, quality, lengths)
        if self.annotate:
            molecules = _Molecules(site_index, methylated, targets, converted)
            reads = replace(reads, truth=_truth(reads, batch, self, molecules, true_bases))
        return reads


def _orient(rows: np.ndarray, reverse: np.ndarray) -> np.ndarray:
    """Template-orientation rows -> sequencing orientation."""
    oriented = rows.copy()
    oriented[reverse] = COMPLEMENT[rows[reverse, ::-1]]
    return oriented


@dataclass(frozen=True)
class _Molecules:
    """Every template base of a batch after methylation and conversion."""

    site_index: np.ndarray  # int64[S], the template index of each site
    methylated: np.ndarray  # bool[T]
    targets: np.ndarray  # bool[T], the cytosines of the informative strand
    converted: np.ndarray  # bool[T]


def _truth(reads: Reads, batch: FragmentBatch, simulator: ReadSimulator,
           molecules: _Molecules, true_bases: np.ndarray) -> ReadTruth:
    """zt and zr of every read, zf (and zx) of every fragment."""
    # Flags shared by zr and zf: haplotype, informative strand, conversion mode.
    bottom = np.isin(batch.strand, (OB, CTOB)).astype(np.uint32)
    if simulator.bisulfite:
        informative = bottom + 1
        read_mode = bottom[reads.fragment] ^ reads.reverse.astype(np.uint32)
        fragment_mode = bottom
    else:
        informative = np.zeros(batch.size, dtype=np.uint32)
        read_mode = np.full(len(reads.fragment), NO_CONVERSION, dtype=np.uint32)
        fragment_mode = np.full(batch.size, NO_CONVERSION, dtype=np.uint32)
    fragment_flags = batch.haplotype.astype(np.uint32) | informative << 2

    states, read_summary = _read_truth(reads, batch, molecules, true_bases,
                                       fragment_flags[reads.fragment] | read_mode << 4)
    fragment_summary = _fragment_truth(
        reads, batch, molecules, fragment_flags | fragment_mode << 4, read_summary)
    return ReadTruth(
        states=states, read_summary=_pack_summary(read_summary),
        fragment_summary=_pack_summary(fragment_summary),
        realization=_realizations(batch, molecules) if simulator.realization else None)


def _read_truth(reads: Reads, batch: FragmentBatch, molecules: _Molecules,
                true_bases: np.ndarray, read_flags: np.ndarray
                ) -> tuple[np.ndarray, np.ndarray]:
    """zt of every read, and its zr summary, unpacked."""
    total = len(batch.bases)
    fragment = reads.fragment
    context = np.zeros(total, dtype=np.uint8)
    context[molecules.site_index] = batch.site_contexts
    asm = np.zeros(total, dtype=bool)
    asm[molecules.site_index] = batch.site_asm
    variant = np.zeros(total, dtype=bool)
    starts = batch.template_offsets[owners(batch.event_offsets)] + batch.event_starts
    lengths = batch.event_ends - batch.event_starts
    variant[np.repeat(starts, lengths) + ranks(lengths)] = True

    # One mate (and read length) at a time, in template orientation. A false
    # methylation call is a converted cytosine that a sequencing error turns back
    # into the original base.
    methylated, targets, converted = molecules.methylated, molecules.targets, molecules.converted
    states = np.zeros(reads.sequence.shape, dtype=np.uint8)
    summary = np.zeros((len(fragment), 12), dtype=np.int64)
    base_flags = np.zeros(len(fragment), dtype=np.int64)
    for rows, length in reads.mates():
        index = (batch.template_offsets[fragment[rows]]
                 + reads.template_start[rows])[:, None] + np.arange(length)
        reverse = reads.reverse[rows]
        observed = _orient(reads.sequence[rows, :length], reverse)
        error = observed != _orient(true_bases[rows, :length], reverse)
        false_methylation = error & converted[index] & (observed == batch.bases[index])
        in_read_context = context[index]
        in_read_methylated = methylated[index] & (in_read_context > 0)
        in_read_converted = converted[index]
        states[rows, :length] = (
            in_read_context
            | in_read_methylated.astype(np.uint8) << 2
            | in_read_converted.astype(np.uint8) << 3
            | variant[index].astype(np.uint8) << 4
            | error.astype(np.uint8) << 5
        )
        base_flags[rows] = (
            np.where(asm[index].any(axis=1), _FLAG_ASM, 0)
            | np.where(in_read_converted.any(axis=1), _FLAG_CONVERTED, 0)
        )
        for column, code in ((1, CG), (3, CHG), (5, CHH)):
            summary[rows, column] = (in_read_context == code).sum(axis=1)
            summary[rows, column + 1] = (in_read_methylated & (in_read_context == code)).sum(axis=1)
        summary[rows, 7] = in_read_converted.sum(axis=1)
        summary[rows, 8] = (targets[index] & ~methylated[index] & ~in_read_converted).sum(axis=1)
        summary[rows, 10] = error.sum(axis=1)
        summary[rows, 11] = false_methylation.sum(axis=1)

    read_events, kind_flags = _events_seen(reads, batch)
    summary[:, 0] = read_flags | kind_flags | base_flags
    summary[:, 9] = read_events
    return states, summary


def _events_seen(reads: Reads, batch: FragmentBatch) -> tuple[np.ndarray, np.ndarray]:
    """The variant events each read sees (it overlaps them, or a deletion lies strictly
    inside it): their number and the flags of their kinds."""
    n_reads = len(reads.fragment)
    per_read = np.diff(batch.event_offsets)[reads.fragment]
    pair_read = np.repeat(np.arange(n_reads), per_read)
    pair_event = batch.event_offsets[reads.fragment][pair_read] + ranks(per_read)
    begin = reads.template_start[pair_read]
    end = begin + reads.length[pair_read]
    event_start = batch.event_starts[pair_event]
    event_end = batch.event_ends[pair_event]
    seen = np.where(
        event_start == event_end,
        (begin < event_start) & (event_start < end),
        np.maximum(event_start, begin) < np.minimum(event_end, end),
    )
    counts = np.bincount(pair_read[seen], minlength=n_reads).astype(np.uint32)
    return counts, _kind_flags(pair_read[seen], batch.event_kinds[pair_event[seen]], n_reads)


def _kind_flags(owner: np.ndarray, kinds: np.ndarray, size: int) -> np.ndarray:
    """The flags of the variant kinds each owner has (``owner[i]`` has an event of ``kinds[i]``)."""
    flags = np.zeros(size, dtype=np.uint32)
    for kind in (SNV, INSERTION, DELETION):
        flags[owner[kinds == kind]] |= _KIND_FLAG[kind]  # repeated owners all write the same value
    return flags


def _fragment_truth(reads: Reads, batch: FragmentBatch, molecules: _Molecules,
                    flags: np.ndarray, read_summary: np.ndarray) -> np.ndarray:
    """The zf summary of every fragment over the complete molecule, unpacked."""
    template_fragment = owners(batch.template_offsets)
    site_fragment = batch.site_fragments()
    site_methylated = molecules.methylated[molecules.site_index]
    methylated, targets, converted = molecules.methylated, molecules.targets, molecules.converted
    summary = np.zeros((batch.size, 12), dtype=np.int64)
    kind_flags = _kind_flags(owners(batch.event_offsets), batch.event_kinds, batch.size)
    converted_per_fragment = np.bincount(template_fragment[converted], minlength=batch.size)
    summary[:, 0] = (
        flags | kind_flags
        | np.where(np.bincount(site_fragment[batch.site_asm], minlength=batch.size) > 0,
                   _FLAG_ASM, 0)
        | np.where(converted_per_fragment > 0, _FLAG_CONVERTED, 0)
    )
    for column, code in ((1, CG), (3, CHG), (5, CHH)):
        selected = batch.site_contexts == code
        summary[:, column] = np.bincount(site_fragment[selected], minlength=batch.size)
        summary[:, column + 1] = np.bincount(
            site_fragment[selected & site_methylated], minlength=batch.size)
    summary[:, 7] = converted_per_fragment
    summary[:, 8] = np.bincount(
        template_fragment[targets & ~methylated & ~converted], minlength=batch.size)
    summary[:, 9] = np.diff(batch.event_offsets)
    # errors and false methylation calls of its reads, which are adjacent
    summary[:, 10:12] = read_summary[:, 10:12].reshape(batch.size, -1, 2).sum(axis=1)
    return summary


def _realizations(batch: FragmentBatch, molecules: _Molecules) -> list[str]:
    """zx of every fragment: its site states, and the conversions of its target cytosines."""
    site_methylated = molecules.methylated[molecules.site_index]
    realizations = []
    for i in range(batch.size):
        bases = slice(batch.template_offsets[i], batch.template_offsets[i + 1])
        realizations.append(_realization(
            site_methylated[batch.site_offsets[i]:batch.site_offsets[i + 1]],
            molecules.converted[bases][molecules.targets[bases]]))
    return realizations


def _pack_summary(values: np.ndarray) -> np.ndarray:
    overflow = (values[:, 1:] > 0xFFFF).any(axis=1)
    values[:, 0] |= np.where(overflow, _FLAG_OVERFLOW, 0)
    values[:, 1:] = np.minimum(values[:, 1:], 0xFFFF)
    return values.astype(np.uint16)


def _realization(states: np.ndarray, conversions: np.ndarray) -> str:
    def bits(values):
        packed = np.packbits(values.astype(np.uint8), bitorder="little").tobytes()
        return base64.urlsafe_b64encode(packed).rstrip(b"=").decode("ascii")

    return f"{len(states):x}.{len(conversions):x}.{bits(states)}.{bits(conversions)}"
