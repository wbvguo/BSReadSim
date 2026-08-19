"""Project decoded protocol columns into scientific processing values."""

from __future__ import annotations

from typing import List, Tuple

import numpy as np

try:
    from ._native import (
        decode_protocol_fragments as _native_decode_fragments,
        pack_protocol_common_columns as _native_pack_common_columns,
    )
except (ImportError, AttributeError):
    _native_decode_fragments = None
    _native_pack_common_columns = None

from .numpy_postprocess import NumpyFragmentBatch
from .postprocess import MethylationModelBatch, PostprocessError
from .model import (
    NO_VARIANT_EVENT,
    CaptureStrand,
    Fragment,
    Mate,
    MethylationAllele,
    MethylationContext,
    MethylationSite,
    MethylationSource,
    SiteReference,
    VariantEvent,
    VariantKind,
)
from .protocol import (
    NO_REFERENCE_POSITION,
    DecodedBatchView,
    Header,
    TruthMode,
)


_U1 = np.dtype("u1")
_U4 = np.dtype("<u4")
_U8 = np.dtype("<u8")
_I8 = np.dtype("<i8")
_F4 = np.dtype("<f4")


def _decode_fragments_python(
    batch: DecodedBatchView,
    header: Header,
) -> Tuple[Fragment, ...]:
    """Reconstruct processing values from validated protocol columns."""

    if not isinstance(batch, DecodedBatchView):
        raise PostprocessError("batch must be a decoded protocol batch")
    if not isinstance(header, Header):
        raise PostprocessError("header must be a protocol Header")
    if (header.truth_columns is TruthMode.FULL) != (batch.truth is not None):
        raise PostprocessError("batch truth columns disagree with the header")

    truth = batch.truth
    fragments = []
    for row in range(batch.fragment_count):
        template_flat_begin = batch.template_offsets[row]
        template_flat_end = batch.template_offsets[row + 1]
        template_length = template_flat_end - template_flat_begin
        template_bases = bytes(
            batch.template_bases[template_flat_begin:template_flat_end]
        )

        base_event_ids = [NO_VARIANT_EVENT] * template_length
        variant_events = []
        event_by_id = {}
        reference_positions = [-1] * template_length
        if truth is not None:
            projection_begin = truth.projection_offsets[row]
            projection_end = truth.projection_offsets[row + 1]
            for index in range(projection_begin, projection_end):
                template_begin = truth.projection_template_begins[index]
                template_end = truth.projection_template_ends[index]
                reference_begin = truth.projection_reference_begins[index]
                for local_offset in range(template_begin, template_end):
                    reference_positions[local_offset] = reference_begin + (
                        local_offset - template_begin
                    )

            event_begin = truth.event_offsets[row]
            event_end = truth.event_offsets[row + 1]
            for index in range(event_begin, event_end):
                ref_begin = truth.event_ref_offsets[index]
                ref_end = truth.event_ref_offsets[index + 1]
                alt_begin = truth.event_alt_offsets[index]
                alt_end = truth.event_alt_offsets[index + 1]
                event = VariantEvent(
                    event_id=truth.event_ids[index],
                    kind=VariantKind(truth.event_kinds[index]),
                    phased_haplotype=truth.event_phased_haplotypes[index],
                    reference_start=truth.event_reference_begins[index],
                    reference_end=truth.event_reference_ends[index],
                    ref_bases=bytes(truth.event_ref_bases[ref_begin:ref_end]),
                    alt_bases=bytes(truth.event_alt_bases[alt_begin:alt_end]),
                )
                variant_events.append(event)
                event_by_id[event.event_id] = event
                for local_offset in range(
                    truth.event_template_begins[index],
                    truth.event_template_ends[index],
                ):
                    base_event_ids[local_offset] = event.event_id

        site_begin = batch.site_offsets[row]
        site_end = batch.site_offsets[row + 1]
        methylation_sites = tuple(
            MethylationSite(
                site_index=local_index,
                template_offset=batch.site_template_offsets[index],
                reference_pos=(
                    -1
                    if truth is None
                    else _reference_position(truth.site_reference_positions[index])
                ),
                context=MethylationContext(batch.site_contexts[index]),
                source=MethylationSource(batch.site_sources[index]),
                allele=MethylationAllele(batch.site_alleles[index]),
                methylation_probability=batch.site_probabilities[index],
            )
            for local_index, index in enumerate(range(site_begin, site_end))
        )

        mate_begin = batch.mate_offsets[row]
        mate_end = batch.mate_offsets[row + 1]
        mates = []
        for index in range(mate_begin, mate_end):
            template_begin = batch.mate_template_begins[index]
            template_end = batch.mate_template_ends[index]
            reverse = bool(batch.mate_reverse_complements[index])
            mapped = [] if truth is None else [
                position
                for position in reference_positions[template_begin:template_end]
                if position != -1
            ]
            if truth is None:
                reference_start = batch.reference_begins[row]
                reference_end = batch.reference_ends[row]
            elif mapped:
                reference_start = min(mapped)
                reference_end = max(mapped) + 1
            else:
                anchors = {
                    event_by_id[event_id].reference_start
                    for event_id in base_event_ids[template_begin:template_end]
                    if event_id != NO_VARIANT_EVENT
                }
                if len(anchors) != 1:
                    raise PostprocessError(
                        "insertion-only v2 mate has ambiguous reference anchors"
                    )
                reference_start = anchors.pop()
                reference_end = reference_start
            site_refs = []
            for site in methylation_sites:
                if not template_begin <= site.template_offset < template_end:
                    continue
                read_offset = (
                    template_end - 1 - site.template_offset
                    if reverse
                    else site.template_offset - template_begin
                )
                site_refs.append(SiteReference(read_offset, site.site_index))
            site_refs.sort(key=lambda item: item.read_offset)
            mates.append(
                Mate(
                    mate_index=batch.mate_indices[index],
                    reverse_complement=reverse,
                    template_start=template_begin,
                    template_end=template_end,
                    reference_start=reference_start,
                    reference_end=reference_end,
                    site_refs=tuple(site_refs),
                )
            )

        fragments.append(
            Fragment(
                fragment_ordinal=batch.first_fragment_ordinal + row,
                contig_index=batch.contig_indices[row],
                haplotype=batch.haplotypes[row],
                capture_strand=CaptureStrand(batch.capture_strands[row]),
                reference_start=batch.reference_begins[row],
                reference_end=batch.reference_ends[row],
                template_bases=template_bases,
                reference_positions=tuple(reference_positions),
                base_event_ids=tuple(base_event_ids),
                variant_events=tuple(variant_events),
                methylation_sites=methylation_sites,
                mates=tuple(mates),
            )
        )
    return tuple(fragments)


def decode_fragments(
    batch: DecodedBatchView,
    header: Header,
) -> Tuple[Fragment, ...]:
    """Reconstruct processing values, preserving truth when it is present."""

    if (
        _native_decode_fragments is None
        or header.truth_columns is TruthMode.NONE
    ):
        return _decode_fragments_python(batch, header)
    if not isinstance(batch, DecodedBatchView):
        raise PostprocessError("batch must be a decoded protocol batch")
    if not isinstance(header, Header):
        raise PostprocessError("header must be a protocol Header")
    if header.truth_columns is not TruthMode.FULL or batch.truth is None:
        raise PostprocessError("typed protocol decoding requires Full Truth")
    truth = batch.truth
    common_columns = (
        batch.contig_indices.raw,
        batch.reference_begins.raw,
        batch.reference_ends.raw,
        batch.template_offsets.raw,
        batch.mate_offsets.raw,
        batch.site_offsets.raw,
        batch.mate_template_begins.raw,
        batch.mate_template_ends.raw,
        batch.site_template_offsets.raw,
        batch.site_probabilities.raw,
        batch.haplotypes.raw,
        batch.capture_strands.raw,
        batch.mate_indices.raw,
        batch.mate_reverse_complements.raw,
        batch.site_contexts.raw,
        batch.site_sources.raw,
        batch.site_alleles.raw,
        batch.template_bases.raw,
    )
    truth_columns = (
        truth.projection_offsets.raw,
        truth.event_offsets.raw,
        truth.original_n_offsets.raw,
        truth.projection_template_begins.raw,
        truth.projection_template_ends.raw,
        truth.projection_reference_begins.raw,
        truth.event_ids.raw,
        truth.event_reference_begins.raw,
        truth.event_reference_ends.raw,
        truth.event_template_begins.raw,
        truth.event_template_ends.raw,
        truth.event_ref_offsets.raw,
        truth.event_alt_offsets.raw,
        truth.site_reference_positions.raw,
        truth.original_n_template_offsets.raw,
        truth.event_kinds.raw,
        truth.event_phased_haplotypes.raw,
        truth.event_ref_bases.raw,
        truth.event_alt_bases.raw,
    )
    try:
        return _native_decode_fragments(
            common_columns,
            truth_columns,
            batch.first_fragment_ordinal,
            tuple(contig.length for contig in header.contigs),
            header.mates_per_fragment,
            header.read_length_r1,
            header.read_length_r2,
        )
    except (BufferError, OverflowError, TypeError, ValueError) as error:
        raise PostprocessError(str(error)) from error


def _decode_common_numpy_batch_python(
    batch: DecodedBatchView,
    header: Header,
) -> NumpyFragmentBatch:
    """Build exact built-in FASTQ inputs directly from v2 common columns."""

    if not isinstance(batch, DecodedBatchView):
        raise PostprocessError("batch must be a decoded protocol batch")
    if not isinstance(header, Header):
        raise PostprocessError("header must be a protocol Header")
    if header.truth_columns is not TruthMode.NONE or batch.truth is not None:
        raise PostprocessError("common-column processing requires TruthMode.NONE")

    ordinals = tuple(
        batch.first_fragment_ordinal + row for row in range(batch.fragment_count)
    )
    contig_names = tuple(
        header.contigs[index].name for index in batch.contig_indices
    )
    template_offsets = _as_bytes(batch.template_offsets, _U8)
    site_offsets = _as_bytes(batch.site_offsets, _U8)
    local_site_indices = np.concatenate(
        tuple(
            np.arange(
                batch.site_offsets[row + 1] - batch.site_offsets[row],
                dtype=_U4,
            )
            for row in range(batch.fragment_count)
        )
    ) if batch.methylation_site_count else np.empty(0, dtype=_U4)
    model = MethylationModelBatch(
        fragment_ordinals=ordinals,
        contig_names=contig_names,
        fragment_ordinal_data=np.asarray(ordinals, dtype=_U8).tobytes(),
        contig_indices=_as_bytes(batch.contig_indices, _U4),
        template_offsets=template_offsets,
        template_bases=_as_bytes(batch.template_bases, _U1),
        site_offsets=site_offsets,
        site_indices=local_site_indices.tobytes(),
        site_template_offsets=_as_bytes(batch.site_template_offsets, _U4),
        site_reference_positions=np.full(
            batch.methylation_site_count, -1, dtype=_I8
        ).tobytes(),
        site_contexts=_as_bytes(batch.site_contexts, _U1),
        site_sources=_as_bytes(batch.site_sources, _U1),
        site_alleles=_as_bytes(batch.site_alleles, _U1),
        site_probabilities=_as_bytes(batch.site_probabilities, _F4),
    )

    site_ref_offsets = [0]
    site_ref_read_offsets = []  # type: List[int]
    site_ref_site_indices = []  # type: List[int]
    mate_reference_starts = []
    mate_reference_ends = []
    for row in range(batch.fragment_count):
        fragment_site_begin = batch.site_offsets[row]
        fragment_site_end = batch.site_offsets[row + 1]
        for mate_index in range(batch.mate_offsets[row], batch.mate_offsets[row + 1]):
            template_begin = batch.mate_template_begins[mate_index]
            template_end = batch.mate_template_ends[mate_index]
            reverse = bool(batch.mate_reverse_complements[mate_index])
            references = []
            for site_index in range(fragment_site_begin, fragment_site_end):
                template_offset = batch.site_template_offsets[site_index]
                if not template_begin <= template_offset < template_end:
                    continue
                read_offset = (
                    template_end - 1 - template_offset
                    if reverse
                    else template_offset - template_begin
                )
                references.append((read_offset, site_index - fragment_site_begin))
            references.sort()
            site_ref_read_offsets.extend(item[0] for item in references)
            site_ref_site_indices.extend(item[1] for item in references)
            site_ref_offsets.append(len(site_ref_read_offsets))
            # NumpyFragmentBatch requires these common coordinate buffers; the
            # FASTQ-only processor never treats them as Truth provenance.
            mate_reference_starts.append(batch.reference_begins[row])
            mate_reference_ends.append(batch.reference_ends[row])

    return NumpyFragmentBatch(
        model=model,
        haplotypes=_as_bytes(batch.haplotypes, _U1),
        capture_strands=_as_bytes(batch.capture_strands, _U1),
        mate_offsets=_as_bytes(batch.mate_offsets, _U8),
        mate_indices=_as_bytes(batch.mate_indices, _U1),
        mate_reverse_complements=_as_bytes(batch.mate_reverse_complements, _U1),
        mate_template_starts=_as_bytes(batch.mate_template_begins, _U4),
        mate_template_ends=_as_bytes(batch.mate_template_ends, _U4),
        mate_reference_starts=np.asarray(
            mate_reference_starts, dtype=_U8
        ).tobytes(),
        mate_reference_ends=np.asarray(mate_reference_ends, dtype=_U8).tobytes(),
        site_ref_offsets=np.asarray(site_ref_offsets, dtype=_U8).tobytes(),
        site_ref_read_offsets=np.asarray(
            site_ref_read_offsets, dtype=_U4
        ).tobytes(),
        site_ref_site_indices=np.asarray(
            site_ref_site_indices, dtype=_U4
        ).tobytes(),
    )


def _decode_common_numpy_batch(
    batch: DecodedBatchView,
    header: Header,
) -> NumpyFragmentBatch:
    """Build exact common-column inputs with a Python reference implementation."""

    if _native_pack_common_columns is None:
        return _decode_common_numpy_batch_python(batch, header)
    if not isinstance(batch, DecodedBatchView):
        raise PostprocessError("batch must be a decoded protocol batch")
    if not isinstance(header, Header):
        raise PostprocessError("header must be a protocol Header")
    if header.truth_columns is not TruthMode.NONE or batch.truth is not None:
        raise PostprocessError("common-column processing requires TruthMode.NONE")

    common_columns = (
        batch.contig_indices.raw,
        batch.reference_begins.raw,
        batch.reference_ends.raw,
        batch.template_offsets.raw,
        batch.mate_offsets.raw,
        batch.site_offsets.raw,
        batch.mate_template_begins.raw,
        batch.mate_template_ends.raw,
        batch.site_template_offsets.raw,
        batch.site_probabilities.raw,
        batch.haplotypes.raw,
        batch.capture_strands.raw,
        batch.mate_indices.raw,
        batch.mate_reverse_complements.raw,
        batch.site_contexts.raw,
        batch.site_sources.raw,
        batch.site_alleles.raw,
        batch.template_bases.raw,
    )
    try:
        ordinals, contig_names, model_columns, fragment_columns = (
            _native_pack_common_columns(
                common_columns,
                batch.first_fragment_ordinal,
                tuple(contig.name for contig in header.contigs),
                tuple(contig.length for contig in header.contigs),
                header.mates_per_fragment,
                header.read_length_r1,
                header.read_length_r2,
            )
        )
    except (BufferError, OverflowError, TypeError, ValueError) as error:
        raise PostprocessError(str(error)) from error
    model = MethylationModelBatch(ordinals, contig_names, *model_columns)
    return NumpyFragmentBatch(model, *fragment_columns)


def _reference_position(value: int) -> int:
    return -1 if value == NO_REFERENCE_POSITION else value


def _as_bytes(values, dtype: np.dtype) -> bytes:
    raw = getattr(values, "raw", None)
    if (
        raw is not None
        and raw.format == "B"
        and raw.nbytes == len(values) * dtype.itemsize
    ):
        if dtype in (_U1, _U4, _F4):
            return raw.tobytes()
    return np.asarray(tuple(values), dtype=dtype).tobytes()


__all__ = [
    "decode_fragments",
]
