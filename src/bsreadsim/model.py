"""In-memory scientific values shared by simulation post-processing stages."""

from dataclasses import dataclass
from enum import IntEnum
from typing import Tuple


NO_VARIANT_EVENT = 0xFFFFFFFF


class CaptureStrand(IntEnum):
    UNKNOWN = 0
    FORWARD = 1
    REVERSE = 2


class VariantKind(IntEnum):
    SNV = 1
    INSERTION = 2
    DELETION = 3


class MethylationContext(IntEnum):
    CG_C = 1
    CHG_C = 3
    CHH_C = 7
    CG_G = 9
    CHG_G = 11
    CHH_G = 15


class MethylationSource(IntEnum):
    CGMAP = 1
    ASM = 2
    BETA = 3
    POOLED_CGMAP = 4


class MethylationAllele(IntEnum):
    SHARED = 0
    REFERENCE_HAPLOTYPE = 1
    ALTERNATE_HAPLOTYPE = 2


@dataclass(frozen=True)
class VariantEvent:
    event_id: int
    kind: VariantKind
    phased_haplotype: int
    reference_start: int
    reference_end: int
    ref_bases: bytes
    alt_bases: bytes


@dataclass(frozen=True)
class MethylationSite:
    site_index: int
    template_offset: int
    reference_pos: int
    context: MethylationContext
    source: MethylationSource
    allele: MethylationAllele
    methylation_probability: float


@dataclass(frozen=True)
class SiteReference:
    read_offset: int
    site_index: int


@dataclass(frozen=True)
class Mate:
    mate_index: int
    reverse_complement: bool
    template_start: int
    template_end: int
    reference_start: int
    reference_end: int
    site_refs: Tuple[SiteReference, ...]


@dataclass(frozen=True)
class Fragment:
    fragment_ordinal: int
    contig_index: int
    haplotype: int
    capture_strand: CaptureStrand
    reference_start: int
    reference_end: int
    template_bases: bytes
    reference_positions: Tuple[int, ...]
    base_event_ids: Tuple[int, ...]
    variant_events: Tuple[VariantEvent, ...]
    methylation_sites: Tuple[MethylationSite, ...]
    mates: Tuple[Mate, ...]


@dataclass(frozen=True)
class FragmentSummary:
    """Fragment aggregates observed independently by Python workers."""

    fragment_count: int
    mate_count: int
    template_base_count: int
    methylation_site_count: int
    per_contig_fragment_counts: Tuple[int, ...]


__all__ = [
    "NO_VARIANT_EVENT",
    "CaptureStrand",
    "Fragment",
    "FragmentSummary",
    "Mate",
    "MethylationAllele",
    "MethylationContext",
    "MethylationSite",
    "MethylationSource",
    "SiteReference",
    "VariantEvent",
    "VariantKind",
]
