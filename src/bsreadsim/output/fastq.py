"""FASTQ encoding for processed fragments."""

from __future__ import annotations


from ..process.batch import ProcessedFragment, ProcessedMate
from ..process.batch import ReadNameError, format_fragment_identifier
from .errors import OutputError


_VALID_SEQUENCE = frozenset("ACGTN")


def _validate_processed_fragment(
    fragment: ProcessedFragment,
    paired_end: bool,
    *,
    require_base_states: bool = True,
) -> None:
    if not isinstance(fragment, ProcessedFragment):
        raise OutputError("fragment must be a ProcessedFragment")
    if not isinstance(paired_end, bool):
        raise OutputError("paired_end must be a boolean")
    if not isinstance(require_base_states, bool):
        raise OutputError("require_base_states must be a boolean")
    try:
        format_fragment_identifier(
            fragment.contig_name,
            fragment.reference_start,
            fragment.reference_end,
            fragment.fragment_ordinal,
        )
    except ReadNameError as error:
        raise OutputError(str(error)) from error
    expected_indices = (0, 1) if paired_end else (0,)
    observed_indices = tuple(sorted(mate.mate_index for mate in fragment.mates))
    if observed_indices != expected_indices:
        raise OutputError("processed mate cardinality disagrees with output mode")
    for mate in fragment.mates:
        _validate_mate(mate, require_base_states=require_base_states)


def _validate_mate(
    mate: ProcessedMate,
    *,
    require_base_states: bool = True,
) -> None:
    if not isinstance(mate, ProcessedMate):
        raise OutputError("fragment mates must be ProcessedMate values")
    if len(mate.sequence) != len(mate.quality):
        raise OutputError("FASTQ sequence and quality lengths differ")
    if not mate.sequence or any(base not in _VALID_SEQUENCE for base in mate.sequence):
        raise OutputError("FASTQ sequence must be non-empty uppercase A/C/G/T/N")
    if any(not 33 <= ord(character) <= 126 for character in mate.quality):
        raise OutputError("FASTQ quality contains a non-Phred+33 character")
    if require_base_states and len(mate.base_states) != len(mate.sequence):
        raise OutputError("mate base states must cover every read base")


def format_fragment_records_trusted(
    fragment: ProcessedFragment,
    *,
    paired_end: bool,
) -> tuple[bytes, bytes | None]:
    """Format a locally produced, already validated process-worker value."""

    records = {}  # type: dict[str, bytes]
    for mate in sorted(fragment.mates, key=lambda item: item.mate_index):
        role = "read1" if mate.mate_index == 0 else "read2"
        records[role] = _fastq_record(fragment, mate).encode("utf-8")
    return records["read1"], records.get("read2")


def _fastq_record(
    fragment: ProcessedFragment,
    mate: ProcessedMate,
) -> str:
    return _fastq_record_fields(
        fragment.contig_name,
        fragment.reference_start,
        fragment.reference_end,
        fragment.fragment_ordinal,
        mate.mate_index,
        mate.sequence,
        mate.quality,
    )


def _fastq_record_fields(
    contig_name: str,
    reference_start: int,
    reference_end: int,
    fragment_ordinal: int,
    mate_index: int,
    sequence: str,
    quality: str,
) -> str:
    """Format trusted FASTQ fields shared by typed and common-column lanes."""

    identifier = format_fragment_identifier(
        contig_name,
        reference_start,
        reference_end,
        fragment_ordinal,
        pair_number=mate_index + 1,
    )
    return "@{}\n{}\n+\n{}\n".format(identifier, sequence, quality)


__all__ = [
    "format_fragment_records_trusted",
]
