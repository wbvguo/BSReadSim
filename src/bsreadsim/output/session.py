"""Transactional FASTQ, annotated BAM, and manifest output.

The output component accepts only format-ready ``ProcessedFragment`` values.
It stages every byte on the destination filesystem, enforces fragment order
and SE/PE cardinality.  Finalization computes data-file identities while they
remain private; an explicit commit publishes the manifest last as the marker
that makes the whole file set valid.
"""

from __future__ import annotations
from contextlib import suppress

from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import tempfile
from typing import BinaryIO

from ..process.batch import ProcessedFragment
from .bam import BamConfig, BamOutput, format_sam_fragment
from .errors import OutputError
from .fastq import (
    format_fragment_records_trusted,
    _validate_processed_fragment,
)


_PREFIX = re.compile(r"^[A-Za-z0-9._-]+$")
BytesLike = bytes | bytearray | memoryview


@dataclass(frozen=True)
class OutputConfig:
    """Validated output policy for one run."""

    directory: Path
    prefix: str
    paired_end: bool
    compression: str = "gzip"
    gzip_level: int = 6
    bam: BamConfig | None = None

    def __post_init__(self) -> None:
        directory = Path(self.directory)
        if not directory.is_absolute():
            raise OutputError("output directory must be an absolute path")
        if (
            not isinstance(self.prefix, str)
            or not 1 <= len(self.prefix) <= 128
            or _PREFIX.fullmatch(self.prefix) is None
        ):
            raise OutputError("output prefix contains unsupported characters")
        if not isinstance(self.paired_end, bool):
            raise OutputError("paired_end must be a boolean")
        if self.compression not in ("none", "gzip"):
            raise OutputError("compression must be 'none' or 'gzip'")
        if (
            isinstance(self.gzip_level, bool)
            or not isinstance(self.gzip_level, int)
            or not 0 <= self.gzip_level <= 9
        ):
            raise OutputError("gzip_level must be an integer in [0, 9]")
        if self.bam is not None and not isinstance(
            self.bam, BamConfig
        ):
            raise OutputError("bam must be a BamConfig or None")
        object.__setattr__(self, "directory", directory)


@dataclass(frozen=True)
class OutputFileSummary:
    """Digest and record count of one finalized data artifact."""

    role: str
    path: Path
    size_bytes: int
    record_count: int
    sha256: str


@dataclass(frozen=True)
class OutputSummary:
    """Accounting returned while finalized data files remain staged."""

    fragment_count: int
    mate_count: int
    files: tuple[OutputFileSummary, ...]

    def file_for_role(self, role: str) -> OutputFileSummary:
        for file_summary in self.files:
            if file_summary.role == role:
                return file_summary
        raise KeyError(role)


class _DigestingSink:
    """Hash and count the exact staged bytes as they reach the raw file."""

    def __init__(self, raw: BinaryIO) -> None:
        self.raw = raw
        self.size_bytes = 0
        self.digest = hashlib.sha256()

    def write(self, value: bytes) -> int:
        view = memoryview(value)
        total = 0
        while view:
            written = self.raw.write(view)
            if written is None or written <= 0:
                raise OSError("staged output write made no progress")
            self.digest.update(view[:written])
            self.size_bytes += written
            total += written
            view = view[written:]
        return total


class _BinaryOutput:
    def __init__(self, path: Path, compression: str, gzip_level: int) -> None:
        self.path = path
        self.raw = path.open("xb")
        self.sink = _DigestingSink(self.raw)
        self.compressed = compression == "gzip"
        if compression == "gzip":
            self.stream = gzip.GzipFile(
                filename="",
                mode="wb",
                compresslevel=gzip_level,
                fileobj=self.sink,
                mtime=0,
            )
        else:
            self.stream = self.sink
        self.closed = False

    def write_text(self, value: str) -> None:
        self.write_bytes(value.encode("utf-8"))

    def write_bytes(self, value: BytesLike) -> None:
        if self.closed:
            raise OutputError("cannot write to a closed staged output")
        if not isinstance(value, (bytes, bytearray, memoryview)):
            raise OutputError("staged output accepts bytes-like values only")
        self.stream.write(value)

    def close(self) -> None:
        if self.closed:
            return
        try:
            if self.compressed:
                self.stream.close()
            self.raw.flush()
            os.fsync(self.raw.fileno())
        finally:
            self.raw.close()
            self.closed = True

    @property
    def size_bytes(self) -> int:
        if not self.closed:
            raise OutputError("staged output identity requires a closed stream")
        return self.sink.size_bytes

    @property
    def sha256(self) -> str:
        if not self.closed:
            raise OutputError("staged output identity requires a closed stream")
        return self.sink.digest.hexdigest()


class OutputSession:
    """Stage, validate, and transactionally expose one set of output files.

    Publication uses hard links with no-replace semantics.  If any link fails,
    already-linked members are removed.  Data links are created first and the
    complete manifest is linked last as the externally visible commit marker.
    """

    def __init__(self, config: OutputConfig) -> None:
        if not isinstance(config, OutputConfig):
            raise OutputError("config must be OutputConfig")
        self.config = config
        self._state = "open"
        self._expected_ordinal = 0
        self._mate_count = 0
        self._published = []  # type: list[Path]
        try:
            config.directory.mkdir(parents=True, exist_ok=True)
            if not config.directory.is_dir():
                raise OutputError("output directory is not a directory")
            self._staging_directory = Path(
                tempfile.mkdtemp(
                    prefix=".{}.staging-".format(config.prefix),
                    dir=str(config.directory),
                )
            )
            self._paths = self._build_paths()
            self._manifest_paths = (
                self._staging_directory
                / "{}.manifest.json".format(config.prefix),
                self.config.directory
                / "{}.manifest.json".format(config.prefix),
            )
            self._require_destinations_absent()
            self._streams = {}  # type: dict[str, _BinaryOutput | BamOutput]
            for role, (staged_path, final_path) in self._paths.items():
                if role == "bam":
                    if config.bam is None:
                        raise OutputError("BAM path has no writer policy")
                    self._streams[role] = BamOutput(
                        staged_path, config.bam
                    )
                else:
                    self._streams[role] = _BinaryOutput(
                        staged_path, config.compression, config.gzip_level
                    )
        except Exception:
            for stream in getattr(self, "_streams", {}).values():
                try:
                    stream.close()
                except Exception:
                    pass
            staging = getattr(self, "_staging_directory", None)
            if staging is not None:
                shutil.rmtree(staging, ignore_errors=True)
            self._state = "aborted"
            raise

    @property
    def staging_directory(self) -> Path:
        return self._staging_directory

    @property
    def committed(self) -> bool:
        return self._state == "committed"

    def __enter__(self) -> "OutputSession":
        self._require_state("open")
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if self._state != "committed":
            self.abort()
        return False

    def write_fragment(self, fragment: ProcessedFragment) -> None:
        """Write one complete SE/PE fragment to every staged artifact."""
        self._require_state("open")
        try:
            self._validate_fragment(fragment)
            if "read1" in self._streams:
                read1, read2 = format_fragment_records_trusted(
                    fragment,
                    paired_end=self.config.paired_end,
                )
                self._streams["read1"].write_bytes(read1)
                if read2 is not None and "read2" in self._streams:
                    self._streams["read2"].write_bytes(read2)
            if self.config.bam is not None:
                bam_records = format_sam_fragment(
                    fragment,
                    paired_end=self.config.paired_end,
                    read_group_id=self.config.bam.read_group_id,
                    contig_length=self.config.bam.reference_length(
                        fragment.contig_name
                    ),
                    include_fragment_summary=(
                        self.config.bam.fragment_summary
                    ),
                )
                self._streams["bam"].write_bytes(b"".join(bam_records))
            self._expected_ordinal += 1
            self._mate_count += len(fragment.mates)
        except Exception:
            self.abort()
            raise

    def write_formatted_batch(
        self,
        first_ordinal: int,
        record_lengths: tuple[tuple[int, int], ...],
        read1: BytesLike | None,
        read2: BytesLike | None,
        *,
        alignment_record_lengths: tuple[int, ...] = (),
        alignment: BytesLike | None = None,
    ) -> None:
        """Write one already validated, consecutive formatted fragment batch."""

        self._require_state("open")
        views = []
        try:
            if (
                isinstance(first_ordinal, bool)
                or not isinstance(first_ordinal, int)
                or first_ordinal < 0
            ):
                raise OutputError("formatted batch first_ordinal must be non-negative")
            if first_ordinal != self._expected_ordinal:
                raise OutputError("formatted batches must be consecutive")
            if not isinstance(record_lengths, tuple) or not record_lengths:
                raise OutputError("formatted batch record lengths must be non-empty")

            fastq_enabled = "read1" in self._streams
            read1_view = None
            if read1 is not None:
                read1_view = _byte_view("formatted read1", read1)
                views.append(read1_view)
            read2_view = None
            if read2 is not None:
                read2_view = _byte_view("formatted read2", read2)
                views.append(read2_view)
            if fastq_enabled:
                if read1_view is None:
                    raise OutputError("formatted FASTQ batch omitted read1")
                if self.config.paired_end != (read2_view is not None):
                    raise OutputError(
                        "formatted batch mate cardinality disagrees with output"
                    )
            elif read1_view is not None or read2_view is not None:
                raise OutputError("BAM-only batch contains unpublished FASTQ bytes")

            alignment_enabled = self.config.bam is not None
            alignment_view = None
            if alignment is not None:
                alignment_view = _byte_view("formatted details alignment", alignment)
                views.append(alignment_view)
            if alignment_enabled != (alignment_view is not None):
                raise OutputError(
                    "formatted details alignment presence disagrees with output policy"
                )
            if not isinstance(alignment_record_lengths, tuple):
                raise OutputError(
                    "formatted details alignment lengths must be an immutable tuple"
                )
            expected_alignment_records = len(record_lengths) * (
                2 if self.config.paired_end else 1
            )
            if alignment_enabled:
                if len(alignment_record_lengths) != expected_alignment_records:
                    raise OutputError(
                        "formatted details alignment count disagrees with mate count"
                    )
                if any(
                    isinstance(value, bool)
                    or not isinstance(value, int)
                    or value <= 0
                    for value in alignment_record_lengths
                ):
                    raise OutputError(
                        "formatted details alignment lengths must be positive"
                    )
            elif alignment_record_lengths:
                raise OutputError(
                    "formatted details alignment lengths disagree with output policy"
                )

            read1_total = 0
            read2_total = 0
            for lengths in record_lengths:
                if (
                    not isinstance(lengths, tuple)
                    or len(lengths) != 2
                    or any(
                        isinstance(value, bool)
                        or not isinstance(value, int)
                        or value < 0
                        for value in lengths
                    )
                ):
                    raise OutputError("formatted batch contains invalid record lengths")
                read1_length, read2_length = lengths
                if fastq_enabled:
                    if read1_length == 0:
                        raise OutputError("formatted read1 records must be non-empty")
                    if self.config.paired_end != (read2_length > 0):
                        raise OutputError(
                            "formatted record mate cardinality disagrees with output"
                        )
                elif lengths != (0, 0):
                    raise OutputError("BAM-only batch has FASTQ record lengths")
                read1_total += read1_length
                read2_total += read2_length
            if read1_view is not None and read1_total != len(read1_view):
                raise OutputError("formatted batch region sizes disagree with records")
            if read2_view is not None and read2_total != len(read2_view):
                raise OutputError("formatted read2 region size disagrees with records")
            if alignment_view is not None and sum(
                alignment_record_lengths
            ) != len(alignment_view):
                raise OutputError(
                    "formatted details alignment region size disagrees with records"
                )

            if read1_view is not None and "read1" in self._streams:
                self._streams["read1"].write_bytes(read1_view)
            if read2_view is not None and "read2" in self._streams:
                self._streams["read2"].write_bytes(read2_view)
            if alignment_view is not None:
                self._streams["bam"].write_bytes(alignment_view)

            fragment_count = len(record_lengths)
            self._expected_ordinal += fragment_count
            self._mate_count += fragment_count * (2 if self.config.paired_end else 1)
        except Exception:
            self.abort()
            raise
        finally:
            for view in views:
                view.release()

    def finalize(self) -> OutputSummary:
        """Close and identify every data file without publishing any of them."""
        self._require_state("open")
        try:
            self._close_streams()
            self._require_destinations_absent()
            self._summary = OutputSummary(
                fragment_count=self._expected_ordinal,
                mate_count=self._mate_count,
                files=self._summaries(),
            )
            self._state = "finalized"
            return self._summary
        except Exception as error:
            self.abort()
            if isinstance(error, OutputError):
                raise
            raise OutputError("cannot finalize staged outputs: {}".format(error)) from error

    def commit(self, manifest_json: str) -> OutputSummary:
        """Publish finalized data and link the complete manifest last."""
        self._require_state("finalized")
        try:
            manifest_bytes = _validate_manifest_json(
                manifest_json, self._summary
            )
            staged_manifest, final_manifest = self._manifest_paths
            with staged_manifest.open("xb") as manifest_file:
                manifest_file.write(manifest_bytes)
                manifest_file.flush()
                os.fsync(manifest_file.fileno())
            self._require_destinations_absent()
            for role in self._ordered_roles():
                staged_path, final_path = self._paths[role]
                os.link(staged_path, final_path)
                self._published.append(final_path)
            os.link(staged_manifest, final_manifest)
            self._published.append(final_manifest)
            for staged_path, _ in self._paths.values():
                staged_path.unlink()
            staged_manifest.unlink()
            self._staging_directory.rmdir()
            self._state = "committed"
            return self._summary
        except Exception as error:
            self.abort()
            if isinstance(error, OutputError):
                raise
            raise OutputError("cannot publish staged outputs: {}".format(error)) from error

    def abort(self) -> None:
        """Best-effort cleanup of this transaction's exact private paths."""
        if self._state in ("aborted", "committed"):
            return
        try:
            self._close_streams()
        except Exception:
            pass
        for final_path in reversed(self._published):
            with suppress(FileNotFoundError):
                final_path.unlink()
        self._published.clear()
        shutil.rmtree(self._staging_directory, ignore_errors=True)
        self._state = "aborted"

    def _build_paths(self) -> dict[str, tuple[Path, Path]]:
        suffix = ".gz" if self.config.compression == "gzip" else ""
        names = {}
        if self.config.bam is None:
            names["read1"] = "{}.R1.fastq{}".format(
                self.config.prefix, suffix
            )
        if self.config.bam is not None:
            names["bam"] = "{}.bam".format(self.config.prefix)
        if self.config.paired_end and self.config.bam is None:
            names["read2"] = "{}.R2.fastq{}".format(self.config.prefix, suffix)
        return {
            role: (
                self._staging_directory / name,
                self.config.directory / name,
            )
            for role, name in names.items()
        }

    def _ordered_roles(self) -> tuple[str, ...]:
        roles = []
        if self.config.bam is None:
            roles.append("read1")
            if self.config.paired_end:
                roles.append("read2")
        if self.config.bam is not None:
            roles.append("bam")
        return tuple(roles)

    def _require_destinations_absent(self) -> None:
        collisions = [
            str(final_path)
            for _, final_path in self._paths.values()
            if os.path.lexists(final_path)
        ]
        if os.path.lexists(self._manifest_paths[1]):
            collisions.append(str(self._manifest_paths[1]))
        if collisions:
            raise OutputError(
                "output destination already exists: {}".format(", ".join(collisions))
            )

    def _validate_fragment(self, fragment: ProcessedFragment) -> None:
        _validate_processed_fragment(
            fragment,
            self.config.paired_end,
            require_base_states=self.config.bam is not None,
        )
        if fragment.fragment_ordinal != self._expected_ordinal:
            raise OutputError("processed fragments must start at zero and be consecutive")

    def _close_streams(self) -> None:
        first_error = None
        for stream in getattr(self, "_streams", {}).values():
            try:
                stream.close()
            except Exception as error:
                if first_error is None:
                    first_error = error
        if first_error is not None:
            raise OutputError(
                "cannot close staged output: {}".format(first_error)
            ) from first_error

    def _summaries(self) -> tuple[OutputFileSummary, ...]:
        summaries = []
        for role in self._ordered_roles():
            _, final_path = self._paths[role]
            stream = self._streams[role]
            record_count = (
                self._mate_count
                if role == "bam"
                else self._expected_ordinal
            )
            summaries.append(
                OutputFileSummary(
                    role=role,
                    path=final_path,
                    size_bytes=stream.size_bytes,
                    record_count=record_count,
                    sha256=stream.sha256,
                )
            )
        return tuple(summaries)

    def _require_state(self, expected: str) -> None:
        if self._state != expected:
            raise OutputError(
                "output transaction is {}, expected {}".format(self._state, expected)
            )


def _byte_view(name: str, value: BytesLike) -> memoryview:
    if not isinstance(value, (bytes, bytearray, memoryview)):
        raise OutputError("{} must be bytes-like".format(name))
    try:
        view = memoryview(value).cast("B")
    except (TypeError, ValueError) as error:
        raise OutputError("{} must be one contiguous byte buffer".format(name)) from error
    return view


def _validate_manifest_json(
    manifest_json: str, summary: OutputSummary
) -> bytes:
    if not isinstance(manifest_json, str):
        raise OutputError("manifest_json must be canonical JSON text")
    try:
        document = json.loads(manifest_json)
    except (json.JSONDecodeError, UnicodeError) as error:
        raise OutputError("manifest_json is not valid JSON") from error
    if not isinstance(document, dict) or document.get("status") != "complete":
        raise OutputError("the commit manifest must have status 'complete'")
    canonical = json.dumps(
        document,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    )
    if manifest_json != canonical:
        raise OutputError("manifest_json must use canonical compact JSON encoding")
    expected_outputs = {
        item.role: {
            "path": str(item.path),
            "record_count": item.record_count,
            "role": item.role,
            "sha256": item.sha256,
            "size_bytes": item.size_bytes,
        }
        for item in summary.files
    }
    observed_entries = document.get("outputs")
    if not isinstance(observed_entries, list):
        raise OutputError("commit manifest outputs section is missing")
    observed_outputs = {
        item.get("role"): item
        for item in observed_entries
        if isinstance(item, dict) and isinstance(item.get("role"), str)
    }
    if (
        len(observed_outputs) != len(observed_entries)
        or observed_outputs != expected_outputs
    ):
        raise OutputError("commit manifest does not identify the staged outputs")

    counts = document.get("counts")
    python_counts = counts.get("python") if isinstance(counts, dict) else None
    expected_counts = {
        "fragment_count": summary.fragment_count,
        "mate_count": summary.mate_count,
        "records_by_role": {
            item.role: item.record_count for item in summary.files
        },
    }
    if python_counts != expected_counts:
        raise OutputError("commit manifest Python counts disagree with staged outputs")
    return manifest_json.encode("utf-8") + b"\n"
