"""Transactional FASTQ, truth, and manifest output for processed fragments.

The output component accepts only format-ready ``ProcessedFragment`` values.
It stages every byte on the destination filesystem, enforces fragment order
and SE/PE cardinality.  Finalization computes data-file identities while they
remain private; an explicit commit publishes the manifest last as the marker
that makes the whole file set valid.
"""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
from typing import BinaryIO, Dict, Optional, Tuple, Union

from .bam import format_truth_sam_fragment
from .read_names import ReadNameError, format_fragment_identifier

from .postprocess import (
    BaseAnnotation,
    ProcessedFragment,
    ProcessedMate,
    _CompactAnnotations,
)
from .model import VariantEvent

try:
    from ._native import (
        canonical_truth_json_bytes as _native_canonical_truth_json_bytes,
    )
except ImportError:
    _native_canonical_truth_json_bytes = None


_PREFIX = re.compile(r"^[A-Za-z0-9._-]+$")
_VALID_SEQUENCE = frozenset("ACGTN")
_BASE_CHARS = "ACGTN"
BytesLike = Union[bytes, bytearray, memoryview]

# Every truth dictionary below is declared in canonical lexicographic key
# order.  Python 3.8+ preserves that insertion order, so re-sorting hundreds of
# fixed-schema per-base dictionaries only burns worker CPU.  The value graph is
# built here from fresh primitive containers and therefore cannot be circular.
_TRUTH_JSON_ENCODER = json.JSONEncoder(
    allow_nan=False,
    check_circular=False,
    ensure_ascii=False,
    separators=(",", ":"),
    sort_keys=False,
)
NATIVE_TRUTH_JSON_AVAILABLE = _native_canonical_truth_json_bytes is not None


class OutputError(RuntimeError):
    """The staged output set cannot be completed safely."""


@dataclass(frozen=True)
class TruthBamConfig:
    """Validated streaming contract for one HTSlib-backed truth BAM."""

    writer_argv: Tuple[str, ...]
    sam_header: bytes
    references: Tuple[Tuple[str, int], ...]
    read_group_id: str

    def __post_init__(self) -> None:
        if (
            not isinstance(self.writer_argv, tuple)
            or not self.writer_argv
            or any(
                not isinstance(value, str) or not value or "\x00" in value
                for value in self.writer_argv
            )
        ):
            raise OutputError("truth BAM writer argv must be non-empty text")
        if type(self.sam_header) is not bytes or not self.sam_header.startswith(
            b"@HD\t"
        ):
            raise OutputError("truth BAM SAM header is invalid")
        if not self.sam_header.endswith(b"\n"):
            raise OutputError("truth BAM SAM header must end with a newline")
        if not isinstance(self.references, tuple) or not self.references:
            raise OutputError("truth BAM references must be non-empty")
        observed = set()
        for reference in self.references:
            if (
                not isinstance(reference, tuple)
                or len(reference) != 2
                or not isinstance(reference[0], str)
                or not reference[0]
                or isinstance(reference[1], bool)
                or not isinstance(reference[1], int)
                or reference[1] <= 0
                or reference[0] in observed
            ):
                raise OutputError("truth BAM references are invalid")
            observed.add(reference[0])
        if not isinstance(self.read_group_id, str) or not self.read_group_id:
            raise OutputError("truth BAM read group identifier is invalid")

    def reference_length(self, name: str) -> int:
        for reference_name, length in self.references:
            if reference_name == name:
                return length
        raise OutputError("truth BAM fragment refers to an unknown contig")


@dataclass(frozen=True)
class OutputConfig:
    """Validated output policy for one run."""

    directory: Path
    prefix: str
    paired_end: bool
    compression: str = "gzip"
    gzip_level: int = 6
    truth: str = "full"
    truth_bam: Optional[TruthBamConfig] = None

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
        if self.truth not in ("full", "none"):
            raise OutputError("truth must be 'full' or 'none'")
        if self.truth_bam is not None and not isinstance(
            self.truth_bam, TruthBamConfig
        ):
            raise OutputError("truth_bam must be a TruthBamConfig or None")
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
    files: Tuple[OutputFileSummary, ...]

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


class _BamOutput:
    """Stream SAM into the bundled HTSlib writer and identify its BAM bytes."""

    _WAIT_SECONDS = 30

    def __init__(self, path: Path, config: TruthBamConfig) -> None:
        self.path = path
        self.closed = False
        self._completed = False
        self._size_bytes = 0
        self._sha256 = ""
        self._stderr_path = path.with_name(path.name + ".writer-stderr")
        self.raw = None  # type: Optional[BinaryIO]
        self.stderr = None  # type: Optional[BinaryIO]
        self.process = None  # type: Optional[subprocess.Popen]
        try:
            self.raw = path.open("xb")
            self.stderr = self._stderr_path.open("xb")
            self.process = subprocess.Popen(
                config.writer_argv,
                stdin=subprocess.PIPE,
                stdout=self.raw,
                stderr=self.stderr,
                close_fds=True,
            )
            if self.process.stdin is None:
                raise OutputError("truth BAM writer has no SAM input stream")
            self.write_bytes(config.sam_header)
        except Exception:
            self._stop_process()
            self._close_files()
            self._unlink_stderr()
            raise

    def write_bytes(self, value: BytesLike) -> None:
        if self.closed:
            raise OutputError("cannot write to a closed truth BAM")
        if not isinstance(value, (bytes, bytearray, memoryview)):
            raise OutputError("truth BAM accepts bytes-like SAM records only")
        if self.process is None or self.process.stdin is None:
            raise OutputError("truth BAM writer is unavailable")
        try:
            self.process.stdin.write(value)
        except (BrokenPipeError, OSError) as error:
            raise OutputError(
                "truth BAM writer closed its SAM input early{}".format(
                    self._stderr_suffix()
                )
            ) from error

    def close(self) -> None:
        if self.closed:
            return
        error = None  # type: Optional[BaseException]
        status = None
        try:
            if self.process is None or self.process.stdin is None:
                raise OutputError("truth BAM writer is unavailable")
            try:
                self.process.stdin.close()
            except (BrokenPipeError, OSError):
                pass
            try:
                status = self.process.wait(timeout=self._WAIT_SECONDS)
            except subprocess.TimeoutExpired as timeout_error:
                self._stop_process()
                raise OutputError(
                    "truth BAM writer did not terminate"
                ) from timeout_error
            if self.raw is None:
                raise OutputError("truth BAM staged file is unavailable")
            self.raw.flush()
            os.fsync(self.raw.fileno())
            if status != 0:
                raise OutputError(
                    "truth BAM writer exited with status {}{}".format(
                        status, self._stderr_suffix()
                    )
                )
        except BaseException as observed:
            error = observed
        finally:
            self._close_files()
            self.closed = True

        if error is None:
            try:
                digest = hashlib.sha256()
                size = 0
                with self.path.open("rb") as staged:
                    while True:
                        block = staged.read(1024 * 1024)
                        if not block:
                            break
                        size += len(block)
                        digest.update(block)
                self._size_bytes = size
                self._sha256 = digest.hexdigest()
                self._completed = True
            except BaseException as observed:
                error = observed
        self._unlink_stderr()
        if error is not None:
            if isinstance(error, OutputError):
                raise error
            raise OutputError("cannot finalize truth BAM: {}".format(error)) from error

    def _stop_process(self) -> None:
        process = self.process
        if process is None:
            return
        if process.stdin is not None:
            try:
                process.stdin.close()
            except OSError:
                pass
        if process.poll() is None:
            process.terminate()
            try:
                process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                process.kill()
                process.wait(timeout=5)

    def _close_files(self) -> None:
        for stream in (self.raw, self.stderr):
            if stream is not None and not stream.closed:
                try:
                    stream.close()
                except OSError:
                    pass

    def _stderr_suffix(self) -> str:
        try:
            if self.stderr is not None and not self.stderr.closed:
                self.stderr.flush()
            with self._stderr_path.open("rb") as stream:
                value = stream.read(8192).decode("utf-8", errors="replace").strip()
        except OSError:
            return ""
        return ": {}".format(value) if value else ""

    def _unlink_stderr(self) -> None:
        try:
            self._stderr_path.unlink()
        except FileNotFoundError:
            pass

    @property
    def size_bytes(self) -> int:
        if not self._completed:
            raise OutputError("staged truth BAM identity requires a closed stream")
        return self._size_bytes

    @property
    def sha256(self) -> str:
        if not self._completed:
            raise OutputError("staged truth BAM identity requires a closed stream")
        return self._sha256


class OutputTransaction:
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
            self._streams = {}  # type: Dict[str, Union[_BinaryOutput, _BamOutput]]
            for role, (staged_path, final_path) in self._paths.items():
                if role == "truth_bam":
                    if config.truth_bam is None:
                        raise OutputError("truth BAM path has no writer policy")
                    self._streams[role] = _BamOutput(
                        staged_path, config.truth_bam
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

    def __enter__(self) -> "OutputTransaction":
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
            read1, read2, truth = _format_fragment_records_trusted(
                fragment,
                paired_end=self.config.paired_end,
                include_truth=self.config.truth == "full",
            )
            self._streams["read1"].write_bytes(read1)
            if read2 is not None:
                self._streams["read2"].write_bytes(read2)
            if truth is not None:
                self._streams["truth"].write_bytes(truth)
            if self.config.truth_bam is not None:
                bam_records = format_truth_sam_fragment(
                    fragment,
                    paired_end=self.config.paired_end,
                    read_group_id=self.config.truth_bam.read_group_id,
                    contig_length=self.config.truth_bam.reference_length(
                        fragment.contig_name
                    ),
                )
                self._streams["truth_bam"].write_bytes(b"".join(bam_records))
            self._expected_ordinal += 1
            self._mate_count += len(fragment.mates)
        except Exception:
            self.abort()
            raise

    def write_formatted_batch(
        self,
        first_ordinal: int,
        record_lengths: Tuple[Tuple[int, int, int], ...],
        read1: BytesLike,
        read2: Optional[BytesLike],
        truth: Optional[BytesLike],
        *,
        alignment_record_lengths: Tuple[int, ...] = (),
        alignment: Optional[BytesLike] = None,
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

            read1_view = _byte_view("formatted read1", read1)
            views.append(read1_view)
            truth_enabled = self.config.truth == "full"
            truth_view = None
            if truth is not None:
                truth_view = _byte_view("formatted truth", truth)
                views.append(truth_view)
            if truth_enabled != (truth_view is not None):
                raise OutputError(
                    "formatted truth presence disagrees with output policy"
                )
            read2_view = None
            if read2 is not None:
                read2_view = _byte_view("formatted read2", read2)
                views.append(read2_view)
            if self.config.paired_end != (read2_view is not None):
                raise OutputError("formatted batch mate cardinality disagrees with output")

            alignment_enabled = self.config.truth_bam is not None
            alignment_view = None
            if alignment is not None:
                alignment_view = _byte_view("formatted truth alignment", alignment)
                views.append(alignment_view)
            if alignment_enabled != (alignment_view is not None):
                raise OutputError(
                    "formatted truth alignment presence disagrees with output policy"
                )
            if not isinstance(alignment_record_lengths, tuple):
                raise OutputError(
                    "formatted truth alignment lengths must be an immutable tuple"
                )
            expected_alignment_records = len(record_lengths) * (
                2 if self.config.paired_end else 1
            )
            if alignment_enabled:
                if len(alignment_record_lengths) != expected_alignment_records:
                    raise OutputError(
                        "formatted truth alignment count disagrees with mate count"
                    )
                if any(
                    isinstance(value, bool)
                    or not isinstance(value, int)
                    or value <= 0
                    for value in alignment_record_lengths
                ):
                    raise OutputError(
                        "formatted truth alignment lengths must be positive"
                    )
            elif alignment_record_lengths:
                raise OutputError(
                    "formatted truth alignment lengths disagree with output policy"
                )

            read1_total = 0
            read2_total = 0
            truth_total = 0
            for lengths in record_lengths:
                if (
                    not isinstance(lengths, tuple)
                    or len(lengths) != 3
                    or any(
                        isinstance(value, bool)
                        or not isinstance(value, int)
                        or value < 0
                        for value in lengths
                    )
                ):
                    raise OutputError("formatted batch contains invalid record lengths")
                read1_length, read2_length, truth_length = lengths
                if read1_length == 0:
                    raise OutputError("formatted read1 records must be non-empty")
                if truth_enabled != (truth_length > 0):
                    raise OutputError(
                        "formatted truth record lengths disagree with output policy"
                    )
                if self.config.paired_end != (read2_length > 0):
                    raise OutputError(
                        "formatted record mate cardinality disagrees with output"
                    )
                read1_total += read1_length
                read2_total += read2_length
                truth_total += truth_length
            if read1_total != len(read1_view):
                raise OutputError("formatted batch region sizes disagree with records")
            if truth_view is not None and truth_total != len(truth_view):
                raise OutputError("formatted truth region size disagrees with records")
            if read2_view is not None and read2_total != len(read2_view):
                raise OutputError("formatted read2 region size disagrees with records")
            if alignment_view is not None and sum(
                alignment_record_lengths
            ) != len(alignment_view):
                raise OutputError(
                    "formatted truth alignment region size disagrees with records"
                )

            self._streams["read1"].write_bytes(read1_view)
            if read2_view is not None:
                self._streams["read2"].write_bytes(read2_view)
            if truth_view is not None:
                self._streams["truth"].write_bytes(truth_view)
            if alignment_view is not None:
                self._streams["truth_bam"].write_bytes(alignment_view)

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
            try:
                final_path.unlink()
            except FileNotFoundError:
                pass
        self._published.clear()
        shutil.rmtree(self._staging_directory, ignore_errors=True)
        self._state = "aborted"

    def _build_paths(self) -> Dict[str, Tuple[Path, Path]]:
        suffix = ".gz" if self.config.compression == "gzip" else ""
        names = {
            "read1": "{}.R1.fastq{}".format(self.config.prefix, suffix),
        }
        if self.config.truth == "full":
            names["truth"] = "{}.truth.jsonl{}".format(
                self.config.prefix, suffix
            )
        if self.config.truth_bam is not None:
            names["truth_bam"] = "{}.truth.bam".format(self.config.prefix)
        if self.config.paired_end:
            names["read2"] = "{}.R2.fastq{}".format(self.config.prefix, suffix)
        return {
            role: (
                self._staging_directory / name,
                self.config.directory / name,
            )
            for role, name in names.items()
        }

    def _ordered_roles(self) -> Tuple[str, ...]:
        roles = ["read1"]
        if self.config.paired_end:
            roles.append("read2")
        if self.config.truth == "full":
            roles.append("truth")
        if self.config.truth_bam is not None:
            roles.append("truth_bam")
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
            require_annotations=(
                self.config.truth == "full" or self.config.truth_bam is not None
            ),
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

    def _summaries(self) -> Tuple[OutputFileSummary, ...]:
        summaries = []
        for role in self._ordered_roles():
            _, final_path = self._paths[role]
            stream = self._streams[role]
            record_count = (
                self._mate_count
                if role == "truth_bam"
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


def _validate_processed_fragment(
    fragment: ProcessedFragment,
    paired_end: bool,
    *,
    require_annotations: bool = True,
) -> None:
    if not isinstance(fragment, ProcessedFragment):
        raise OutputError("fragment must be a ProcessedFragment")
    if not isinstance(paired_end, bool):
        raise OutputError("paired_end must be a boolean")
    if not isinstance(require_annotations, bool):
        raise OutputError("require_annotations must be a boolean")
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
        _validate_mate(mate, require_annotations=require_annotations)


def _validate_mate(
    mate: ProcessedMate,
    *,
    require_annotations: bool = True,
) -> None:
    if not isinstance(mate, ProcessedMate):
        raise OutputError("fragment mates must be ProcessedMate values")
    if len(mate.sequence) != len(mate.quality):
        raise OutputError("FASTQ sequence and quality lengths differ")
    if not mate.sequence or any(base not in _VALID_SEQUENCE for base in mate.sequence):
        raise OutputError("FASTQ sequence must be non-empty uppercase A/C/G/T/N")
    if any(not 33 <= ord(character) <= 126 for character in mate.quality):
        raise OutputError("FASTQ quality contains a non-Phred+33 character")
    if require_annotations and len(mate.annotations) != len(mate.sequence):
        raise OutputError("mate annotations must cover every read base")


def _format_fragment_records(
    fragment: ProcessedFragment,
    *,
    paired_end: bool,
    include_truth: bool = True,
) -> Tuple[bytes, Optional[bytes], Optional[bytes]]:
    """Validate and format one fragment without owning publication order."""

    _validate_processed_fragment(
        fragment,
        paired_end,
        require_annotations=include_truth,
    )
    return _format_fragment_records_trusted(
        fragment,
        paired_end=paired_end,
        include_truth=include_truth,
    )


def _format_fragment_records_trusted(
    fragment: ProcessedFragment,
    *,
    paired_end: bool,
    include_truth: bool = True,
) -> Tuple[bytes, Optional[bytes], Optional[bytes]]:
    """Format a locally produced, already validated process-worker value."""

    records = {}  # type: Dict[str, bytes]
    for mate in sorted(fragment.mates, key=lambda item: item.mate_index):
        role = "read1" if mate.mate_index == 0 else "read2"
        records[role] = _fastq_record(fragment, mate).encode("utf-8")
    truth = (
        _canonical_truth_json_bytes(fragment, newline=True)
        if include_truth
        else None
    )
    return records["read1"], records.get("read2"), truth


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


def _canonical_truth_json(fragment: ProcessedFragment) -> str:
    if _native_canonical_truth_json_bytes is not None:
        return _native_canonical_truth_json_bytes(fragment).decode("utf-8")
    return _canonical_truth_json_python(fragment)


def _canonical_truth_json_bytes(
    fragment: ProcessedFragment,
    *,
    newline: bool = False,
) -> bytes:
    if not isinstance(newline, bool):
        raise OutputError("truth JSON newline flag must be a boolean")
    if _native_canonical_truth_json_bytes is not None:
        return _native_canonical_truth_json_bytes(fragment, newline)
    suffix = "\n" if newline else ""
    return (_canonical_truth_json_python(fragment) + suffix).encode("utf-8")


def _canonical_truth_json_python(fragment: ProcessedFragment) -> str:
    """Reference encoder retained for optional-native equivalence checks."""

    value = {
        "contig": fragment.contig_name,
        "fragment_conversion_mode": fragment.fragment_conversion_mode.name,
        "fragment_ordinal": fragment.fragment_ordinal,
        "haplotype": fragment.haplotype,
        "mates": [_mate_truth(mate) for mate in fragment.mates],
        "site_states": [
            {
                "allele": site.allele.name,
                "context": site.context.name,
                "methylated": site.methylated,
                "probability": site.probability,
                "reference_pos": site.reference_pos,
                "site_index": site.site_index,
                "source": site.source.name,
                "template_offset": site.template_offset,
            }
            for site in fragment.site_states
        ],
        "variant_events": [
            _variant_event_truth(event) for event in fragment.variant_events
        ],
    }
    return _TRUTH_JSON_ENCODER.encode(value)


def _mate_truth(mate: ProcessedMate) -> dict:
    if isinstance(mate.annotations, _CompactAnnotations):
        annotations = [
            _compact_annotation_truth(mate.annotations, offset)
            for offset in range(len(mate.annotations))
        ]
    else:
        annotations = [_annotation_truth(item) for item in mate.annotations]
    return {
        "annotations": annotations,
        "conversion_mode": mate.conversion_mode.name,
        "mate_index": mate.mate_index,
        "reference_end": mate.reference_end,
        "reference_start": mate.reference_start,
        "reverse_complement": mate.reverse_complement,
        "sequence": mate.sequence,
    }


def _annotation_truth(annotation: BaseAnnotation) -> dict:
    return {
        "conversion_attempted": annotation.conversion_attempted,
        "conversion_succeeded": annotation.conversion_succeeded,
        "final_base": annotation.final_base,
        "methylated": annotation.methylated,
        "oriented_base": annotation.oriented_base,
        "post_conversion_base": annotation.post_conversion_base,
        "quality_phred": annotation.quality_phred,
        "read_offset": annotation.read_offset,
        "reference_pos": annotation.reference_pos,
        "sequencing_error": annotation.sequencing_error,
        "site_index": annotation.site_index,
        "variant_event_id": annotation.variant_event_id,
    }


def _compact_annotation_truth(
    annotations: _CompactAnnotations,
    offset: int,
) -> dict:
    return {
        "conversion_attempted": annotations.attempted[offset],
        "conversion_succeeded": annotations.succeeded[offset],
        "final_base": annotations.final_bases[offset],
        "methylated": annotations.methylated[offset],
        "oriented_base": annotations.oriented_bases[offset],
        "post_conversion_base": annotations.post_conversion_bases[offset],
        "quality_phred": annotations.quality_phreds[offset],
        "read_offset": offset,
        "reference_pos": annotations.reference_positions[offset],
        "sequencing_error": annotations.error_flags[offset],
        "site_index": annotations.site_indices[offset],
        "variant_event_id": annotations.variant_event_ids[offset],
    }


def _variant_event_truth(event: VariantEvent) -> dict:
    if not isinstance(event, VariantEvent):
        raise OutputError("variant events must be protocol VariantEvent values")
    try:
        ref_bases = "".join(_BASE_CHARS[base] for base in event.ref_bases)
        alt_bases = "".join(_BASE_CHARS[base] for base in event.alt_bases)
    except (IndexError, TypeError) as error:
        raise OutputError("variant event bases are outside protocol encoding") from error
    return {
        "alt_bases": alt_bases,
        "event_id": event.event_id,
        "kind": event.kind.name,
        "phased_haplotype": event.phased_haplotype,
        "ref_bases": ref_bases,
        "reference_end": event.reference_end,
        "reference_start": event.reference_start,
    }


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
    try:
        # Delayed import avoids a module cycle: manifest construction consumes
        # OutputSummary, while commit verifies the completed manifest value.
        from .manifest import (  # pylint: disable=import-outside-toplevel
            MANIFEST_SCHEMA_VERSION,
            ManifestError,
            verify_complete_manifest,
        )

        if document.get("manifest_schema_version") != MANIFEST_SCHEMA_VERSION:
            raise OutputError("manifest schema version is unsupported")
        verify_complete_manifest(document)
    except ManifestError as error:
        raise OutputError("commit manifest verification failed: {}".format(error)) from error

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
