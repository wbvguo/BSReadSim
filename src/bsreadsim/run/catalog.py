"""Immutable catalog orchestration through the htsim core."""

from __future__ import annotations
from contextlib import suppress

import os
from pathlib import Path
import struct
import subprocess
import tempfile
from collections.abc import Mapping
import uuid
import zlib

from ..htsim.launch import CoreExecutableError, build_core_argv, resolve_core_executable
from .config import normalize_run_config
from .prepare import prepare_run


PathLike = str | os.PathLike


class CatalogError(RuntimeError):
    """An htsim candidate catalog could not be exported safely."""


class _BgzfWriter:
    """Small deterministic BGZF stream used for text exchange output."""

    _MAXIMUM_BLOCK_INPUT = 32768
    _EOF = bytes.fromhex(
        "1f8b08040000000000ff0600424302001b0003000000000000000000"
    )

    def __init__(self, output) -> None:
        self._output = output
        self._buffer = bytearray()
        self._closed = False

    def __enter__(self) -> _BgzfWriter:
        return self

    def __exit__(self, exception_type, exception, traceback) -> None:
        self.close()

    def write(self, data: bytes) -> int:
        if self._closed:
            raise ValueError("BGZF stream is closed")
        view = memoryview(data)
        total = len(view)
        if self._buffer:
            count = min(
                len(view), self._MAXIMUM_BLOCK_INPUT - len(self._buffer)
            )
            self._buffer.extend(view[:count])
            view = view[count:]
            if len(self._buffer) == self._MAXIMUM_BLOCK_INPUT:
                self._write_block(self._buffer)
                self._buffer.clear()
        while len(view) >= self._MAXIMUM_BLOCK_INPUT:
            self._write_block(view[: self._MAXIMUM_BLOCK_INPUT])
            view = view[self._MAXIMUM_BLOCK_INPUT :]
        self._buffer.extend(view)
        return total

    def close(self) -> None:
        if self._closed:
            return
        if self._buffer:
            self._write_block(self._buffer)
            self._buffer.clear()
        self._output.write(self._EOF)
        self._closed = True

    def _write_block(self, data) -> None:
        raw = bytes(data)
        compressor = zlib.compressobj(6, zlib.DEFLATED, -15)
        payload = compressor.compress(raw) + compressor.flush()
        block_size = 18 + len(payload) + 8
        if block_size > 65536:
            raise CatalogError("BGZF block exceeds 64 KiB")
        header = (
            b"\x1f\x8b\x08\x04"
            + struct.pack("<I", 0)
            + b"\x00\xff"
            + struct.pack("<H", 6)
            + b"BC"
            + struct.pack("<HH", 2, block_size - 1)
        )
        footer = struct.pack(
            "<II", zlib.crc32(raw) & 0xFFFFFFFF, len(raw)
        )
        self._output.write(header)
        self._output.write(payload)
        self._output.write(footer)


def export_rrbs_catalog(
    document: Mapping[str, object],
    output_path: PathLike,
    *,
    base_directory: PathLike = ".",
    core_executable: PathLike | None = None,
) -> Path:
    """Export the exact htsim RRBS candidate table to a new plain BED file."""

    loaded = normalize_run_config(document, base_directory)
    normalized = loaded.normalized
    if normalized["technology"] != "RRBS":
        raise CatalogError("candidate export requires RRBS CLI arguments")
    rrbs = normalized["rrbs"]
    coverage = normalized["coverage"]
    if (
        not isinstance(rrbs, Mapping)
        or "candidate_bed" in rrbs
        or not isinstance(coverage, Mapping)
        or coverage["kind"] != "uniform"
    ):
        raise CatalogError(
            "candidate export requires uniform RRBS arguments without a candidate BED"
        )

    effective = loaded
    if effective.master_seed is None:
        effective = effective.with_master_seed(0)

    prepared = prepare_run(effective)
    return _export_catalog(
        prepared,
        output_path,
        subcommand="rrbs-catalog",
        label="RRBS candidate",
        core_executable=core_executable,
    )


def export_methdb_catalog(
    document: Mapping[str, object],
    output_path: PathLike,
    *,
    base_directory: PathLike = ".",
    core_executable: PathLike | None = None,
) -> Path:
    """Export the exact fixed site-probability catalog used by a run."""

    loaded = normalize_run_config(document, base_directory)
    if "methdb" in loaded.normalized["inputs"]:
        raise CatalogError("cannot save a MethDB while loading another MethDB")
    effective = loaded if loaded.master_seed is not None else loaded.with_master_seed(0)
    return _export_catalog(
        prepare_run(effective),
        output_path,
        subcommand="methdb-catalog",
        label="MethDB",
        core_executable=core_executable,
    )


def export_methdb_bed(
    input_path: PathLike,
    output_path: PathLike,
    *,
    compressed: bool = True,
    core_executable: PathLike | None = None,
) -> Path:
    """Decode one MethDB snapshot into a new human-readable extended BED."""

    try:
        executable = resolve_core_executable(core_executable)
    except CoreExecutableError as error:
        raise CatalogError(str(error)) from error
    try:
        source = Path(input_path).expanduser().resolve(strict=True)
    except (TypeError, OSError) as error:
        raise CatalogError("cannot resolve input MethDB: {}".format(input_path)) from error
    if not source.is_file():
        raise CatalogError("input MethDB is not a regular file: {}".format(source))
    try:
        destination = Path(output_path).expanduser().resolve(strict=False)
    except (TypeError, OSError) as error:
        raise CatalogError("MethDB BED output path is invalid") from error
    expected_suffix = ".bed.gz" if compressed else ".bed"
    if not str(destination).lower().endswith(expected_suffix):
        raise CatalogError(
            "MethDB BED output path must end in {}".format(expected_suffix)
        )

    created = False
    try:
        destination.parent.mkdir(parents=True, exist_ok=True)
        with destination.open("xb") as output:
            created = True
            argv = (str(executable), "methdb-export", str(source))
            if compressed:
                with tempfile.TemporaryFile() as errors:
                    process = subprocess.Popen(
                        argv,
                        stdout=subprocess.PIPE,
                        stderr=errors,
                    )
                    if process.stdout is None:
                        process.kill()
                        process.wait()
                        raise CatalogError(
                            "htsim MethDB BED export did not open stdout"
                        )
                    try:
                        with _BgzfWriter(output) as encoded:
                            while True:
                                block = process.stdout.read(1024 * 1024)
                                if not block:
                                    break
                                encoded.write(block)
                    except BaseException:
                        process.kill()
                        process.wait()
                        raise
                    finally:
                        process.stdout.close()
                    returncode = process.wait()
                    errors.seek(0)
                    stderr = errors.read()
            else:
                completed = subprocess.run(
                    argv,
                    stdout=output,
                    stderr=subprocess.PIPE,
                    check=False,
                )
                returncode = completed.returncode
                stderr = completed.stderr
        if returncode != 0:
            detail = stderr.decode("utf-8", errors="replace").strip()
            raise CatalogError(
                "htsim MethDB BED export failed{}".format(
                    ": " + detail if detail else ""
                )
            )
    except CatalogError:
        if created:
            with suppress(OSError):
                destination.unlink()
        raise
    except FileExistsError as error:
        raise CatalogError(
            "MethDB BED output already exists: {}".format(destination)
        ) from error
    except OSError as error:
        if created:
            with suppress(OSError):
                destination.unlink()
        raise CatalogError(
            "cannot write MethDB BED output {}: {}".format(destination, error)
        ) from error

    return destination


def export_variant_catalog(
    document: Mapping[str, object],
    output_path: PathLike,
    *,
    base_directory: PathLike = ".",
    core_executable: PathLike | None = None,
) -> Path:
    """Export the normalized phased variant truth as a new ``.vcf.gz`` file."""

    loaded = normalize_run_config(document, base_directory)
    if not str(output_path).lower().endswith(".vcf.gz"):
        raise CatalogError("variant VCF output path must end in .vcf.gz")

    effective = loaded if loaded.master_seed is not None else loaded.with_master_seed(0)
    return _export_catalog(
        prepare_run(effective),
        output_path,
        subcommand="variant-catalog",
        label="variant VCF",
        core_executable=core_executable,
        gzip_output=True,
    )


def _export_catalog(
    prepared,
    output_path: PathLike,
    *,
    subcommand: str,
    label: str,
    core_executable: PathLike | None,
    gzip_output: bool = False,
) -> Path:
    try:
        executable = resolve_core_executable(core_executable)
    except CoreExecutableError as error:
        raise CatalogError(str(error)) from error
    base_argv = build_core_argv(
        prepared,
        str(uuid.uuid4()),
        executable,
        emit_details=False,
    )
    argv = (base_argv[0], subcommand, *base_argv[1:])

    try:
        destination = Path(output_path).expanduser().resolve(strict=False)
    except (TypeError, OSError) as error:
        raise CatalogError("{} output path is invalid".format(label)) from error
    if not destination.parent.is_dir():
        raise CatalogError(
            "{} output parent is not a directory: {}".format(label, destination.parent)
        )

    created = False
    try:
        with destination.open("xb") as output:
            created = True
            if gzip_output:
                with tempfile.TemporaryFile() as errors:
                    process = subprocess.Popen(
                        argv,
                        stdout=subprocess.PIPE,
                        stderr=errors,
                    )
                    if process.stdout is None:
                        process.kill()
                        process.wait()
                        raise CatalogError(
                            "htsim {} export did not open stdout".format(label)
                        )
                    try:
                        with _BgzfWriter(output) as compressed:
                            while True:
                                block = process.stdout.read(1024 * 1024)
                                if not block:
                                    break
                                compressed.write(block)
                    except BaseException:
                        process.kill()
                        process.wait()
                        raise
                    finally:
                        process.stdout.close()
                    returncode = process.wait()
                    errors.seek(0)
                    stderr = errors.read()
            else:
                completed = subprocess.run(
                    argv,
                    stdout=output,
                    stderr=subprocess.PIPE,
                    check=False,
                )
                returncode = completed.returncode
                stderr = completed.stderr
        if returncode != 0:
            detail = stderr.decode("utf-8", errors="replace").strip()
            raise CatalogError(
                "htsim {} export failed{}".format(
                    label,
                    ": " + detail if detail else ""
                )
            )
    except CatalogError:
        if created:
            with suppress(OSError):
                destination.unlink()
        raise
    except FileExistsError as error:
        raise CatalogError(
            "{} output already exists: {}".format(label, destination)
        ) from error
    except OSError as error:
        if created:
            with suppress(OSError):
                destination.unlink()
        raise CatalogError(
            "cannot write {} output {}: {}".format(label, destination, error)
        ) from error

    return destination


__all__ = [
    "CatalogError",
    "export_methdb_bed",
    "export_methdb_catalog",
    "export_rrbs_catalog",
    "export_variant_catalog",
]
