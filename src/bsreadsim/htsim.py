"""htsim, the C++ core (``HtsimCore``), and the fragment stream of a simulation
run (``HtsimProcess``). No other module knows its command line. The package
bundles the executable as ``bin/htsim``; ``main`` is the ``htsim`` command.

The stream (format in htsim/src/stream.h) is

    header   one line of JSON (htsim's version, which must be the package's)
    batch*   u64 payload size, then the payload (``FragmentBatch.decode``)
    end      u64 zero
    summary  one line of JSON
"""

from __future__ import annotations

from collections import deque
from collections.abc import Iterator, Sequence
from contextlib import suppress
from dataclasses import dataclass
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import threading
from typing import IO, NoReturn

from . import __version__
from .errors import BSReadSimError
from .settings import Settings

CORE_FILENAME = "htsim.exe" if os.name == "nt" else "htsim"
BUNDLED_CORE = Path(__file__).resolve().parent / "bin" / CORE_FILENAME
_SIZE = struct.Struct("<Q")

# Settings passed to htsim as --kebab-case options (plus a few computed ones)
CORE_FIELDS = (
    "seed", "seed_mut", "seed_phase", "seed_meth", "reference", "vcf", "cgmap",
    "bed_methyl", "methbg", "methbed", "methdb", "asm", "asm_bed", "technology",
    "directional", "paired_end", "read_length", "insert_min", "insert_mean", "insert_max",
    "insert_sd", "depth", "max_ambiguous_fraction", "mutation_rate", "indel_fraction",
    "indel_extension_probability", "homozygous_only", "cpg_only", "pool_meth", "beta_cg",
    "beta_chg", "beta_chh", "sampling", "gc_profile", "cut_sites", "rrbs_candidates",
    "targets", "center_sd",
)
# What ``build`` makes -> the htsim command that makes it
BUILD_COMMANDS = {"rrbs": "rrbs-catalog", "variants": "variant-catalog", "methdb": "methdb-build"}


# ---------------------------------------------------------------- the executable

@dataclass(frozen=True)
class HtsimCore:
    """The htsim executable: its command line, its commands, and simulation runs.

    The commands that write a file (``build``: RRBS candidates, a variant VCF,
    or a MethDB; ``export_methdb``) never replace an existing file and remove
    what they wrote when htsim fails.
    """

    path: Path

    @classmethod
    def find(cls, path: os.PathLike | str | None = None) -> HtsimCore:
        """The htsim bundled with the package, or the executable at ``path``."""
        candidate = BUNDLED_CORE if path is None else Path(path).expanduser()
        if not candidate.is_file() or not os.access(candidate, os.X_OK):
            raise BSReadSimError(f"htsim is not an executable file: {candidate}")
        return cls(candidate.resolve())

    def argv(self, settings: Settings, command: str | None = None, **extra) -> list[str]:
        """The command line of ``settings``; ``extra`` adds options like threads."""
        options = {name: getattr(settings, name) for name in CORE_FIELDS}
        options["fragments"] = settings.fragments
        options.update(details=False, batch_size=1024, threads=1)
        options.update(extra)
        argv = [str(self.path)] + ([command] if command else [])
        for name, value in options.items():
            if value is None or value == ():
                continue
            if isinstance(value, bool):
                value = "true" if value else "false"
            elif isinstance(value, (tuple, list)):
                value = ",".join(map(str, value))
            argv += ["--" + name.replace("_", "-"), str(value)]
        return argv

    def stream(self, settings: Settings, **extra) -> HtsimProcess:
        """Start the simulation run of ``settings``; ``extra`` as for ``argv``.

        htsim must be of this package's version: then its stream, its random
        numbers, and everything else it does are the ones this package expects.
        """
        process = HtsimProcess(self.argv(settings, **extra))
        if process.header.core_version != __version__:
            process.kill()
            raise BSReadSimError(
                f"htsim {process.header.core_version} does not match bsreadsim {__version__}")
        return process

    def validate(self, settings: Settings) -> dict:
        """htsim's check of the inputs of ``settings`` against the reference."""
        return json.loads(_run(self.argv(settings, "validate-inputs")).stdout)

    def build(self, settings: Settings, kind: str, output: str | os.PathLike) -> Path:
        """Build RRBS candidates (``kind`` rrbs), a variant VCF (variants), or a MethDB (methdb)."""
        destination = _new_file(output)
        _run_into(self.argv(settings, BUILD_COMMANDS[kind], output=destination), destination)
        return destination

    def export_methdb(self, source: str | os.PathLike, output: str | os.PathLike) -> Path:
        """Decode a MethDB into extended BED; ``output`` ending in .gz is compressed."""
        source = Path(source).expanduser()
        if not source.is_file():
            raise BSReadSimError(f"input MethDB is not a file: {source}")
        destination = _new_file(output)
        _run_into([str(self.path), "methdb-export", str(source.resolve()), str(destination)],
                  destination)
        return destination

    def sam_to_bam(self, level: int, threads: int) -> list[str]:
        """The command line that turns SAM text on stdin into BAM on stdout."""
        return [str(self.path), "--sam-to-bam", str(level), str(threads)]


def _run(argv: Sequence[str]) -> subprocess.CompletedProcess:
    """Run an htsim command that does not stream fragments."""
    try:
        completed = subprocess.run(list(argv), capture_output=True, check=False)
    except OSError as error:
        raise BSReadSimError(f"cannot run htsim: {error}") from error
    if completed.returncode != 0:
        raise BSReadSimError(_core_message(completed.stderr.decode("utf-8", errors="replace"))
                             or f"htsim exited with status {completed.returncode}")
    return completed


def _core_message(stderr: str) -> str:
    """htsim's last error line, without its program prefix."""
    lines = [line for line in stderr.strip().splitlines() if line.strip()]
    return lines[-1].removeprefix("htsim: ") if lines else ""


def _new_file(output: str | os.PathLike) -> Path:
    destination = Path(output).expanduser().resolve()
    if destination.exists():
        raise BSReadSimError(f"output already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    return destination


def _run_into(argv: Sequence[str], destination: Path) -> None:
    """Run a command that writes ``destination``; remove it if the command fails."""
    try:
        _run(argv)
    except BSReadSimError:
        with suppress(OSError):
            destination.unlink()
        raise


# ---------------------------------------------------------------- the fragment stream

@dataclass(frozen=True)
class Contig:
    name: str
    length: int
    md5: str  # hex digest of the uppercase sequence


@dataclass(frozen=True)
class Header:
    core_version: str
    contigs: tuple[Contig, ...]


@dataclass(frozen=True)
class Summary:
    fragment_count: int
    mate_count: int
    template_base_count: int
    methylation_site_count: int
    skipped_fragment_count: int
    per_contig_fragment_counts: tuple[int, ...]


def _parse_header(line: bytes) -> Header:
    try:
        document = json.loads(line)
        return Header(
            core_version=document["version"],
            contigs=tuple(Contig(name, length, md5)
                          for name, length, md5 in document["contigs"]),
        )
    except (ValueError, KeyError, TypeError) as error:
        raise BSReadSimError("htsim did not begin its output with a stream header") from error


def _parse_summary(line: bytes) -> Summary:
    document = json.loads(line)
    document["per_contig_fragment_counts"] = tuple(document["per_contig_fragment_counts"])
    return Summary(**document)


class HtsimProcess:
    """One htsim run whose stdout is a fragment stream.

    Use as a context manager: read ``header``, iterate ``batches()`` to the
    end, then read ``summary``. Leaving the block early kills htsim.
    """

    def __init__(self, argv: Sequence[str | os.PathLike]) -> None:
        self.argv = tuple(str(value) for value in argv)
        self.summary: Summary | None = None
        self._stderr = deque(maxlen=200)
        try:
            self._process = subprocess.Popen(
                self.argv, stdin=subprocess.DEVNULL, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, bufsize=1024 * 1024)
        except OSError as error:
            raise BSReadSimError(f"cannot run htsim: {error}") from error
        self._stderr_thread = threading.Thread(
            target=self._collect_stderr, name="htsim-stderr", daemon=True)
        self._stderr_thread.start()
        line = self._process.stdout.readline()
        if not line:
            self._fail("htsim exited before writing its stream header")
        self.header = _parse_header(line)

    def __enter__(self) -> HtsimProcess:
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if exc_type is not None or self.summary is None:
            self.kill()
        self._process.wait()
        self._stderr_thread.join()
        if exc_type is None:
            if self.summary is None:
                self._fail("the fragment stream was not read to its end")
            if self._process.returncode != 0:
                self._fail(f"htsim exited with status {self._process.returncode}")
        return False

    def kill(self) -> None:
        """Stop htsim before the end of its stream."""
        self._process.kill()
        self._process.wait()
        self._stderr_thread.join()

    def batches(self) -> Iterator[bytes]:
        """Yield every batch payload, then read the summary."""
        stdout = self._process.stdout
        while True:
            size = _SIZE.unpack(self._read(stdout, _SIZE.size))[0]
            if size == 0:
                break
            yield self._read(stdout, size)
        line = stdout.readline()
        if not line:
            self._fail("htsim stream ended without its summary")
        self.summary = _parse_summary(line)

    def _read(self, stdout: IO[bytes], size: int) -> bytes:
        data = stdout.read(size)
        if len(data) != size:
            self._process.wait()
            self._fail("htsim stream was truncated")
        return data

    def _collect_stderr(self) -> None:
        for line in self._process.stderr:
            self._stderr.append(line.decode("utf-8", errors="replace").rstrip())

    def _fail(self, message: str) -> NoReturn:
        self._process.wait()
        self._stderr_thread.join(timeout=5)
        raise BSReadSimError(_core_message("\n".join(self._stderr)) or message)


# ---------------------------------------------------------------- the htsim command

def main() -> None:
    """Run the bundled htsim with the command's arguments (it replaces this process)."""
    try:
        core = HtsimCore.find()
    except BSReadSimError as error:
        sys.exit(f"htsim: {error}")
    os.execv(core.path, ["htsim", *sys.argv[1:]])
