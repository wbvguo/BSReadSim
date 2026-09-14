"""Everything a run writes: FASTQ and SAM records, BAM, the manifest, and the files.

A run's output format (``FastqFormat`` or ``SamFormat``) turns the reads of a
batch into the records of each output file.

* FASTQ read names are ``@<contig>:<start>-<end>:<ordinal-hex>/<mate>`` with
  the one-based inclusive fragment envelope.
* SAM records place each read at its simulated origin with an indel-aware
  CIGAR. SEQ, QUAL, and zt are reference-forward. Tags: RG, AS, MQ/MC
  (paired), Bismark XG/XR/YS (bisulfite), zt per-base truth, zr read summary,
  and optionally zf fragment summary and zx fragment realization (see
  docs/outputs/index.md). A SAM file gets the text; htsim turns it into BAM.
* The files are written to a hidden staging directory inside the output
  directory. ``RunFiles.publish`` hard-links them into place (never replacing
  an existing file) and writes the manifest last: a dataset is complete only
  once its manifest exists. Leaving the context without publishing removes
  everything that was staged.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from contextlib import suppress
from dataclasses import dataclass
import gzip
import hashlib
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import tempfile
import zlib

import numpy as np

from . import __version__
from .batch import BASE_LETTERS, COMPLEMENT, CTOB, OB, STRAND_NAMES, FragmentBatch
from .errors import BSReadSimError
from .htsim import Header, Summary
from .input import sha256_file
from .reads import G_TO_A, Reads
from .settings import Settings

# The manifest's own format version: it changes only when a released reader could not
# read new manifests (htsim and the package share one version, the package's).
MANIFEST_VERSION = 1


# ---------------------------------------------------------------- FASTQ

def fragment_names(batch: FragmentBatch, contig_names: Sequence[str]) -> list[str]:
    return [
        f"{contig_names[contig]}:{start + 1}-{end}:{ordinal:x}"
        for contig, start, end, ordinal in zip(
            batch.contig.tolist(), batch.start.tolist(), batch.end.tolist(),
            range(batch.first_ordinal, batch.first_ordinal + batch.size),
            strict=True,
        )
    ]


@dataclass(frozen=True)
class FastqFormat:
    """FASTQ records, a file per mate; with a gzip level, each batch is one gzip member."""

    contig_names: tuple[str, ...]
    gzip_level: int | None = None

    def records(self, batch: FragmentBatch, reads: Reads) -> dict[str, bytes]:
        """The records of a batch by file (read1, read2 if paired), in fragment order."""
        names = fragment_names(batch, self.contig_names)
        records = {}
        for mate, (rows, length) in enumerate(reads.mates(), 1):
            sequence = reads.sequence[rows, :length]
            body = np.empty((len(sequence), 2 * length + 5), dtype=np.uint8)
            body[:, 0] = ord("\n")
            body[:, 1:length + 1] = BASE_LETTERS[sequence]
            body[:, length + 1:length + 4] = np.frombuffer(b"\n+\n", dtype=np.uint8)
            body[:, length + 4:2 * length + 4] = reads.quality[rows, :length] + 33
            body[:, -1] = ord("\n")
            width = body.shape[1]
            data = body.tobytes()
            text = b"".join([
                f"@{names[fragment]}/{mate}".encode("ascii")
                + data[row * width:(row + 1) * width]
                for row, fragment in enumerate(reads.fragment[rows].tolist())
            ])
            records[f"read{mate}"] = text if self.gzip_level is None \
                else _gzip_member(text, self.gzip_level)
        return records


def _gzip_member(data: bytes, level: int) -> bytes:
    """One deterministic gzip member (no name, zero mtime)."""
    compressor = zlib.compressobj(level, zlib.DEFLATED, 31)
    return compressor.compress(data) + compressor.flush()


# ---------------------------------------------------------------- SAM

MAPQ = 60
STATE_ALPHABET = b"ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_"  # zt codes
_STATE_LETTERS = np.frombuffer(STATE_ALPHABET, dtype=np.uint8)


@dataclass(frozen=True)
class SamFormat:
    """SAM records of the reads of a run, for the SAM file or for htsim to make the BAM file."""

    role: str  # the output file: bam or sam
    contig_names: tuple[str, ...]
    contig_lengths: tuple[int, ...]
    read_group: str  # the run ID
    sample: str
    paired_end: bool
    bisulfite: bool
    fragment_summary: bool = False
    fragment_realization: bool = False

    @classmethod
    def of_run(cls, header: Header, settings: Settings, run_id: str) -> SamFormat:
        return cls(
            role=settings.format,
            contig_names=tuple(contig.name for contig in header.contigs),
            contig_lengths=tuple(contig.length for contig in header.contigs),
            read_group=run_id, sample=settings.prefix, paired_end=settings.paired_end,
            bisulfite=settings.bisulfite, fragment_summary=settings.fragment_summary,
            fragment_realization=settings.fragment_realization)

    def header(self) -> bytes:
        """The SAM header, with @CO lines that declare the BSReadSim tags."""
        def enabled(flag: bool) -> int:
            return 1 if flag else 0

        lines = ["@HD\tVN:1.6\tSO:unsorted"]
        lines += [f"@SQ\tSN:{name}\tLN:{length}"
                  for name, length in zip(self.contig_names, self.contig_lengths, strict=True)]
        lines += [
            f"@RG\tID:{self.read_group}\tSM:{self.sample}",
            f"@PG\tID:bsreadsim\tPN:bsreadsim\tVN:{__version__}",
            "@CO\tDetails alignments; MAPQ 60 denotes simulated origin, "
            "not calibrated mapping confidence",
            "@CO\tAS_SCHEME=details-max;AS_MAX=query_length",
            f"@CO\tBSREADSIM_ZT=state64;ALPHABET={STATE_ALPHABET.decode()}",
            "@CO\tBSREADSIM_XG=bismark-genome-conversion;ENABLED="
            f"{enabled(self.bisulfite)};VALUES=CT|GA",
            "@CO\tBSREADSIM_XR=bismark-read-conversion;ENABLED="
            f"{enabled(self.bisulfite)};VALUES=CT|GA",
            "@CO\tBSREADSIM_YS=bismark-strand-id;ENABLED="
            f"{enabled(self.bisulfite)};VALUES=OT|OB|CTOT|CTOB",
            "@CO\tBSREADSIM_ZR=u16x12;REQUIRED=1",
            f"@CO\tBSREADSIM_ZF=u16x12;ENABLED={enabled(self.fragment_summary)}",
            "@CO\tBSREADSIM_ZX=packed-b64url;ENABLED="
            f"{enabled(self.fragment_realization)};BIT_ORDER=LSB0",
        ]
        return ("\n".join(lines) + "\n").encode("ascii")

    def records(self, batch: FragmentBatch, reads: Reads) -> dict[str, bytes]:
        """The SAM text of a batch, in read order (mates adjacent, R1 first)."""
        return {self.role: b"".join(self.lines(batch, reads))}

    def lines(self, batch: FragmentBatch, reads: Reads) -> list[bytes]:
        placed = self._place(batch, reads)
        names = [name.encode() for name in fragment_names(batch, self.contig_names)]
        names = [names[i] for i in reads.fragment.tolist()]
        contigs = [self.contig_names[i].encode() for i in batch.contig[reads.fragment].tolist()]
        flags, mate_fields = _mate_fields(reads, placed, self.paired_end)
        tags = self._tags(batch, reads, placed)
        return [
            b"%s\t%d\t%s\t%d\t60\t%s\t%s\t%s\t%s\t%s\n" % values
            for values in zip(names, flags.tolist(), contigs, placed.positions.tolist(),
                              placed.cigars, mate_fields, placed.sequences, placed.qualities,
                              tags, strict=True)
        ]

    def _place(self, batch: FragmentBatch, reads: Reads) -> _Placed:
        """The alignment and forward strings of every read, one mate (and length) at a time."""
        count = len(reads.fragment)
        placed = _Placed(cigars=[b""] * count, positions=np.zeros(count, dtype=np.int64),
                         starts=np.zeros(count, dtype=np.int64),
                         ends=np.zeros(count, dtype=np.int64),
                         sequences=[b""] * count, qualities=[b""] * count, states=[b""] * count)
        for rows, length in reads.mates():
            fragment = reads.fragment[rows]
            index = (batch.template_offsets[fragment]
                     + reads.template_start[rows])[:, None] + np.arange(length)
            placed.cigars[rows], placed.positions[rows], placed.starts[rows], \
                placed.ends[rows] = _alignments(batch, fragment, reads.template_start[rows],
                                                index, self.contig_lengths)
            flip = reads.reverse[rows]
            sequence = reads.sequence[rows, :length]
            quality = reads.quality[rows, :length]
            forward_sequence = sequence.copy()
            forward_quality = quality.copy()
            forward_sequence[flip] = COMPLEMENT[sequence[flip, ::-1]]
            forward_quality[flip] = quality[flip, ::-1]
            placed.sequences[rows] = _rows(BASE_LETTERS[forward_sequence])
            placed.qualities[rows] = _rows(forward_quality + 33)
            placed.states[rows] = _rows(_STATE_LETTERS[reads.truth.states[rows, :length]])
        return placed

    def _tags(self, batch: FragmentBatch, reads: Reads, placed: _Placed) -> list[bytes]:
        """RG, AS, MQ/MC (paired), XG/XR/YS (bisulfite), zt, zr, and zf/zx when asked for."""
        fragment = reads.fragment
        labels = [b"RG:Z:" + self.read_group.encode() + b"\tAS:i:%d" % length
                  for length in reads.lengths]
        tags = [labels[mate] for mate in (reads.mate - 1).tolist()]
        if self.paired_end:
            mate = np.arange(len(tags)) ^ 1
            tags = [
                tag + b"\tMQ:i:%d\tMC:Z:%s" % (MAPQ, placed.cigars[other])
                for tag, other in zip(tags, mate.tolist(), strict=True)
            ]
        if self.bisulfite:
            bottom = np.isin(batch.strand, (OB, CTOB))[fragment]
            genome = np.where(bottom, b"GA", b"CT")
            read = np.where(bottom ^ reads.reverse == G_TO_A, b"GA", b"CT")
            strand = np.array([name.encode() for name in STRAND_NAMES])[batch.strand[fragment]]
            tags = [
                tag + b"\tXG:Z:%s\tXR:Z:%s\tYS:Z:%s" % (g, r, s)
                for tag, g, r, s in zip(tags, genome.tolist(), read.tolist(), strand.tolist(),
                                        strict=True)
            ]
        tags = [
            tag + b"\tzt:Z:" + state + b"\tzr:B:S," + summary
            for tag, state, summary in zip(
                tags, placed.states, _summaries(reads.truth.read_summary), strict=True)
        ]
        if self.fragment_summary:
            fragment_summaries = _summaries(reads.truth.fragment_summary)
            tags = [tag + b"\tzf:B:S," + fragment_summaries[i]
                    for tag, i in zip(tags, fragment.tolist(), strict=True)]
        if self.fragment_realization:
            realizations = [value.encode() for value in reads.truth.realization]
            tags = [tag + b"\tzx:Z:" + realizations[i]
                    for tag, i in zip(tags, fragment.tolist(), strict=True)]
        return tags


@dataclass
class _Placed:
    """Where each read of a batch aligns, and its reference-forward SEQ, QUAL, and zt."""

    cigars: list[bytes]
    positions: np.ndarray  # int64[R], 1-based POS
    starts: np.ndarray  # int64[R], the aligned reference interval [start, end)
    ends: np.ndarray  # int64[R]
    sequences: list[bytes]
    qualities: list[bytes]
    states: list[bytes]


def _mate_fields(reads: Reads, placed: _Placed, paired: bool) -> tuple[np.ndarray, list[bytes]]:
    """FLAG, and RNEXT/PNEXT/TLEN, of every read."""
    reverse = reads.reverse
    if not paired:
        return np.where(reverse, 0x10, 0), [b"*\t0\t0"] * len(reverse)
    mate = np.arange(len(reverse)) ^ 1  # mates are adjacent: R1 then R2
    flags = 0x1 | 0x2 | np.where(reads.mate == 1, 0x40, 0x80) \
        | np.where(reverse, 0x10, 0) | np.where(reverse[mate], 0x20, 0)
    starts, ends = placed.starts, placed.ends
    span = np.maximum(ends, ends[mate]) - np.minimum(starts, starts[mate])
    leftmost = (starts < starts[mate]) | ((starts == starts[mate]) & (reads.mate == 1))
    fields = [
        b"=\t%d\t%d" % values for values in zip(
            placed.positions[mate].tolist(), np.where(leftmost, span, -span).tolist(),
            strict=True)
    ]
    return flags, fields


def _rows(characters: np.ndarray) -> list[bytes]:
    """Split a 2-D uint8 array into one bytes object per row."""
    data = np.ascontiguousarray(characters, dtype=np.uint8).tobytes()
    width = characters.shape[1]
    return [data[start:start + width] for start in range(0, len(data), width)]


def _summaries(values: np.ndarray) -> list[bytes]:
    return [",".join(map(str, row)).encode() for row in values.tolist()]


def _alignments(batch: FragmentBatch, fragments: np.ndarray, template_starts: np.ndarray,
                index: np.ndarray, contig_lengths: Sequence[int]
                ) -> tuple[list[bytes], np.ndarray, np.ndarray, np.ndarray]:
    """CIGARs (bytes), 1-based POS, and reference [start, end) of reads of one length."""
    positions = batch.reference_positions[index]
    length = index.shape[1]
    starts = positions[:, 0].copy()
    ends = starts + length
    cigars = [b"%dM" % length] * len(starts)
    simple = (positions >= 0).all(axis=1) & (
        positions[:, -1] - positions[:, 0] == length - 1)
    for row in np.flatnonzero(~simple).tolist():
        fragment = fragments[row]
        template = batch.reference_positions[
            batch.template_offsets[fragment]:batch.template_offsets[fragment + 1]]
        cigar, _position, starts[row], ends[row] = _indel_alignment(
            positions[row], template, int(template_starts[row]),
            contig_lengths[batch.contig[fragment]])
        cigars[row] = cigar.encode()
    return cigars, starts + 1, starts, ends


def _indel_alignment(positions: np.ndarray, template: np.ndarray, template_start: int,
                     contig_length: int) -> tuple[str, int, int, int]:
    """The CIGAR, POS, and reference [start, end) of a read with an indel or a gap."""
    operations = []

    def add(operation: str, count: int = 1) -> None:
        if operations and operations[-1][0] == operation:
            operations[-1][1] += count
        else:
            operations.append([operation, count])

    previous = None
    for position in positions.tolist():
        if position < 0:
            add("I")
            continue
        if previous is not None and position > previous + 1:
            add("D", position - previous - 1)
        add("M")
        previous = position
    cigar = "".join(f"{count}{operation}" for operation, count in operations)
    mapped = positions[positions >= 0]
    if mapped.size:
        start, end = int(mapped[0]), int(mapped[-1]) + 1
        return cigar, start + 1, start, end
    # The read lies entirely inside an insertion: anchor it at the next base.
    following = template[template_start + len(positions):]
    following = following[following >= 0]
    if following.size:
        anchor = int(following[0])
    else:
        preceding = template[:template_start]
        anchor = int(preceding[preceding >= 0][-1]) + 1
    start = min(anchor, contig_length - 1)
    return cigar, start + 1, start, start + 1


# ---------------------------------------------------------------- files

@dataclass(frozen=True)
class OutputFile:
    role: str  # read1, read2, bam, sam, truth.methdb, truth.vcf
    path: Path  # published location
    size_bytes: int
    sha256: str
    record_count: int | None  # None when not counted


class BamWriter:
    """Feeds SAM text to a command (htsim's SAM to BAM) that writes BAM to ``path``."""

    def __init__(self, path: Path, argv: Sequence[str], sam_header: bytes) -> None:
        self.path = path
        with path.open("xb") as output:
            self._process = subprocess.Popen(
                list(argv), stdin=subprocess.PIPE, stdout=output, stderr=subprocess.PIPE)
        self.write(sam_header)

    def write(self, data: bytes) -> None:
        try:
            self._process.stdin.write(data)
        except OSError as error:
            self._process.kill()
            raise BSReadSimError(f"the BAM writer stopped early: {self._error()}") from error

    def close(self) -> None:
        try:
            self._process.stdin.close()
        except OSError:
            pass
        if self._process.wait() != 0:
            raise BSReadSimError(f"the BAM writer failed: {self._error()}")

    def kill(self) -> None:
        self._process.kill()
        self._process.wait()

    def _error(self) -> str:
        self._process.wait()
        return self._process.stderr.read().decode("utf-8", errors="replace").strip()


class RunFiles:
    """The files of one run: reads (FASTQ, BAM, or SAM) and truth files, staged until
    published."""

    def __init__(self, directory: Path, prefix: str, *, paired_end: bool,
                 format: str) -> None:
        self.directory = directory
        self.prefix = prefix
        self.fragments = 0
        self._mates = 2 if paired_end else 1
        self._counts: dict[str, int] = {}
        self._digests: dict = {}
        self._streams: dict = {}
        if format in ("bam", "sam"):
            names = {format: f"{prefix}.{format}"}
        else:
            suffix = ".gz" if format == "fastq.gz" else ""
            names = {"read1": f"{prefix}.R1.fastq{suffix}"}
            if paired_end:
                names["read2"] = f"{prefix}.R2.fastq{suffix}"
        self._names = dict(names)
        self.manifest_path = directory / (prefix + ".manifest.json")
        directory.mkdir(parents=True, exist_ok=True)
        self._require_absent(list(names.values()))
        self.staging = Path(tempfile.mkdtemp(prefix=f".{prefix}.staging-",
                                             dir=directory))
        try:
            for role, name in names.items():
                if role != "bam":
                    self._streams[role] = (self.staging / name).open("xb")
                    self._digests[role] = hashlib.sha256()
        except BaseException:
            self.abort()
            raise

    def __enter__(self) -> "RunFiles":
        return self

    def __exit__(self, exc_type, exc_value, traceback) -> bool:
        if self.staging.exists():
            self.abort()
        return False

    def start_alignments(self, sam_header: bytes, sam_to_bam: Sequence[str]) -> None:
        """Begin the SAM or BAM file with its header; ``sam_to_bam`` is the command that
        makes BAM of the SAM text."""
        if "bam" in self._names:
            self._streams["bam"] = BamWriter(self.staging / self._names["bam"], sam_to_bam,
                                             sam_header)
        else:
            self.write(0, {"sam": sam_header})

    def truth_path(self, name: str) -> Path:
        """Where a truth file named ``name`` is staged."""
        self._require_absent([str(Path("truth") / name)])
        path = self.staging / "truth" / name
        path.parent.mkdir(exist_ok=True)
        return path

    def write(self, fragments: int, records: Mapping[str, bytes]) -> None:
        """Append the records of ``fragments`` consecutive fragments, by file role."""
        for role, data in records.items():
            self._streams[role].write(data)
            if role in self._digests:
                self._digests[role].update(data)
        self.fragments += fragments

    def add_truth(self, role: str, name: str) -> None:
        """Include the truth file staged at ``truth_path(name)`` (truth.methdb or truth.vcf)."""
        self._names[role] = str(Path("truth") / name)

    def finish(self) -> list[OutputFile]:
        """Close every file and identify it; nothing is published yet."""
        for role, stream in self._streams.items():
            stream.close()
            self._counts[role] = self.fragments * (self._mates if role in ("bam", "sam") else 1)
        self._streams = {}
        files = []
        for role, name in self._names.items():
            staged = self.staging / name
            if role == "truth.vcf":
                self._counts[role] = _vcf_record_count(staged)
            if role in self._digests:
                sha256 = self._digests[role].hexdigest()
            else:
                sha256 = sha256_file(staged)
            files.append(OutputFile(role, self.directory / name, staged.stat().st_size,
                                    sha256, self._counts.get(role)))
        return files

    def publish(self, manifest: dict) -> None:
        """Link the staged files into place, then write the manifest."""
        staged_manifest = self.staging / self.manifest_path.name
        staged_manifest.write_text(
            json.dumps(manifest, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
            encoding="utf-8")
        published = []
        try:
            for name in [*self._names.values(), self.manifest_path.name]:
                final = self.directory / name
                final.parent.mkdir(exist_ok=True)
                os.link(self.staging / name, final)
                published.append(final)
        except OSError as error:
            for path in published:
                with suppress(OSError):
                    path.unlink()
            self.abort()
            raise BSReadSimError(f"cannot publish the outputs: {error}") from error
        shutil.rmtree(self.staging)

    def abort(self) -> None:
        for stream in self._streams.values():
            with suppress(Exception):
                stream.kill() if isinstance(stream, BamWriter) else stream.close()
        self._streams = {}
        shutil.rmtree(self.staging, ignore_errors=True)

    def _require_absent(self, names: list[str]) -> None:
        existing = [str(self.directory / name) for name in [*names, self.manifest_path.name]
                    if os.path.lexists(self.directory / name)]
        if existing:
            raise BSReadSimError(f"output destination already exists: {', '.join(existing)}")


def _vcf_record_count(path: Path) -> int:
    with gzip.open(path, "rb") as source:
        return sum(1 for line in source if line.strip() and not line.startswith(b"#"))


# ---------------------------------------------------------------- manifest

def build_manifest(settings: Settings, argv: Sequence[str], header: Header, summary: Summary,
                   inputs: list[dict], outputs: list[OutputFile], *, run_id: str,
                   methylation_model: str) -> dict:
    """``inputs`` holds role, format, path, size_bytes, and sha256 of each input;
    ``methylation_model`` names the model the reads were simulated with."""
    return {
        "version": MANIFEST_VERSION,
        "status": "complete",
        "run_id": run_id,
        "summary": {
            "technology": settings.technology,
            "output_format": settings.format,
            "paired_end": settings.paired_end,
            "fragment_count": summary.fragment_count,
            "skipped_fragment_count": summary.skipped_fragment_count,
            "read_count": summary.mate_count,
            "read_base_count": summary.fragment_count * sum(
                settings.read_lengths[:2 if settings.paired_end else 1]),
            "template_base_count": summary.template_base_count,
            "methylation_site_count": summary.methylation_site_count,
            "output_file_count": len(outputs),
            "output_size_bytes": sum(item.size_bytes for item in outputs),
        },
        "inputs": inputs,
        "outputs": [
            {"role": item.role, "path": str(item.path), "record_count": item.record_count,
             "size_bytes": item.size_bytes, "sha256": item.sha256}
            for item in outputs
        ],
        "command": {
            "interface": "cli",
            "user_command": shlex.join(argv),
            "full_command": full_command(settings, argv),
        },
        "details": {
            "configuration": settings.as_dict(),
            "configuration_sha256": settings.sha256,
            "randomness": {
                "master_seed": str(settings.seed),
                "mutation_seed": str(settings.seed_mut),
                "phasing_seed": str(settings.seed_phase),
                "methylation_seed": str(settings.seed_meth),
            },
            "models": {"methylation_state": (
                {"requested": settings.meth_model, "effective": methylation_model}
                if settings.bisulfite else {"requested": "disabled", "effective": "disabled"})},
            "contigs": [
                {"index": index, "name": contig.name, "length": contig.length,
                 "md5": contig.md5,
                 "fragment_count": summary.per_contig_fragment_counts[index]}
                for index, contig in enumerate(header.contigs)
            ],
            "software_versions": {"core": header.core_version, "python": __version__},
        },
    }


def full_command(settings: Settings, argv: Sequence[str]) -> str:
    """The received command with every omitted seed made explicit."""
    full = list(argv)
    seeds = [("--seed", settings.seed), ("--seed-mut", settings.seed_mut),
             ("--seed-phase", settings.seed_phase)]
    if settings.bisulfite:
        seeds.append(("--seed-meth", settings.seed_meth))
    for option, value in seeds:
        if not any(arg == option or arg.startswith(option + "=") for arg in argv):
            full += [option, str(value)]
    return shlex.join(full)
