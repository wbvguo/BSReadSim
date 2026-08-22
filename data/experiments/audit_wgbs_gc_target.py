#!/usr/bin/env python3
"""Audit one completed WGBS target run by replaying its core stream.

The FASTQ sequence is not used for GC auditing because bisulfite conversion
changes C to T.  Instead, the script reconstructs the manifest's exact
normalized configuration, verifies all current input bytes, and replays the
same deterministic protocol stream with the same run ID.  Matching stream
SHA-256 proves that the audited fragment coordinates are those used for the
completed FASTQ run.
"""

from __future__ import annotations

import argparse
import copy
import csv
from dataclasses import dataclass
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
import time
from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bsreadsim.run.config import normalize_run_config  # noqa: E402
from bsreadsim.native.launch import build_core_argv  # noqa: E402
from bsreadsim.native.subprocess import CoreProcess  # noqa: E402
from bsreadsim.run.manifest import (  # noqa: E402
    validate_header_projection,
    verify_complete_manifest,
)
from bsreadsim.run.prepare import PreparedRun, prepare_run  # noqa: E402
from bsreadsim.native.launch import resolve_core_executable  # noqa: E402


_NUMBER = re.compile(
    rb"^-?(?:(?:[0-9]+(?:\.[0-9]*)?)|(?:\.[0-9]+))(?:[eE][+-]?[0-9]+)?$"
)


class AuditError(RuntimeError):
    """The run cannot support a reference-only GC audit."""


@dataclass(frozen=True)
class CandidateSpace:
    contig_name: str
    contig_length: int
    calibration_length: int
    gc_prefix: np.ndarray
    bin_by_start: np.ndarray
    candidate_counts: np.ndarray

    @property
    def valid_start_count(self) -> int:
        return int(self.candidate_counts.sum())


def _mapping(value: object, name: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise AuditError("{} must be an object".format(name))
    return value


def _load_manifest(path: Path) -> dict[str, Any]:
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError) as error:
        raise AuditError("cannot read manifest {}: {}".format(path, error)) from error
    verify_complete_manifest(document)
    return document


def _prepare_manifest_run(manifest: Mapping[str, Any]) -> PreparedRun:
    config_section = _mapping(manifest.get("config"), "manifest.config")
    normalized = copy.deepcopy(
        dict(_mapping(config_section.get("normalized"), "manifest.config.normalized"))
    )
    output = normalized.get("output")
    if not isinstance(output, dict) or output.get("truth") != "none":
        raise AuditError("the audit requires a completed production run")
    output.pop("truth")
    loaded = normalize_run_config(normalized, Path("/"))
    if loaded.sha256 != config_section.get("sha256"):
        raise AuditError("reconstructed normalized config SHA-256 changed")
    prepared = prepare_run(loaded)

    manifested_inputs = {
        item["role"]: item
        for item in manifest.get("inputs", [])
        if isinstance(item, Mapping) and isinstance(item.get("role"), str)
    }
    if len(manifested_inputs) != len(prepared.files):
        raise AuditError("manifest input roles disagree with prepared inputs")
    for identity in prepared.files:
        item = manifested_inputs.get(identity.role)
        if item is None:
            raise AuditError("manifest omitted input role {}".format(identity.role))
        if (
            item.get("path") != str(identity.path)
            or item.get("size_bytes") != identity.size_bytes
            or item.get("sha256") != identity.sha256
        ):
            raise AuditError(
                "current input bytes disagree with manifest role {}".format(
                    identity.role
                )
            )
    return prepared


def _decoded_bytes(path: Path) -> tuple[bytes, bytes]:
    try:
        raw = path.read_bytes()
        decoded = gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw
    except (OSError, EOFError, gzip.BadGzipFile) as error:
        raise AuditError("cannot decode {}: {}".format(path, error)) from error
    return raw, decoded


def _load_profile(path: Path, expected_sha256: str) -> np.ndarray:
    raw, decoded = _decoded_bytes(path)
    if hashlib.sha256(raw).hexdigest() != expected_sha256:
        raise AuditError("coverage profile SHA-256 changed")
    lines = decoded.splitlines()
    if len(lines) < 2:
        raise AuditError("coverage profile requires at least two probabilities")
    probabilities = np.empty(len(lines), dtype=np.float64)
    for index, line in enumerate(lines):
        if not line or _NUMBER.fullmatch(line) is None:
            raise AuditError(
                "coverage profile line {} is not one strict probability".format(
                    index + 1
                )
            )
        value = float(line)
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise AuditError(
                "coverage profile line {} is outside [0,1]".format(index + 1)
            )
        probabilities[index] = value
    if not math.isclose(
        float(math.fsum(probabilities)), 1.0, rel_tol=0.0, abs_tol=1e-9
    ):
        raise AuditError("target profile probabilities must sum to one")
    return probabilities


def _load_single_fasta(path: Path) -> tuple[str, bytes]:
    _, decoded = _decoded_bytes(path)
    name = None
    sequence_parts = []  # type: list[bytes]
    for line_number, line in enumerate(decoded.splitlines(), 1):
        if line.startswith(b">"):
            if name is not None:
                raise AuditError("GC audit requires exactly one FASTA contig")
            fields = line[1:].split()
            if not fields:
                raise AuditError("FASTA header is empty")
            try:
                name = fields[0].decode("ascii")
            except UnicodeDecodeError as error:
                raise AuditError("FASTA contig name is not ASCII") from error
            continue
        if name is None:
            raise AuditError("FASTA sequence precedes its header")
        upper = line.upper()
        if not upper or any(base not in b"ACGTN" for base in upper):
            raise AuditError(
                "FASTA line {} contains an empty or invalid sequence".format(
                    line_number
                )
            )
        sequence_parts.append(upper)
    if name is None or not sequence_parts:
        raise AuditError("FASTA contains no sequence")
    return name, b"".join(sequence_parts)


def _bin_for_gc_counts(
    gc_counts: np.ndarray,
    fragment_length: Any,
    bin_count: int,
) -> np.ndarray:
    product = gc_counts.astype(np.uint64, copy=False) * np.uint64(bin_count - 1)
    result = product // np.uint64(fragment_length)
    remainder = product % np.uint64(fragment_length)
    result += (np.uint64(2) * remainder >= fragment_length)
    return result


def _candidate_space(
    reference_path: Path,
    fragments: Mapping[str, Any],
    bin_count: int,
    *,
    chunk_starts: int,
) -> CandidateSpace:
    if fragments.get("paired_end") is not True:
        raise AuditError("the current GC audit requires paired-end fragments")
    read_length_1 = int(fragments["read_length_1"])
    read_length_2 = int(fragments["read_length_2"])
    if read_length_1 != read_length_2:
        raise AuditError("the current GC audit requires equal mate lengths")
    insert_min = int(fragments["insert_min"])
    insert_mean = int(fragments["insert_mean"])
    insert_max = int(fragments["insert_max"])
    if not insert_min <= insert_mean <= insert_max:
        raise AuditError("insert length parameters are inconsistent")
    calibration_length = insert_mean

    contig_name, sequence = _load_single_fasta(reference_path)
    contig_length = len(sequence)
    if calibration_length > contig_length:
        raise AuditError("mean insert length exceeds the reference")
    encoded = np.frombuffer(sequence, dtype=np.uint8)
    gc_mask = (encoded == ord("G")) | (encoded == ord("C"))
    n_mask = encoded == ord("N")
    gc_prefix = np.empty(contig_length + 1, dtype=np.uint32)
    n_prefix = np.empty(contig_length + 1, dtype=np.uint32)
    gc_prefix[0] = 0
    n_prefix[0] = 0
    np.cumsum(gc_mask, dtype=np.uint32, out=gc_prefix[1:])
    np.cumsum(n_mask, dtype=np.uint32, out=n_prefix[1:])

    possible = contig_length - calibration_length + 1
    bin_dtype = (
        np.uint8
        if bin_count <= 256
        else np.uint16
        if bin_count <= 65536
        else np.uint32
    )
    bin_by_start = np.empty(possible, dtype=bin_dtype)
    candidate_counts = np.zeros(bin_count, dtype=np.uint64)
    maximum_n = math.floor(
        float(fragments["max_ambiguous_fraction"]) * read_length_1
    )
    second_offset = calibration_length - read_length_2

    for begin in range(0, possible, chunk_starts):
        end = min(begin + chunk_starts, possible)
        starts = np.arange(begin, end, dtype=np.int64)
        first_n = n_prefix[starts + read_length_1] - n_prefix[starts]
        second_starts = starts + second_offset
        second_n = (
            n_prefix[second_starts + read_length_2] - n_prefix[second_starts]
        )
        valid = (first_n <= maximum_n) & (second_n <= maximum_n)
        gc_counts = gc_prefix[starts + calibration_length] - gc_prefix[starts]
        bins = _bin_for_gc_counts(gc_counts, calibration_length, bin_count)
        bin_by_start[begin:end] = bins.astype(bin_dtype, copy=False)
        candidate_counts += np.bincount(
            bins[valid].astype(np.int64, copy=False), minlength=bin_count
        ).astype(np.uint64, copy=False)

    if int(candidate_counts.sum()) == 0:
        raise AuditError("reference has no eligible mean-insert starts")
    return CandidateSpace(
        contig_name=contig_name,
        contig_length=contig_length,
        calibration_length=calibration_length,
        gc_prefix=gc_prefix,
        bin_by_start=bin_by_start,
        candidate_counts=candidate_counts,
    )


def _replay_core_histogram(
    manifest: Mapping[str, Any],
    prepared: PreparedRun,
    core: Path,
    candidates: CandidateSpace,
    bin_count: int,
) -> tuple[np.ndarray, np.ndarray, Mapping[str, Any]]:
    config = prepared.config.normalized
    fragments = _mapping(config["fragments"], "config.fragments")
    execution = _mapping(config["execution"], "config.execution")
    insert_min = int(fragments["insert_min"])
    insert_max = int(fragments["insert_max"])
    protocol_batch_fragments = min(
        64, int(execution["max_in_flight_fragments"])
    )
    argv = build_core_argv(
        prepared,
        str(manifest["run_id"]),
        core,
        truth_columns="none",
        protocol_batch_fragments=protocol_batch_fragments,
    )

    observed = np.zeros(bin_count, dtype=np.uint64)
    observed_inserts = np.zeros(insert_max - insert_min + 1, dtype=np.uint64)
    start_buffer = np.empty(1024 * 1024, dtype=np.uint32)
    length_buffer = np.empty(len(start_buffer), dtype=np.uint32)
    buffered = 0

    def flush() -> None:
        nonlocal buffered
        if buffered == 0:
            return
        starts = start_buffer[:buffered].astype(np.uint64, copy=False)
        lengths = length_buffer[:buffered].astype(np.uint64, copy=False)
        ends = starts + lengths
        gc_counts = candidates.gc_prefix[ends] - candidates.gc_prefix[starts]
        bins = _bin_for_gc_counts(gc_counts, lengths, bin_count)
        observed[:] += np.bincount(
            bins.astype(np.int64, copy=False), minlength=bin_count
        ).astype(np.uint64, copy=False)
        observed_inserts[:] += np.bincount(
            (lengths - insert_min).astype(np.int64, copy=False),
            minlength=len(observed_inserts),
        ).astype(np.uint64, copy=False)
        buffered = 0

    process = CoreProcess(
        argv,
        read_length=int(fragments["read_length_1"]),
        paired_end=True,
        protocol_idle_timeout_seconds=300.0,
        exit_timeout_seconds=30.0,
    )
    with process as running:
        validate_header_projection(prepared, running.header)
        if len(running.header.contigs) != 1:
            raise AuditError("replayed core stream does not contain one contig")
        contig = running.header.contigs[0]
        if (
            contig.name != candidates.contig_name
            or contig.length != candidates.contig_length
        ):
            raise AuditError("replayed core contig disagrees with FASTA audit")

        for batch in running.iter_batches():
            starts = np.frombuffer(batch.reference_begins.raw, dtype="<u4")
            ends = np.frombuffer(batch.reference_ends.raw, dtype="<u4")
            contigs = np.frombuffer(batch.contig_indices.raw, dtype="<u4")
            offsets = np.frombuffer(batch.template_offsets.raw, dtype="<u4")
            reference_lengths = (
                ends.astype(np.uint64) - starts.astype(np.uint64)
            )
            template_lengths = np.diff(offsets.astype(np.uint64))
            if (
                np.any(contigs != 0)
                or np.any(ends < starts)
                or np.any(reference_lengths < insert_min)
                or np.any(reference_lengths > insert_max)
                or np.any(template_lengths != reference_lengths)
                or np.any(ends > candidates.contig_length)
            ):
                raise AuditError("replayed fragment batch violates insert audit")
            copied = 0
            while copied < len(starts):
                available = len(start_buffer) - buffered
                take = min(available, len(starts) - copied)
                start_buffer[buffered : buffered + take] = starts[copied : copied + take]
                length_buffer[buffered : buffered + take] = reference_lengths[
                    copied : copied + take
                ]
                buffered += take
                copied += take
                if buffered == len(start_buffer):
                    flush()
    flush()

    trailer = process.trailer
    manifested_counts = _mapping(
        _mapping(manifest["counts"], "manifest.counts")["core"],
        "manifest.counts.core",
    )
    observed_stream = trailer.stream_sha256.hex()
    if observed_stream != manifest.get("stream_sha256"):
        raise AuditError("replayed protocol stream SHA-256 differs from FASTQ run")
    if (
        trailer.fragment_count != manifested_counts.get("fragment_count")
        or trailer.skipped_fragment_count
        != manifested_counts.get("skipped_fragment_count")
        or int(observed.sum()) != trailer.fragment_count
    ):
        raise AuditError("replayed fragment accounting differs from manifest")
    return observed, observed_inserts, {
        "stream_sha256": observed_stream,
        "fragment_count": trailer.fragment_count,
        "skipped_fragment_count": trailer.skipped_fragment_count,
        "protocol_batch_fragments": protocol_batch_fragments,
    }


def _overlap(first: np.ndarray, second: np.ndarray) -> float:
    return float(np.minimum(first, second).sum())


def _distribution_rows(
    probabilities: np.ndarray,
    candidates: np.ndarray,
    observed: np.ndarray,
    uniform_observed: np.ndarray | None = None,
    *,
    project_unreachable: bool = False,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    candidate_total = int(candidates.sum())
    observed_total = int(observed.sum())
    candidate_distribution = candidates.astype(np.float64) / candidate_total
    target = probabilities / float(math.fsum(probabilities))
    unreachable = (target > 0.0) & (candidates == 0)
    if np.any(unreachable) and not project_unreachable:
        raise AuditError(
            "positive target bin {} has no eligible fragment start".format(
                int(np.flatnonzero(unreachable)[0])
            )
        )
    projected_target = target.copy()
    projected_target[unreachable] = 0.0
    projected_total = float(math.fsum(projected_target))
    if projected_total <= 0.0:
        raise AuditError("target has no reachable positive probability")
    projected_target /= projected_total
    ratios = np.divide(
        projected_target,
        candidates,
        out=np.zeros_like(target),
        where=candidates > 0,
    )
    calibrated_acceptance = ratios / float(ratios.max())
    observed_distribution = observed.astype(np.float64) / observed_total
    if uniform_observed is not None:
        uniform_total = int(uniform_observed.sum())
        if uniform_total == 0:
            raise AuditError("uniform baseline contains no fragments")
        uniform_distribution = (
            uniform_observed.astype(np.float64) / uniform_total
        )
    else:
        uniform_distribution = None
    overall_coverage = observed_total / candidate_total

    rows = []  # type: list[dict[str, Any]]
    last_bin = len(probabilities) - 1
    for index in range(len(probabilities)):
        candidate_count = int(candidates[index])
        observed_count = int(observed[index])
        bin_coverage = (
            observed_count / candidate_count if candidate_count else None
        )
        row = {
            "bin": index,
            "gc_ratio_center": index / last_bin,
            "input_target_probability": float(probabilities[index]),
            "normalized_target_probability": float(target[index]),
            "eligible_start_count": candidate_count,
            "eligible_start_probability": float(candidate_distribution[index]),
            "calibrated_acceptance_probability": float(
                calibrated_acceptance[index]
            ),
            "expected_fragment_probability": float(projected_target[index]),
            "observed_fragment_count": observed_count,
            "observed_fragment_probability": float(
                observed_distribution[index]
            ),
            "mean_fragments_per_eligible_start": bin_coverage,
            "relative_bin_coverage": (
                None if bin_coverage is None else bin_coverage / overall_coverage
            ),
            "observed_minus_expected": float(
                observed_distribution[index] - projected_target[index]
            ),
        }
        if uniform_distribution is not None:
            row["uniform_fragment_count"] = int(uniform_observed[index])
            row["uniform_fragment_probability"] = float(
                uniform_distribution[index]
            )
        rows.append(row)

    target_overlap = _overlap(target, observed_distribution)
    metrics = {
        "target_vs_observed_overlap": target_overlap,
        "target_vs_observed_total_variation": 1.0 - target_overlap,
        "candidate_vs_observed_overlap": _overlap(
            candidate_distribution, observed_distribution
        ),
        "target_vs_observed_max_abs_error": float(
            np.max(np.abs(target - observed_distribution))
        ),
        "expected_proposal_acceptance": float(
            np.dot(candidate_distribution, calibrated_acceptance)
        ),
    }
    if uniform_distribution is not None:
        metrics.update(
            {
                "target_vs_uniform_overlap": _overlap(
                    target, uniform_distribution
                ),
                "profile_vs_uniform_overlap": _overlap(
                    observed_distribution, uniform_distribution
                ),
            }
        )
    if project_unreachable:
        projected_overlap = _overlap(projected_target, observed_distribution)
        metrics.update(
            {
                "dropped_target_probability": float(target[unreachable].sum()),
                "projected_target_vs_observed_overlap": projected_overlap,
                "projected_target_vs_observed_total_variation": (
                    1.0 - projected_overlap
                ),
                "projected_target_vs_observed_max_abs_error": float(
                    np.max(np.abs(projected_target - observed_distribution))
                ),
            }
        )
    return rows, metrics


def _insert_metrics(
    insert_min: int,
    profile_counts: np.ndarray,
    uniform_counts: np.ndarray,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    profile_total = int(profile_counts.sum())
    uniform_total = int(uniform_counts.sum())
    if profile_total == 0 or uniform_total == 0:
        raise AuditError("insert audit requires non-empty profile and uniform runs")
    profile = profile_counts.astype(np.float64) / profile_total
    uniform = uniform_counts.astype(np.float64) / uniform_total
    lengths = np.arange(
        insert_min,
        insert_min + len(profile_counts),
        dtype=np.float64,
    )
    if len(profile) != len(uniform):
        raise AuditError("profile and uniform insert domains disagree")
    profile_mean = float(np.dot(lengths, profile))
    uniform_mean = float(np.dot(lengths, uniform))
    profile_stddev = float(
        math.sqrt(np.dot((lengths - profile_mean) ** 2, profile))
    )
    uniform_stddev = float(
        math.sqrt(np.dot((lengths - uniform_mean) ** 2, uniform))
    )
    rows = [
        {
            "insert_length": int(length),
            "profile_count": int(profile_counts[index]),
            "profile_probability": float(profile[index]),
            "uniform_count": int(uniform_counts[index]),
            "uniform_probability": float(uniform[index]),
        }
        for index, length in enumerate(lengths)
        if profile_counts[index] != 0 or uniform_counts[index] != 0
    ]
    overlap = _overlap(profile, uniform)
    return rows, {
        "profile_vs_uniform_insert_overlap": overlap,
        "profile_vs_uniform_insert_total_variation": 1.0 - overlap,
        "profile_insert_mean": profile_mean,
        "uniform_insert_mean": uniform_mean,
        "profile_minus_uniform_insert_mean": profile_mean - uniform_mean,
        "profile_insert_stddev": profile_stddev,
        "uniform_insert_stddev": uniform_stddev,
        "profile_minus_uniform_insert_stddev": (
            profile_stddev - uniform_stddev
        ),
    }


def _write_tsv(
    path: Path,
    rows: Sequence[Mapping[str, Any]],
) -> None:
    fields = tuple(rows[0])
    with path.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(
            output,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    field: (
                        "NA"
                        if row[field] is None
                        else "{:.12g}".format(row[field])
                        if isinstance(row[field], float)
                        else row[field]
                    )
                    for field in fields
                }
            )


def _write_plot(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as error:
        raise AuditError("--plot requires matplotlib") from error

    x = [100.0 * float(row["gc_ratio_center"]) for row in rows]
    figure, axis = plt.subplots(figsize=(10, 6))
    axis.plot(
        x,
        [row["normalized_target_probability"] for row in rows],
        label="input target distribution",
        linewidth=2,
    )
    axis.plot(
        x,
        [row["eligible_start_probability"] for row in rows],
        label="uniform opportunity distribution",
        linewidth=2,
    )
    axis.plot(
        x,
        [row["observed_fragment_probability"] for row in rows],
        label="observed target-profile output",
        linewidth=1.5,
        linestyle="--",
    )
    if "uniform_fragment_probability" in rows[0]:
        axis.plot(
            x,
            [row["uniform_fragment_probability"] for row in rows],
            label="observed variable-insert uniform output",
            linewidth=1.25,
            linestyle=":",
        )
    axis.set_xlabel("fragment GC ratio (%)")
    axis.set_ylabel("probability mass per bin")
    axis.set_title("WGBS target GC-distribution audit")
    axis.grid(alpha=0.2)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=160)
    plt.close(figure)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument(
        "--baseline-manifest",
        type=Path,
        help="matching uniform run used to quantify insert-distribution drift",
    )
    parser.add_argument("--core", type=Path)
    parser.add_argument("--output-tsv", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--candidate-chunk-starts", type=int, default=1000000)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    arguments = build_parser().parse_args(argv)
    if arguments.candidate_chunk_starts <= 0:
        raise AuditError("--candidate-chunk-starts must be positive")
    for output in (arguments.output_tsv, arguments.output_json, arguments.plot):
        if output is not None and output.exists():
            raise AuditError("refusing to overwrite {}".format(output))

    started = time.perf_counter()
    manifest_path = arguments.manifest.expanduser().resolve(strict=True)
    manifest = _load_manifest(manifest_path)
    prepared = _prepare_manifest_run(manifest)
    config = prepared.config.normalized
    if config.get("technology") != "WGBS":
        raise AuditError("manifest is not a WGBS run")
    mutation = _mapping(config["mutation"], "config.mutation")
    inputs = _mapping(config["inputs"], "config.inputs")
    coverage = _mapping(config["coverage"], "config.coverage")
    if mutation.get("rate") != 0 or "vcf" in inputs:
        raise AuditError("exact reference-start audit requires no variants")
    if coverage.get("kind") != "profile":
        raise AuditError("manifest did not use profile coverage")
    artifact = _mapping(coverage["artifact"], "config.coverage.artifact")
    if artifact.get("version") != "wgbs-gc-target-v1":
        raise AuditError("manifest did not use the target-profile contract")

    profile = _load_profile(Path(artifact["path"]), str(artifact["sha256"]))
    fragments = _mapping(config["fragments"], "config.fragments")
    variable_insert = (
        int(fragments["insert_min"]) != int(fragments["insert_mean"])
        or int(fragments["insert_min"]) != int(fragments["insert_max"])
        or float(fragments["insert_stddev"]) != 0.0
    )
    candidate_started = time.perf_counter()
    candidates = _candidate_space(
        Path(config["reference"]),
        fragments,
        len(profile),
        chunk_starts=arguments.candidate_chunk_starts,
    )
    candidate_seconds = time.perf_counter() - candidate_started

    core = resolve_core_executable(arguments.core)
    replay_started = time.perf_counter()
    observed, observed_inserts, replay = _replay_core_histogram(
        manifest, prepared, core, candidates, len(profile)
    )
    replay_seconds = time.perf_counter() - replay_started
    baseline_manifest_path = None
    baseline_replay = None
    uniform_observed = None
    uniform_inserts = None
    baseline_seconds = 0.0
    if arguments.baseline_manifest is not None:
        baseline_manifest_path = (
            arguments.baseline_manifest.expanduser().resolve(strict=True)
        )
        baseline_manifest = _load_manifest(baseline_manifest_path)
        baseline_prepared = _prepare_manifest_run(baseline_manifest)
        baseline_config = baseline_prepared.config.normalized
        baseline_coverage = _mapping(
            baseline_config["coverage"], "baseline.config.coverage"
        )
        if baseline_config.get("technology") != "WGBS" or (
            baseline_coverage.get("kind") != "uniform"
        ):
            raise AuditError("baseline manifest is not uniform WGBS")
        if (
            baseline_config["fragments"] != config["fragments"]
            or baseline_config["seed"] != config["seed"]
            or baseline_config["mutation"] != config["mutation"]
            or baseline_prepared.file_for_role("reference").sha256
            != prepared.file_for_role("reference").sha256
        ):
            raise AuditError(
                "uniform baseline does not match target run fragments, seed, mutation, and reference"
            )
        baseline_started = time.perf_counter()
        uniform_observed, uniform_inserts, baseline_replay = (
            _replay_core_histogram(
                baseline_manifest,
                baseline_prepared,
                core,
                candidates,
                len(profile),
            )
        )
        baseline_seconds = time.perf_counter() - baseline_started
    rows, metrics = _distribution_rows(
        profile,
        candidates.candidate_counts,
        observed,
        uniform_observed,
        project_unreachable=variable_insert,
    )
    attempts = replay["fragment_count"] + replay["skipped_fragment_count"]
    metrics["observed_proposal_acceptance"] = (
        replay["fragment_count"] / attempts
    )
    insert_rows = None
    if uniform_inserts is not None and baseline_replay is not None:
        insert_rows, insert_summary = _insert_metrics(
            int(fragments["insert_min"]),
            observed_inserts,
            uniform_inserts,
        )
        metrics.update(insert_summary)
        baseline_attempts = (
            baseline_replay["fragment_count"]
            + baseline_replay["skipped_fragment_count"]
        )
        metrics["uniform_observed_proposal_acceptance"] = (
            baseline_replay["fragment_count"] / baseline_attempts
        )

    output_tsv = arguments.output_tsv.expanduser().resolve(strict=False)
    output_json = arguments.output_json.expanduser().resolve(strict=False)
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    if output_json.parent != output_tsv.parent:
        output_json.parent.mkdir(parents=True, exist_ok=True)
    _write_tsv(output_tsv, rows)
    if arguments.plot is not None:
        plot = arguments.plot.expanduser().resolve(strict=False)
        if plot.parent != output_tsv.parent and plot.parent != output_json.parent:
            plot.parent.mkdir(parents=True, exist_ok=True)
        _write_plot(plot, rows)
    else:
        plot = None

    summary = {
        "experiment_schema_version": "1.1",
        "source_manifest": str(manifest_path),
        "source_baseline_manifest": (
            None
            if baseline_manifest_path is None
            else str(baseline_manifest_path)
        ),
        "run_id": manifest["run_id"],
        "reference": {
            "path": config["reference"],
            "contig": candidates.contig_name,
            "length": candidates.contig_length,
            "sha256": prepared.file_for_role("reference").sha256,
        },
        "profile": {
            "path": artifact["path"],
            "version": artifact["version"],
            "sha256": artifact["sha256"],
            "bin_count": len(profile),
        },
        "simulation": {
            "depth": fragments.get("depth"),
            "read_length_1": fragments["read_length_1"],
            "read_length_2": fragments["read_length_2"],
            "insert_min": fragments["insert_min"],
            "insert_mean": fragments["insert_mean"],
            "insert_max": fragments["insert_max"],
            "insert_stddev": fragments["insert_stddev"],
            "calibration_insert_length": candidates.calibration_length,
            "max_ambiguous_fraction": fragments["max_ambiguous_fraction"],
            "master_seed": config["seed"],
            "workers": _mapping(config["execution"], "config.execution")[
                "workers"
            ],
            "core_workers": _mapping(config["execution"], "config.execution")[
                "core_workers"
            ],
        },
        "accounting": {
            "eligible_start_count": candidates.valid_start_count,
            "fragment_count": replay["fragment_count"],
            "skipped_fragment_count": replay["skipped_fragment_count"],
            "protocol_batch_fragments": replay["protocol_batch_fragments"],
            "stream_sha256": replay["stream_sha256"],
            "stream_matches_fastq_run": True,
            "uniform_baseline": baseline_replay,
        },
        "metrics": metrics,
        "insert_distribution": insert_rows,
        "timing_seconds": {
            "candidate_histogram": candidate_seconds,
            "core_replay": replay_seconds,
            "uniform_core_replay": baseline_seconds,
            "total": time.perf_counter() - started,
        },
        "outputs": {
            "per_bin_tsv": str(output_tsv),
            "plot": None if plot is None else str(plot),
        },
    }
    output_json.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(output_json)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AuditError as error:
        print("audit_wgbs_gc_target: error: {}".format(error), file=sys.stderr)
        raise SystemExit(1)
