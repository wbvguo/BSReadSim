#!/usr/bin/env python3
"""Paired one-CPU benchmark for uniform and target-GC WGBS coverage."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from collections.abc import Mapping, Sequence
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bsreadsim.cli import build_parser, build_run_document  # noqa: E402
from bsreadsim.run.config import normalize_run_config  # noqa: E402
from bsreadsim.native.launch import (  # noqa: E402
    build_core_argv,
    resolve_core_executable,
)
from bsreadsim.run.manifest import verify_complete_manifest  # noqa: E402
from bsreadsim.run.execute import run_document  # noqa: E402
from bsreadsim.run.prepare import prepare_run  # noqa: E402


RUN_IDS = {
    "uniform": "00000000-0000-4000-8000-0000000000b4",
    "profile": "00000000-0000-4000-8000-0000000000b5",
}
PAIRED_ORDER = (
    (0, "uniform"),
    (0, "profile"),
    (1, "profile"),
    (1, "uniform"),
    (2, "uniform"),
    (2, "profile"),
)


def _arguments(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--fragments", type=int, default=500000)
    parser.add_argument("--warmup-fragments", type=int, default=50000)
    parser.add_argument("--read-length", type=int, default=150)
    parser.add_argument("--insert-min", type=int, default=150)
    parser.add_argument("--insert-mean", type=int, default=400)
    parser.add_argument("--insert-max", type=int, default=1000)
    parser.add_argument("--insert-stddev", type=float, default=25.0)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=Path("/tmp/bsreadsim-coverage-benchmark"),
    )
    return parser.parse_args(argv)


def _sha256(path: Path) -> str:
    with path.open("rb") as source:
        return hashlib.file_digest(source, "sha256").hexdigest()


def _command(*arguments: object) -> str | None:
    try:
        completed = subprocess.run(
            [str(value) for value in arguments],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    return completed.stdout.strip()


def _native_extension_identity() -> Mapping[str, str]:
    try:
        module = importlib.import_module("bsreadsim._native")
    except ImportError as error:
        raise RuntimeError(
            "formal FASTQ benchmarking requires the bsreadsim native extension; "
            "build it with 'python3 setup.py build_ext --inplace'"
        ) from error
    origin = getattr(module, "__file__", None)
    if not isinstance(origin, str):
        raise RuntimeError("bsreadsim native extension has no filesystem identity")
    path = Path(origin).resolve(strict=True)
    return {"path": str(path), "sha256": _sha256(path)}


def _usage(kind: int) -> tuple[float, float]:
    observed = resource.getrusage(kind)
    return observed.ru_utime, observed.ru_stime


def _usage_delta(
    self_before: tuple[float, float],
    children_before: tuple[float, float],
) -> Mapping[str, float]:
    self_after = _usage(resource.RUSAGE_SELF)
    children_after = _usage(resource.RUSAGE_CHILDREN)
    return {
        "user_seconds": (
            self_after[0] - self_before[0]
            + children_after[0] - children_before[0]
        ),
        "system_seconds": (
            self_after[1] - self_before[1]
            + children_after[1] - children_before[1]
        ),
    }


def _direct_arguments(
    reference: Path,
    profile: Path,
    output: Path,
    fragments: int,
    coverage: str,
    read_length: int,
    insert_min: int,
    insert_mean: int,
    insert_max: int,
    insert_stddev: float,
) -> list[str]:
    result = [
        "run",
        "--reference", str(reference),
        "--output", str(output),
        "--read-pairs", str(fragments),
        "--seed", "20260815",
        "--mutation-rate", "0",
        "--read-length", str(read_length),
        "--insert-min", str(insert_min),
        "--insert-mean", str(insert_mean),
        "--insert-max", str(insert_max),
        "--insert-stddev", str(insert_stddev),
        "--max-ambiguous-fraction", "0.05",
        "--error-rate", "0",
        "--workers", "1",
        "--core-workers", "1",
        "--chunk-size", "4096",
        "--max-in-flight-fragments", "4096",
        "--compression", "gzip",
        "--gzip-level", "1",
        "--prefix", "bench",
    ]
    if coverage == "profile":
        result.extend(("--coverage-profile", str(profile)))
    return result


def _document(
    reference: Path,
    profile: Path,
    output: Path,
    fragments: int,
    coverage: str,
    read_length: int,
    insert_min: int,
    insert_mean: int,
    insert_max: int,
    insert_stddev: float,
) -> Mapping[str, object]:
    parsed = build_parser().parse_args(
        _direct_arguments(
            reference,
            profile,
            output,
            fragments,
            coverage,
            read_length,
            insert_min,
            insert_mean,
            insert_max,
            insert_stddev,
        )
    )
    return build_run_document(parsed, REPOSITORY_ROOT)


def _measure_core(
    document: Mapping[str, object],
    core: Path,
    coverage: str,
) -> Mapping[str, Any]:
    loaded = normalize_run_config(document, REPOSITORY_ROOT)
    prepared = prepare_run(loaded)
    execution = prepared.config.normalized["execution"]
    if not isinstance(execution, Mapping):
        raise RuntimeError("normalized execution section is invalid")
    protocol_batch_fragments = min(
        64, int(execution["max_in_flight_fragments"])
    )
    argv = build_core_argv(
        prepared,
        RUN_IDS[coverage],
        core,
        truth_columns="none",
        protocol_batch_fragments=protocol_batch_fragments,
    )
    self_before = _usage(resource.RUSAGE_SELF)
    children_before = _usage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    completed = subprocess.run(
        argv,
        cwd=str(REPOSITORY_ROOT),
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        check=False,
    )
    wall_seconds = time.perf_counter() - started
    cpu = _usage_delta(self_before, children_before)
    if completed.returncode != 0 or completed.stderr:
        raise RuntimeError(
            "core benchmark failed for {}: status={} stderr={!r}".format(
                coverage, completed.returncode, completed.stderr[-4096:]
            )
        )
    fragments = int(prepared.config.normalized["fragments"]["read_pairs"])
    return {
        "wall_seconds": wall_seconds,
        "fragments_per_second": fragments / wall_seconds,
        "cpu": cpu,
        "stdout_sink": "/dev/null",
        "protocol_batch_fragments": protocol_batch_fragments,
    }


def _measure_e2e(
    document: Mapping[str, object],
    core: Path,
    coverage: str,
) -> Mapping[str, Any]:
    self_before = _usage(resource.RUSAGE_SELF)
    children_before = _usage(resource.RUSAGE_CHILDREN)
    started = time.perf_counter()
    result = run_document(
        document,
        base_directory=REPOSITORY_ROOT,
        core_executable=core,
        run_id=RUN_IDS[coverage],
        mode="production",
    )
    wall_seconds = time.perf_counter() - started
    cpu = _usage_delta(self_before, children_before)
    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    verify_complete_manifest(manifest)
    fragments = int(manifest["counts"]["core"]["fragment_count"])
    outputs = [
        {
            "role": item["role"],
            "record_count": item["record_count"],
            "size_bytes": item["size_bytes"],
            "sha256": item["sha256"],
        }
        for item in manifest["outputs"]
    ]
    return {
        "wall_seconds": wall_seconds,
        "fragments_per_second": fragments / wall_seconds,
        "cpu": cpu,
        "stream_sha256": manifest["stream_sha256"],
        "core_counts": manifest["counts"]["core"],
        "python_counts": manifest["counts"]["python"],
        "outputs": outputs,
        "output_bytes": sum(int(item["size_bytes"]) for item in outputs),
        "manifest": str(result.manifest_path),
    }


def _e2e_signature(measurement: Mapping[str, Any]) -> object:
    return (
        measurement["core_counts"],
        measurement["python_counts"],
        tuple(
            (
                item["role"],
                item["record_count"],
                item["sha256"],
            )
            for item in measurement["outputs"]
        ),
    )


def _summarize(measurements: Sequence[Mapping[str, Any]]) -> Mapping[str, Any]:
    by_coverage = {}  # type: dict[str, list[Mapping[str, Any]]]
    for measurement in measurements:
        by_coverage.setdefault(str(measurement["coverage"]), []).append(measurement)
    if set(by_coverage) != {"uniform", "profile"}:
        raise RuntimeError("benchmark did not produce both coverage modes")

    summaries = {}  # type: dict[str, Mapping[str, float]]
    for coverage, values in by_coverage.items():
        walls = [float(item["wall_seconds"]) for item in values]
        rates = [float(item["fragments_per_second"]) for item in values]
        summaries[coverage] = {
            "median_wall_seconds": statistics.median(walls),
            "median_fragments_per_second": statistics.median(rates),
            "minimum_fragments_per_second": min(rates),
            "maximum_fragments_per_second": max(rates),
        }

    paired = []  # type: list[Mapping[str, float]]
    for pair_index in range(3):
        pair = [
            item for item in measurements if int(item["pair_index"]) == pair_index
        ]
        if len(pair) != 2:
            raise RuntimeError("paired benchmark is incomplete")
        rates = {str(item["coverage"]): float(item["fragments_per_second"]) for item in pair}
        paired.append(
            {
                "pair_index": pair_index,
                "profile_to_uniform_throughput_ratio": (
                    rates["profile"] / rates["uniform"]
                ),
            }
        )

    uniform_rate = summaries["uniform"]["median_fragments_per_second"]
    profile_rate = summaries["profile"]["median_fragments_per_second"]
    ratio = profile_rate / uniform_rate
    return {
        "uniform": summaries["uniform"],
        "profile": summaries["profile"],
        "profile_to_uniform_throughput_ratio": ratio,
        "profile_throughput_change_percent": 100.0 * (ratio - 1.0),
        "profile_wall_overhead_percent": 100.0 * (
            summaries["profile"]["median_wall_seconds"]
            / summaries["uniform"]["median_wall_seconds"]
            - 1.0
        ),
        "paired_ratios": paired,
    }


def _run_lane(
    lane: str,
    reference: Path,
    profile: Path,
    core: Path,
    workspace: Path,
    fragments: int,
    warmup_fragments: int,
    read_length: int,
    insert_min: int,
    insert_mean: int,
    insert_max: int,
    insert_stddev: float,
) -> tuple[list[Mapping[str, Any]], Mapping[str, Any]]:
    measure = _measure_core if lane == "core-producer" else _measure_e2e
    for coverage in ("uniform", "profile"):
        output = workspace / "{}-warmup-{}".format(lane, coverage)
        document = _document(
            reference,
            profile,
            output,
            warmup_fragments,
            coverage,
            read_length,
            insert_min,
            insert_mean,
            insert_max,
            insert_stddev,
        )
        print(
            "warmup lane={} coverage={}".format(lane, coverage),
            file=sys.stderr,
            flush=True,
        )
        measure(document, core, coverage)

    measurements = []  # type: list[Mapping[str, Any]]
    for order_index, (pair_index, coverage) in enumerate(PAIRED_ORDER):
        output = workspace / "{}-{}-{}".format(lane, order_index, coverage)
        document = _document(
            reference,
            profile,
            output,
            fragments,
            coverage,
            read_length,
            insert_min,
            insert_mean,
            insert_max,
            insert_stddev,
        )
        print(
            "measure lane={} order={} pair={} coverage={}".format(
                lane, order_index, pair_index, coverage
            ),
            file=sys.stderr,
            flush=True,
        )
        observed = dict(measure(document, core, coverage))
        observed.update(
            {
                "order_index": order_index,
                "pair_index": pair_index,
                "coverage": coverage,
            }
        )
        measurements.append(observed)

    summary = dict(_summarize(measurements))
    if lane == "FASTQ-E2E":
        for coverage in ("uniform", "profile"):
            values = [
                item for item in measurements if item["coverage"] == coverage
            ]
            expected = _e2e_signature(values[0])
            if any(_e2e_signature(item) != expected for item in values[1:]):
                raise RuntimeError(
                    "{} repetitions changed counts or FASTQ bytes".format(coverage)
                )
            skipped = int(values[0]["core_counts"]["skipped_fragment_count"])
            fixed_insert = (
                insert_min == insert_mean == insert_max
                and insert_stddev == 0.0
            )
            if coverage == "uniform" and fixed_insert and skipped != 0:
                raise RuntimeError("uniform sampling reported a rejection")
            if coverage == "profile" and skipped <= 0:
                raise RuntimeError("profile sampling reported no rejections")
            summary[coverage]["output_bytes"] = int(values[0]["output_bytes"])
        profile_counts = next(
            item["core_counts"]
            for item in measurements
            if item["coverage"] == "profile"
        )
        profile_fragments = int(profile_counts["fragment_count"])
        profile_skips = int(profile_counts["skipped_fragment_count"])
        summary["profile_observed_proposal_acceptance"] = (
            profile_fragments / (profile_fragments + profile_skips)
        )
    return measurements, summary


def main(argv: Sequence[str] | None = None) -> int:
    arguments = _arguments(argv)
    if arguments.fragments <= 0 or arguments.warmup_fragments <= 0:
        raise SystemExit("fragment counts must be positive")
    if arguments.read_length <= 0:
        raise SystemExit("--read-length must be positive")
    if not (
        arguments.read_length <= arguments.insert_min
        <= arguments.insert_mean
        <= arguments.insert_max
    ):
        raise SystemExit(
            "read-length <= insert-min <= insert-mean <= insert-max must hold"
        )
    if not math.isfinite(arguments.insert_stddev) or arguments.insert_stddev < 0:
        raise SystemExit("--insert-stddev must be finite and non-negative")
    if arguments.output_json.exists():
        raise SystemExit("refusing to overwrite {}".format(arguments.output_json))
    affinity = sorted(os.sched_getaffinity(0))
    if len(affinity) != 1:
        raise SystemExit(
            "benchmark requires exactly one logical CPU affinity, observed {}".format(
                affinity
            )
        )

    reference = arguments.reference.expanduser().resolve(strict=True)
    profile_path = arguments.profile.expanduser().resolve(strict=True)
    core = resolve_core_executable(arguments.core).resolve(strict=True)
    native_extension = _native_extension_identity()
    workspace_root = arguments.workspace_root.expanduser().resolve(strict=False)
    workspace_root.mkdir(parents=True, exist_ok=True)
    workspace = Path(
        tempfile.mkdtemp(prefix="paired-", dir=str(workspace_root))
    ).resolve()

    lanes = {}
    for lane in ("core-producer", "FASTQ-E2E"):
        measurements, summary = _run_lane(
            lane,
            reference,
            profile_path,
            core,
            workspace,
            arguments.fragments,
            arguments.warmup_fragments,
            arguments.read_length,
            arguments.insert_min,
            arguments.insert_mean,
            arguments.insert_max,
            arguments.insert_stddev,
        )
        lanes[lane] = {
            "measurements": measurements,
            "summary": summary,
        }

    report = {
        "schema": "bsreadsim-wgbs-coverage-paired-benchmark-1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "workload": {
            "fragments": arguments.fragments,
            "reads": 2 * arguments.fragments,
            "warmup_fragments": arguments.warmup_fragments,
            "reference": str(reference),
            "reference_sha256": _sha256(reference),
            "profile": str(profile_path),
            "profile_sha256": _sha256(profile_path),
            "paired_end": True,
            "read_length": arguments.read_length,
            "insert_min": arguments.insert_min,
            "insert_mean": arguments.insert_mean,
            "insert_max": arguments.insert_max,
            "insert_stddev": arguments.insert_stddev,
            "mutation_rate": 0,
            "error_rate": 0,
            "compression": "gzip",
            "gzip_level": 1,
            "workers": 1,
            "core_workers": 1,
            "chunk_size": 4096,
            "max_in_flight_fragments": 4096,
        },
        "execution_lane": {
            "logical_cpu": affinity[0],
            "affinity": affinity,
            "paired_order": [coverage for _, coverage in PAIRED_ORDER],
            "workspace": str(workspace),
        },
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "processor": platform.processor(),
            "lscpu_json": _command("lscpu", "-J"),
            "core": str(core),
            "core_sha256": _sha256(core),
            "core_version": _command(core, "--version"),
            "python_native_extension": native_extension,
            "git_commit": _command(
                "git", "-C", REPOSITORY_ROOT, "rev-parse", "HEAD"
            ),
            "git_status": _command(
                "git", "-C", REPOSITORY_ROOT, "status", "--short"
            ),
        },
        "lanes": lanes,
    }
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    output = arguments.output_json.expanduser().resolve(strict=False)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(encoded, encoding="utf-8")
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
