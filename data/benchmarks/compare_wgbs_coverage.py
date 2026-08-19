#!/usr/bin/env python3
"""Paired one-CPU benchmark for uniform and target-GC WGBS coverage."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import statistics
import subprocess
import sys
import tempfile
import time
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if str(SOURCE_ROOT) not in sys.path:
    sys.path.insert(0, str(SOURCE_ROOT))

from bsreadsim.cli import build_parser, build_run_document  # noqa: E402
from bsreadsim.config import normalize_run_config  # noqa: E402
from bsreadsim.core_argv import build_core_argv  # noqa: E402
from bsreadsim.manifest import verify_complete_manifest  # noqa: E402
from bsreadsim.pipeline import resolve_core_executable, run_document  # noqa: E402
from bsreadsim.preparation import prepare_run  # noqa: E402


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


def _arguments(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--profile", type=Path, required=True)
    parser.add_argument("--core", type=Path, required=True)
    parser.add_argument("--fragments", type=int, default=500000)
    parser.add_argument("--warmup-fragments", type=int, default=50000)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=Path("/tmp/bsreadsim-coverage-benchmark"),
    )
    return parser.parse_args(argv)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _command(*arguments: object) -> Optional[str]:
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


def _usage(kind: int) -> Tuple[float, float]:
    observed = resource.getrusage(kind)
    return observed.ru_utime, observed.ru_stime


def _usage_delta(
    self_before: Tuple[float, float],
    children_before: Tuple[float, float],
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
) -> List[str]:
    result = [
        "run",
        "--reference", str(reference),
        "--output", str(output),
        "--read-pairs", str(fragments),
        "--seed", "20260815",
        "--mutation-rate", "0",
        "--read-length", "150",
        "--insert-size", "300",
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
) -> Mapping[str, object]:
    parsed = build_parser().parse_args(
        _direct_arguments(reference, profile, output, fragments, coverage)
    )
    return build_run_document(parsed, REPOSITORY_ROOT)


def _measure_core(
    document: Mapping[str, object],
    core: Path,
    coverage: str,
) -> Mapping[str, Any]:
    loaded = normalize_run_config(document, REPOSITORY_ROOT, mode="production")
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
    by_coverage = {}  # type: Dict[str, List[Mapping[str, Any]]]
    for measurement in measurements:
        by_coverage.setdefault(str(measurement["coverage"]), []).append(measurement)
    if set(by_coverage) != {"uniform", "profile"}:
        raise RuntimeError("benchmark did not produce both coverage modes")

    summaries = {}  # type: Dict[str, Mapping[str, float]]
    for coverage, values in by_coverage.items():
        walls = [float(item["wall_seconds"]) for item in values]
        rates = [float(item["fragments_per_second"]) for item in values]
        summaries[coverage] = {
            "median_wall_seconds": statistics.median(walls),
            "median_fragments_per_second": statistics.median(rates),
            "minimum_fragments_per_second": min(rates),
            "maximum_fragments_per_second": max(rates),
        }

    paired = []  # type: List[Mapping[str, float]]
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
) -> Tuple[List[Mapping[str, Any]], Mapping[str, Any]]:
    measure = _measure_core if lane == "core-producer" else _measure_e2e
    for coverage in ("uniform", "profile"):
        output = workspace / "{}-warmup-{}".format(lane, coverage)
        document = _document(
            reference, profile, output, warmup_fragments, coverage
        )
        print(
            "warmup lane={} coverage={}".format(lane, coverage),
            file=sys.stderr,
            flush=True,
        )
        measure(document, core, coverage)

    measurements = []  # type: List[Mapping[str, Any]]
    for order_index, (pair_index, coverage) in enumerate(PAIRED_ORDER):
        output = workspace / "{}-{}-{}".format(lane, order_index, coverage)
        document = _document(reference, profile, output, fragments, coverage)
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
            if coverage == "uniform" and skipped != 0:
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


def main(argv: Optional[Sequence[str]] = None) -> int:
    arguments = _arguments(argv)
    if arguments.fragments <= 0 or arguments.warmup_fragments <= 0:
        raise SystemExit("fragment counts must be positive")
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
            "read_length": 150,
            "insert_length": 300,
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
