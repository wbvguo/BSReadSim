#!/usr/bin/env python3
"""Measure the current BSReadSim pipeline with repeatable output gates."""

from __future__ import annotations

import argparse
import copy
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
import threading
import time
from typing import Dict, Mapping, Optional, Tuple


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = REPOSITORY_ROOT / "src"
if (
    os.environ.get("BSREADSIM_BENCHMARK_USE_INSTALLED") != "1"
    and str(SOURCE_ROOT) not in sys.path
):
    sys.path.insert(0, str(SOURCE_ROOT))

from bsreadsim.cli import build_parser, build_run_document  # noqa: E402
from bsreadsim.config import normalize_run_config  # noqa: E402
from bsreadsim.manifest import verify_complete_manifest  # noqa: E402
from bsreadsim.pipeline import resolve_core_executable, run_document  # noqa: E402
from bsreadsim.preparation import materialize_master_seed  # noqa: E402


RUN_ID = "00000000-0000-4000-8000-0000000000b1"


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--core", type=Path)
    parser.add_argument("--mode", choices=("production", "debug"), default="production")
    parser.add_argument("--repetitions", type=int, default=5)
    parser.add_argument("--warmup", action="store_true")
    parser.add_argument(
        "--warmup-fragments",
        type=int,
        help="use this many read pairs for warmup; implies --warmup",
    )
    parser.add_argument("--memory-sample-ms", type=int, default=25)
    parser.add_argument(
        "--workspace-root",
        type=Path,
        default=Path("/tmp/bsreadsim-benchmarks"),
    )
    parser.add_argument("--json-out", type=Path)
    parser.add_argument(
        "run_arguments",
        nargs=argparse.REMAINDER,
        help="BSReadSim run arguments after '--'; omit -o/--output, --core, and --mode",
    )
    return parser.parse_args()


class _MemoryMonitor:
    """Sample aggregate Linux process-tree RSS and PSS."""

    def __init__(self, interval_ms: int) -> None:
        self.interval_seconds = interval_ms / 1000.0
        self.stop_event = threading.Event()
        self.thread = threading.Thread(target=self._run, daemon=True)
        self.samples = 0
        self.peak_processes = 0
        self.peak_rss_kib = 0
        self.peak_pss_kib = 0

    def start(self) -> None:
        self.thread.start()

    def stop(self) -> Mapping[str, int]:
        self.stop_event.set()
        self.thread.join()
        return {
            "sample_interval_ms": round(self.interval_seconds * 1000),
            "samples": self.samples,
            "peak_processes": self.peak_processes,
            "peak_rss_kib": self.peak_rss_kib,
            "peak_pss_kib": self.peak_pss_kib,
        }

    def _run(self) -> None:
        root_pid = os.getpid()
        while True:
            memory = tuple(
                sample
                for pid in _process_tree(root_pid)
                for sample in [_process_memory(pid)]
                if sample is not None
            )
            self.samples += 1
            self.peak_processes = max(self.peak_processes, len(memory))
            self.peak_rss_kib = max(
                self.peak_rss_kib, sum(item[0] for item in memory)
            )
            self.peak_pss_kib = max(
                self.peak_pss_kib, sum(item[1] for item in memory)
            )
            if self.stop_event.wait(self.interval_seconds):
                return


def _process_tree(root_pid: int):
    pending = [root_pid]
    observed = set()
    while pending:
        pid = pending.pop()
        if pid in observed:
            continue
        observed.add(pid)
        yield pid
        try:
            children = Path(
                "/proc/{}/task/{}/children".format(pid, pid)
            ).read_text(encoding="ascii")
        except (FileNotFoundError, PermissionError, ProcessLookupError):
            continue
        pending.extend(int(value) for value in children.split())


def _process_memory(pid: int) -> Optional[Tuple[int, int]]:
    rss_kib = None
    pss_kib = None
    try:
        with Path("/proc/{}/smaps_rollup".format(pid)).open(
            "r", encoding="ascii"
        ) as source:
            for line in source:
                if line.startswith("Rss:"):
                    rss_kib = int(line.split()[1])
                elif line.startswith("Pss:"):
                    pss_kib = int(line.split()[1])
    except (FileNotFoundError, PermissionError, ProcessLookupError):
        return None
    if rss_kib is None or pss_kib is None:
        return None
    return rss_kib, pss_kib


def _resource_usage(kind: int) -> Dict[str, float]:
    usage = resource.getrusage(kind)
    return {
        "user_seconds": usage.ru_utime,
        "system_seconds": usage.ru_stime,
        "input_blocks": float(usage.ru_inblock),
        "output_blocks": float(usage.ru_oublock),
    }


def _process_io() -> Dict[str, int]:
    try:
        lines = Path("/proc/self/io").read_text(encoding="ascii").splitlines()
    except (FileNotFoundError, PermissionError):
        return {}
    return {
        key: int(value.strip())
        for key, value in (line.split(":", 1) for line in lines)
    }


def _difference(after: Mapping[str, float], before: Mapping[str, float]):
    return {
        key: after[key] - before.get(key, 0)
        for key in sorted(after)
    }


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


def _document_for_output(
    template: Mapping[str, object], output: Path
) -> Mapping[str, object]:
    document = copy.deepcopy(dict(template))
    document["output"]["directory"] = str(output)  # type: ignore[index]
    return document


def _document_for_warmup(
    template: Mapping[str, object],
    output: Path,
    fragments: Optional[int],
) -> Mapping[str, object]:
    document = _document_for_output(template, output)
    if fragments is None:
        return document
    fragment_config = document.get("fragments")
    if not isinstance(fragment_config, dict) or "read_pairs" not in fragment_config:
        raise ValueError(
            "--warmup-fragments requires a run configured with --read-pairs"
        )
    if fragment_config.get("depth") is not None:
        raise ValueError(
            "--warmup-fragments cannot replace a depth-based workload"
        )
    fragment_config["read_pairs"] = fragments
    return document


def _expected_output_roles(
    template: Mapping[str, object], mode: str
) -> set[str]:
    fragments = template.get("fragments")
    output = template.get("output")
    if not isinstance(fragments, Mapping) or not isinstance(output, Mapping):
        raise ValueError("normalized benchmark template is incomplete")
    roles = {"read1"}
    if fragments.get("paired_end"):
        roles.add("read2")
    if mode == "debug":
        roles.add("truth")
    if output.get("truth_bam"):
        roles.add("truth_bam")
    return roles


def _measure(
    document: Mapping[str, object],
    core: Path,
    mode: str,
    memory_sample_ms: int,
) -> Mapping[str, object]:
    self_before = _resource_usage(resource.RUSAGE_SELF)
    children_before = _resource_usage(resource.RUSAGE_CHILDREN)
    io_before = _process_io()
    monitor = None if memory_sample_ms == 0 else _MemoryMonitor(memory_sample_ms)
    if monitor is not None:
        monitor.start()
    started = time.perf_counter()
    try:
        result = run_document(
            document,
            base_directory=Path.cwd(),
            core_executable=core,
            run_id=RUN_ID,
            mode=mode,
        )
    finally:
        memory = None if monitor is None else monitor.stop()
    wall_seconds = time.perf_counter() - started

    manifest = json.loads(result.manifest_path.read_text(encoding="utf-8"))
    verify_complete_manifest(manifest)
    fragments = manifest["counts"]["core"]["fragment_count"]
    reads = manifest["counts"]["core"]["mate_count"]
    outputs = [
        {
            "role": item.role,
            "record_count": item.record_count,
            "size_bytes": item.size_bytes,
            "sha256": item.sha256,
        }
        for item in result.outputs.files
    ]
    return {
        "wall_seconds": wall_seconds,
        "fragments_per_second": fragments / wall_seconds,
        "reads_per_second": reads / wall_seconds,
        "fragments": fragments,
        "reads": reads,
        "outputs": outputs,
        "output_bytes": sum(item["size_bytes"] for item in outputs),
        "core_counts": manifest["counts"]["core"],
        "python_counts": manifest["counts"]["python"],
        "resource_usage": {
            "self": _difference(
                _resource_usage(resource.RUSAGE_SELF), self_before
            ),
            "children": _difference(
                _resource_usage(resource.RUSAGE_CHILDREN), children_before
            ),
        },
        "process_io": _difference(_process_io(), io_before),
        "process_tree_memory": memory,
    }


def _output_signature(measurement: Mapping[str, object]) -> object:
    return (
        measurement["core_counts"],
        tuple(
            (item["role"], item["record_count"], item["sha256"])
            for item in measurement["outputs"]  # type: ignore[union-attr]
        ),
    )


def _direct_run_arguments(arguments: argparse.Namespace) -> list[str]:
    values = list(arguments.run_arguments)
    if values[:1] == ["--"]:
        values.pop(0)
    if not values:
        raise SystemExit("provide BSReadSim run arguments after '--'")
    for value in values:
        if (
            value in {"-o", "--output", "--core", "--mode"}
            or (value.startswith("-o") and not value.startswith("--"))
            or value.startswith("--output=")
            or value.startswith("--core=")
            or value.startswith("--mode=")
        ):
            raise SystemExit(
                "benchmark manages -o/--output, --core, and --mode"
            )
    return values


def main() -> int:
    arguments = _arguments()
    if arguments.repetitions < 1:
        raise SystemExit("repetitions must be positive")
    if arguments.warmup_fragments is not None and arguments.warmup_fragments < 1:
        raise SystemExit("warmup-fragments must be positive")
    if arguments.memory_sample_ms < 0:
        raise SystemExit("memory-sample-ms must be non-negative")

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")
    workspace = (
        arguments.workspace_root / "{}-{}".format(stamp, os.getpid())
    ).resolve()
    workspace.mkdir(parents=True)

    run_arguments = _direct_run_arguments(arguments)
    parsed_run = build_parser().parse_args(
        [
            "run",
            "-o", str(workspace / "output-placeholder"),
            *run_arguments,
        ]
    )
    invocation_directory = Path.cwd().resolve()
    loaded = materialize_master_seed(
        normalize_run_config(
            build_run_document(parsed_run, invocation_directory),
            invocation_directory,
            mode=arguments.mode,
        )
    )
    template = loaded.as_dict()
    template["output"].pop("truth")
    core = resolve_core_executable(arguments.core)

    warmup = arguments.warmup or arguments.warmup_fragments is not None
    if warmup:
        try:
            warmup_document = _document_for_warmup(
                template,
                workspace / "warmup-output",
                arguments.warmup_fragments,
            )
        except ValueError as error:
            raise SystemExit(str(error)) from error
        _measure(
            warmup_document,
            core,
            arguments.mode,
            arguments.memory_sample_ms,
        )

    measurements = []
    for index in range(arguments.repetitions):
        measurements.append(
            _measure(
                _document_for_output(
                    template,
                    workspace / "run-{}-output".format(index + 1),
                ),
                core,
                arguments.mode,
                arguments.memory_sample_ms,
            )
        )

    expected = _output_signature(measurements[0])
    if any(_output_signature(item) != expected for item in measurements[1:]):
        raise SystemExit("repetitions produced different counts or output bytes")
    expected_roles = _expected_output_roles(template, arguments.mode)
    observed_roles = {item["role"] for item in measurements[0]["outputs"]}
    if observed_roles != expected_roles:
        raise SystemExit("output roles disagree with the selected mode")

    walls = [item["wall_seconds"] for item in measurements]
    throughputs = [item["fragments_per_second"] for item in measurements]
    read_throughputs = [item["reads_per_second"] for item in measurements]
    memory_values = [
        item["process_tree_memory"]["peak_pss_kib"]
        for item in measurements
        if item["process_tree_memory"] is not None
    ]
    report = {
        "schema": "bsreadsim-current-benchmark-1",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "workspace": str(workspace),
        "run_arguments": run_arguments,
        "effective_seed": str(loaded.master_seed),
        "mode": arguments.mode,
        "warmup": warmup,
        "warmup_fragments": arguments.warmup_fragments,
        "repetitions": arguments.repetitions,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "cpu_count": os.cpu_count(),
            "core": str(core),
            "core_sha256": _sha256(core),
            "core_version": _command(core, "--version"),
            "git_commit": _command("git", "-C", REPOSITORY_ROOT, "rev-parse", "HEAD"),
            "git_status": _command("git", "-C", REPOSITORY_ROOT, "status", "--short"),
        },
        "measurements": measurements,
        "summary": {
            "median_wall_seconds": statistics.median(walls),
            "median_fragments_per_second": statistics.median(throughputs),
            "minimum_fragments_per_second": min(throughputs),
            "maximum_fragments_per_second": max(throughputs),
            "median_reads_per_second": statistics.median(read_throughputs),
            "minimum_reads_per_second": min(read_throughputs),
            "maximum_reads_per_second": max(read_throughputs),
            "peak_pss_kib": max(memory_values) if memory_values else None,
            "output_bytes": measurements[0]["output_bytes"],
        },
    }
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if arguments.json_out is not None:
        output_path = arguments.json_out.expanduser().resolve(strict=False)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(encoded, encoding="utf-8")
    sys.stdout.write(encoded)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
