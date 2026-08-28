"""Strict projection of a prepared run onto the C++ core argv contract."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import json
import math
import os
from pathlib import Path
import re
import uuid

from ..run.config import LoadedRunConfig, UINT64_MAX
from ..run.prepare import FileDigest, PreparedRun


PathLike = str | os.PathLike
UINT32_MAX = (1 << 32) - 1
_HEX_DIGEST = re.compile(r"^[0-9a-f]{64}$")


class CoreArgvError(ValueError):
    """A prepared run cannot be represented by the core CLI contract."""


def build_core_argv(
    prepared: PreparedRun,
    run_id: str,
    core_executable: PathLike,
    *,
    emit_details: bool = False,
    protocol_batch_fragments: int = 64,
    methdb_output_path: PathLike | None = None,
) -> tuple[str, ...]:
    """Return one complete, deterministic ``htsim-core`` argv tuple.

    Only fields owned by the C++ generator cross this boundary. Input paths
    come from ``PreparedRun.files`` rather than being copied directly from the
    normalized config; file digests remain manifest metadata and are not core
    command-line inputs.
    """

    executable = _path_argument("core_executable", core_executable)
    canonical_run_id = _canonical_uuid(run_id)
    if not isinstance(emit_details, bool):
        raise CoreArgvError("emit_details must be a boolean")
    if (
        isinstance(protocol_batch_fragments, bool)
        or not isinstance(protocol_batch_fragments, int)
        or not 1 <= protocol_batch_fragments <= 64
    ):
        raise CoreArgvError("protocol_batch_fragments must be in [1, 64]")
    config, roles = _validate_prepared_run(prepared)

    fragments = _mapping(config, "fragments")
    execution = _mapping(config, "execution")
    mutation = _mapping(config, "mutation")
    sequencing = _mapping(config, "sequencing")
    seeds = _mapping(config, "seeds")
    methylation = _mapping(config, "methylation")
    beta = _mapping(methylation, "beta")
    coverage = _mapping(config, "coverage")
    inputs = _mapping(config, "inputs")

    arguments = [
        executable,
        "--emit-details",
        str(emit_details).lower(),
        "--run-id",
        canonical_run_id,
        "--config-sha256",
        _digest("config.sha256", prepared.config.sha256),
        "--seed",
        _unsigned("config.seed", prepared.config.master_seed, UINT64_MAX),
        "--seed-mut",
        _unsigned(
            "seeds.mutation",
            int(_text("seeds.mutation", seeds["mutation"]), 10),
            UINT64_MAX,
        ),
        "--seed-phase",
        _unsigned(
            "seeds.phasing",
            int(_text("seeds.phasing", seeds["phasing"]), 10),
            UINT64_MAX,
        ),
        "--seed-meth",
        _unsigned(
            "seeds.methylation",
            int(_text("seeds.methylation", seeds["methylation"]), 10),
            UINT64_MAX,
        ),
    ]
    arguments.extend(
        (
            "--protocol-batch-fragments",
            str(protocol_batch_fragments),
        )
    )

    reference = roles["reference"]
    _append_file(arguments, "--reference", reference)

    for input_name in ("vcf", "cgmap", "bed_methyl", "methdb", "asm", "asm_bed"):
        if input_name in inputs:
            option_name = input_name.replace("_", "-")
            _append_file(
                arguments,
                "--" + option_name,
                roles["input." + input_name],
            )
    if methdb_output_path is not None:
        if "methdb" in inputs:
            raise CoreArgvError(
                "methdb_output_path cannot be combined with a MethDB input"
            )
        arguments.extend(
            (
                "--methdb-output",
                _path_argument("methdb_output_path", methdb_output_path),
            )
        )

    technology = _text("technology", config["technology"])
    arguments.extend(("--technology", technology))
    arguments.extend(
        (
            "--directional",
            _boolean("sequencing.directional", sequencing["directional"]),
        )
    )

    paired_end = _boolean("fragments.paired_end", fragments["paired_end"])
    arguments.extend(("--paired-end", paired_end))
    arguments.extend(
        (
            "--read-length-1",
            _unsigned(
                "fragments.read_length_1", fragments["read_length_1"], UINT32_MAX
            ),
        )
    )
    if fragments["paired_end"]:
        arguments.extend(
            (
                "--read-length-2",
                _unsigned(
                    "fragments.read_length_2",
                    fragments["read_length_2"],
                    UINT32_MAX,
                ),
            )
        )

    for option, field in (
        ("--insert-min", "insert_min"),
        ("--insert-mean", "insert_mean"),
        ("--insert-max", "insert_max"),
    ):
        arguments.extend(
            (option, _unsigned("fragments." + field, fragments[field], UINT32_MAX))
        )
    arguments.extend(
        (
            "--insert-sd",
            _number("fragments.insert_sd", fragments["insert_sd"]),
        )
    )
    if "depth" in fragments:
        arguments.extend(("--depth", _number("fragments.depth", fragments["depth"])))
    else:
        arguments.extend(
            (
                "--fragments",
                _unsigned(
                    "fragments.count", fragments["count"], UINT32_MAX
                ),
            )
        )
    arguments.extend(
        (
            "--max-ambiguous-fraction",
            _number(
                "fragments.max_ambiguous_fraction",
                fragments["max_ambiguous_fraction"],
            ),
            "--chunk-size",
            _unsigned("execution.chunk_size", execution["chunk_size"], UINT32_MAX),
            "--core-workers",
            _unsigned("execution.core_workers", execution["core_workers"], 64),
        )
    )

    if technology == "RRBS":
        rrbs = _mapping(config, "rrbs")
        for cut_site in rrbs["cut_sites"]:
            arguments.extend(("--rrbs-cut-site", _text("rrbs.cut_site", cut_site)))
        if "candidate_bed" in rrbs:
            candidate_bed = Path(
                _text("rrbs.candidate_bed", rrbs["candidate_bed"])
            )
            if not candidate_bed.is_absolute():
                raise CoreArgvError(
                    "normalized rrbs.candidate_bed path must be absolute"
                )
            arguments.extend(
                (
                    "--rrbs-candidate-bed",
                    _path_text("rrbs.candidate_bed", candidate_bed),
                )
            )
    elif technology in ("TBS", "WES", "TS"):
        tbs = _mapping(config, "tbs")
        _append_file(arguments, "--tbs-bed", roles["input.tbs-bed"])
        arguments.extend(
            (
                "--tbs-center-stddev",
                _number(
                    "tbs.fragment_center_stddev", tbs["fragment_center_stddev"]
                ),
            )
        )
    elif technology not in ("WGBS", "WGS"):
        raise CoreArgvError(
            "technology must be WGBS, RRBS, TBS, WGS, WES, or TS"
        )

    coverage_kind = _text("coverage.kind", coverage["kind"])
    arguments.extend(("--coverage", coverage_kind))
    if coverage_kind == "profile" and "artifact" in coverage:
        profile = roles["model.coverage"]
        arguments.extend(
            (
                "--coverage-profile",
                _path_text("model.coverage.path", profile.path),
            )
        )
    elif coverage_kind == "profile" and technology == "RRBS":
        if "candidate_bed" not in _mapping(config, "rrbs"):
            raise CoreArgvError(
                "RRBS profile coverage requires rrbs.candidate_bed"
            )
    elif coverage_kind not in ("uniform", "target-score"):
        raise CoreArgvError(
            "coverage.kind must be uniform, profile, or target-score"
        )

    arguments.extend(
        (
            "--mutation-rate",
            _number("mutation.rate", mutation["rate"]),
            "--indel-fraction",
            _number("mutation.indel_fraction", mutation["indel_fraction"]),
            "--indel-extension-probability",
            _number(
                "mutation.indel_extension_probability",
                mutation["indel_extension_probability"],
            ),
            "--homozygous-only",
            _boolean("mutation.homozygous_only", mutation["homozygous_only"]),
            "--collect-non-cpg",
            _boolean("methylation.collect_non_cpg", methylation["collect_non_cpg"]),
            "--cgmap-pool",
            _boolean("methylation.cgmap_pool", methylation["cgmap_pool"]),
            "--update-variant-boundaries",
            _boolean(
                "methylation.update_variant_boundaries",
                methylation["update_variant_boundaries"],
            ),
            "--beta-cg",
            _beta_pair("methylation.beta.CG", beta["CG"]),
            "--beta-chg",
            _beta_pair("methylation.beta.CHG", beta["CHG"]),
            "--beta-chh",
            _beta_pair("methylation.beta.CHH", beta["CHH"]),
        )
    )

    for index, argument in enumerate(arguments):
        _text("argv[{}]".format(index), argument)
    return tuple(arguments)


def _validate_prepared_run(
    prepared: PreparedRun,
) -> tuple[Mapping[str, Any], dict[str, FileDigest]]:
    if not isinstance(prepared, PreparedRun):
        raise CoreArgvError("prepared must be a PreparedRun")
    if not isinstance(prepared.config, LoadedRunConfig):
        raise CoreArgvError("prepared.config must be a LoadedRunConfig")
    if prepared.config.master_seed is None:
        raise CoreArgvError("prepared config must contain a materialized seed")

    if not isinstance(prepared.files, tuple):
        raise CoreArgvError("prepared files must be a tuple")
    config = prepared.config.normalized
    expected = _expected_roles(config)
    observed = {}  # type: dict[str, FileDigest]
    for index, file_digest in enumerate(prepared.files):
        if not isinstance(file_digest, FileDigest):
            raise CoreArgvError(
                "prepared.files[{}] must be a FileDigest".format(index)
            )
        role = _text("prepared.files[{}].role".format(index), file_digest.role)
        if role in observed:
            raise CoreArgvError("duplicate prepared file role: {}".format(role))
        _validate_file_digest(file_digest, index)
        observed[role] = file_digest

    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    if missing:
        raise CoreArgvError("missing prepared file roles: {}".format(", ".join(missing)))
    if unexpected:
        raise CoreArgvError(
            "unexpected prepared file roles: {}".format(", ".join(unexpected))
        )

    for role, (path, declared_sha256) in expected.items():
        file_digest = observed[role]
        if file_digest.path != path:
            raise CoreArgvError(
                "prepared file path for {} disagrees with normalized config".format(role)
            )
        if declared_sha256 is None:
            if file_digest.declared_sha256 is not None:
                raise CoreArgvError(
                    "prepared file {} has an unexpected declared digest".format(role)
                )
        elif (
            file_digest.declared_sha256 != declared_sha256
            or file_digest.sha256 != declared_sha256
        ):
            raise CoreArgvError(
                "prepared file digest for {} disagrees with its declaration".format(role)
            )
    return config, observed


def _expected_roles(
    config: Mapping[str, Any],
) -> dict[str, tuple[Path, str | None]]:
    expected = {}  # type: dict[str, tuple[Path, str | None]]

    def add(role: str, path_value: Any, declared_sha256: str | None = None) -> None:
        if role in expected:
            raise CoreArgvError("duplicate expected file role: {}".format(role))
        path = Path(_text(role + ".path", path_value))
        if not path.is_absolute():
            raise CoreArgvError("normalized path for {} must be absolute".format(role))
        if declared_sha256 is not None:
            declared_sha256 = _digest(role + ".declared_sha256", declared_sha256)
        expected[role] = (path, declared_sha256)

    add("reference", config["reference"])
    inputs = _mapping(config, "inputs")
    for name in ("vcf", "cgmap", "bed_methyl", "methdb", "asm", "asm_bed"):
        if name in inputs:
            add("input." + name, inputs[name])
    if config["technology"] in ("TBS", "WES", "TS"):
        add("input.tbs-bed", _mapping(config, "tbs")["bed"])

    model_containers = (
        ("model.coverage", config["coverage"]),
        ("model.quality", _mapping(config, "sequencing")["quality"]),
        ("model.error", _mapping(config, "sequencing")["error"]),
    )
    for role, container_value in model_containers:
        container = _mapping_value(role, container_value)
        artifact_value = container.get("artifact")
        if artifact_value is None:
            continue
        artifact = _mapping_value(role + ".artifact", artifact_value)
        add(role, artifact["path"], artifact["sha256"])
    return expected


def _validate_file_digest(file_digest: FileDigest, index: int) -> None:
    label = "prepared.files[{}]".format(index)
    if not isinstance(file_digest.path, Path):
        raise CoreArgvError(label + ".path must be a Path")
    if not file_digest.path.is_absolute():
        raise CoreArgvError(label + ".path must be absolute")
    _path_text(label + ".path", file_digest.path)
    if isinstance(file_digest.size_bytes, bool) or not isinstance(
        file_digest.size_bytes, int
    ):
        raise CoreArgvError(label + ".size_bytes must be a non-negative integer")
    if file_digest.size_bytes < 0:
        raise CoreArgvError(label + ".size_bytes must be a non-negative integer")
    _digest(label + ".sha256", file_digest.sha256)
    if file_digest.declared_sha256 is not None:
        _digest(label + ".declared_sha256", file_digest.declared_sha256)


def _mapping(parent: Mapping[str, Any], key: str) -> Mapping[str, Any]:
    try:
        value = parent[key]
    except (KeyError, TypeError) as error:
        raise CoreArgvError("missing normalized config field: {}".format(key)) from error
    return _mapping_value(key, value)


def _mapping_value(name: str, value: Any) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise CoreArgvError("{} must be an object".format(name))
    return value


def _canonical_uuid(value: Any) -> str:
    if not isinstance(value, str):
        raise CoreArgvError("run_id must be canonical lowercase UUID text")
    try:
        parsed = uuid.UUID(value)
    except (ValueError, AttributeError) as error:
        raise CoreArgvError("run_id must be canonical lowercase UUID text") from error
    if str(parsed) != value:
        raise CoreArgvError("run_id must be canonical lowercase UUID text")
    return value


def _path_argument(name: str, value: PathLike) -> str:
    try:
        path = os.fspath(value)
    except TypeError as error:
        raise CoreArgvError("{} must be a text path".format(name)) from error
    if not isinstance(path, str):
        raise CoreArgvError("{} must resolve to text, not bytes".format(name))
    return _text(name, path)


def _path_text(name: str, value: Path) -> str:
    return _text(name, str(value))


def _text(name: str, value: Any) -> str:
    if not isinstance(value, str) or not value or "\x00" in value:
        raise CoreArgvError("{} must be non-empty text without NUL".format(name))
    return value


def _digest(name: str, value: Any) -> str:
    if not isinstance(value, str) or _HEX_DIGEST.fullmatch(value) is None:
        raise CoreArgvError("{} must contain 64 lowercase hex digits".format(name))
    return value


def _unsigned(name: str, value: Any, maximum: int) -> str:
    if (
        isinstance(value, bool)
        or not isinstance(value, int)
        or value < 0
        or value > maximum
    ):
        raise CoreArgvError("{} is outside its unsigned integer range".format(name))
    return str(value)


def _number(name: str, value: Any) -> str:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CoreArgvError("{} must be a finite JSON number".format(name))
    try:
        converted = float(value)
    except (OverflowError, ValueError) as error:
        raise CoreArgvError("{} must be a finite JSON number".format(name)) from error
    if not math.isfinite(converted):
        raise CoreArgvError("{} must be a finite JSON number".format(name))
    try:
        return json.dumps(value, allow_nan=False, separators=(",", ":"))
    except (TypeError, ValueError, OverflowError) as error:
        raise CoreArgvError("{} must be a finite JSON number".format(name)) from error


def _boolean(name: str, value: Any) -> str:
    if not isinstance(value, bool):
        raise CoreArgvError("{} must be a boolean".format(name))
    return "true" if value else "false"


def _beta_pair(name: str, value: Any) -> str:
    if not isinstance(value, list) or len(value) != 2:
        raise CoreArgvError("{} must contain two numbers".format(name))
    return "{},{}".format(
        _number(name + "[0]", value[0]),
        _number(name + "[1]", value[1]),
    )


def _append_file(arguments: list, option: str, file_digest: FileDigest) -> None:
    arguments.extend(
        (
            option,
            _path_text(file_digest.role + ".path", file_digest.path),
        )
    )


CORE_FILENAME = "htsim-core.exe" if os.name == "nt" else "htsim-core"


class CoreExecutableError(RuntimeError):
    """The requested or installed generation core is unavailable."""


def packaged_core_candidate() -> Path:
    """Return the expected core location without asserting that it exists."""
    return Path(__file__).resolve().parent / CORE_FILENAME


def resolve_core_executable(value: PathLike | None = None) -> Path:
    """Resolve the bundled core or one explicit filesystem override."""

    if value is None:
        candidate = packaged_core_candidate()
    else:
        try:
            text = os.fspath(value)
        except TypeError as error:
            raise CoreExecutableError("core executable must be a text path") from error
        if not isinstance(text, str) or not text or "\x00" in text:
            raise CoreExecutableError(
                "core executable must be a non-empty text path"
            )
        candidate = Path(text).expanduser()

    try:
        resolved = candidate.resolve(strict=True)
    except OSError as error:
        raise CoreExecutableError(
            "cannot resolve core executable {}: {}".format(candidate, error)
        ) from error
    if not resolved.is_file():
        raise CoreExecutableError(
            "core executable is not a regular file: {}".format(resolved)
        )
    if not os.access(str(resolved), os.X_OK):
        raise CoreExecutableError(
            "core executable is not executable: {}".format(resolved)
        )
    return resolved


__all__ = [
    "CORE_FILENAME",
    "CoreArgvError",
    "CoreExecutableError",
    "build_core_argv",
    "packaged_core_candidate",
    "resolve_core_executable",
]
