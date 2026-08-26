"""Validate and identify one CLI-derived simulation specification."""

from __future__ import annotations

import copy
from dataclasses import dataclass
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from collections.abc import Iterable, Mapping
from typing import Any

from jsonschema import Draft202012Validator, validators
from jsonschema.exceptions import SchemaError, ValidationError


PathLike = str | Path
UINT64_MAX = (1 << 64) - 1
RUN_CONFIG_SCHEMA_FILENAME = "run-config.schema.json"


class ConfigError(ValueError):
    """A simulation configuration cannot be used."""


class ConfigValidationError(ConfigError):
    """A CLI-derived run document violates the current internal format."""


@dataclass(frozen=True)
class LoadedRunConfig:
    """Canonical configuration identity passed from normalization to launch.

    The normalized document is reconstructed from canonical JSON on access, so
    callers cannot silently mutate a frozen identity and force later layers to
    repeat the complete schema validation.
    """

    master_seed: int | None
    canonical_json: str
    sha256: str

    def __post_init__(self) -> None:
        try:
            normalized = json.loads(self.canonical_json)
        except (TypeError, json.JSONDecodeError) as error:
            raise ConfigValidationError("canonical config JSON is invalid") from error
        if not isinstance(normalized, dict):
            raise ConfigValidationError("canonical config must be a JSON object")
        if _canonical_json(normalized) != self.canonical_json:
            raise ConfigValidationError("config JSON is not canonical")
        observed_sha256 = hashlib.sha256(
            self.canonical_json.encode("utf-8")
        ).hexdigest()
        if observed_sha256 != self.sha256:
            raise ConfigValidationError("config SHA-256 disagrees with canonical JSON")
        if _parse_seed(normalized.get("seed")) != self.master_seed:
            raise ConfigValidationError("config seed identity is inconsistent")

    @property
    def normalized(self) -> dict[str, Any]:
        """Return an independent mutable copy of the normalized document."""
        value = json.loads(self.canonical_json)
        if not isinstance(value, dict):  # Guaranteed by ``__post_init__``.
            raise ConfigValidationError("canonical config must be a JSON object")
        return value

    def as_dict(self) -> dict[str, Any]:
        return self.normalized

    def with_master_seed(self, seed: int) -> "LoadedRunConfig":
        """Insert a generated seed without revalidating an unchanged config."""
        if (
            isinstance(seed, bool)
            or not isinstance(seed, int)
            or not 0 <= seed <= UINT64_MAX
        ):
            raise ConfigValidationError("master seed must be an unsigned 64-bit integer")
        normalized = self.normalized
        normalized["seed"] = str(seed)
        return _freeze_config(normalized)


def normalize_run_config(
    document: Mapping[str, Any],
    base_directory: PathLike,
) -> LoadedRunConfig:
    """Validate values once, apply defaults, resolve paths, and freeze identity.

    Referenced files are not opened here. Preparation hashes them immediately
    before launch, after an omitted seed has been materialized.
    """
    if not isinstance(document, Mapping):
        raise ConfigValidationError("run specification must be a mapping")
    normalized = copy.deepcopy(dict(document))
    errors = sorted(
        _DefaultingDraft202012Validator(_load_schema()).iter_errors(normalized),
        key=_validation_error_sort_key,
    )
    if errors:
        raise ConfigValidationError(_format_schema_errors(errors))

    master_seed = _parse_seed(normalized.get("seed"))
    for name, seed in normalized["seeds"].items():
        if int(seed, 10) > UINT64_MAX:
            raise ConfigValidationError(
                "$.seeds.{}: value exceeds unsigned 64-bit maximum".format(name)
            )
    _validate_cross_field_rules(normalized)

    try:
        base_path = Path(base_directory).expanduser().resolve(strict=False)
    except (TypeError, OSError) as error:
        raise ConfigValidationError("base directory is invalid") from error
    _resolve_paths(normalized, base_path)
    _validate_model_hash_consistency(normalized)
    return _freeze_config(normalized, master_seed=master_seed)


@lru_cache(maxsize=1)
def _load_schema() -> dict[str, Any]:
    schema_path = Path(__file__).parent.parent / RUN_CONFIG_SCHEMA_FILENAME
    try:
        schema = json.loads(schema_path.read_text(encoding="utf-8"))
        Draft202012Validator.check_schema(schema)
    except (OSError, UnicodeError, json.JSONDecodeError, SchemaError) as error:
        raise ConfigError(
            "cannot load packaged normalization schema {}: {}".format(
                schema_path, error
            )
        ) from error
    return schema


def _set_defaults(validator, properties, instance, schema):
    if isinstance(instance, dict):
        for property_name, subschema in properties.items():
            if (
                isinstance(subschema, Mapping)
                and "default" in subschema
                and property_name not in instance
            ):
                instance[property_name] = copy.deepcopy(subschema["default"])
    yield from Draft202012Validator.VALIDATORS["properties"](
        validator, properties, instance, schema
    )


_DefaultingDraft202012Validator = validators.extend(
    Draft202012Validator,
    {"properties": _set_defaults},
)


def _validation_error_sort_key(error: ValidationError) -> tuple[str, str]:
    return (_json_path(error.absolute_path), error.message)


def _format_schema_errors(errors: Iterable[ValidationError]) -> str:
    return "run specification validation failed: " + "; ".join(
        "{}: {}".format(_json_path(error.absolute_path), error.message)
        for error in errors
    )


def _json_path(parts: Iterable[Any]) -> str:
    path = "$"
    for part in parts:
        path += "[{}]".format(part) if isinstance(part, int) else ".{}".format(part)
    return path


def _parse_seed(seed: Any) -> int | None:
    if seed is None:
        return None
    if not isinstance(seed, str) or not seed.isdecimal() or (
        len(seed) > 1 and seed.startswith("0")
    ):
        raise ConfigValidationError("$.seed must be a canonical decimal string")
    value = int(seed, 10)
    if value > UINT64_MAX:
        raise ConfigValidationError("$.seed: value exceeds unsigned 64-bit maximum")
    return value


def _validate_cross_field_rules(config: Mapping[str, Any]) -> None:
    fragments = config["fragments"]
    insert_min = fragments["insert_min"]
    insert_mean = fragments["insert_mean"]
    insert_max = fragments["insert_max"]
    if not insert_min <= insert_mean <= insert_max:
        raise ConfigValidationError(
            "$.fragments: insert_min <= insert_mean <= insert_max must hold"
        )

    read_lengths = [fragments["read_length_1"]]
    if fragments["paired_end"]:
        read_lengths.append(fragments["read_length_2"])
    fixed_mean_insert = (
        config["technology"] in ("WGBS", "TBS", "WGS", "WES", "TS")
        and fragments["insert_sd"] == 0
    )
    read_boundary = insert_mean if fixed_mean_insert else insert_min
    if any(read_length > read_boundary for read_length in read_lengths):
        raise ConfigValidationError(
            "$.fragments: every read length must be <= {}".format(
                "insert_mean" if fixed_mean_insert else "insert_min"
            )
        )

    expected_sections = {
        "WGBS": (False, False),
        "RRBS": (True, False),
        "TBS": (False, True),
        "WGS": (False, False),
        "WES": (False, True),
        "TS": (False, True),
    }
    observed_sections = ("rrbs" in config, "tbs" in config)
    if observed_sections != expected_sections[config["technology"]]:
        raise ConfigValidationError(
            "$.technology: {} requires exactly its matching technology section".format(
                config["technology"]
            )
        )

    output = config["output"]
    if output["fragment_summary"] and output["format"] != "bam":
        raise ConfigValidationError(
            "$.output.fragment_summary: requires format='bam'"
        )
    if output["fragment_realization"] and output["format"] != "bam":
        raise ConfigValidationError(
            "$.output.fragment_realization: requires format='bam'"
        )
    if output["fragment_realization"] and not output["fragment_summary"]:
        raise ConfigValidationError(
            "$.output.fragment_realization: requires fragment_summary=true"
        )

    technology = config["technology"]
    coverage = config["coverage"]
    coverage_kind = coverage["kind"]
    rrbs = config.get("rrbs", {})
    has_rrbs_candidates = "candidate_bed" in rrbs
    if coverage_kind == "profile":
        # Artifact-backed profiles retain their generic normalized-model
        # shape; the released capability gate decides which technology can use
        # them.  The no-artifact form is specifically the RRBS scored-BED
        # exchange and therefore needs its candidate path here.
        if "artifact" not in coverage:
            if technology != "RRBS":
                raise ConfigValidationError(
                    "$.coverage: profile coverage without an artifact requires RRBS"
                )
            if not has_rrbs_candidates:
                raise ConfigValidationError(
                    "$.rrbs.candidate_bed: RRBS profile coverage requires a candidate BED"
                )
    inputs = config["inputs"]
    if technology in ("WGS", "WES", "TS"):
        unsupported = sorted(set(inputs) - {"vcf"})
        if unsupported:
            raise ConfigValidationError(
                "$.inputs: standard sequencing forbids methylation input(s): {}".format(
                    ", ".join(unsupported)
                )
            )
        if config["output"]["save_methdb"]:
            raise ConfigValidationError(
                "$.output.save_methdb: standard sequencing has no MethDB truth"
            )
        if config["output"]["fragment_realization"]:
            raise ConfigValidationError(
                "$.output.fragment_realization: standard sequencing has no "
                "methylation/conversion realization"
            )
        sequencing = config["sequencing"]
        if sequencing["conversion_rate"] != 0:
            raise ConfigValidationError(
                "$.sequencing.conversion_rate: standard sequencing requires 0"
            )
        if not sequencing["directional"]:
            raise ConfigValidationError(
                "$.sequencing.directional: standard sequencing requires the "
                "inert value true"
            )
    if "methdb" in inputs:
        conflicts = ("cgmap", "bed_methyl", "asm", "asm_bed")
        if any(name in inputs for name in conflicts):
            raise ConfigValidationError(
                "$.inputs.methdb: cannot be combined with methylation overlays"
            )
        if config["methylation"]["cgmap_pool"]:
            raise ConfigValidationError(
                "$.inputs.methdb: cannot be combined with cgmap_pool=true"
            )
def _resolve_paths(config: dict[str, Any], base_directory: Path) -> None:
    for path in (
        ("reference",),
        ("output", "directory"),
        ("inputs", "vcf"),
        ("inputs", "cgmap"),
        ("inputs", "bed_methyl"),
        ("inputs", "methdb"),
        ("inputs", "asm"),
        ("inputs", "asm_bed"),
        ("rrbs", "candidate_bed"),
        ("tbs", "bed"),
        ("coverage", "artifact", "path"),
        ("sequencing", "quality", "artifact", "path"),
        ("sequencing", "error", "artifact", "path"),
    ):
        _resolve_path_at(config, path, base_directory)


def _resolve_path_at(
    config: dict[str, Any], path: tuple[str, ...], base_directory: Path
) -> None:
    parent: Any = config
    for part in path[:-1]:
        if not isinstance(parent, dict) or part not in parent:
            return
        parent = parent[part]
    leaf = path[-1]
    if not isinstance(parent, dict) or leaf not in parent:
        return
    value = Path(parent[leaf])
    if not value.is_absolute():
        value = base_directory / value
    parent[leaf] = str(value.resolve(strict=False))


def _validate_model_hash_consistency(config: Mapping[str, Any]) -> None:
    seen_hashes: dict[str, str] = {}
    for label, artifact in _iter_model_artifacts(config):
        path = artifact["path"]
        sha256 = artifact["sha256"]
        previous = seen_hashes.get(path)
        if previous is not None and previous != sha256:
            raise ConfigValidationError(
                "{}: model path {} has conflicting sha256 declarations".format(
                    label, path
                )
            )
        seen_hashes[path] = sha256


def _iter_model_artifacts(
    config: Mapping[str, Any]
) -> Iterable[tuple[str, Mapping[str, Any]]]:
    candidates = (
        ("$.coverage.artifact", config["coverage"]),
        ("$.sequencing.quality.artifact", config["sequencing"]["quality"]),
        ("$.sequencing.error.artifact", config["sequencing"]["error"]),
    )
    for label, container in candidates:
        artifact = container.get("artifact")
        if artifact is not None:
            yield label, artifact


def _canonical_json(document: Mapping[str, Any]) -> str:
    try:
        return json.dumps(
            document,
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        )
    except (TypeError, ValueError, UnicodeError) as error:
        raise ConfigValidationError("normalized config is not JSON-compatible") from error


def _freeze_config(
    normalized: Mapping[str, Any],
    *,
    master_seed: int | None = None,
) -> LoadedRunConfig:
    canonical_json = _canonical_json(normalized)
    if master_seed is None:
        master_seed = _parse_seed(normalized.get("seed"))
    return LoadedRunConfig(
        master_seed=master_seed,
        canonical_json=canonical_json,
        sha256=hashlib.sha256(canonical_json.encode("utf-8")).hexdigest(),
    )


__all__ = [
    "ConfigError",
    "ConfigValidationError",
    "LoadedRunConfig",
    "PathLike",
    "UINT64_MAX",
    "normalize_run_config",
]
