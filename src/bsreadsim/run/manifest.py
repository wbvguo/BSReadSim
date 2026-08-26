"""Construction and verification of the final BSReadSim run manifest."""

from __future__ import annotations

import copy
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path
import shlex
from collections.abc import Mapping, Sequence
from typing import Any

from .. import __version__
from ..output.bam import (
    ANNOTATION_FRAGMENT_REALIZATION_SCHEMA,
    ANNOTATION_FRAGMENT_SUMMARY_SCHEMA,
    ANNOTATION_READ_SUMMARY_SCHEMA,
    ANNOTATION_STATE_SCHEMA,
    BAM_CONTRACT,
    BAM_MAPQ,
)
from ..output import OutputFileSummary, OutputSummary
from .prepare import FileDigest, PreparedRun
from ..htsim.protocol import (
    AmbiguityPolicy,
    BaseEncoding,
    PROTOCOL_VERSION,
    Header,
    Technology,
    Trailer,
)
from ..process.batch import READ_NAME_CONTRACT
from ..rng import RNG_CONTRACT, STAGE_NAMES
from .invocation import FullCommandError, build_full_run_argv


MANIFEST_VERSION = 2
COORDINATE_CONVENTION = "0-based-half-open"
METHDB_MAGIC = b"methdb"
METHDB_VERSION = 1


class ManifestError(ValueError):
    """Run evidence is inconsistent and cannot form a complete manifest."""


@dataclass(frozen=True)
class CompleteManifest:
    """Canonical complete-manifest value ready for transactional commit."""

    document: dict[str, Any]
    canonical_json: str
    reproducibility_sha256: str

    def as_dict(self) -> dict[str, Any]:
        return copy.deepcopy(self.document)


def build_complete_manifest(
    prepared: PreparedRun,
    header: Header,
    trailer: Trailer,
    outputs: OutputSummary,
    *,
    invocation_argv: Sequence[str] | None = None,
) -> CompleteManifest:
    """Cross-check core/Python evidence and build the final commit marker."""
    if not isinstance(prepared, PreparedRun):
        raise ManifestError("prepared must be a PreparedRun")
    config = prepared.config
    if config.master_seed is None:
        raise ManifestError("the effective config must contain a materialized seed")
    if not isinstance(header, Header) or not isinstance(trailer, Trailer):
        raise ManifestError("header and trailer must use the protocol contract")
    validate_header_projection(prepared, header)
    if not isinstance(outputs, OutputSummary):
        raise ManifestError("outputs must be an OutputSummary")

    normalized = config.normalized
    _validate_counts(normalized, header, trailer, outputs)
    _validate_output_paths(normalized, outputs)

    bisulfite = normalized["technology"] in ("WGBS", "RRBS", "TBS")
    methylation_state_model = (
        {
            "contract": "bernoulli-site",
            "effective": "bernoulli",
            "requested": normalized["methylation"]["state_model"],
        }
        if bisulfite
        else {
            "effective": "disabled",
            "requested": "disabled",
        }
    )
    fragments = normalized["fragments"]
    read_base_count = trailer.fragment_count * (
        fragments["read_length_1"]
        + (fragments["read_length_2"] if fragments["paired_end"] else 0)
    )
    details = {
        "configuration": _manifest_effective_config(normalized),
        "configuration_sha256": config.sha256,
        "contigs": [
            {
                "fragment_count": trailer.per_contig_fragment_counts[index],
                "index": index,
                "length": contig.length,
                "name": contig.name,
                "reference_sequence_sha256": contig.reference_sha256.hex(),
            }
            for index, contig in enumerate(header.contigs)
        ],
        "contracts": {
            "read_name": READ_NAME_CONTRACT,
            "rng": RNG_CONTRACT,
        },
        "models": {"methylation_state": methylation_state_model},
        "protocol_version": PROTOCOL_VERSION,
        "randomness": {
            "master_seed": str(config.master_seed),
            "methylation_seed": normalized["seeds"]["methylation"],
            "mutation_seed": normalized["seeds"]["mutation"],
            "phasing_seed": normalized["seeds"]["phasing"],
            "stages": list(STAGE_NAMES),
        },
        "reproducibility": {
            "numerical_tolerance_exceptions": [],
            "scope": (
                "same released core/Python versions, effective config, "
                "input/model digests, and master seed"
            ),
        },
        "software_versions": {
            "core": header.core_version,
            "python": __version__,
        },
        "stream_sha256": trailer.stream_sha256.hex(),
    }

    document = {
        "command": _command_manifest_entry(invocation_argv, normalized),
        "details": details,
        "inputs": [
            _input_manifest_entry(file_digest, normalized)
            for file_digest in prepared.files
        ],
        "outputs": [
            {
                "path": str(item.path),
                "record_count": item.record_count,
                "role": item.role,
                "sha256": item.sha256,
                "size_bytes": item.size_bytes,
            }
            for item in outputs.files
        ],
        "run_id": header.run_id,
        "status": "complete",
        "summary": {
            "fragment_count": trailer.fragment_count,
            "methylation_site_count": trailer.methylation_site_count,
            "output_file_count": len(outputs.files),
            "output_format": normalized["output"]["format"],
            "output_size_bytes": sum(item.size_bytes for item in outputs.files),
            "paired_end": fragments["paired_end"],
            "read_base_count": read_base_count,
            "read_count": trailer.mate_count,
            "skipped_fragment_count": trailer.skipped_fragment_count,
            "technology": normalized["technology"],
            "template_base_count": trailer.template_base_count,
        },
        "version": MANIFEST_VERSION,
    }
    if normalized["output"]["format"] == "bam":
        fragment_summary = bool(normalized["output"]["fragment_summary"])
        fragment_realization = bool(normalized["output"]["fragment_realization"])
        details["alignment"] = {
            "coordinate_convention": COORDINATE_CONVENTION,
            "format": "BAM",
            "mapq": BAM_MAPQ,
            "mapq_semantics": "simulated-origin-not-calibrated-confidence",
            "sam_version": "1.6",
            "sort_order": "unsorted",
            "fastq_sidecars": False,
            "fastq_recovery": "samtools fastq",
            "tags": {
                "zx": {
                    "required": fragment_realization,
                    "schema": ANNOTATION_FRAGMENT_REALIZATION_SCHEMA,
                    "scope": "complete-physical-fragment-realization",
                },
                "zf": {
                    "required": fragment_summary,
                    "schema": ANNOTATION_FRAGMENT_SUMMARY_SCHEMA,
                    "scope": "complete-physical-fragment",
                },
                "zr": {
                    "required": True,
                    "schema": ANNOTATION_READ_SUMMARY_SCHEMA,
                    "scope": "single-read",
                },
                "zt": {
                    "alphabet": "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_",
                    "required": True,
                    "schema": ANNOTATION_STATE_SCHEMA,
                    "scope": "one-character-per-bam-seq-base",
                },
            },
        }
        details["contracts"]["bam"] = BAM_CONTRACT

    identity_json = _canonical_json(document)
    reproducibility_digest = hashlib.sha256(
        identity_json.encode("utf-8")
    ).hexdigest()
    details["reproducibility"]["sha256"] = reproducibility_digest
    verify_complete_manifest(document)
    canonical_json = _canonical_json(document)
    return CompleteManifest(
        document=document,
        canonical_json=canonical_json,
        reproducibility_sha256=reproducibility_digest,
    )


def validate_header_projection(prepared: PreparedRun, header: Header) -> None:
    """Validate all Python-owned identity and execution header fields."""
    if not isinstance(prepared, PreparedRun):
        raise ManifestError("prepared must be a PreparedRun")
    if not isinstance(header, Header):
        raise ManifestError("header must be a protocol Header")
    normalized_output = prepared.config.normalized["output"]
    expected_has_details = normalized_output["format"] == "bam"
    _validate_protocol_projection(prepared, header)

    normalized = prepared.config.normalized
    technology = {
        "WGBS": Technology.WGBS,
        "RRBS": Technology.RRBS,
        "TBS": Technology.TBS,
        "WGS": Technology.WGS,
        "WES": Technology.WES,
        "TS": Technology.TS,
    }.get(normalized["technology"])
    if header.technology is not technology:
        raise ManifestError("core technology disagrees with the effective config")
    fragments = normalized["fragments"]
    if not isinstance(fragments, Mapping):
        raise ManifestError("normalized fragments section is invalid")
    paired_end = fragments["paired_end"]
    expected_mates = 2 if paired_end else 1
    if header.mates_per_fragment != expected_mates:
        raise ManifestError("core mate count disagrees with the effective config")
    if header.read_length_r1 != fragments["read_length_1"]:
        raise ManifestError("core R1 length disagrees with the effective config")
    expected_r2 = fragments["read_length_2"] if paired_end else 0
    if header.read_length_r2 != expected_r2:
        raise ManifestError("core R2 length disagrees with the effective config")
    if header.has_details is not expected_has_details:
        raise ManifestError("core detail policy disagrees with the requested output")
    if header.base_encoding is not BaseEncoding.ACGTN_U8:
        raise ManifestError("core base encoding is unsupported")
    if header.ambiguity_policy is not AmbiguityPolicy.PRESERVE_N:
        raise ManifestError("core ambiguity policy is unsupported")


def verify_complete_manifest(document: Mapping[str, Any]) -> None:
    """Fail unless *document* is complete and its stable digest verifies."""
    if not isinstance(document, Mapping):
        raise ManifestError("manifest must be an object")
    if document.get("status") != "complete":
        raise ManifestError("manifest status must be complete")
    manifest_version = document.get("version")
    if type(manifest_version) is not int or manifest_version != MANIFEST_VERSION:
        raise ManifestError("manifest version is unsupported")
    details = document.get("details")
    if not isinstance(details, Mapping):
        raise ManifestError("manifest details section is missing")
    normalized_config = _manifest_run_config(
        details.get("configuration"),
        details.get("configuration_sha256"),
        details.get("randomness"),
    )
    _validate_command(document.get("command"), normalized_config)
    reproducibility = details.get("reproducibility")
    if not isinstance(reproducibility, Mapping):
        raise ManifestError("manifest reproducibility section is missing")
    observed_digest = reproducibility.get("sha256")
    if (
        not isinstance(observed_digest, str)
        or len(observed_digest) != 64
        or any(character not in "0123456789abcdef" for character in observed_digest)
    ):
        raise ManifestError("manifest reproducibility SHA-256 is invalid")

    identity = copy.deepcopy(dict(document))
    identity_details = identity.get("details")
    if not isinstance(identity_details, dict):
        raise ManifestError("manifest details section is invalid")
    identity_reproducibility = identity_details.get("reproducibility")
    if not isinstance(identity_reproducibility, dict):
        raise ManifestError("manifest reproducibility section is invalid")
    identity_reproducibility.pop("sha256", None)
    expected_digest = hashlib.sha256(
        _canonical_json(identity).encode("utf-8")
    ).hexdigest()
    if observed_digest != expected_digest:
        raise ManifestError("manifest reproducibility SHA-256 mismatch")


def _validate_protocol_projection(
    prepared: PreparedRun, header: Header
) -> None:
    config = prepared.config
    if header.rng_contract != RNG_CONTRACT:
        raise ManifestError("core RNG contract disagrees with Python")
    if header.master_seed != config.master_seed:
        raise ManifestError("core master seed disagrees with the effective config")
    try:
        expected_config_digest = bytes.fromhex(config.sha256)
    except ValueError as error:
        raise ManifestError("effective config SHA-256 is invalid") from error
    if header.normalized_config_sha256 != expected_config_digest:
        raise ManifestError("core normalized-config digest disagrees with Python")
    if not header.contigs:
        raise ManifestError("core header contains no contigs")


def _validate_counts(
    config: Mapping[str, Any],
    header: Header,
    trailer: Trailer,
    outputs: OutputSummary,
) -> None:
    if (
        config["technology"] in ("WGS", "WES", "TS")
        and trailer.methylation_site_count != 0
    ):
        raise ManifestError(
            "standard sequencing cannot report methylation sites"
        )
    paired_end = bool(config["fragments"]["paired_end"])
    expected_mates = trailer.fragment_count * (2 if paired_end else 1)
    if trailer.mate_count != expected_mates:
        raise ManifestError("core mate count disagrees with SE/PE configuration")
    if outputs.fragment_count != trailer.fragment_count:
        raise ManifestError("Python/core fragment counts disagree")
    if outputs.mate_count != trailer.mate_count:
        raise ManifestError("Python/core mate counts disagree")
    if len(trailer.per_contig_fragment_counts) != len(header.contigs):
        raise ManifestError("per-contig count cardinality disagrees with header")
    if sum(trailer.per_contig_fragment_counts) != trailer.fragment_count:
        raise ManifestError("per-contig counts do not sum to the fragment count")

    bam = config["output"]["format"] == "bam"
    expected_roles = set()
    if not bam:
        expected_roles.add("read1")
        if paired_end:
            expected_roles.add("read2")
    if bam:
        expected_roles.add("bam")
    if config["output"]["save_methdb"]:
        expected_roles.add("truth.methdb")
    if config["output"]["save_vcf"]:
        expected_roles.add("truth.vcf")
    observed_roles = {item.role for item in outputs.files}
    if observed_roles != expected_roles or len(outputs.files) != len(expected_roles):
        raise ManifestError("output roles disagree with SE/PE configuration")
    for item in outputs.files:
        if not isinstance(item, OutputFileSummary):
            raise ManifestError("output files contain an invalid summary")
        if item.role == "bam":
            expected_record_count = trailer.mate_count
        elif item.role in ("read1", "read2"):
            expected_record_count = trailer.fragment_count
        elif item.role == "truth.methdb":
            expected_record_count = trailer.methylation_site_count
        else:
            expected_record_count = None
        if expected_record_count is not None and item.record_count != expected_record_count:
            raise ManifestError("output record count disagrees with core evidence")
        if item.record_count < 0:
            raise ManifestError("output record count must be non-negative")
        if item.size_bytes < 0:
            raise ManifestError("output byte size must be non-negative")
        _require_hex_digest("output SHA-256", item.sha256)


def _validate_output_paths(
    config: Mapping[str, Any], outputs: OutputSummary
) -> None:
    output_config = config["output"]
    directory = Path(output_config["directory"])
    prefix = output_config["prefix"]
    for item in outputs.files:
        expected_parent = (
            directory / "truth" if item.role.startswith("truth.") else directory
        )
        if item.path.parent != expected_parent:
            raise ManifestError("output path escaped its normalized output location")
        if not item.path.name.startswith(prefix + "."):
            raise ManifestError("output path disagrees with the normalized prefix")


def _input_manifest_entry(
    file_digest: FileDigest,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    if not isinstance(file_digest, FileDigest):
        raise ManifestError("prepared files contain an invalid digest")
    _require_hex_digest("input SHA-256", file_digest.sha256)
    entry = {
        "format": (
            _methdb_input_format(file_digest.path)
            if file_digest.role == "input.methdb"
            else _input_format(file_digest.role)
        ),
        "path": str(file_digest.path),
        "role": file_digest.role,
        "sha256": file_digest.sha256,
        "size_bytes": file_digest.size_bytes,
    }
    if file_digest.role == "input.methdb":
        entry["format_version"] = METHDB_VERSION
    artifact = _artifact_for_role(config, file_digest.role)
    if artifact is not None:
        entry["declared_sha256"] = artifact["sha256"]
    return entry


def _command_manifest_entry(
    invocation_argv: Sequence[str] | None,
    config: Mapping[str, Any],
) -> dict[str, Any]:
    if invocation_argv is None:
        return {"interface": "python-api"}
    if isinstance(invocation_argv, (str, bytes)):
        raise ManifestError("invocation argv must be a sequence of strings")
    argv = list(invocation_argv)
    if not argv or any(not isinstance(token, str) for token in argv):
        raise ManifestError("CLI invocation argv must contain strings")
    full_argv = _full_run_argv(config, argv)
    return {
        "full_command": shlex.join(full_argv),
        "interface": "cli",
        "user_command": shlex.join(argv),
    }


def _validate_command(
    value: object, config: Mapping[str, Any]
) -> None:
    if not isinstance(value, Mapping):
        raise ManifestError("manifest command section is missing")
    interface = value.get("interface")
    if interface == "python-api":
        if set(value) != {"interface"}:
            raise ManifestError("Python API command section is invalid")
        return
    if interface != "cli":
        raise ManifestError("manifest command interface is invalid")
    command = value.get("user_command")
    full_command = value.get("full_command")
    if not isinstance(command, str):
        raise ManifestError("CLI command section is invalid")
    try:
        argv = shlex.split(command)
    except ValueError as error:
        raise ManifestError("CLI command section is invalid") from error
    if (
        not argv
        or any(not isinstance(token, str) for token in argv)
        or command != shlex.join(argv)
        or not isinstance(full_command, str)
        or full_command != shlex.join(_full_run_argv(config, argv))
        or set(value)
        != {
            "full_command",
            "interface",
            "user_command",
        }
    ):
        raise ManifestError("CLI command section is invalid")


def _full_run_argv(
    config: Mapping[str, Any], argv: Sequence[str]
) -> list[str]:
    try:
        return build_full_run_argv(config, argv)
    except FullCommandError as error:
        raise ManifestError(
            "CLI command cannot form a full command: {}".format(error)
        ) from error


_MODEL_DECLARATION_PATHS = (
    ("coverage",),
    ("sequencing", "quality"),
    ("sequencing", "error"),
)


def _manifest_effective_config(
    normalized: Mapping[str, Any],
) -> dict[str, Any]:
    effective = copy.deepcopy(dict(normalized))
    if "seed" not in effective or "seeds" not in effective:
        raise ManifestError("effective config randomness is incomplete")
    effective.pop("seed")
    effective.pop("seeds")
    _rename_model_discriminators(effective, "kind", "type")
    _serialize_beta_pairs(effective)
    return effective


def _manifest_run_config(
    value: object, sha256_value: object, randomness_value: object
) -> dict[str, Any]:
    if not isinstance(value, Mapping):
        raise ManifestError("manifest configuration section is missing")
    normalized = copy.deepcopy(dict(value))
    if "seed" in normalized or "seeds" in normalized:
        raise ManifestError("manifest configuration duplicates randomness")
    _rename_model_discriminators(normalized, "type", "kind")
    _deserialize_beta_pairs(normalized)

    if not isinstance(randomness_value, Mapping):
        raise ManifestError("manifest randomness section is missing")
    if set(randomness_value) != {
        "master_seed",
        "methylation_seed",
        "mutation_seed",
        "phasing_seed",
        "stages",
    }:
        raise ManifestError("manifest randomness section is invalid")
    master_seed = randomness_value.get("master_seed")
    stage_seeds = {
        "methylation": randomness_value.get("methylation_seed"),
        "mutation": randomness_value.get("mutation_seed"),
        "phasing": randomness_value.get("phasing_seed"),
    }
    stages = randomness_value.get("stages")
    if (
        not _is_seed_text(master_seed)
        or any(not _is_seed_text(seed) for seed in stage_seeds.values())
        or stages != list(STAGE_NAMES)
    ):
        raise ManifestError("manifest randomness section is invalid")
    normalized["seed"] = master_seed
    normalized["seeds"] = copy.deepcopy(dict(stage_seeds))

    _require_hex_digest("configuration SHA-256", sha256_value)
    try:
        config_json = json.dumps(
            normalized,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError, UnicodeError) as error:
        raise ManifestError(
            "manifest effective configuration is not JSON-compatible"
        ) from error
    observed_sha256 = hashlib.sha256(config_json.encode("utf-8")).hexdigest()
    if observed_sha256 != sha256_value:
        raise ManifestError("manifest configuration SHA-256 mismatch")
    return normalized


def _rename_model_discriminators(
    config: dict[str, Any], source: str, destination: str
) -> None:
    for path in _MODEL_DECLARATION_PATHS:
        declaration: object = config
        for component in path:
            if not isinstance(declaration, dict):
                raise ManifestError(
                    "manifest model declaration is invalid: {}".format(
                        ".".join(path)
                    )
                )
            declaration = declaration.get(component)
        if (
            not isinstance(declaration, dict)
            or source not in declaration
            or destination in declaration
        ):
            raise ManifestError(
                "manifest model discriminator is invalid: {}".format(
                    ".".join(path)
                )
            )
        declaration[destination] = declaration.pop(source)


def _serialize_beta_pairs(config: dict[str, Any]) -> None:
    beta = _beta_section(config)
    for context in ("CG", "CHG", "CHH"):
        pair = beta.get(context)
        if not isinstance(pair, list) or len(pair) != 2:
            raise ManifestError("manifest beta declaration is invalid")
        beta[context] = ",".join(_number_text(value) for value in pair)


def _deserialize_beta_pairs(config: dict[str, Any]) -> None:
    beta = _beta_section(config)
    for context in ("CG", "CHG", "CHH"):
        value = beta.get(context)
        if not isinstance(value, str):
            raise ManifestError("manifest beta declaration is invalid")
        tokens = value.split(",")
        if len(tokens) != 2:
            raise ManifestError("manifest beta declaration is invalid")
        try:
            pair = [json.loads(token) for token in tokens]
        except json.JSONDecodeError as error:
            raise ManifestError("manifest beta declaration is invalid") from error
        if (
            any(
                isinstance(number, bool)
                or not isinstance(number, (int, float))
                or not math.isfinite(number)
                for number in pair
            )
            or ",".join(_number_text(number) for number in pair) != value
        ):
            raise ManifestError("manifest beta declaration is invalid")
        beta[context] = pair


def _beta_section(config: dict[str, Any]) -> dict[str, Any]:
    methylation = config.get("methylation")
    beta = methylation.get("beta") if isinstance(methylation, dict) else None
    if not isinstance(beta, dict):
        raise ManifestError("manifest beta section is missing")
    return beta


def _number_text(value: object) -> str:
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        raise ManifestError("manifest beta value is invalid")
    return json.dumps(value, allow_nan=False, separators=(",", ":"))


def _is_seed_text(value: object) -> bool:
    if not isinstance(value, str) or not value or not value.isdigit():
        return False
    if len(value) > 1 and value.startswith("0"):
        return False
    return int(value, 10) <= (1 << 64) - 1


def _input_format(role: str) -> str:
    return {
        "reference": "fasta",
        "input.vcf": "vcf",
        "input.cgmap": "cgmap",
        "input.bed_methyl": "bedMethyl",
        "input.asm": "asm",
        "input.asm_bed": "asm-bed",
        "input.tbs-bed": "bed",
        "model.coverage": "tsv",
        "model.quality": "json",
        "model.error": "json",
    }.get(role, "model")


def _methdb_input_format(path: Path) -> str:
    try:
        with path.open("rb") as input_file:
            prefix = input_file.read(len(METHDB_MAGIC) + 1)
    except OSError as error:
        raise ManifestError("cannot inspect MethDB format: {}".format(error)) from error
    if prefix[: len(METHDB_MAGIC)] != METHDB_MAGIC:
        raise ManifestError("input MethDB has an unknown format magic")
    if len(prefix) != len(METHDB_MAGIC) + 1 or prefix[-1] != METHDB_VERSION:
        raise ManifestError("input MethDB has an unsupported format version")
    return "methdb"


def _artifact_for_role(
    config: Mapping[str, Any], role: str
) -> Mapping[str, Any] | None:
    containers = {
        "model.coverage": config["coverage"],
        "model.quality": config["sequencing"]["quality"],
        "model.error": config["sequencing"]["error"],
    }
    container = containers.get(role)
    if not isinstance(container, Mapping):
        return None
    artifact = container.get("artifact")
    return artifact if isinstance(artifact, Mapping) else None


def _require_hex_digest(name: str, value: str) -> None:
    if (
        not isinstance(value, str)
        or len(value) != 64
        or any(character not in "0123456789abcdef" for character in value)
    ):
        raise ManifestError("{} is invalid".format(name))


def _canonical_json(document: Mapping[str, Any]) -> str:
    try:
        return json.dumps(
            document,
            allow_nan=False,
            ensure_ascii=False,
            indent=2,
            sort_keys=True,
        )
    except (TypeError, ValueError, UnicodeError) as error:
        raise ManifestError("manifest is not canonical JSON data") from error
