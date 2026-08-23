"""Construction and verification of the final BSReadSim run manifest."""

from __future__ import annotations

import copy
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from collections.abc import Mapping
from typing import Any

from .. import __version__
from ..output.bam import BAM_CONTRACT, BAM_MAPQ
from .config import RUN_CONFIG_SCHEMA_VERSION
from ..output import OutputFileSummary, OutputSummary
from .prepare import FileDigest, PreparedRun
from ..htsim.protocol import (
    AmbiguityPolicy,
    BaseEncoding,
    PROTOCOL_MAJOR,
    PROTOCOL_MINOR,
    Header,
    Technology,
    Trailer,
)
from ..process.batch import READ_NAME_CONTRACT
from ..rng import RNG_CONTRACT, STAGE_NAMES


MANIFEST_SCHEMA_VERSION = "1.1"
COORDINATE_CONVENTION = "0-based-half-open"


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
) -> CompleteManifest:
    """Cross-check core/Python evidence and build the final commit marker."""
    if not isinstance(prepared, PreparedRun):
        raise ManifestError("prepared must be a PreparedRun")
    config = prepared.config
    if config.master_seed is None:
        raise ManifestError("the effective config must contain a materialized seed")
    if not isinstance(header, Header) or not isinstance(trailer, Trailer):
        raise ManifestError("header and trailer must use the protocol contract")
    protocol_version = (PROTOCOL_MAJOR, PROTOCOL_MINOR)
    validate_header_projection(prepared, header)
    if not isinstance(outputs, OutputSummary):
        raise ManifestError("outputs must be an OutputSummary")

    _validate_counts(config.normalized, header, trailer, outputs)
    _validate_output_paths(config.normalized, outputs)

    document = {
        "config": {
            "normalized": config.as_dict(),
            "sha256": config.sha256,
        },
        "contigs": [
            {
                "index": index,
                "length": contig.length,
                "name": contig.name,
                "reference_sequence_sha256": contig.reference_sha256.hex(),
            }
            for index, contig in enumerate(header.contigs)
        ],
        "counts": {
            "core": {
                "fragment_count": trailer.fragment_count,
                "mate_count": trailer.mate_count,
                "methylation_site_count": trailer.methylation_site_count,
                "per_contig_fragment_counts": list(
                    trailer.per_contig_fragment_counts
                ),
                "skipped_fragment_count": trailer.skipped_fragment_count,
                "template_base_count": trailer.template_base_count,
            },
            "python": {
                "fragment_count": outputs.fragment_count,
                "mate_count": outputs.mate_count,
                "records_by_role": {
                    item.role: item.record_count for item in outputs.files
                },
            },
        },
        "inputs": [
            _input_manifest_entry(file_digest, config.normalized)
            for file_digest in prepared.files
        ],
        "manifest_schema_version": MANIFEST_SCHEMA_VERSION,
        "models": {
            "methylation_state": {
                "contract": "bernoulli-site-v1",
                "effective": "bernoulli",
                "requested": config.normalized["methylation"]["state_model"],
            }
        },
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
        "randomness": {
            "catalog_seed": config.normalized["methylation"]["catalog_seed"],
            "contract": RNG_CONTRACT,
            "master_seed": str(config.master_seed),
            "stages": list(STAGE_NAMES),
        },
        "run_id": header.run_id,
        "status": "complete",
        "stream_sha256": trailer.stream_sha256.hex(),
        "versions": {
            "config_schema": RUN_CONFIG_SCHEMA_VERSION,
            "core": header.core_version,
            "manifest": MANIFEST_SCHEMA_VERSION,
            "protocol": "{}.{}".format(*protocol_version),
            "python": __version__,
            "read_name": READ_NAME_CONTRACT,
            "rng": RNG_CONTRACT,
        },
        "reproducibility": {
            "numerical_tolerance_exceptions": [],
            "scope": (
                "same released core/Python versions, normalized config, "
                "input/model digests, and master seed"
            ),
        },
    }
    if config.normalized["output"]["bam"]:
        fragment_summary = bool(
            config.normalized["output"]["fragment_summary"]
        )
        fragment_realization = bool(
            config.normalized["output"]["fragment_realization"]
        )
        document["annotation_alignment"] = {
            "contract": BAM_CONTRACT,
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
                    "schema": "packed-b64url-v1",
                    "scope": "complete-physical-fragment-realization",
                },
                "zf": {
                    "required": fragment_summary,
                    "schema": "u16x12-v1",
                    "scope": "complete-physical-fragment",
                },
                "zr": {
                    "required": True,
                    "schema": "u16x12-v1",
                    "scope": "single-read",
                },
                "zt": {
                    "alphabet": "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789-_",
                    "required": True,
                    "schema": "state64-v1",
                    "scope": "one-character-per-bam-seq-base",
                },
            },
        }
        document["versions"]["bam"] = BAM_CONTRACT

    identity_json = _canonical_json(document)
    reproducibility_digest = hashlib.sha256(
        identity_json.encode("utf-8")
    ).hexdigest()
    document["reproducibility"]["sha256"] = reproducibility_digest
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
    expected_has_details = bool(normalized_output["bam"])
    _validate_protocol_projection(prepared, header)

    normalized = prepared.config.normalized
    technology = {
        "WGBS": Technology.WGBS,
        "RRBS": Technology.RRBS,
        "TBS": Technology.TBS,
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
    reproducibility = document.get("reproducibility")
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
    identity_reproducibility = identity.get("reproducibility")
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
    if header.config_schema_version != RUN_CONFIG_SCHEMA_VERSION:
        raise ManifestError("core config-schema version disagrees with Python")
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

    bam = bool(config["output"]["bam"])
    expected_roles = set()
    if not bam:
        expected_roles.add("read1")
        if paired_end:
            expected_roles.add("read2")
    if bam:
        expected_roles.add("bam")
    observed_roles = {item.role for item in outputs.files}
    if observed_roles != expected_roles or len(outputs.files) != len(expected_roles):
        raise ManifestError("output roles disagree with SE/PE configuration")
    for item in outputs.files:
        if not isinstance(item, OutputFileSummary):
            raise ManifestError("output files contain an invalid summary")
        expected_record_count = (
            trailer.mate_count
            if item.role == "bam"
            else trailer.fragment_count
        )
        if item.record_count != expected_record_count:
            raise ManifestError("output record count disagrees with core fragments")
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
        if item.path.parent != directory:
            raise ManifestError("output path escaped the normalized output directory")
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
        "format": _input_format(file_digest.role),
        "path": str(file_digest.path),
        "role": file_digest.role,
        "sha256": file_digest.sha256,
        "size_bytes": file_digest.size_bytes,
    }
    artifact = _artifact_for_role(config, file_digest.role)
    if artifact is not None:
        entry["format"] = artifact["format"]
        entry["model_version"] = artifact["version"]
        entry["declared_sha256"] = artifact["sha256"]
    return entry


def _input_format(role: str) -> str:
    return {
        "reference": "fasta",
        "input.vcf": "vcf",
        "input.cgmap": "cgmap",
        "input.bed_methyl": "bedMethyl",
        "input.methdb": "bsreadsim-methdb-v1",
        "input.asm": "asm",
        "input.asm_bed": "asm-bed",
        "input.tbs-bed": "bed",
    }.get(role, "model")


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
            separators=(",", ":"),
            sort_keys=True,
        )
    except (TypeError, ValueError, UnicodeError) as error:
        raise ManifestError("manifest is not canonical JSON data") from error
