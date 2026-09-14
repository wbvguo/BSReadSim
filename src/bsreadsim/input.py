"""The input files of a command: which settings name them, how the manifest
records them (format, size, SHA-256), and their validation by htsim."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path

from .htsim import HtsimCore
from .settings import Settings

# Settings fields that name input files -> their format in the manifest
INPUT_FORMATS = {
    "reference": "fasta",
    "vcf": "vcf",
    "cgmap": "cgmap",
    "bed_methyl": "bedMethyl",
    "methbg": "MethBG",
    "methbed": "MethBED",
    "methdb": "methdb",
    "asm": "asm",
    "asm_bed": "asm-bed",
    "gc_profile": "gc-profile",
    "rrbs_candidates": "rrbs-candidates",
    "targets": "bed",
    "quality_model": "json",
    "error_model": "json",
}
# Methylation profile names in the validation summary -> settings fields
_PROFILE_FIELDS = {"cgmap": "cgmap", "bedmethyl": "bed_methyl", "methbg": "methbg",
                   "methbed": "methbed"}


def input_files(settings: Settings) -> dict[str, Path]:
    """The input files of ``settings`` by field, in a fixed order."""
    return {name: getattr(settings, name) for name in INPUT_FORMATS
            if getattr(settings, name) is not None}


def describe_inputs(settings: Settings) -> list[dict]:
    """The manifest's record of each input file: role, format, path, size, and SHA-256."""
    return [{"role": name, "format": INPUT_FORMATS[name], "path": str(path),
             "size_bytes": path.stat().st_size, "sha256": sha256_file(path)}
            for name, path in input_files(settings).items()]


def validate_inputs(settings: Settings, core: HtsimCore) -> dict:
    """htsim's check of the inputs against the reference, with each file's SHA-256."""
    summary = core.validate(settings)
    files = {"reference": settings.reference, "vcf": settings.vcf}
    if summary["methylation"]:
        files["methylation"] = getattr(settings, _PROFILE_FIELDS[summary["methylation"]["format"]])
    if summary["asm"]:
        files["asm"] = settings.asm_bed or settings.asm
    for key, path in files.items():
        if summary.get(key) is not None:
            summary[key]["sha256"] = sha256_file(path)
    return summary


def sha256_file(path: os.PathLike | str) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()
