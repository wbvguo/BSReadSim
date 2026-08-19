"""RRBS candidate-catalog exchange driven by the native generator."""

from __future__ import annotations

from pathlib import Path
import subprocess
from typing import Mapping, Optional, Union
import os
import uuid

from .config import normalize_run_config
from .core_argv import build_core_argv
from .preparation import prepare_run
from .runtime import CoreExecutableError, resolve_core_executable


PathLike = Union[str, os.PathLike]


class CatalogError(RuntimeError):
    """A native candidate catalog could not be exported safely."""


def export_rrbs_catalog(
    document: Mapping[str, object],
    output_path: PathLike,
    *,
    base_directory: PathLike = ".",
    core_executable: Optional[PathLike] = None,
) -> Path:
    """Export the exact native RRBS candidate table to a new plain BED file.

    ``document`` is an in-memory projection built from direct CLI arguments;
    this path neither reads nor writes a JSON configuration file.  It uses the
    same normalization and native argv boundary as simulation, and it never
    generates or embeds a candidate hash.
    """

    loaded = normalize_run_config(
        document, base_directory, mode="production"
    )
    normalized = loaded.normalized
    if normalized["technology"] != "RRBS":
        raise CatalogError("candidate export requires RRBS CLI arguments")
    rrbs = normalized["rrbs"]
    coverage = normalized["coverage"]
    if (
        not isinstance(rrbs, Mapping)
        or "candidate_bed" in rrbs
        or not isinstance(coverage, Mapping)
        or coverage["kind"] != "uniform"
    ):
        raise CatalogError(
            "candidate export requires uniform RRBS arguments without a candidate BED"
        )

    inputs = normalized["inputs"]
    mutation = normalized["mutation"]
    if loaded.master_seed is None and (
        "vcf" in inputs or mutation["rate"] > 0
    ):
        raise CatalogError(
            "RRBS candidate export with VCF or de novo variants requires "
            "an explicit --seed"
        )
    effective = loaded
    if effective.master_seed is None:
        # Reference-only candidate geometry is seed independent.  A fixed
        # internal seed avoids manufacturing an unrecorded random identity for
        # this catalog-only process.
        effective = effective.with_master_seed(0)

    prepared = prepare_run(effective)
    try:
        executable = resolve_core_executable(core_executable)
    except CoreExecutableError as error:
        raise CatalogError(str(error)) from error
    base_argv = build_core_argv(
        prepared,
        str(uuid.uuid4()),
        executable,
        truth_columns="none",
    )
    argv = (base_argv[0], "rrbs-catalog", *base_argv[1:])

    try:
        destination = Path(output_path).expanduser().resolve(strict=False)
    except (TypeError, OSError) as error:
        raise CatalogError("RRBS candidate output path is invalid") from error
    if not destination.parent.is_dir():
        raise CatalogError(
            "RRBS candidate output parent is not a directory: {}".format(
                destination.parent
            )
        )

    created = False
    try:
        with destination.open("xb") as output:
            created = True
            completed = subprocess.run(
                argv,
                stdout=output,
                stderr=subprocess.PIPE,
                check=False,
            )
        if completed.returncode != 0:
            detail = completed.stderr.decode("utf-8", errors="replace").strip()
            raise CatalogError(
                "native RRBS candidate export failed{}".format(
                    ": " + detail if detail else ""
                )
            )
    except CatalogError:
        if created:
            try:
                destination.unlink()
            except OSError:
                pass
        raise
    except FileExistsError as error:
        raise CatalogError(
            "RRBS candidate output already exists: {}".format(destination)
        ) from error
    except OSError as error:
        if created:
            try:
                destination.unlink()
            except OSError:
                pass
        raise CatalogError(
            "cannot write RRBS candidate output {}: {}".format(
                destination, error
            )
        ) from error

    return destination


__all__ = ["CatalogError", "export_rrbs_catalog"]
