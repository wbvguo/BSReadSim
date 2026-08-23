"""Immutable catalog orchestration through the htsim core."""

from __future__ import annotations
from contextlib import suppress

import os
from pathlib import Path
import subprocess
from collections.abc import Mapping
import uuid

from ..htsim.launch import CoreExecutableError, build_core_argv, resolve_core_executable
from .config import normalize_run_config
from .prepare import prepare_run


PathLike = str | os.PathLike


class CatalogError(RuntimeError):
    """An htsim candidate catalog could not be exported safely."""


def export_rrbs_catalog(
    document: Mapping[str, object],
    output_path: PathLike,
    *,
    base_directory: PathLike = ".",
    core_executable: PathLike | None = None,
) -> Path:
    """Export the exact htsim RRBS candidate table to a new plain BED file."""

    loaded = normalize_run_config(document, base_directory)
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

    effective = loaded
    if effective.master_seed is None:
        effective = effective.with_master_seed(0)

    prepared = prepare_run(effective)
    return _export_catalog(
        prepared,
        output_path,
        subcommand="rrbs-catalog",
        label="RRBS candidate",
        core_executable=core_executable,
    )


def export_methdb_catalog(
    document: Mapping[str, object],
    output_path: PathLike,
    *,
    base_directory: PathLike = ".",
    core_executable: PathLike | None = None,
) -> Path:
    """Export the exact fixed site-probability catalog used by a run."""

    loaded = normalize_run_config(document, base_directory)
    if "methdb" in loaded.normalized["inputs"]:
        raise CatalogError("cannot save a MethDB while loading another MethDB")
    effective = loaded if loaded.master_seed is not None else loaded.with_master_seed(0)
    return _export_catalog(
        prepare_run(effective),
        output_path,
        subcommand="methdb-catalog",
        label="MethDB",
        core_executable=core_executable,
    )


def _export_catalog(
    prepared,
    output_path: PathLike,
    *,
    subcommand: str,
    label: str,
    core_executable: PathLike | None,
) -> Path:
    try:
        executable = resolve_core_executable(core_executable)
    except CoreExecutableError as error:
        raise CatalogError(str(error)) from error
    base_argv = build_core_argv(
        prepared,
        str(uuid.uuid4()),
        executable,
        emit_details=False,
    )
    argv = (base_argv[0], subcommand, *base_argv[1:])

    try:
        destination = Path(output_path).expanduser().resolve(strict=False)
    except (TypeError, OSError) as error:
        raise CatalogError("{} output path is invalid".format(label)) from error
    if not destination.parent.is_dir():
        raise CatalogError(
            "{} output parent is not a directory: {}".format(label, destination.parent)
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
                "htsim {} export failed{}".format(
                    label,
                    ": " + detail if detail else ""
                )
            )
    except CatalogError:
        if created:
            with suppress(OSError):
                destination.unlink()
        raise
    except FileExistsError as error:
        raise CatalogError(
            "{} output already exists: {}".format(label, destination)
        ) from error
    except OSError as error:
        if created:
            with suppress(OSError):
                destination.unlink()
        raise CatalogError(
            "cannot write {} output {}: {}".format(label, destination, error)
        ) from error

    return destination


__all__ = ["CatalogError", "export_methdb_catalog", "export_rrbs_catalog"]
