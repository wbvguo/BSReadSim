"""Locate the private C++ core without importing the simulation pipeline."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
from typing import Optional, Union


PathLike = Union[str, os.PathLike]
CORE_FILENAME = "htsim-core.exe" if os.name == "nt" else "htsim-core"


class CoreExecutableError(RuntimeError):
    """The requested or installed generation core is unavailable."""


def packaged_core_candidate() -> Path:
    """Return the expected core location without asserting that it exists."""
    return Path(__file__).resolve().parent / "_bin" / CORE_FILENAME


def resolve_core_executable(value: Optional[PathLike] = None) -> Path:
    """Resolve the packaged core, a PATH command, or an explicit override."""
    if value is None:
        packaged = packaged_core_candidate()
        if packaged.exists():
            candidate = packaged
        else:
            discovered = shutil.which("htsim-core")
            if discovered is None:
                raise CoreExecutableError(
                    "cannot find the packaged htsim-core or one on PATH; "
                    "pass --core explicitly"
                )
            candidate = Path(discovered)
    else:
        try:
            text = os.fspath(value)
        except TypeError as error:
            raise CoreExecutableError("core executable must be a text path") from error
        if not isinstance(text, str) or not text or "\x00" in text:
            raise CoreExecutableError(
                "core executable must be a non-empty text path"
            )
        path_candidate = Path(text).expanduser()
        if os.path.dirname(text) or path_candidate.exists():
            candidate = path_candidate
        else:
            discovered = shutil.which(text)
            candidate = Path(discovered) if discovered is not None else Path(text)

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
    "CoreExecutableError",
    "packaged_core_candidate",
    "resolve_core_executable",
]
