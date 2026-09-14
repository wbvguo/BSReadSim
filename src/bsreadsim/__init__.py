"""BSReadSim simulates bisulfite and standard sequencing reads.

Use the ``bsreadsim`` command, or ``simulate`` with the same arguments.
"""

from __future__ import annotations

from collections.abc import Sequence
import os

from ._version import __version__


def simulate(arguments: Sequence[str], core: str | os.PathLike | None = None):
    """Run ``bsreadsim run ARGUMENTS...``; returns its ``pipeline.RunResult``.

    With more than one thread the reads are made in worker processes, which import the
    calling script again: call this under ``if __name__ == "__main__":``.
    """
    from .cli import build_parser, from_arguments   # here, so that importing bsreadsim stays light
    from .pipeline import run

    argv = ["run", *arguments]
    return run(from_arguments(build_parser().parse_args(argv)), ["bsreadsim", *argv], core)


__all__ = ["__version__", "simulate"]
