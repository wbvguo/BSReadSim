"""Validate the publishable contents and metadata of built distributions."""

from __future__ import annotations

from email.parser import BytesParser
from email.policy import default
from pathlib import Path, PurePosixPath
import sys
import tarfile
import zipfile


def _only(directory: Path, pattern: str, label: str) -> Path:
    matches = sorted(directory.glob(pattern))
    if len(matches) != 1:
        raise SystemExit(
            "expected exactly one {} in {}, found {}".format(
                label, directory, len(matches)
            )
        )
    return matches[0]


def _require_members(
    archive: Path,
    members: set[str],
    expected: set[str],
) -> None:
    missing = sorted(expected - members)
    if missing:
        raise SystemExit(
            "{} is missing required members: {}".format(
                archive.name, ", ".join(missing)
            )
        )


def _check_wheel(wheel: Path) -> None:
    with zipfile.ZipFile(wheel) as archive:
        names = set(archive.namelist())
        _require_members(
            wheel,
            names,
            {
                "bsreadsim/_licenses/HTSlib.txt",
                "bsreadsim/_licenses/htscodecs.txt",
                "bsreadsim/data/example/test.fa",
                "bsreadsim/data/registry.json",
                "bsreadsim/htsim/htsim-core",
                "bsreadsim/run-config.schema.json",
            },
        )

        extension_members = sorted(
            name
            for name in names
            if name.startswith("bsreadsim/_cext") and name.endswith(".so")
        )
        if len(extension_members) != 1:
            raise SystemExit(
                "{} must contain exactly one compiled _cext shared object".format(
                    wheel.name
                )
            )

        leaked_sources = sorted(
            name
            for name in names
            if name.startswith("bsreadsim/cext/")
            and PurePosixPath(name).suffix in {".c", ".h"}
        )
        if leaked_sources:
            raise SystemExit(
                "{} leaked C extension sources: {}".format(
                    wheel.name, ", ".join(leaked_sources)
                )
            )

        metadata_members = sorted(
            name for name in names if name.endswith(".dist-info/METADATA")
        )
        if len(metadata_members) != 1:
            raise SystemExit(
                "{} must contain exactly one METADATA file".format(wheel.name)
            )
        metadata = BytesParser(policy=default).parsebytes(
            archive.read(metadata_members[0])
        )
        if metadata.get("License-Expression") != "MIT":
            raise SystemExit(
                "{} does not declare License-Expression: MIT".format(wheel.name)
            )
        license_files = metadata.get_all("License-File", [])
        if "LICENSE" not in license_files:
            raise SystemExit(
                "{} does not declare LICENSE in core metadata".format(wheel.name)
            )
        installed_licenses = [
            name for name in names if name.endswith(".dist-info/licenses/LICENSE")
        ]
        if len(installed_licenses) != 1:
            raise SystemExit(
                "{} must install exactly one project LICENSE".format(wheel.name)
            )


def _check_sdist(sdist: Path) -> None:
    with tarfile.open(sdist, mode="r:gz") as archive:
        raw_names = archive.getnames()

    roots = {
        PurePosixPath(name).parts[0]
        for name in raw_names
        if PurePosixPath(name).parts
    }
    if len(roots) != 1:
        raise SystemExit(
            "{} must have one top-level directory".format(sdist.name)
        )
    root = next(iter(roots))
    names = {
        PurePosixPath(*PurePosixPath(name).parts[1:]).as_posix()
        for name in raw_names
        if PurePosixPath(name).parts[0] == root
        and len(PurePosixPath(name).parts) > 1
    }
    _require_members(
        sdist,
        names,
        {
            ".gitmodules",
            "CMakeLists.txt",
            "LICENSE",
            "MANIFEST.in",
            "README.md",
            "data/example/test.fa",
            "data/registry.json",
            "htsim/src/core.cpp",
            "htslib/LICENSE",
            "htslib/htscodecs/LICENSE.md",
            "pyproject.toml",
            "setup.py",
            "src/bsreadsim/cext/api.h",
            "src/bsreadsim/cext/module.c",
            "src/bsreadsim/cext/protocol.c",
            "src/bsreadsim/cext/sam.c",
        },
    )


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit(
            "usage: check_distribution_contents.py DISTRIBUTION_DIRECTORY"
        )
    directory = Path(argv[1]).resolve()
    if not directory.is_dir():
        raise SystemExit("distribution directory does not exist: {}".format(directory))
    sdist = _only(directory, "*.tar.gz", "source distribution")
    wheel = _only(directory, "*.whl", "wheel")
    _check_sdist(sdist)
    _check_wheel(wheel)
    print("validated {} and {}".format(sdist.name, wheel.name))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
