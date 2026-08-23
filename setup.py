"""Setuptools build hook for the private C++ generation core.

The public package remains Python-owned, but a regular wheel is only usable
when it carries the matching ``htsim-core`` executable.  The executable is
built in setuptools' temporary build tree and copied into the wheel package;
the source checkout is never modified by this hook.
"""

from pathlib import Path
import os
import shutil
import stat
import subprocess

from setuptools import Extension, setup
from setuptools.command.build_py import build_py
from setuptools.dist import Distribution
from setuptools.errors import CompileError


SOURCE_ROOT = Path(__file__).resolve().parent
CORE_TARGET = "htsim_core"
CORE_FILENAME = "htsim-core.exe" if os.name == "nt" else "htsim-core"


class BinaryDistribution(Distribution):
    """Mark wheels as platform-specific because they contain compiled code."""

    def has_ext_modules(self):
        return True


class BuildPythonWithCore(build_py):
    """Build and stage the C++ core after copying the Python package."""

    def run(self):
        self._reset_package_build_directory()
        super().run()
        cmake = shutil.which("cmake")
        if cmake is None:
            raise CompileError(
                "building bsreadsim requires CMake 3.20+, a C++17 compiler, "
                "zlib development files, Autoconf, Make, and initialized "
                "recursive Git submodules"
            )

        build_command = self.get_finalized_command("build")
        core_build_root = Path(build_command.build_temp).resolve() / "htsim-core"
        build_directory = core_build_root / "build"
        install_directory = core_build_root / "install"
        build_directory.mkdir(parents=True, exist_ok=True)
        install_directory.mkdir(parents=True, exist_ok=True)

        self._run_cmake(
            (
                cmake,
                "-S",
                str(SOURCE_ROOT),
                "-B",
                str(build_directory),
                "-DCMAKE_BUILD_TYPE=Release",
                "-DBUILD_TESTING=OFF",
                "-DCMAKE_INSTALL_PREFIX={}".format(install_directory),
            ),
            "configure",
        )
        self._run_cmake(
            (
                cmake,
                "--build",
                str(build_directory),
                "--config",
                "Release",
                "--target",
                CORE_TARGET,
                "--parallel",
            ),
            "build",
        )
        self._run_cmake(
            (
                cmake,
                "--install",
                str(build_directory),
                "--config",
                "Release",
                "--prefix",
                str(install_directory),
            ),
            "install",
        )

        built_core = install_directory / "bin" / CORE_FILENAME
        if not built_core.is_file() or built_core.stat().st_size == 0:
            raise CompileError(
                "CMake did not install the expected core: {}".format(built_core)
            )
        package_core = (
            Path(self.build_lib) / "bsreadsim" / "htsim" / CORE_FILENAME
        )
        package_core.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(built_core), str(package_core))
        if os.name != "nt":
            package_core.chmod(
                package_core.stat().st_mode
                | stat.S_IXUSR
                | stat.S_IXGRP
                | stat.S_IXOTH
            )

        license_directory = Path(self.build_lib) / "bsreadsim" / "_licenses"
        license_directory.mkdir(parents=True, exist_ok=True)
        third_party_licenses = (
            (SOURCE_ROOT / "htslib" / "LICENSE", "HTSlib.txt"),
            (
                SOURCE_ROOT
                / "htslib"
                / "htscodecs"
                / "LICENSE.md",
                "htscodecs.txt",
            ),
        )
        for source, filename in third_party_licenses:
            if not source.is_file():
                raise CompileError(
                    "missing third-party license; initialize recursive Git "
                    "submodules: {}".format(source)
                )
            shutil.copy2(str(source), str(license_directory / filename))

    def _reset_package_build_directory(self):
        """Prevent removed package files from leaking out of a reused build."""
        build_lib = Path(self.build_lib).resolve()
        package_directory = Path(self.build_lib) / "bsreadsim"
        source_package = (SOURCE_ROOT / "src" / "bsreadsim").resolve()
        if package_directory.is_symlink():
            raise CompileError(
                "refusing to replace symlinked package staging directory: {}".format(
                    package_directory
                )
            )
        resolved_package = package_directory.resolve(strict=False)
        if (
            resolved_package == source_package
            or resolved_package.parent != build_lib
            or resolved_package.name != "bsreadsim"
        ):
            raise CompileError(
                "refusing to clean unsafe package staging directory: {}".format(
                    resolved_package
                )
            )
        if package_directory.exists():
            if not package_directory.is_dir():
                raise CompileError(
                    "package staging path is not a directory: {}".format(
                        package_directory
                    )
                )
            self.announce(
                "resetting package staging directory: {}".format(
                    package_directory
                ),
                level=2,
            )
            shutil.rmtree(str(package_directory))

    def _run_cmake(self, argv, phase):
        self.announce("CMake {}: {}".format(phase, " ".join(argv)), level=2)
        try:
            subprocess.run(argv, cwd=str(SOURCE_ROOT), check=True)
        except (OSError, subprocess.CalledProcessError) as error:
            raise CompileError(
                "CMake {} failed while building htsim-core: {}".format(
                    phase, error
                )
            ) from error


setup(
    cmdclass={"build_py": BuildPythonWithCore},
    distclass=BinaryDistribution,
    ext_modules=[
        Extension(
            "bsreadsim._cext",
            depends=["src/bsreadsim/cext/api.h"],
            sources=[
                "src/bsreadsim/cext/module.c",
                "src/bsreadsim/cext/protocol.c",
                "src/bsreadsim/cext/sam.c",
            ],
        )
    ],
)
