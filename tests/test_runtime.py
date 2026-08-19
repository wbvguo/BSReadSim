"""Tests for installed and explicitly overridden core discovery."""

from __future__ import annotations

import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPOSITORY_ROOT / "src"))

from bsreadsim.runtime import (  # noqa: E402
    CoreExecutableError,
    resolve_core_executable,
)


class RuntimeTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.directory = Path(self.temporary_directory.name).resolve()

    def executable(self, name: str) -> Path:
        path = self.directory / name
        path.write_bytes(b"#!/bin/sh\nexit 0\n")
        path.chmod(0o755)
        return path

    def test_packaged_core_precedes_path_discovery(self) -> None:
        packaged = self.executable("packaged-core")
        with mock.patch(
            "bsreadsim.runtime.packaged_core_candidate", return_value=packaged
        ), mock.patch("bsreadsim.runtime.shutil.which") as which:
            self.assertEqual(resolve_core_executable(), packaged)
        which.assert_not_called()

    def test_missing_packaged_core_uses_path(self) -> None:
        discovered = self.executable("path-core")
        with mock.patch(
            "bsreadsim.runtime.packaged_core_candidate",
            return_value=self.directory / "missing-packaged-core",
        ), mock.patch(
            "bsreadsim.runtime.shutil.which", return_value=str(discovered)
        ):
            self.assertEqual(resolve_core_executable(), discovered)

    def test_explicit_override_and_invalid_core_fail_closed(self) -> None:
        self.assertEqual(
            resolve_core_executable(sys.executable),
            Path(sys.executable).resolve(),
        )
        with self.assertRaisesRegex(CoreExecutableError, "cannot resolve"):
            resolve_core_executable(self.directory / "missing-core")
        non_executable = self.directory / "not-executable"
        non_executable.write_bytes(b"not executable")
        if os.name != "nt":
            # Some WSL test environments inherit a Windows TEMP directory on
            # drvfs, where chmod cannot clear the synthetic executable bits.
            # Mock the permission probe so this remains a unit test of the
            # fail-closed resolver boundary rather than the host mount mode.
            with mock.patch(
                "bsreadsim.runtime.os.access", return_value=False
            ) as access:
                with self.assertRaisesRegex(
                    CoreExecutableError, "not executable"
                ):
                    resolve_core_executable(non_executable)
            access.assert_called_once_with(
                str(non_executable.resolve()), os.X_OK
            )

    def test_absent_default_reports_both_discovery_boundaries(self) -> None:
        with mock.patch(
            "bsreadsim.runtime.packaged_core_candidate",
            return_value=self.directory / "missing-packaged-core",
        ), mock.patch("bsreadsim.runtime.shutil.which", return_value=None):
            with self.assertRaisesRegex(
                CoreExecutableError, "packaged htsim-core or one on PATH"
            ):
                resolve_core_executable()
