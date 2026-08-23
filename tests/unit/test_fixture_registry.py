"""Integrity checks for the flat test-fixture registry."""

from __future__ import annotations

import json
from pathlib import Path
import unittest


FIXTURE_ROOT = Path(__file__).resolve().parents[1] / "fixtures"
REGISTRY_FILES = {"README.md", "fixtures.json"}


class FixtureRegistryTests(unittest.TestCase):
    def test_registry_covers_every_fixture_exactly_once(self) -> None:
        document = json.loads(
            (FIXTURE_ROOT / "fixtures.json").read_text(encoding="utf-8")
        )
        entries = document["fixtures"]
        registered = [entry["file"] for entry in entries]
        self.assertEqual(len(registered), len(set(registered)))

        files = {
            path.name
            for path in FIXTURE_ROOT.iterdir()
            if path.is_file() and path.name not in REGISTRY_FILES
        }
        self.assertEqual(set(registered), files)

        for entry in entries:
            with self.subTest(file=entry["file"]):
                self.assertIs(type(entry.get("valid")), bool)
                self.assertTrue(entry.get("purpose"))
                if not entry["valid"]:
                    self.assertTrue(entry.get("expected_result"))


if __name__ == "__main__":
    unittest.main()
