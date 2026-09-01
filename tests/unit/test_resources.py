"""Contracts for immutable data bundled with BSReadSim releases."""

from __future__ import annotations

import hashlib
from pathlib import Path
import tempfile
import unittest

from bsreadsim.resources import (
    ResourceError,
    copy_resource,
    get_resource,
    list_resources,
    read_resource,
)


class BundledResourceTests(unittest.TestCase):
    def test_test_fasta_is_registered_and_verified(self) -> None:
        resources = list_resources()
        names = [resource.name for resource in resources]

        self.assertEqual(names, sorted(names))
        self.assertIn("test-fasta", names)
        metadata = get_resource("test-fasta")
        payload = read_resource(metadata.name)
        sequence = b"".join(
            line for line in payload.splitlines() if not line.startswith(b">")
        )

        self.assertEqual(metadata.kind, "example.reference")
        self.assertEqual(metadata.format, "fasta")
        self.assertEqual(len(payload), metadata.size_bytes)
        self.assertEqual(hashlib.sha256(payload).hexdigest(), metadata.sha256)
        self.assertEqual(len(sequence), 1_961_600)
        self.assertIn(b"CG", sequence)

    def test_copy_is_exact_and_refuses_to_overwrite(self) -> None:
        with tempfile.TemporaryDirectory() as temporary:
            output = Path(temporary) / "nested" / "test.fa"
            copied = copy_resource("test-fasta", output)

            self.assertEqual(copied, output.resolve())
            self.assertEqual(copied.read_bytes(), read_resource("test-fasta"))
            with self.assertRaisesRegex(ResourceError, "already exists"):
                copy_resource("test-fasta", output)

    def test_unknown_resource_is_rejected(self) -> None:
        with self.assertRaisesRegex(ResourceError, "unknown bundled resource"):
            get_resource("missing-resource")


if __name__ == "__main__":
    unittest.main()
