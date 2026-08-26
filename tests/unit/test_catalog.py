"""Catalog exchange container tests."""

import gzip
import io
import unittest

from bsreadsim.run.catalog import _BgzfWriter


class BgzfWriterTests(unittest.TestCase):
    def test_multiblock_stream_is_deterministic_and_gzip_compatible(self) -> None:
        payload = bytes(range(256)) * 500

        def encode() -> bytes:
            output = io.BytesIO()
            with _BgzfWriter(output) as writer:
                writer.write(payload[:17001])
                writer.write(payload[17001:90003])
                writer.write(payload[90003:])
            return output.getvalue()

        first = encode()
        self.assertEqual(first, encode())
        self.assertEqual(gzip.decompress(first), payload)

        offset = 0
        blocks = 0
        while offset < len(first):
            self.assertEqual(first[offset : offset + 2], b"\x1f\x8b")
            block_size = int.from_bytes(
                first[offset + 16 : offset + 18], "little"
            ) + 1
            self.assertGreaterEqual(block_size, 28)
            self.assertLessEqual(block_size, 65536)
            offset += block_size
            blocks += 1
        self.assertEqual(offset, len(first))
        self.assertGreaterEqual(blocks, 5)  # four data blocks plus BGZF EOF


if __name__ == "__main__":
    unittest.main()
