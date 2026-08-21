"""Cross-check the C++ protocol writer against frozen Python fixtures."""

from pathlib import Path
import subprocess
import sys

from bsreadsim.native.protocol import CoreReportedError
from stream_support import read_stream


PREAMBLE_SIZE = 16
FRAME_ENVELOPE_SIZE = 16
CRC_SIZE = 4


def fixture(root: Path, name: str) -> bytes:
    return bytes.fromhex((root / (name + ".hex")).read_text(encoding="ascii"))


def frame_bounds(stream: bytes):
    bounds = []
    cursor = PREAMBLE_SIZE
    while cursor < len(stream):
        if len(stream) - cursor < FRAME_ENVELOPE_SIZE:
            raise SystemExit("C++ stream ends inside a frame envelope")
        payload_length = int.from_bytes(stream[cursor : cursor + 4], "little")
        payload_start = cursor + FRAME_ENVELOPE_SIZE
        payload_end = payload_start + payload_length
        end = payload_end + CRC_SIZE
        if end > len(stream):
            raise SystemExit("C++ stream ends inside a frame")
        bounds.append((cursor, payload_start, payload_end, end))
        cursor = end
    return bounds


def main(argv):
    if len(argv) != 4:
        raise SystemExit(
            "usage: check_cpp_protocol.py WRITER_TEST OUTPUT_PREFIX FIXTURE_ROOT"
        )
    writer = Path(argv[1]).resolve(strict=True)
    output_prefix = Path(argv[2]).resolve(strict=False)
    fixture_root = Path(argv[3]).resolve(strict=True)
    subprocess.run(
        [str(writer), "--write-fixtures", str(output_prefix)],
        check=True,
    )
    none = Path(str(output_prefix) + "-none.bin").read_bytes()
    full = Path(str(output_prefix) + "-full.bin").read_bytes()
    error = Path(str(output_prefix) + "-error.bin").read_bytes()

    expected_none = b"".join(
        fixture(fixture_root, name)
        for name in ("header-none", "batch-none", "trailer-none")
    )
    if none != expected_none:
        raise SystemExit("C++ no-Details stream differs from frozen Python bytes")
    if error != fixture(fixture_root, "header-none") + fixture(
        fixture_root, "error"
    ):
        raise SystemExit("C++ error stream differs from frozen Python bytes")

    full_frames = frame_bounds(full)
    if len(full_frames) != 3:
        raise SystemExit("C++ Full-Details stream has the wrong frame count")
    batch_start, _, _, batch_end = full_frames[1]
    if full[batch_start:batch_end] != fixture(fixture_root, "batch-full"):
        raise SystemExit("C++ Full-Details batch differs from frozen Python bytes")

    decoded_none = read_stream(
        none,
        expected_skipped_fragment_count=7,
    )
    if (
        decoded_none.trailer.fragment_count != 2
        or decoded_none.trailer.fragment_batch_count != 1
        or tuple(decoded_none.batches[0].template_offsets) != (0, 4, 8)
        or not decoded_none.batches[0].raw_payload.readonly
    ):
        raise SystemExit("Python decoded incorrect C++ no-Details values")
    decoded_full = read_stream(full)
    details = decoded_full.batches[0].details
    if details is None or tuple(details.variant_kinds) != (1, 2, 3):
        raise SystemExit("Python decoded incorrect C++ Full-Details events")
    try:
        read_stream(error)
    except CoreReportedError as reported:
        if reported.error_code != 1204 or reported.message != "batch exceeds limit":
            raise SystemExit("Python decoded incorrect C++ error fields")
    else:
        raise SystemExit("Python accepted the C++ terminal error stream")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
