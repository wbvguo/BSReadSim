"""Mutation gate for the required installed C-extension validator."""

from dataclasses import replace
from pathlib import Path
import struct
import sys

import bsreadsim.htsim.protocol as protocol
import bsreadsim.process.fragment as protocol_adapter


FRAME_ENVELOPE = struct.Struct("<IBBHQ")
PREAMBLE_SIZE = struct.calcsize("<8sHHI")


def fixture(root: Path, name: str) -> bytes:
    path = root / "tests" / "fixtures" / (name + ".hex")
    return bytes.fromhex(path.read_text(encoding="ascii"))


def frame_payload(frame: bytes, *, offset: int = 0) -> tuple[int, bytes]:
    payload_length, _frame_type, flags, reserved, _sequence = (
        FRAME_ENVELOPE.unpack_from(frame, offset)
    )
    if reserved != 0:
        raise AssertionError("fixture frame has a nonzero reserved field")
    begin = offset + FRAME_ENVELOPE.size
    end = begin + payload_length
    return flags, frame[begin:end]


def accepts(payload: bytes, flags: int, header) -> bool:
    try:
        protocol._decode_batch_payload(
            payload,
            flags,
            header,
            expected_first_ordinal=0,
        )
    except (MemoryError, protocol.ProtocolError):
        return False
    return True


def check_fixture(
    root: Path,
    name: str,
    header,
) -> tuple[int, int]:
    flags, payload = frame_payload(fixture(root, name))
    candidates = [payload]
    candidates.extend(payload[:length] for length in range(len(payload)))
    candidates.extend(
        payload + suffix
        for suffix in (b"\x00", b"\x01", b"\x00\x00", b"\x00\x00\x00\x00")
    )
    for offset in range(len(payload)):
        for mask in (0x01, 0x80, 0xFF):
            changed = bytearray(payload)
            changed[offset] ^= mask
            candidates.append(bytes(changed))

    accepted = sum(accepts(candidate, flags, header) for candidate in candidates)
    return len(candidates), accepted


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit(
            "usage: check_installed_cext.py REPOSITORY_ROOT"
        )
    root = Path(argv[1]).resolve()
    header_frame = fixture(root, "header-none")
    _flags, header_payload = frame_payload(header_frame, offset=PREAMBLE_SIZE)
    without_annotations_header = protocol._decode_header_payload(header_payload)
    full_annotations_header = replace(
        without_annotations_header,
        has_details=True,
    )

    totals = []
    totals.append(
        check_fixture(root, "batch-none", without_annotations_header)
    )
    totals.append(
        check_fixture(root, "batch-full", full_annotations_header)
    )
    common_flags, common_payload = frame_payload(fixture(root, "batch-none"))
    common_batch = protocol._decode_batch_payload(
        common_payload,
        common_flags,
        without_annotations_header,
        expected_first_ordinal=0,
    )
    cext_common = protocol_adapter.decode_common_numpy_batch(
        common_batch,
        without_annotations_header,
    )
    if (
        cext_common.fragment_count != common_batch.fragment_count
        or cext_common.mate_count != common_batch.mate_count
    ):
        raise AssertionError(
            "C-extension common-column adapter changed fixture cardinality"
        )
    full_flags, full_payload = frame_payload(fixture(root, "batch-full"))
    full_batch = protocol._decode_batch_payload(
        full_payload,
        full_flags,
        full_annotations_header,
        expected_first_ordinal=0,
    )
    cext_fragments = protocol_adapter.decode_fragments(
        full_batch,
        full_annotations_header,
    )
    if (
        len(cext_fragments) != full_batch.fragment_count
        or sum(len(fragment.mates) for fragment in cext_fragments)
        != full_batch.mate_count
    ):
        raise AssertionError(
            "C-extension Full-Details adapter changed fixture cardinality"
        )
    if tuple(cext_fragments[1].reference_positions) != (20, 21, -1, 22, 23):
        raise AssertionError("C-extension Full-Details adapter changed reference projection")
    print(
        "protocol C-extension differential gate: {} candidates, {} accepted; "
        "common and Full-Details adapters fixture-validated".format(
            sum(total for total, _accepted in totals),
            sum(accepted for _total, accepted in totals),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
