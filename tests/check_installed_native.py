"""Differential gate for the installed native protocol validator.

Run this script outside the repository import path with an installed wheel and
pass the repository root as its only argument.  The frozen wire fixtures are
decoded once through the native validator and once through the retained Python
reference validator after deterministic corruption of every payload byte.
"""

from dataclasses import replace
from pathlib import Path
import struct
import sys

import bsreadsim.protocol as protocol
import bsreadsim.protocol_adapter as protocol_adapter


FRAME_ENVELOPE = struct.Struct("<IBBHQ")
PREAMBLE_SIZE = struct.calcsize("<8sHHI")


def fixture(root: Path, name: str) -> bytes:
    path = root / "tests" / "fixtures" / "protocol-v2" / (name + ".hex")
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


def accepts(
    payload: bytes,
    flags: int,
    header,
    validator,
) -> bool:
    previous = protocol._native_validate_batch_columns
    protocol._native_validate_batch_columns = validator
    try:
        protocol._decode_batch_payload(
            payload,
            flags,
            header,
            expected_first_ordinal=0,
        )
    except (MemoryError, protocol.ProtocolError):
        return False
    finally:
        protocol._native_validate_batch_columns = previous
    return True


def check_fixture(
    root: Path,
    name: str,
    header,
    native_validator,
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

    accepted = 0
    for ordinal, candidate in enumerate(candidates):
        reference_result = accepts(candidate, flags, header, None)
        native_result = accepts(candidate, flags, header, native_validator)
        if native_result != reference_result:
            raise AssertionError(
                "{} candidate {} disagrees: native={} reference={}".format(
                    name,
                    ordinal,
                    native_result,
                    reference_result,
                )
            )
        accepted += int(native_result)
    return len(candidates), accepted


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit(
            "usage: check_installed_native.py REPOSITORY_ROOT"
        )
    root = Path(argv[1]).resolve()
    native_validator = protocol._native_validate_batch_columns
    if native_validator is None:
        raise AssertionError("installed wheel has no protocol native validator")

    header_frame = fixture(root, "header-none")
    _flags, header_payload = frame_payload(header_frame, offset=PREAMBLE_SIZE)
    no_truth_header = protocol._decode_header_payload(header_payload)
    full_truth_header = replace(
        no_truth_header,
        truth_columns=protocol.TruthMode.FULL,
    )

    totals = []
    totals.append(
        check_fixture(root, "batch-none", no_truth_header, native_validator)
    )
    totals.append(
        check_fixture(root, "batch-full", full_truth_header, native_validator)
    )
    common_flags, common_payload = frame_payload(fixture(root, "batch-none"))
    common_batch = protocol._decode_batch_payload(
        common_payload,
        common_flags,
        no_truth_header,
        expected_first_ordinal=0,
    )
    reference_common = (
        protocol_adapter._decode_common_numpy_batch_python(
            common_batch,
            no_truth_header,
        )
    )
    native_common = protocol_adapter._decode_common_numpy_batch(
        common_batch,
        no_truth_header,
    )
    if native_common != reference_common:
        raise AssertionError(
            "native and reference common-column adapters disagree"
        )
    full_flags, full_payload = frame_payload(fixture(root, "batch-full"))
    full_batch = protocol._decode_batch_payload(
        full_payload,
        full_flags,
        full_truth_header,
        expected_first_ordinal=0,
    )
    reference_fragments = protocol_adapter._decode_fragments_python(
        full_batch,
        full_truth_header,
    )
    native_fragments = protocol_adapter.decode_fragments(
        full_batch,
        full_truth_header,
    )
    if native_fragments != reference_fragments:
        raise AssertionError(
            "native and reference Full-Truth adapters disagree"
        )
    print(
        "protocol native differential gate: {} candidates, {} accepted; "
        "common and Full-Truth adapters exact".format(
            sum(total for total, _accepted in totals),
            sum(accepted for _total, accepted in totals),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
