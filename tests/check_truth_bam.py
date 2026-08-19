"""Generate truth BAM and independently inspect its BAM binary contract."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

from bsreadsim.manifest import verify_complete_manifest


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)
CORE = struct.Struct("<iiIIiiii")
INT32 = struct.Struct("<i")
UINT32 = struct.Struct("<I")
CIGAR_OPERATIONS = "MIDNSHP=XB"
SEQUENCE_CODES = "=ACMGRSVTWYHKDBN"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def _environment() -> dict:
    environment = os.environ.copy()
    python_paths = []
    if environment.get("BSREADSIM_TEST_INSTALLED") != "1":
        python_paths.append(str(REPOSITORY_ROOT / "src"))
    if environment.get("PYTHONPATH"):
        python_paths.append(environment["PYTHONPATH"])
    if python_paths:
        environment["PYTHONPATH"] = os.pathsep.join(python_paths)
    else:
        environment.pop("PYTHONPATH", None)
    return environment


def _take(data: bytes, cursor: int, length: int):
    end = cursor + length
    if length < 0 or end > len(data):
        raise SystemExit("BAM structure is truncated")
    return data[cursor:end], end


def _i32(data: bytes, cursor: int):
    value, end = _take(data, cursor, 4)
    return INT32.unpack(value)[0], end


def _parse_aux(data: bytes) -> dict:
    values = {}
    cursor = 0
    while cursor < len(data):
        prefix, cursor = _take(data, cursor, 3)
        tag = prefix[:2].decode("ascii")
        kind = chr(prefix[2])
        if kind != "Z":
            raise SystemExit("truth BAM emitted an unexpected auxiliary type")
        end = data.find(b"\x00", cursor)
        if end < 0:
            raise SystemExit("truth BAM contains an unterminated auxiliary value")
        if tag in values:
            raise SystemExit("truth BAM contains a duplicate auxiliary tag")
        values[tag] = data[cursor:end].decode("ascii")
        cursor = end + 1
    return values


def _parse_record(block: bytes) -> dict:
    if len(block) < CORE.size:
        raise SystemExit("BAM alignment core is truncated")
    (
        reference_id,
        position,
        bin_mq_nl,
        flag_nc,
        sequence_length,
        next_reference_id,
        next_position,
        template_length,
    ) = CORE.unpack_from(block)
    query_name_length = bin_mq_nl & 0xFF
    mapq = (bin_mq_nl >> 8) & 0xFF
    cigar_count = flag_nc & 0xFFFF
    flag = flag_nc >> 16
    cursor = CORE.size
    query_name_data, cursor = _take(block, cursor, query_name_length)
    if not query_name_data.endswith(b"\x00"):
        raise SystemExit("BAM query name is not NUL terminated")
    query_name = query_name_data[:-1].decode("ascii")
    cigar_data, cursor = _take(block, cursor, cigar_count * 4)
    cigar_parts = []
    query_consumed = 0
    reference_consumed = 0
    for offset in range(0, len(cigar_data), 4):
        encoded = UINT32.unpack_from(cigar_data, offset)[0]
        length = encoded >> 4
        operation_index = encoded & 0xF
        if not length or operation_index >= len(CIGAR_OPERATIONS):
            raise SystemExit("BAM CIGAR operation is invalid")
        operation = CIGAR_OPERATIONS[operation_index]
        cigar_parts.append("{}{}".format(length, operation))
        if operation in "MIS=X":
            query_consumed += length
        if operation in "MDN=X":
            reference_consumed += length
    packed_sequence, cursor = _take(block, cursor, (sequence_length + 1) // 2)
    sequence = "".join(
        SEQUENCE_CODES[
            (packed_sequence[index // 2] >> (4 if index % 2 == 0 else 0)) & 0xF
        ]
        for index in range(sequence_length)
    )
    quality_data, cursor = _take(block, cursor, sequence_length)
    if any(value == 0xFF for value in quality_data):
        raise SystemExit("truth BAM unexpectedly omitted base qualities")
    quality = "".join(chr(value + 33) for value in quality_data)
    if query_consumed != sequence_length:
        raise SystemExit("truth BAM CIGAR does not consume its complete query")
    return {
        "aux": _parse_aux(block[cursor:]),
        "cigar": "".join(cigar_parts),
        "flag": flag,
        "mapq": mapq,
        "next_position": next_position,
        "next_reference_id": next_reference_id,
        "position": position,
        "query_name": query_name,
        "quality": quality,
        "reference_consumed": max(reference_consumed, 1),
        "reference_id": reference_id,
        "sequence": sequence,
        "template_length": template_length,
    }


def _parse_bam(path: Path):
    raw = path.read_bytes()
    if not raw.endswith(BGZF_EOF):
        raise SystemExit(
            "truth BAM is missing the canonical BGZF EOF marker: {}".format(
                raw[-28:].hex()
            )
        )
    try:
        data = gzip.decompress(raw)
    except (OSError, EOFError) as error:
        raise SystemExit("truth BAM BGZF stream cannot be decompressed") from error
    if not data.startswith(b"BAM\x01"):
        raise SystemExit("truth BAM magic is invalid")
    cursor = 4
    header_length, cursor = _i32(data, cursor)
    header_data, cursor = _take(data, cursor, header_length)
    header = header_data.decode("ascii")
    reference_count, cursor = _i32(data, cursor)
    references = []
    for _ in range(reference_count):
        name_length, cursor = _i32(data, cursor)
        name_data, cursor = _take(data, cursor, name_length)
        length, cursor = _i32(data, cursor)
        if not name_data.endswith(b"\x00"):
            raise SystemExit("BAM reference name is not NUL terminated")
        references.append((name_data[:-1].decode("ascii"), length))
    records = []
    while cursor < len(data):
        block_size, cursor = _i32(data, cursor)
        block, cursor = _take(data, cursor, block_size)
        records.append(_parse_record(block))
    return header, tuple(references), tuple(records), raw


def _fastq(path: Path) -> dict:
    with gzip.open(path, "rt", encoding="ascii", newline="") as stream:
        lines = stream.read().splitlines()
    if len(lines) % 4:
        raise SystemExit("FASTQ output is truncated")
    values = {}
    for offset in range(0, len(lines), 4):
        if not lines[offset].startswith("@") or lines[offset + 2] != "+":
            raise SystemExit("FASTQ record structure is invalid")
        values[lines[offset][1:]] = (lines[offset + 1], lines[offset + 3])
    return values


def _check_htslib_indel_serialization(root: Path, core: Path) -> None:
    sam = (
        b"@HD\tVN:1.6\tSO:unsorted\n"
        b"@SQ\tSN:chr1\tLN:100\n"
        b"indel\t0\tchr1\t11\t60\t2M1I1M2D1M\t*\t0\t0\tACGTN\tABCDE\n"
        b"pure-insertion\t0\tchr1\t21\t60\t3I\t*\t0\t0\tACG\tABC\n"
    )
    completed = subprocess.run(
        [str(core), "--sam-to-bam", "0"],
        input=sam,
        check=False,
        capture_output=True,
    )
    if completed.returncode != 0:
        raise SystemExit(
            "HTSlib rejected truth indel CIGARs: {!r}".format(completed.stderr)
        )
    path = root / "indel-contract.bam"
    path.write_bytes(completed.stdout)
    _, references, records, _ = _parse_bam(path)
    if references != (("chr1", 100),):
        raise SystemExit("HTSlib indel BAM reference dictionary changed")
    if tuple(record["cigar"] for record in records) != (
        "2M1I1M2D1M",
        "3I",
    ):
        raise SystemExit("HTSlib changed a truth indel CIGAR")

    malformed = subprocess.run(
        [str(core), "--sam-to-bam", "6"],
        input=b"not a SAM stream\n",
        check=False,
        capture_output=True,
    )
    if malformed.returncode == 0:
        raise SystemExit("HTSlib helper accepted a malformed SAM stream")


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_truth_bam.py CORE_EXECUTABLE")
    core = Path(sys.argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-truth-bam-") as temporary:
        root = Path(temporary).resolve()
        _check_htslib_indel_serialization(root, core)
        (root / "tiny.fa").write_bytes(b">chr1\nACGTCGTAA\n")
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "bsreadsim",
                "run",
                "-r",
                "tiny.fa",
                "-o",
                "output",
                "-n",
                "4",
                "--seed",
                "17",
                "--mutation-rate",
                "0",
                "--read-length",
                "3",
                "--insert-size",
                "5",
                "--max-ambiguous-fraction",
                "0",
                "--conversion-rate",
                "1",
                "--phred",
                "37",
                "--error-rate",
                "0",
                "--workers",
                "2",
                "--core-workers",
                "1",
                "--chunk-size",
                "4",
                "--max-in-flight-fragments",
                "4",
                "--prefix",
                "sample",
                "--truth-bam",
                "--core",
                str(core),
            ],
            cwd=str(root),
            env=_environment(),
            check=False,
            capture_output=True,
            text=True,
        )
        if completed.returncode != 0 or completed.stderr:
            raise SystemExit(
                "truth BAM CLI failed: status={} stderr={!r}".format(
                    completed.returncode, completed.stderr
                )
            )
        manifest_path = Path(completed.stdout.strip())
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        verify_complete_manifest(manifest)
        roles = {item["role"]: item for item in manifest["outputs"]}
        if set(roles) != {"read1", "read2", "truth_bam"}:
            raise SystemExit("truth BAM run emitted the wrong output roles")
        if roles["truth_bam"]["record_count"] != 8:
            raise SystemExit("truth BAM manifest recorded the wrong alignment count")
        bam_path = root / "output" / "sample.truth.bam"
        if roles["truth_bam"]["sha256"] != hashlib.sha256(
            bam_path.read_bytes()
        ).hexdigest():
            raise SystemExit("truth BAM manifest digest is wrong")
        if manifest["config"]["normalized"]["output"]["truth"] != "none":
            raise SystemExit("truth BAM unexpectedly enabled JSON Full Truth")
        if manifest.get("truth_alignment", {}).get("sam_version") != "1.6":
            raise SystemExit("truth BAM manifest omitted its SAM contract")
        if manifest["versions"].get("read_name") != "bsreadsim-read-name-v1":
            raise SystemExit("truth BAM manifest omitted read-name v1")

        header, references, records, _ = _parse_bam(bam_path)
        if not header.startswith("@HD\tVN:1.6\tSO:unsorted\n"):
            raise SystemExit("truth BAM header has the wrong sort contract")
        if "@RG\t" not in header or "@PG\tID:bsreadsim" not in header:
            raise SystemExit("truth BAM header omitted provenance records")
        if references != (("chr1", 9),):
            raise SystemExit("truth BAM reference dictionary is wrong")
        if len(records) != 8:
            raise SystemExit("truth BAM contains the wrong record count")

        fastq = {}
        fastq.update(_fastq(root / "output" / "sample.R1.fastq.gz"))
        fastq.update(_fastq(root / "output" / "sample.R2.fastq.gz"))
        for pair_offset in range(0, len(records), 2):
            first, second = records[pair_offset : pair_offset + 2]
            if first["query_name"] != second["query_name"]:
                raise SystemExit("truth BAM paired records do not share QNAME")
            try:
                locus, ordinal_text = first["query_name"].rsplit(":", 1)
                contig, interval = locus.rsplit(":", 1)
                left_text, right_text = interval.split("-", 1)
                left = int(left_text)
                right = int(right_text)
                ordinal = int(ordinal_text)
            except (ValueError, TypeError) as error:
                raise SystemExit("truth BAM QNAME violates read-name v1") from error
            if (
                contig != "chr1"
                or not 1 <= left <= right <= 9
                or ordinal != pair_offset // 2
            ):
                raise SystemExit("truth BAM QNAME fields violate read-name v1")
            if first["flag"] != 99 or second["flag"] != 147:
                raise SystemExit("truth BAM emitted incorrect paired-end flags")
            if first["mapq"] != 60 or second["mapq"] != 60:
                raise SystemExit("truth BAM emitted the wrong synthetic MAPQ")
            if first["template_length"] != -second["template_length"]:
                raise SystemExit("truth BAM TLEN signs are inconsistent")
            if first["next_position"] != second["position"] or second[
                "next_position"
            ] != first["position"]:
                raise SystemExit("truth BAM mate positions are inconsistent")
            if first["aux"].get("MC") != second["cigar"] or second[
                "aux"
            ].get("MC") != first["cigar"]:
                raise SystemExit("truth BAM MC tags are inconsistent")

        for record in records:
            if record["reference_id"] != 0 or record["next_reference_id"] != 0:
                raise SystemExit("truth BAM record refers to the wrong reference")
            if record["position"] < 0 or (
                record["position"] + record["reference_consumed"] > 9
            ):
                raise SystemExit("truth BAM alignment lies outside the reference")
            mate = 1 if record["flag"] & 0x40 else 2
            identifier = "{}/{}".format(record["query_name"], mate)
            sequence = record["sequence"]
            quality = record["quality"]
            if record["flag"] & 0x10:
                sequence = sequence.translate(COMPLEMENT)[::-1]
                quality = quality[::-1]
            if fastq.get(identifier) != (sequence, quality):
                raise SystemExit(
                    "truth BAM does not round-trip to the emitted FASTQ record"
                )
            if (
                record["aux"].get("RG") is None
                or record["aux"].get("PG") != "bsreadsim"
            ):
                raise SystemExit("truth BAM record omitted RG/PG provenance tags")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
