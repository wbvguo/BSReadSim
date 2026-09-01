"""Generate BAM and independently inspect its BAM binary contract."""

from __future__ import annotations

import base64
import gzip
import hashlib
import json
import os
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

from bsreadsim.run.manifest import verify_complete_manifest


BGZF_EOF = bytes.fromhex(
    "1f8b08040000000000ff0600424302001b0003000000000000000000"
)
CORE = struct.Struct("<iiIIiiii")
INT32 = struct.Struct("<i")
UINT32 = struct.Struct("<I")
CIGAR_OPERATIONS = "MIDNSHP=XB"
SEQUENCE_CODES = "=ACMGRSVTWYHKDBN"
INTEGER_AUX_TYPES = {
    "c": ("<b", 1),
    "C": ("<B", 1),
    "s": ("<h", 2),
    "S": ("<H", 2),
    "i": ("<i", 4),
    "I": ("<I", 4),
}


def _environment() -> dict:
    return os.environ.copy()


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
        if tag in values:
            raise SystemExit("BAM contains a duplicate auxiliary tag")
        if kind == "Z":
            end = data.find(b"\x00", cursor)
            if end < 0:
                raise SystemExit("BAM contains an unterminated string tag")
            value = data[cursor:end].decode("ascii")
            cursor = end + 1
        elif kind in INTEGER_AUX_TYPES:
            value_format, size = INTEGER_AUX_TYPES[kind]
            encoded, cursor = _take(data, cursor, size)
            value = struct.unpack(value_format, encoded)[0]
        elif kind == "B":
            subtype_data, cursor = _take(data, cursor, 1)
            subtype = subtype_data.decode("ascii")
            if subtype not in INTEGER_AUX_TYPES:
                raise SystemExit("BAM emitted an unsupported array subtype")
            count_data, cursor = _take(data, cursor, 4)
            count = INT32.unpack(count_data)[0]
            if count < 0:
                raise SystemExit("BAM emitted a negative array count")
            value_format, size = INTEGER_AUX_TYPES[subtype]
            payload, cursor = _take(data, cursor, count * size)
            scalar_format = value_format[-1]
            value = (
                struct.unpack("<{}{}".format(count, scalar_format), payload)
                if count
                else ()
            )
            kind = "B:" + subtype
        else:
            raise SystemExit("BAM emitted an unexpected auxiliary type")
        values[tag] = (kind, value)
    return values


def _decode_zx(value: str) -> tuple[int, int]:
    parts = value.split(".")
    if len(parts) != 4:
        raise SystemExit("BAM zx does not have four compact fields")
    try:
        site_count = int(parts[0], 16)
        convertible_count = int(parts[1], 16)
    except ValueError as error:
        raise SystemExit("BAM zx counts are not lowercase hexadecimal") from error
    if parts[0] != format(site_count, "x") or parts[1] != format(
        convertible_count, "x"
    ):
        raise SystemExit("BAM zx counts are not canonical hexadecimal")

    payloads = []
    for text, count in zip(
        parts[2:], (site_count, convertible_count), strict=True
    ):
        try:
            payload = base64.b64decode(
                text + "=" * ((-len(text)) % 4),
                altchars=b"-_",
                validate=True,
            )
        except (ValueError, base64.binascii.Error) as error:
            raise SystemExit("BAM zx bitset is not base64url") from error
        if len(payload) != (count + 7) // 8:
            raise SystemExit("BAM zx bitset length disagrees with its count")
        if count % 8 and payload and payload[-1] >> (count % 8):
            raise SystemExit("BAM zx has nonzero padding bits")
        payloads.append(payload)
    return site_count, convertible_count


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
        raise SystemExit("BAM unexpectedly omitted base qualities")
    quality = "".join(chr(value + 33) for value in quality_data)
    if query_consumed != sequence_length:
        raise SystemExit("BAM CIGAR does not consume its complete query")
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
            "BAM is missing the canonical BGZF EOF marker: {}".format(
                raw[-28:].hex()
            )
        )
    try:
        data = gzip.decompress(raw)
    except (OSError, EOFError) as error:
        raise SystemExit("BAM BGZF stream cannot be decompressed") from error
    if not data.startswith(b"BAM\x01"):
        raise SystemExit("BAM magic is invalid")
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


def _check_htslib_indel_serialization(root: Path, core: Path) -> None:
    sam = (
        b"@HD\tVN:1.6\tSO:unsorted\n"
        b"@SQ\tSN:chr1\tLN:100\n"
        b"indel\t0\tchr1\t11\t60\t2M1I1M2D1M\t*\t0\t0\tACGTN\tABCDE\n"
        b"pure-insertion\t0\tchr1\t21\t60\t3I\t*\t0\t0\tACG\tABC\n"
    )
    completed = subprocess.run(
        [str(core), "--sam-to-bam", "0", "0"],
        input=sam,
        check=False,
        capture_output=True,
    )
    parallel = subprocess.run(
        [str(core), "--sam-to-bam", "0", "4"],
        input=sam,
        check=False,
        capture_output=True,
    )
    if parallel.returncode != 0 or parallel.stderr:
        raise SystemExit("parallel HTSlib BAM writer rejected valid SAM")
    if parallel.stdout != completed.stdout:
        raise SystemExit("BAM bytes changed with compression thread count")
    if completed.returncode != 0:
        raise SystemExit(
            "HTSlib rejected details indel CIGARs: {!r}".format(completed.stderr)
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
        raise SystemExit("HTSlib changed a details indel CIGAR")

    malformed = subprocess.run(
        [str(core), "--sam-to-bam", "6", "0"],
        input=b"not a SAM stream\n",
        check=False,
        capture_output=True,
    )
    if malformed.returncode == 0:
        raise SystemExit("HTSlib helper accepted a malformed SAM stream")


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_bam.py CORE_EXECUTABLE")
    core = Path(sys.argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-bam-") as temporary:
        root = Path(temporary).resolve()
        _check_htslib_indel_serialization(root, core)
        (root / "tiny.fa").write_bytes(b">chr1\nACGTCGTAA\n")
        completed = subprocess.run(
            [
                sys.executable,
                "-m",
                "bsreadsim",
                "run",
                "wgbs",
                "-r",
                "tiny.fa",
                "-o",
                "output",
                "-n",
                "256",
                "--seed",
                "17",
                "--mutation-rate",
                "0",
                "--read-length",
                "3",
                "--insert-mean",
                "5",
                "--insert-sd",
                "0",
                "--max-ambiguous-fraction",
                "0",
                "--conversion-rate",
                "1",
                "--phred",
                "37",
                "--error-rate",
                "0",
                "--threads",
                "16",
                "--prefix",
                "sample",
                "--format",
                "bam",
                "--fragment-realization",
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
                "BAM CLI failed: status={} stderr={!r}".format(
                    completed.returncode, completed.stderr
                )
            )
        manifest_path = Path(completed.stdout.strip())
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        verify_complete_manifest(manifest)
        roles = {item["role"]: item for item in manifest["outputs"]}
        if set(roles) != {"bam"}:
            raise SystemExit("BAM run emitted the wrong output roles")
        if roles["bam"]["record_count"] != 256:
            raise SystemExit("BAM manifest recorded the wrong alignment count")
        bam_path = root / "output" / "sample.bam"
        if roles["bam"]["sha256"] != hashlib.sha256(
            bam_path.read_bytes()
        ).hexdigest():
            raise SystemExit("BAM manifest digest is wrong")
        details = manifest["details"]
        if "details" in details["configuration"]["output"]:
            raise SystemExit("BAM manifest retained the deleted JSON policy")
        if details.get("alignment", {}).get("sam_version") != "1.6":
            raise SystemExit("BAM manifest omitted its SAM contract")
        if details["contracts"].get("read_name") != "bsreadsim-read-name":
            raise SystemExit("BAM manifest omitted the read-name contract")
        if details["models"].get("library_orientation") != {
            "effective": "directional-ot-ob-equal",
            "rng_stage": "library-orientation",
        }:
            raise SystemExit("BAM manifest omitted its library orientation model")
        if details["alignment"]["tags"]["zf"]["required"] is not True:
            raise SystemExit("BAM manifest omitted the requested zf policy")
        if details["alignment"]["tags"]["zx"]["required"] is not True:
            raise SystemExit("BAM manifest omitted the requested zx policy")
        expected_bisulfite_tags = {
            "XG": {
                "required": True,
                "schema": "bismark-genome-conversion",
                "scope": "bisulfite-genome-conversion",
                "values": ["CT", "GA"],
            },
            "XR": {
                "required": True,
                "schema": "bismark-read-conversion",
                "scope": "bisulfite-read-conversion",
                "values": ["CT", "GA"],
            },
            "YS": {
                "required": True,
                "schema": "bismark-strand-id",
                "scope": "bisulfite-library-strand",
                "values": ["OT", "OB", "CTOT", "CTOB"],
            },
        }
        if {
            name: details["alignment"]["tags"].get(name)
            for name in expected_bisulfite_tags
        } != expected_bisulfite_tags:
            raise SystemExit("BAM manifest omitted the Bismark tag policy")

        header, references, records, _ = _parse_bam(bam_path)
        if not header.startswith("@HD\tVN:1.6\tSO:unsorted\n"):
            raise SystemExit("BAM header has the wrong sort contract")
        if "@RG\t" not in header or "@PG\tID:bsreadsim" not in header:
            raise SystemExit("BAM header omitted provenance records")
        if "BSREADSIM_ZX=packed-b64url;ENABLED=1;BIT_ORDER=LSB0" not in header:
            raise SystemExit("BAM header omitted the zx schema")
        for declaration in (
            "BSREADSIM_XG=bismark-genome-conversion;ENABLED=1;VALUES=CT|GA",
            "BSREADSIM_XR=bismark-read-conversion;ENABLED=1;VALUES=CT|GA",
            "BSREADSIM_YS=bismark-strand-id;ENABLED=1;"
            "VALUES=OT|OB|CTOT|CTOB",
        ):
            if declaration not in header:
                raise SystemExit("BAM header omitted a Bismark tag schema")
        if references != (("chr1", 9),):
            raise SystemExit("BAM reference dictionary is wrong")
        if len(records) != 256:
            raise SystemExit("BAM contains the wrong record count")
        if (root / "output" / "sample.R1.fastq.gz").exists() or (
            root / "output" / "sample.R2.fastq.gz"
        ).exists():
            raise SystemExit("BAM run emitted forbidden FASTQ sidecars")

        observed_informative_strands = set()
        for pair_offset in range(0, len(records), 2):
            first, second = records[pair_offset : pair_offset + 2]
            if first["query_name"] != second["query_name"]:
                raise SystemExit("BAM paired records do not share QNAME")
            try:
                locus, ordinal_text = first["query_name"].rsplit(":", 1)
                contig, interval = locus.rsplit(":", 1)
                left_text, right_text = interval.split("-", 1)
                left = int(left_text)
                right = int(right_text)
                ordinal = int(ordinal_text, 16)
            except (ValueError, TypeError) as error:
                raise SystemExit("BAM QNAME violates the read-name contract") from error
            if (
                contig != "chr1"
                or not 1 <= left <= right <= 9
                or ordinal != pair_offset // 2
                or ordinal_text != format(ordinal, "x")
            ):
                raise SystemExit("BAM QNAME fields violate the read-name contract")
            if first["mapq"] != 60 or second["mapq"] != 60:
                raise SystemExit("BAM emitted the wrong synthetic MAPQ")
            if first["template_length"] != -second["template_length"]:
                raise SystemExit("BAM TLEN signs are inconsistent")
            if first["next_position"] != second["position"] or second[
                "next_position"
            ] != first["position"]:
                raise SystemExit("BAM mate positions are inconsistent")
            if first["aux"].get("MC") != ("Z", second["cigar"]) or second[
                "aux"
            ].get("MC") != ("Z", first["cigar"]):
                raise SystemExit("BAM MC tags are inconsistent")
            if first["aux"].get("zx") != second["aux"].get("zx"):
                raise SystemExit("BAM paired records do not share fragment zx")
            conversion_tags = (
                tuple(first["aux"].get(name) for name in ("XG", "XR", "YS")),
                tuple(second["aux"].get(name) for name in ("XG", "XR", "YS")),
            )
            if conversion_tags == (
                (("Z", "CT"), ("Z", "CT"), ("Z", "OT")),
                (("Z", "CT"), ("Z", "GA"), ("Z", "OT")),
            ):
                observed_informative_strands.add("W")
                expected_flags = (99, 147)
            elif conversion_tags == (
                (("Z", "GA"), ("Z", "CT"), ("Z", "OB")),
                (("Z", "GA"), ("Z", "GA"), ("Z", "OB")),
            ):
                observed_informative_strands.add("C")
                expected_flags = (83, 163)
            else:
                raise SystemExit(
                    "BAM paired records have incorrect XG/XR/YS tags: {!r}".format(
                        conversion_tags
                    )
                )
            if (first["flag"], second["flag"]) != expected_flags:
                raise SystemExit("BAM flags disagree with its library strand")

        if observed_informative_strands != {"W", "C"}:
            raise SystemExit(
                "directional BAM omitted an informative strand: {!r}".format(
                    observed_informative_strands
                )
            )

        for record in records:
            if record["reference_id"] != 0 or record["next_reference_id"] != 0:
                raise SystemExit("BAM record refers to the wrong reference")
            if record["position"] < 0 or (
                record["position"] + record["reference_consumed"] > 9
            ):
                raise SystemExit("BAM alignment lies outside the reference")
            if len(record["sequence"]) != len(record["quality"]):
                raise SystemExit("BAM cannot recover complete FASTQ fields")
            aux = record["aux"]
            if aux.get("RG") != ("Z", manifest["run_id"]):
                raise SystemExit("BAM record omitted RG provenance")
            if aux.get("AS")[1] != len(record["sequence"]):
                raise SystemExit("BAM AS is not the maximum origin score")
            if aux.get("MQ")[1] != 60:
                raise SystemExit("BAM MQ tag disagrees with MAPQ")
            xg = aux.get("XG")
            xr = aux.get("XR")
            ys = aux.get("YS")
            if xg is None or xg[0] != "Z" or xg[1] not in {"CT", "GA"}:
                raise SystemExit("BAM XG is not a recognized genome conversion")
            if xr is None or xr[0] != "Z" or xr[1] not in {"CT", "GA"}:
                raise SystemExit("BAM XR is not a recognized read conversion")
            if ys is None or ys[0] != "Z" or ys[1] not in {
                "OT", "OB", "CTOT", "CTOB",
            }:
                raise SystemExit("BAM YS is not a recognized library strand")
            zt = aux.get("zt")
            if zt is None or zt[0] != "Z" or len(zt[1]) != len(record["sequence"]):
                raise SystemExit("BAM zt does not cover BAM SEQ")
            for name in ("zr", "zf"):
                value = aux.get(name)
                if value is None or value[0] != "B:S" or len(value[1]) != 12:
                    raise SystemExit("BAM {} violates u16x12".format(name))
            flags = aux["zr"][1][0]
            expected_xg = {1: "CT", 2: "GA"}.get((flags >> 2) & 0x3)
            expected_xr = {0: "CT", 1: "GA"}.get((flags >> 4) & 0x7)
            if xg[1] != expected_xg or xr[1] != expected_xr:
                raise SystemExit("BAM XG/XR disagree with packed zr flags")
            zx = aux.get("zx")
            if zx is None or zx[0] != "Z":
                raise SystemExit("BAM record omitted fragment realization zx")
            site_count, convertible_count = _decode_zx(zx[1])
            if site_count == 0 or convertible_count == 0:
                raise SystemExit("BAM zx unexpectedly contains empty details domains")

        undirectional = subprocess.run(
            [
                sys.executable,
                "-m",
                "bsreadsim",
                "run",
                "wgbs",
                "--reference",
                "tiny.fa",
                "--output",
                "undirectional",
                "--reads",
                "128",
                "--seed",
                "29",
                "--mutation-rate",
                "0",
                "--read-length",
                "3",
                "--insert-mean",
                "5",
                "--insert-sd",
                "0",
                "--max-ambiguous-fraction",
                "0",
                "--conversion-rate",
                "1",
                "--phred",
                "37",
                "--error-rate",
                "0",
                "--undirectional",
                "--format",
                "bam",
                "--prefix",
                "sample",
                "--core",
                str(core),
            ],
            cwd=str(root),
            env=_environment(),
            check=False,
            capture_output=True,
            text=True,
        )
        if undirectional.returncode != 0 or undirectional.stderr:
            raise SystemExit(
                "undirectional BAM CLI failed: status={} stderr={!r}".format(
                    undirectional.returncode,
                    undirectional.stderr,
                )
            )
        _, _, undirectional_records, _ = _parse_bam(
            root / "undirectional" / "sample.bam"
        )
        observed_orientations = set()
        for offset in range(0, len(undirectional_records), 2):
            first, second = undirectional_records[offset : offset + 2]
            pair = tuple(
                tuple(record["aux"][name][1] for name in ("XG", "XR", "YS"))
                for record in (first, second)
            )
            observed_orientations.add(pair)
        expected_orientations = {
            (("CT", "CT", "OT"), ("CT", "GA", "OT")),
            (("GA", "CT", "OB"), ("GA", "GA", "OB")),
            (("CT", "GA", "CTOT"), ("CT", "CT", "CTOT")),
            (("GA", "GA", "CTOB"), ("GA", "CT", "CTOB")),
        }
        if observed_orientations != expected_orientations:
            raise SystemExit(
                "undirectional BAM has incorrect XG/XR/YS orientations: {!r}".format(
                    observed_orientations
                )
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
