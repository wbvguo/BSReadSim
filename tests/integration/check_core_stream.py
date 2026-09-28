"""Check invariants of the fragment stream written by the real htsim."""

from __future__ import annotations

from collections import Counter
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile
from collections.abc import Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from bsreadsim.batch import DELETION, INSERTION, SNV, STRAND_NAMES, FragmentBatch
from bsreadsim.htsim import HtsimCore, _parse_header, _parse_summary
from bsreadsim.settings import Settings

COMMON_COLUMNS = (
    "contig", "start", "end", "haplotype", "strand", "template_offsets", "bases",
    "site_offsets", "site_positions", "site_probabilities", "site_contexts", "site_asm",
)


def _arguments(core: Path, reference: Path, vcf: Path) -> list[str]:
    settings = Settings(
        reference=reference, vcf=vcf, seed=81985529216486895, seed_mut=17, seed_phase=19,
        seed_meth=23, read_length=(12,), insert_min=24, insert_mean=24, insert_max=24,
        insert_sd=0, reads=514, max_ambiguous_fraction=1, indel_extension_probability=0.3,
        beta_cg=(2, 5), beta_chg=(3, 4), beta_chh=(5, 2))
    return HtsimCore(core).argv(settings, batch_size=64)


def _run(arguments: Sequence[str]) -> bytes:
    result = subprocess.run(
        arguments, stdin=subprocess.DEVNULL, capture_output=True, check=False)
    if result.returncode != 0 or result.stderr:
        raise SystemExit("htsim failed: status={} stderr={!r}".format(
            result.returncode, result.stderr))
    return result.stdout


def _parse(data: bytes):
    """(header, batches, summary) of a complete stream."""
    newline = data.index(b"\n")
    header = _parse_header(data[:newline + 1])
    cursor = newline + 1
    batches = []
    while True:
        size = struct.unpack_from("<Q", data, cursor)[0]
        cursor += 8
        if size == 0:
            break
        batches.append(FragmentBatch.decode(data[cursor:cursor + size]))
        cursor += size
    summary = _parse_summary(data[cursor:])
    return header, batches, summary


def _replace(arguments: Sequence[str], option: str, value: str) -> list[str]:
    changed = list(arguments)
    changed[changed.index(option) + 1] = value
    return changed


def _with_details(arguments: Sequence[str], enabled: bool) -> list[str]:
    return _replace(arguments, "--details", str(enabled).lower())


def _common(batches) -> list:
    return [
        tuple(getattr(batch, name).tolist() for name in COMMON_COLUMNS)
        for batch in batches
    ]


def _strands(batches) -> Counter:
    return Counter(STRAND_NAMES[s] for batch in batches for s in batch.strand.tolist())


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit("usage: check_core_stream.py HTSIM_EXECUTABLE")
    core = Path(argv[1]).resolve(strict=True)

    with tempfile.TemporaryDirectory(prefix="htsim-stream-") as directory:
        root = Path(directory)
        sequence = bytearray(b"ACGT" * 40)
        sequence[90] = ord("N")
        reference = root / "reference.fa"
        reference.write_bytes(b">chrParity\n" + bytes(sequence) + b"\n")
        vcf = root / "variants.vcf"
        vcf.write_bytes(
            b"##fileformat=VCFv4.3\n"
            b"#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            b"chrParity\t25\t.\tA\tT\t.\tPASS\t.\tGT\t1|1\n"
            b"chrParity\t41\t.\tA\tAGG\t.\tPASS\t.\tGT\t1|1\n"
            b"chrParity\t61\t.\tAC\tA\t.\tPASS\t.\tGT\t1|1\n"
        )
        arguments = _arguments(core, reference, vcf)

        compact_bytes = _run(_with_details(arguments, False))
        header, compact, summary = _parse(compact_bytes)
        details = json.loads(compact_bytes.split(b"\n", 1)[0])["details"]
        if details or any(len(b.event_kinds) for b in compact):
            raise SystemExit("a stream without details carried variant events")

        full_arguments = _with_details(arguments, True)
        full_bytes = _run(full_arguments)
        full_header, full, full_summary = _parse(full_bytes)
        if [batch.size for batch in full] != [64, 64, 64, 64, 1]:
            raise SystemExit("batch boundaries changed")
        if [b.first_ordinal for b in full] != [0, 64, 128, 192, 256]:
            raise SystemExit("batch ordinals changed")
        if _common(compact) != _common(full) or summary != full_summary:
            raise SystemExit("requesting details changed the simulation")
        if summary.fragment_count != 257 or summary.mate_count != 514:
            raise SystemExit("summary counts are wrong")

        kinds = {int(k) for batch in full for k in batch.event_kinds}
        if kinds != {SNV, INSERTION, DELETION}:
            raise SystemExit("not every variant kind reached the stream")

        directional = _strands(full)
        if set(directional) != {"OT", "OB"}:
            raise SystemExit("directional WGBS strands: {!r}".format(directional))
        nondirectional_arguments = _replace(full_arguments, "--directional", "false")
        nondirectional_bytes = _run(nondirectional_arguments)
        if set(_strands(_parse(nondirectional_bytes)[1])) != set(STRAND_NAMES):
            raise SystemExit("non-directional WGBS omitted a library strand")

        # The thread count must not change the stream.
        for threads in ("2", "4", "8"):
            if _run(_replace(full_arguments, "--threads", threads)) != full_bytes:
                raise SystemExit("{} threads changed the stream".format(threads))
            changed = _run(_replace(nondirectional_arguments, "--threads", threads))
            if changed != nondirectional_bytes:
                raise SystemExit("{} threads changed the non-directional stream".format(threads))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
