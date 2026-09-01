"""Validate the sole columnar protocol on the real C++ core."""

from __future__ import annotations

import hashlib
from pathlib import Path
import subprocess
import sys
import tempfile
from collections.abc import Sequence

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from bsreadsim.htsim.subprocess import CoreProcess
from bsreadsim.process.batch import VariantKind
from bsreadsim.htsim.protocol import ProtocolStream
from tests.helpers.stream_support import read_stream
from bsreadsim.process.fragment import decode_fragments


def _arguments(core: Path, reference: Path, vcf: Path) -> list[str]:
    return [
        str(core),
        "--run-id", "00000000-0000-4000-8000-000000000002",
        "--config-sha256", hashlib.sha256(b"core-protocol-contract").hexdigest(),
        "--seed", "81985529216486895",
        "--seed-mut", "17",
        "--seed-phase", "19",
        "--seed-meth", "23",
        "--reference", str(reference),
        "--vcf", str(vcf),
        "--technology", "WGBS",
        "--directional", "true",
        "--paired-end", "true",
        "--read-length-1", "12",
        "--read-length-2", "12",
        "--insert-min", "24",
        "--insert-mean", "24",
        "--insert-max", "24",
        "--insert-sd", "0",
        "--fragments", "257",
        "--max-ambiguous-fraction", "1",
        "--chunk-size", "7",
        "--core-workers", "1",
        "--protocol-batch-fragments", "64",
        "--mutation-rate", "0",
        "--indel-fraction", "0.15",
        "--indel-extension-probability", "0.3",
        "--homozygous-only", "false",
        "--collect-non-cpg", "true",
        "--pool-meth", "false",
        "--update-variant-boundaries", "true",
        "--beta-cg", "2,5",
        "--beta-chg", "3,4",
        "--beta-chh", "5,2",
        "--coverage", "uniform",
    ]


def _run(arguments: Sequence[str]) -> subprocess.CompletedProcess[bytes]:
    return subprocess.run(
        arguments,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        check=False,
    )


def _require_success(result: subprocess.CompletedProcess[bytes], name: str) -> bytes:
    if result.returncode != 0 or result.stderr:
        raise SystemExit(
            "{} failed: status={} stderr={!r}".format(
                name, result.returncode, result.stderr
            )
        )
    return result.stdout


def _replace(arguments: Sequence[str], option: str, value: str) -> list[str]:
    changed = list(arguments)
    index = changed.index(option)
    changed[index + 1] = value
    return changed


def _with_details(arguments: Sequence[str], enabled: bool) -> list[str]:
    return list(arguments) + ["--emit-details", str(enabled).lower()]


def _common_signature(stream: ProtocolStream) -> tuple:
    return tuple(
        (
            batch.first_fragment_ordinal,
            tuple(batch.contig_indices),
            tuple(batch.reference_starts),
            tuple(batch.reference_ends),
            tuple(batch.template_offsets),
            tuple(batch.mate_offsets),
            tuple(batch.site_offsets),
            tuple(batch.mate_template_starts),
            tuple(batch.mate_template_ends),
            tuple(batch.site_template_offsets),
            tuple(batch.site_probabilities),
            tuple(batch.haplotypes),
            tuple(batch.capture_strands),
            tuple(batch.mate_indices),
            tuple(batch.mate_reverse_complements),
            tuple(batch.site_contexts),
            tuple(batch.methylation_sources),
            tuple(batch.site_alleles),
            bytes(batch.template_bases),
        )
        for batch in stream.batches
    )


def _counts(stream: ProtocolStream) -> tuple:
    trailer = stream.trailer
    return (
        trailer.fragment_count,
        trailer.mate_count,
        trailer.template_base_count,
        trailer.methylation_site_count,
        trailer.skipped_fragment_count,
        tuple(trailer.per_contig_fragment_counts),
    )


def _library_orientation_counts(stream: ProtocolStream) -> dict[str, int]:
    counts = {"OT": 0, "OB": 0, "CTOT": 0, "CTOB": 0}
    signatures = {
        (1, (0, 1)): "OT",
        (2, (1, 0)): "OB",
        (1, (1, 0)): "CTOT",
        (2, (0, 1)): "CTOB",
    }
    for batch in stream.batches:
        for fragment_index, capture_strand in enumerate(batch.capture_strands):
            mate_begin = batch.mate_offsets[fragment_index]
            mate_end = batch.mate_offsets[fragment_index + 1]
            reverse = tuple(batch.mate_reverse_complements[mate_begin:mate_end])
            orientation = signatures.get((capture_strand, reverse))
            if orientation is None:
                raise SystemExit(
                    "core emitted an invalid library orientation: {!r}".format(
                        (capture_strand, reverse)
                    )
                )
            counts[orientation] += 1
    return counts


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        raise SystemExit("usage: check_core_stream.py CORE_EXECUTABLE")
    core = Path(argv[1]).resolve(strict=True)

    with tempfile.TemporaryDirectory(prefix="htsim-core-protocol-") as directory:
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

        default_bytes = _require_success(_run(arguments), "default invocation")
        without_annotations_bytes = _require_success(
            _run(_with_details(arguments, False)), "explicit no-Details invocation"
        )
        if default_bytes != without_annotations_bytes:
            raise SystemExit("explicit no-Details changed default protocol bytes")
        without_details = read_stream(without_annotations_bytes, core_exit_status=0)
        if without_details.header.has_details:
            raise SystemExit("default core invocation did not select no-Details")
        if any(batch.details is not None for batch in without_details.batches):
            raise SystemExit("no-Details stream emitted provenance columns")

        full_arguments = _with_details(arguments, True)
        full_bytes = _require_success(_run(full_arguments), "Full-Details invocation")
        full = read_stream(full_bytes, core_exit_status=0)
        if [batch.fragment_count for batch in full.batches] != [64, 64, 64, 64, 1]:
            raise SystemExit("canonical protocol batch boundaries changed")
        if _common_signature(without_details) != _common_signature(full):
            raise SystemExit("Details selection changed common simulation columns")
        if _counts(without_details) != _counts(full):
            raise SystemExit("Details selection changed simulation counts")
        if len(without_annotations_bytes) >= len(full_bytes):
            raise SystemExit("no-Details stream did not remove provenance bytes")

        directional_counts = _library_orientation_counts(full)
        if (
            directional_counts["OT"] == 0
            or directional_counts["OB"] == 0
            or directional_counts["CTOT"] != 0
            or directional_counts["CTOB"] != 0
        ):
            raise SystemExit(
                "directional WGBS did not contain only independent OT/OB fragments: "
                "{!r}".format(directional_counts)
            )

        nondirectional_arguments = _replace(
            full_arguments, "--directional", "false"
        )
        nondirectional_bytes = _require_success(
            _run(nondirectional_arguments), "non-directional invocation"
        )
        nondirectional = read_stream(
            nondirectional_bytes, core_exit_status=0
        )
        nondirectional_counts = _library_orientation_counts(nondirectional)
        if any(count == 0 for count in nondirectional_counts.values()):
            raise SystemExit(
                "non-directional WGBS omitted a library orientation: {!r}".format(
                    nondirectional_counts
                )
            )
        nondirectional_rechunked = _require_success(
            _run(_replace(nondirectional_arguments, "--chunk-size", "31")),
            "non-directional rechunked invocation",
        )
        if nondirectional_rechunked != nondirectional_bytes:
            raise SystemExit(
                "non-directional library orientation changed with chunk size"
            )

        variant_kinds = {
            int(event.kind)
            for batch in full.batches
            for fragment in decode_fragments(batch, full.header)
            for event in fragment.variants
        }
        if variant_kinds != {
            int(VariantKind.SNV),
            int(VariantKind.INSERTION),
            int(VariantKind.DELETION),
        }:
            raise SystemExit("fixture did not preserve all variant event kinds")

        for worker_count in (2, 4):
            parallel = _require_success(
                _run(_replace(full_arguments, "--core-workers", str(worker_count))),
                "core worker {}".format(worker_count),
            )
            if parallel != full_bytes:
                raise SystemExit("core worker count changed protocol bytes")
        large_chunk_arguments = _replace(
            full_arguments, "--chunk-size", "257"
        )
        large_chunk_serial = _require_success(
            _run(large_chunk_arguments), "large-chunk serial invocation"
        )
        large_chunk_parallel = _require_success(
            _run(_replace(large_chunk_arguments, "--core-workers", "4")),
            "large-chunk parallel invocation",
        )
        if large_chunk_parallel != large_chunk_serial:
            raise SystemExit("parallel fragment construction changed protocol bytes")
        rechunked = _require_success(
            _run(_replace(full_arguments, "--chunk-size", "31")),
            "rechunked invocation",
        )
        if rechunked != full_bytes:
            raise SystemExit("generation chunk size changed protocol bytes")

        retired = _run(list(arguments) + ["--protocol-major", "1"])
        if retired.returncode == 0 or retired.stdout:
            raise SystemExit("retired protocol selector did not fail at zero stdout")

        supervised = CoreProcess(
            full_arguments,
            expected_header=full.header,
            read_length=12,
            paired_end=True,
            expected_skipped_fragment_count=0,
            stdout_buffer_bytes=128,
        )
        with supervised as running:
            supervised_batches = tuple(running.iter_batches())
        observed = ProtocolStream(supervised.header, supervised_batches, supervised.trailer)
        if _common_signature(observed) != _common_signature(full):
            raise SystemExit("supervised batch transport changed protocol columns")
        if not all(batch.raw_payload.readonly for batch in supervised_batches):
            raise SystemExit("supervised payload is mutable")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
