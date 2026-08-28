"""Exercise MethDB v2 sidecar generation, reload, and model fallback."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile

from bsreadsim.run.manifest import verify_complete_manifest


def _run(root: Path, core: Path, *arguments: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "bsreadsim", "run", "wgbs", *arguments, "--core", str(core)],
        cwd=str(root),
        env=os.environ.copy(),
        check=False,
        capture_output=True,
        text=True,
    )


def _manifest(completed: subprocess.CompletedProcess) -> dict:
    if completed.returncode != 0:
        raise SystemExit(
            "MethDB run failed: status={} stderr={!r}".format(
                completed.returncode, completed.stderr
            )
        )
    path = Path(completed.stdout.strip())
    document = json.loads(path.read_text(encoding="utf-8"))
    verify_complete_manifest(document)
    return document


def _export(
    root: Path, core: Path, input_path: Path, output_path: Path, *arguments: str
) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            sys.executable,
            "-m",
            "bsreadsim",
            "export",
            "methdb",
            "-i",
            str(input_path),
            "-o",
            str(output_path),
            *arguments,
            "--core",
            str(core),
        ],
        cwd=str(root),
        env=os.environ.copy(),
        check=False,
        capture_output=True,
        text=True,
    )


def main() -> int:
    if len(sys.argv) != 2:
        raise SystemExit("usage: check_methdb.py CORE_EXECUTABLE")
    core = Path(sys.argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-methdb-") as temporary:
        root = Path(temporary).resolve()
        (root / "reference.fa").write_bytes(
            b">chr1\n" + b"ACGT" * 300 + b"\n"
        )
        (root / "variants.vcf").write_text(
            "##fileformat=VCFv4.3\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            "chr1\t25\tlinked\tA\tT\t.\tPASS\t.\tGT\t1|0\n",
            encoding="ascii",
        )
        common = (
            "-r", "reference.fa",
            "-n", "16",
            "--seed-meth", "23",
            "--mutation-rate", "0",
            "--read-length", "50",
            "--insert-mean", "150",
            "--insert-sd", "0",
            "--max-ambiguous-fraction", "0",
            "--fragment-realization",
            "--format", "bam",
        )
        exported = _run(
            root,
            core,
            *common,
            "-o", "exported",
            "--seed", "17",
            "--save-methdb",
        )
        if exported.stderr:
            raise SystemExit(f"MethDB export emitted unexpected stderr: {exported.stderr!r}")
        first = _manifest(exported)
        snapshot = root / "exported" / "truth" / "sim.methdb"
        if not snapshot.is_file() or snapshot.stat().st_size == 0:
            raise SystemExit("MethDB export did not publish a snapshot")
        snapshot_bytes = snapshot.read_bytes()
        if not snapshot_bytes.startswith(b"methdb\x02"):
            raise SystemExit("MethDB snapshot omitted its format version")
        snapshot_sha256 = hashlib.sha256(snapshot_bytes).hexdigest()
        if {item["role"] for item in first["outputs"]} != {"bam", "truth.methdb"}:
            raise SystemExit("rich BAM MethDB run emitted the wrong output set")
        if first["details"]["randomness"]["methylation_seed"] != "23":
            raise SystemExit("MethDB manifest lost its methylation seed")
        first_methdb = [
            item for item in first["outputs"] if item["role"] == "truth.methdb"
        ]
        if (
            len(first_methdb) != 1
            or first_methdb[0]["sha256"] != snapshot_sha256
        ):
            raise SystemExit("MethDB export run did not record its exact snapshot")

        compressed_bed = root / "decoded.bed.gz"
        decoded = _export(root, core, snapshot, compressed_bed)
        if decoded.returncode != 0 or decoded.stderr:
            raise SystemExit(
                "compressed MethDB BED export failed: {!r}".format(decoded.stderr)
            )
        bed_bytes = gzip.decompress(compressed_bed.read_bytes())
        bed_text = bed_bytes.decode("utf-8")
        if not bed_text.startswith("#format\tmethdb-bed-v2\n"):
            raise SystemExit("MethDB BED export lost its format marker")
        if "#columns\tchrom\tchromStart\tchromEnd\tname\tscore\tstrand\t" not in bed_text:
            raise SystemExit("MethDB BED export lost its column contract")
        data_rows = [line.split("\t") for line in bed_text.splitlines() if line and not line.startswith("#")]
        if not data_rows or any(len(row) != 16 for row in data_rows):
            raise SystemExit("MethDB BED export emitted invalid extended BED rows")
        if any(int(row[2]) != int(row[1]) + 1 for row in data_rows):
            raise SystemExit("MethDB BED export did not emit one-base intervals")

        repeated_bed = root / "decoded-repeated.bed.gz"
        repeated = _export(root, core, snapshot, repeated_bed)
        if repeated.returncode != 0 or repeated_bed.read_bytes() != compressed_bed.read_bytes():
            raise SystemExit("compressed MethDB BED export is not deterministic")

        plain_bed = root / "decoded.bed"
        plain = _export(
            root, core, snapshot, plain_bed, "--no-compression"
        )
        if plain.returncode != 0 or plain.stderr:
            raise SystemExit(
                "plain MethDB BED export failed: {!r}".format(plain.stderr)
            )
        if plain_bed.read_bytes() != bed_bytes:
            raise SystemExit("plain and compressed MethDB BED contents disagree")

        wrong_suffix = root / "wrong.bed"
        rejected_export = _export(root, core, snapshot, wrong_suffix)
        if (
            rejected_export.returncode == 0
            or "must end in .bed.gz" not in rejected_export.stderr
            or wrong_suffix.exists()
        ):
            raise SystemExit("compressed MethDB export accepted a .bed output path")

        loaded = _run(
            root,
            core,
            *common,
            "-o", "loaded",
            "--seed", "999",
            "--methdb", str(snapshot),
            "--methylation-model", "bilstm",
        )
        if "falling back to BernoulliStateModel" not in loaded.stderr:
            raise SystemExit("BiLSTM request omitted its explicit fallback warning")
        second = _manifest(loaded)
        state_model = second["details"]["models"]["methylation_state"]
        if state_model != {
            "contract": "bernoulli-site",
            "effective": "bernoulli",
            "requested": "bilstm",
        }:
            raise SystemExit("BiLSTM fallback is not auditable in the manifest")
        if second["details"]["randomness"]["master_seed"] != "999":
            raise SystemExit("loaded MethDB incorrectly froze the simulation seed")
        second_methdb = [
            item for item in second["inputs"] if item["role"] == "input.methdb"
        ]
        if (
            len(second_methdb) != 1
            or second_methdb[0]["sha256"] != snapshot_sha256
            or second_methdb[0].get("format_version") != 2
        ):
            raise SystemExit("loaded run did not authenticate the fixed MethDB")

        rebound = _run(
            root,
            core,
            "-r", "reference.fa",
            "-o", "rebound",
            "-n", "2",
            "--seed", "999",
            "--seed-meth", "24",
            "--mutation-rate", "0",
            "--read-length", "50",
            "--insert-mean", "150",
            "--insert-sd", "0",
            "--format", "bam",
            "--methdb", str(snapshot),
        )
        rebound_manifest = _manifest(rebound)
        if rebound_manifest["details"]["randomness"]["methylation_seed"] != "24":
            raise SystemExit("MethDB v2 incorrectly froze an inert build seed")

        variant_common = (
            "-r", "reference.fa",
            "-n", "16",
            "--seed-meth", "23",
            "--read-length", "50",
            "--insert-mean", "150",
            "--insert-sd", "0",
            "--max-ambiguous-fraction", "0",
            "--fragment-realization",
            "--format", "bam",
        )
        variant_exported = _run(
            root,
            core,
            *variant_common,
            "-o", "variant-exported",
            "--seed", "17",
            "--vcf", "variants.vcf",
            "--save-methdb",
        )
        _manifest(variant_exported)
        variant_snapshot = (
            root / "variant-exported" / "truth" / "sim.methdb"
        )
        variant_bed = root / "variant.methdb.bed"
        variant_decoded = _export(
            root,
            core,
            variant_snapshot,
            variant_bed,
            "--no-compression",
        )
        if variant_decoded.returncode != 0 or (
            "#variant\tchr1\t0\t24\t25\tSNV\tA\tT\t1\tlinked\tvcf\n"
            not in variant_bed.read_text(encoding="utf-8")
        ):
            raise SystemExit("MethDB v2 did not preserve embedded VCF authority")

        variant_loaded = _run(
            root,
            core,
            *variant_common,
            "-o", "variant-loaded",
            "--seed", "18",
            "--methdb", str(variant_snapshot),
        )
        _manifest(variant_loaded)

        mixed_authority = _run(
            root,
            core,
            *variant_common,
            "-o", "mixed-authority",
            "--seed", "18",
            "--methdb", str(variant_snapshot),
            "--vcf", "variants.vcf",
        )
        if (
            mixed_authority.returncode == 0
            or "embeds variants" not in mixed_authority.stderr
        ):
            raise SystemExit("MethDB v2 accepted two variant authorities")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
