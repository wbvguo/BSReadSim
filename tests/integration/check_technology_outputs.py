"""Exercise bisulfite and standard technologies through annotated BAM."""

from __future__ import annotations

import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPOSITORY_ROOT))

from bsreadsim.run.manifest import verify_complete_manifest
from tests.integration.check_bam import _parse_bam


FIXTURE_ROOT = REPOSITORY_ROOT / "tests" / "fixtures"


def _environment() -> dict:
    return os.environ.copy()


def _qname_envelope(query_name: str):
    locus, ordinal_text = query_name.rsplit(":", 1)
    contig, interval = locus.rsplit(":", 1)
    left_text, right_text = interval.split("-", 1)
    return contig, int(left_text), int(right_text), int(ordinal_text, 16)


def _run(core: Path, root: Path, name: str, technology_arguments):
    output = root / name
    command = [
        sys.executable,
        "-m",
        "bsreadsim",
        "run",
        technology_arguments[0],
        "-r",
        str(root / "mock-reference.fa"),
        "-o",
        str(output),
        "-n",
        "512",
        "--seed",
        "20260813",
        "--read-length",
        "4",
        "--max-ambiguous-fraction",
        "0",
        "--mutation-rate",
        "0",
        "--phred",
        "35",
        "--error-rate",
        "0",
        "--threads",
        "4",
        "--prefix",
        "sample",
        "--format",
        "bam",
        "--fragment-summary",
        "--core",
        str(core),
        *technology_arguments[1:],
    ]
    if technology_arguments[0] in ("wgbs", "rrbs", "tbs"):
        command.extend(("--conversion-rate", "1"))
    completed = subprocess.run(
        command,
        cwd=str(root),
        env=_environment(),
        check=False,
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0 or completed.stderr:
        raise SystemExit(
            "{} pipeline failed: status={} stderr={!r}".format(
                name, completed.returncode, completed.stderr
            )
        )
    manifest = json.loads(
        Path(completed.stdout.strip()).read_text(encoding="utf-8")
    )
    verify_complete_manifest(manifest)
    if {item["role"] for item in manifest["outputs"]} != {"bam"}:
        raise SystemExit("{} emitted the wrong artifact roles".format(name))
    _, _, records, _ = _parse_bam(output / "sample.bam")
    if len(records) != 512:
        raise SystemExit("{} emitted the wrong BAM record count".format(name))
    return manifest, records


def _validate_standard(
    name: str, manifest: dict, records: list, references: tuple[str, ...]
) -> None:
    if manifest["summary"]["methylation_site_count"] != 0:
        raise SystemExit("{} emitted methylation sites".format(name))
    effective = manifest["details"]["configuration"]
    if effective["technology"] != name.upper():
        raise SystemExit("{} lost its technology identity".format(name))
    if manifest["details"]["models"]["methylation_state"] != {
        "effective": "disabled",
        "requested": "disabled",
    }:
        raise SystemExit("{} manifest enabled methylation".format(name))
    if effective["sequencing"]["conversion_rate"] != 0:
        raise SystemExit("{} retained bisulfite chemistry".format(name))
    tag_policy = manifest["details"]["alignment"]["tags"]
    if any(tag_policy[tag]["required"] for tag in ("XG", "XR", "YS")):
        raise SystemExit("{} required bisulfite-only BAM tags".format(name))
    for record in records:
        if any(tag in record["aux"] for tag in ("XG", "XR", "YS")):
            raise SystemExit("{} emitted bisulfite-only BAM tags".format(name))
        summary = record["aux"]["zf"][1]
        if ((summary[0] >> 4) & 0x7) != 2 or any(summary[1:9]):
            raise SystemExit("{} emitted bisulfite annotations".format(name))
        start = record["position"]
        end = start + record["reference_consumed"]
        expected = references[record["reference_id"]][start:end]
        if record["sequence"] != expected:
            raise SystemExit(
                "{} changed a reference read sequence: qname={} flag={} "
                "position={} observed={} expected={}".format(
                    name,
                    record["query_name"],
                    record["flag"],
                    record["position"],
                    record["sequence"],
                    expected,
                )
            )


def main(argv) -> int:
    if len(argv) != 2:
        raise SystemExit("usage: check_technology_outputs.py CORE_EXECUTABLE")
    core = Path(argv[1]).resolve(strict=True)
    with tempfile.TemporaryDirectory(prefix="bsreadsim-technologies-") as value:
        root = Path(value).resolve()
        shutil.copy2(FIXTURE_ROOT / "mock-reference.fa", root)
        shutil.copy2(FIXTURE_ROOT / "mock-targets.bed", root)

        _, rrbs = _run(
            core,
            root,
            "rrbs",
            (
                "rrbs",
                "--cut-site",
                "C|CGG",
                "--insert-min",
                "8",
                "--insert-mean",
                "8",
                "--insert-max",
                "8",
                "--insert-sd",
                "0",
            ),
        )
        rrbs_modes = set()
        for record in rrbs[::2]:
            contig, left, right, _ = _qname_envelope(record["query_name"])
            if contig != "chrMock" or right - left + 1 != 8:
                raise SystemExit("RRBS emitted a non-MspI fragment envelope")
            rrbs_modes.add((record["aux"]["zf"][1][0] >> 4) & 0x7)
        if rrbs_modes != {0, 1}:
            raise SystemExit(
                "directional RRBS omitted Watson or Crick fragments: {!r}".format(
                    rrbs_modes
                )
            )

        _, tbs = _run(
            core,
            root,
            "tbs",
            (
                "tbs",
                "--targets",
                str(root / "mock-targets.bed"),
                "--sampling",
                "score",
                "--fragment-center-stddev",
                "0",
                "--insert-mean",
                "12",
                "--insert-sd",
                "0",
            ),
        )
        expected_modes = {(11, 22): 0, (27, 38): 1}
        observed = set()
        for record in tbs[::2]:
            contig, left, right, _ = _qname_envelope(record["query_name"])
            envelope = (left, right)
            if contig != "chrMock" or envelope not in expected_modes:
                raise SystemExit("TBS selected an invalid or zero-weight target")
            mode = (record["aux"]["zf"][1][0] >> 4) & 0x7
            if mode != expected_modes[envelope]:
                raise SystemExit("TBS strand did not control conversion mode")
            observed.add(envelope)
        if observed != set(expected_modes):
            raise SystemExit("TBS weighted fixture did not exercise both targets")

        references = []
        for line in (root / "mock-reference.fa").read_text(
            encoding="ascii"
        ).splitlines():
            if line.startswith(">"):
                references.append("")
            else:
                references[-1] += line.strip()
        standard_cases = (
            (
                "wgs",
                "--insert-mean",
                "8",
                "--insert-sd",
                "0",
            ),
            (
                "wes",
                "--targets",
                str(root / "mock-targets.bed"),
                "--fragment-center-stddev",
                "0",
                "--insert-mean",
                "12",
                "--insert-sd",
                "0",
            ),
            (
                "ts",
                "--targets",
                str(root / "mock-targets.bed"),
                "--fragment-center-stddev",
                "0",
                "--insert-mean",
                "12",
                "--insert-sd",
                "0",
            ),
        )
        for arguments in standard_cases:
            name = arguments[0]
            manifest, records = _run(core, root, name, arguments)
            _validate_standard(name, manifest, records, tuple(references))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
