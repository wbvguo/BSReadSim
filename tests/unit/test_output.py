"""The output of a run: staging and publication of its files, and its manifest."""

import gzip
import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from bsreadsim import __version__
from bsreadsim.errors import BSReadSimError
from bsreadsim.htsim import Contig, Header, Summary
from bsreadsim.output import MANIFEST_VERSION, OutputFile, RunFiles, build_manifest
from bsreadsim.settings import Settings


class RunFilesTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.directory = Path(self.temporary.name) / "out"

    def tearDown(self):
        self.temporary.cleanup()

    def visible(self):
        return sorted(str(path.relative_to(self.directory))
                      for path in self.directory.rglob("*") if path.is_file()
                      and not path.relative_to(self.directory).parts[0].startswith("."))

    def test_files_appear_only_when_published(self):
        with RunFiles(self.directory, "sample", paired_end=True, format="fastq") as files:
            files.write(2, {"read1": b"@r1\n", "read2": b"@r2\n"})
            files.write(1, {"read1": b"@r3\n", "read2": b"@r4\n"})
            truth = files.truth_path("sample.variants.vcf.gz")
            truth.write_bytes(gzip.compress(b"##fileformat=VCFv4.2\n#CHROM\n"
                                            b"chr1\t5\t.\tA\tG\nchr1\t9\t.\tC\tT\n"))
            files.add_truth("truth.vcf", truth.name)
            self.assertEqual(self.visible(), [])
            outputs = files.finish()
            files.publish({"status": "complete"})
        self.assertEqual(self.visible(), ["sample.R1.fastq", "sample.R2.fastq",
                                          "sample.manifest.json",
                                          "truth/sample.variants.vcf.gz"])
        self.assertEqual([(item.role, item.record_count) for item in outputs],
                         [("read1", 3), ("read2", 3), ("truth.vcf", 2)])
        self.assertEqual(outputs[0].sha256, hashlib.sha256(b"@r1\n@r3\n").hexdigest())
        self.assertEqual(outputs[0].path, self.directory / "sample.R1.fastq")
        manifest = json.loads((self.directory / "sample.manifest.json").read_text())
        self.assertEqual(manifest, {"status": "complete"})
        self.assertFalse(any(path.name.startswith(".") for path in self.directory.iterdir()))

    def test_failure_leaves_nothing_behind(self):
        with self.assertRaises(RuntimeError):
            with RunFiles(self.directory, "sample", paired_end=False, format="fastq.gz") as files:
                files.write(1, {"read1": b"x"})
                raise RuntimeError("stop")
        self.assertEqual(list(self.directory.iterdir()), [])

    def test_existing_outputs_are_never_replaced(self):
        self.directory.mkdir()
        (self.directory / "sample.R1.fastq.gz").write_bytes(b"old")
        with self.assertRaises(BSReadSimError):
            RunFiles(self.directory, "sample", paired_end=False, format="fastq.gz")
        self.assertEqual((self.directory / "sample.R1.fastq.gz").read_bytes(), b"old")
        self.assertEqual(len(list(self.directory.iterdir())), 1)


def manifest(settings: Settings):
    header = Header(core_version=__version__,
                    contigs=(Contig("chr1", 50, "0" * 32), Contig("chr2", 30, "0" * 32)))
    summary = Summary(fragment_count=3, mate_count=6, template_base_count=900,
                      methylation_site_count=40, skipped_fragment_count=1,
                      per_contig_fragment_counts=(2, 1))
    outputs = [OutputFile("read1", Path("/out/sim.R1.fastq.gz"), 10, "a" * 64, 3),
               OutputFile("read2", Path("/out/sim.R2.fastq.gz"), 12, "b" * 64, 3)]
    inputs = [{"role": "reference", "format": "fasta", "path": "/ref.fa", "size_bytes": 80,
               "sha256": "c" * 64}]
    argv = ["bsreadsim", "run", "wgbs", "-r", "/ref.fa", "-o", "/out"]
    return build_manifest(settings, argv, header, summary, inputs, outputs, run_id="run",
                          methylation_model="bernoulli")


class ManifestTests(unittest.TestCase):
    def test_sections(self):
        settings = Settings(reference=Path("/ref.fa"), output=Path("/out"), seed=1 << 63,
                            meth_model="bilstm")
        document = manifest(settings)
        self.assertEqual((document["version"], document["status"], document["run_id"]),
                         (MANIFEST_VERSION, "complete", "run"))
        self.assertEqual(document["summary"]["read_count"], 6)
        self.assertEqual(document["summary"]["read_base_count"], 600)
        self.assertEqual(document["summary"]["output_size_bytes"], 22)
        self.assertEqual([item["role"] for item in document["outputs"]], ["read1", "read2"])
        details = document["details"]
        self.assertEqual(details["configuration"], settings.as_dict())
        self.assertEqual(details["configuration_sha256"], settings.sha256)
        self.assertEqual(details["randomness"]["master_seed"], str(1 << 63))
        self.assertEqual(details["models"]["methylation_state"],
                         {"requested": "bilstm", "effective": "bernoulli"})
        self.assertEqual([contig["fragment_count"] for contig in details["contigs"]], [2, 1])
        self.assertIn("--seed {}".format(1 << 63), document["command"]["full_command"])

    def test_standard_assays_disable_methylation(self):
        document = manifest(Settings(technology="WGS", reference=Path("/ref.fa")))
        self.assertEqual(document["details"]["models"]["methylation_state"],
                         {"requested": "disabled", "effective": "disabled"})
        self.assertNotIn("--seed-meth", document["command"]["full_command"])


if __name__ == "__main__":
    unittest.main()
