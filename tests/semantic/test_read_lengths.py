"""Paired mates may have different read lengths (--read-length READ1,READ2)."""

import json
import unittest

from tests.semantic import support

READ1, READ2 = 120, 70
READS = 8_000


def bam_run():
    # error-free, all four strands, and indels, so placement and CIGARs are exercised
    return support.simulate(
        "mate-lengths", "wgbs", "-n", READS, "-l", "{},{}".format(READ1, READ2),
        "--threads", 2, "--insert-mean", 300, "--insert-sd", 30,
        "--insert-min", 150, "--insert-max", 500, "--mutation-rate", 0.004,
        "--indel-fraction", 0.5, "--undirectional", "-e", 0,
        "--format", "bam", "--save-vcf",
    )


def query_length(record) -> int:
    return sum(length for length, operation in record.cigar if operation in "MI")


class MateLengthTests(unittest.TestCase):
    def test_each_mate_has_its_length(self):
        run = bam_run()
        self.assertEqual(len(run.bam.records), READS)
        for record in run.bam.records:
            length = READ1 if record.mate_number == 1 else READ2
            self.assertEqual(len(record.sequence), length, record.query_name)
            self.assertEqual(len(record.quality), length)
            self.assertEqual(len(record.tags["zt"]), length)
            self.assertEqual(query_length(record), length)
            self.assertEqual(record.tags["AS"], length)

    def test_reads_start_at_the_fragment_ends(self):
        run = bam_run()
        for name, records in run.fragments.items():
            _contig, start, end = support.envelope(name)
            forward = [r for r in records if not r.is_reverse]
            reverse = [r for r in records if r.is_reverse]
            self.assertEqual((len(forward), len(reverse)), (1, 1), name)
            self.assertEqual(forward[0].position, start, name)
            self.assertEqual(reverse[0].reference_end, end, name)
            self.assertEqual({abs(r.template_length) for r in records}, {end - start}, name)

    def test_reads_reconstruct_their_haplotype(self):
        run = bam_run()
        problems, operations = support.reconstruction_problems(run, run.haplotypes)
        self.assertEqual(problems, [])
        self.assertGreater(operations["I"] + operations["D"], 0)

    def test_fastq_mates_and_read_bases(self):
        run = support.simulate(
            "mate-lengths-fastq", "wgbs", "-n", 2_000, "-l", "{},{}".format(READ1, READ2),
            "--insert-min", 150, "--format", "fastq")
        for suffix, length in ((".R1.fastq", READ1), (".R2.fastq", READ2)):
            lines = run.file(suffix).read_text().splitlines()
            self.assertEqual(len(lines), 4 * 1_000)
            self.assertEqual({len(line) for line in lines[1::4]}, {length})
            self.assertEqual({len(line) for line in lines[3::4]}, {length})
        summary = json.loads(run.file(".manifest.json").read_text())["summary"]
        self.assertEqual(summary["read_base_count"], 1_000 * (READ1 + READ2))

    def test_single_end_takes_one_length(self):
        completed = support.bsreadsim(
            "run", "wgbs", "-r", support.REFERENCE, "-o", support.WORK / "single-two-lengths",
            "-n", 100, "--single-end", "-l", "100,80", check=False)
        self.assertNotEqual(completed.returncode, 0)
        self.assertIn("one length for single-end", completed.stderr)

    def test_read_lengths_are_at_most_10000(self):
        for value in ("0", "10001", "100,10001"):
            completed = support.bsreadsim(
                "run", "wgbs", "-r", support.REFERENCE, "-o", support.WORK / "too-long",
                "-n", 100, "-l", value, "--insert-min", 20_000, "--insert-mean", 20_000,
                "--insert-max", 20_000, check=False)
            self.assertNotEqual(completed.returncode, 0, value)
            self.assertIn("--read-length must be in [1, 10000]", completed.stderr, value)


if __name__ == "__main__":
    unittest.main()
