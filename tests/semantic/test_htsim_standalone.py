"""htsim on its own: only --reference is required, and --output writes WGS/WES/TS reads
(--format: fastq.gz, fastq, bam, sam)."""

import gzip
import json
from pathlib import Path
import subprocess
import unittest

from bsreadsim import __version__
from tests.helpers.bam_reader import read_bam
from tests.semantic import support

FRAGMENTS = 3_000


def htsim(*arguments, check: bool = True) -> subprocess.CompletedProcess:
    completed = subprocess.run(
        [str(support.HTSIM), *(str(argument) for argument in arguments)],
        capture_output=True, cwd=str(support.WORK), check=False)
    if check and completed.returncode != 0:
        raise AssertionError("htsim failed:\n" + completed.stderr.decode())
    return completed


def fastq_run(name, *options):
    prefix = support.WORK / name
    htsim("--reference", support.REFERENCE, "--technology", "WGS", "--fragments", FRAGMENTS,
          "--output", prefix, *options)
    return prefix


def records(path):
    with gzip.open(path, "rt") as handle:
        lines = handle.read().splitlines()
    return [lines[i:i + 4] for i in range(0, len(lines), 4)]


def pairs(prefix):
    return list(zip(records(str(prefix) + ".R1.fastq.gz"), records(str(prefix) + ".R2.fastq.gz"),
                    strict=True))


def fragment_ends(name, read1, read2):
    """The reference bases the forward and the reverse read cover, read-oriented."""
    contig, start, end = support.envelope(name[1:-2])
    sequence = support.GENOME[contig]
    return sequence[start:start + read1], support.reverse_complement(sequence[end - read2:end])


def core_fields(record):
    return (record.query_name, record.flag, record.reference_id, record.position, record.mapq,
            record.cigar, record.next_reference_id, record.next_position, record.template_length,
            record.sequence, record.quality)


class StandaloneReadsTests(unittest.TestCase):
    def test_error_free_reads_are_the_fragment_ends(self):
        prefix = fastq_run("htsim-exact", "--error-rate", 0, "--phred", 30)
        read1_forward = 0
        for ordinal, (first, second) in enumerate(pairs(prefix)):
            self.assertEqual((first[0][-2:], second[0][-2:]), ("/1", "/2"))
            self.assertEqual(first[0][:-2], second[0][:-2])
            self.assertTrue(first[0].endswith(":{:x}/1".format(ordinal)), first[0])
            _contig, start, end = support.envelope(first[0][1:-2])
            self.assertTrue(100 <= end - start <= 1000)
            forward, reverse = fragment_ends(first[0], 100, 100)
            self.assertIn((first[1], second[1]), ((forward, reverse), (reverse, forward)))
            read1_forward += first[1] == forward
            for record in (first, second):
                self.assertEqual(record[3], "?" * 100)
        self.assertTrue(0.4 < read1_forward / FRAGMENTS < 0.6, read1_forward)

    def test_threads_do_not_change_the_reads(self):
        one = fastq_run("htsim-t1", "--threads", 1)
        four = fastq_run("htsim-t4", "--threads", 4)
        self.assertEqual(pairs(one), pairs(four))

    def test_substitutions_follow_the_error_rate(self):
        prefix = fastq_run("htsim-errors", "--error-rate", 0.02, "--read-length", 150,
                           "--insert-min", 150)
        errors = bases = 0
        for first, second in pairs(prefix):
            forward, reverse = fragment_ends(first[0], 150, 150)
            truth = (forward, reverse) if sum(a != b for a, b in zip(first[1], forward, strict=True)) \
                < sum(a != b for a, b in zip(first[1], reverse, strict=True)) else (reverse, forward)
            for read, expected in zip((first[1], second[1]), truth, strict=True):
                for observed, base in zip(read, expected, strict=True):
                    if base != "N":
                        bases += 1
                        errors += observed != base
        rate, sd = errors / bases, (0.02 * 0.98 / bases) ** 0.5
        self.assertLess(abs(rate - 0.02), 4 * sd, rate)

    def test_single_end_and_mate_lengths(self):
        single = fastq_run("htsim-single", "--paired-end", "false")
        self.assertEqual(len(records(str(single) + ".R1.fastq.gz")), FRAGMENTS)
        self.assertFalse((support.WORK / "htsim-single.R2.fastq.gz").exists())
        mates = fastq_run("htsim-mates", "--read-length", "120,80", "--insert-min", 150)
        for first, second in pairs(mates):
            self.assertEqual((len(first[1]), len(second[1])), (120, 80))
            self.assertEqual((len(first[3]), len(second[3])), (120, 80))

    def test_misuse_is_rejected(self):
        reference = ("--reference", support.REFERENCE)
        cases = (
            ((*reference, "--fragments", 10, "--output", "x"), "bisulfite reads come from"),
            ((*reference, "--technology", "WGS", "--fragments", 10, "--phred", 30), "need --output"),
            ((*reference, "--technology", "WGS", "--fragments", 10, "--format", "bam"),
             "need --output"),
            ((*reference, "--technology", "WGS", "--output", "x"), "one of --fragments and --depth"),
            ((*reference, "--technology", "WGS", "--output", "x", "--format", "cram"),
             "must be one of fastq.gz, fastq, bam, sam"),
            ((*reference, "--technology", "WGS", "--output", "x", "--format", "bam,bam"), "repeats"),
            (("--technology", "WGS", "--fragments", 10), "--reference is required"),
            (("variant-catalog", *reference, "--technology", "WGS", "--output", "x.vcf.gz",
              "--format", "bam"), "does not write reads"),
        )
        for arguments, message in cases:
            completed = htsim(*arguments, check=False)
            self.assertNotEqual(completed.returncode, 0, arguments)
            self.assertIn(message, completed.stderr.decode(), arguments)

    def test_bam_records_are_the_packages(self):
        # the same fragments (same options and seeds) without errors: the package's BAM
        # and htsim's agree in every field but the package's RG and truth tags
        variants = ("--seed-mut", 11, "--seed-phase", 12, "--mutation-rate", 0.01,
                    "--indel-fraction", 0.5)
        package = support.simulate("htsim-bam-package", "wgs", "-n", 2 * FRAGMENTS, *variants,
                                   "-e", 0, "-f", "bam", seed="77")
        htsim("--reference", support.REFERENCE, "--technology", "WGS", "--fragments", FRAGMENTS,
              "--seed", 77, *variants, "--error-rate", 0, "--output", support.WORK / "htsim",
              "--format", "bam")
        own = read_bam(support.WORK / "htsim.bam")
        self.assertEqual(own.references, package.bam.references)
        self.assertEqual(len(own.records), 2 * FRAGMENTS)
        indels = 0
        for mine, theirs in zip(own.records, package.bam.records, strict=True):
            self.assertEqual(core_fields(mine), core_fields(theirs), mine.query_name)
            self.assertEqual({tag: mine.tags[tag] for tag in ("AS", "MQ", "MC")},
                             {tag: theirs.tags[tag] for tag in ("AS", "MQ", "MC")})
            indels += any(op in "ID" for _length, op in mine.cigar)
        self.assertGreater(indels, 0)

    def test_every_format_holds_the_same_reads(self):
        prefix = fastq_run("htsim-all", "--error-rate", 0.02, "--format", "fastq.gz,fastq,bam,sam")
        for mate in ("R1", "R2"):
            with gzip.open(f"{prefix}.{mate}.fastq.gz", "rb") as compressed:
                self.assertEqual(compressed.read(), Path(f"{prefix}.{mate}.fastq").read_bytes())
        records = read_bam(f"{prefix}.bam").records
        reads = [read for pair in pairs(prefix) for read in pair]
        self.assertEqual(len(records), len(reads))
        for record, read in zip(records, reads, strict=True):
            self.assertEqual(record.query_name, read[0][1:-2])
            sequence = support.reverse_complement(record.sequence) if record.is_reverse \
                else record.sequence
            self.assertEqual(sequence, read[1])
        lines = [line.split("\t") for line in Path(f"{prefix}.sam").read_text().splitlines()
                 if not line.startswith("@")]
        self.assertEqual(
            [(line[0], int(line[1]), int(line[3]) - 1, line[5], line[9]) for line in lines],
            [(record.query_name, record.flag, record.position,
              "".join(f"{length}{op}" for length, op in record.cigar), record.sequence)
             for record in records])

    def test_fastq_gz_is_the_default(self):
        prefix = fastq_run("htsim-default")
        self.assertEqual(sorted(path.name for path in support.WORK.glob("htsim-default.*")),
                         ["htsim-default.R1.fastq.gz", "htsim-default.R2.fastq.gz"])
        self.assertEqual(len(pairs(prefix)), FRAGMENTS)

    def test_the_stream_needs_only_a_reference_and_an_amount(self):
        completed = htsim("--reference", support.REFERENCE, "--fragments", 10)
        header = json.loads(completed.stdout.split(b"\n", 1)[0])
        self.assertEqual((header["version"], header["technology"], header["read_lengths"]),
                         (__version__, "WGBS", [100, 100]))


if __name__ == "__main__":
    unittest.main()
