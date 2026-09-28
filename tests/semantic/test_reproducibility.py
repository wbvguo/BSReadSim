"""Fixed seeds reproduce outputs; the sample layer is independent of sequencing; SAM and
BAM hold the same records."""

import unittest

from tests.semantic import support

SAMPLE_SEEDS = ("--seed-mut", 101, "--seed-phase", 102, "--seed-meth", 103)


def fastq_run(name, *options, seed=support.SEED, read_length=100, threads=1):
    return support.simulate(
        name, "wgbs", "-n", 6_000, "-l", read_length, "--threads", threads,
        "--insert-mean", 300, "--insert-sd", 30, "--insert-min", 150,
        "--insert-max", 500, "--mutation-rate", 0.003, "--indel-fraction", 0.2,
        "--format", "fastq", "--save-truth", *SAMPLE_SEEDS, *options, seed=seed,
    )


def truth(run):
    sites = [
        (s.contig, s.position, s.strand, s.set, s.origin, s.context, s.probability_u16)
        for s in run.methdb_sites
    ]
    return run.variants, sites


def read_fastq(path):
    lines = path.read_text().splitlines()
    return [lines[i:i + 4] for i in range(0, len(lines), 4)]


def sam_tag(text):
    """A SAM tag as the BAM reader decodes it: i -> int, B -> tuple of ints, others text."""
    tag, kind, value = text.split(":", 2)
    if kind == "i":
        return tag, int(value)
    if kind == "B":
        return tag, tuple(int(item) for item in value.split(",")[1:])
    return tag, value


class ReproducibilityTests(unittest.TestCase):
    def test_thread_count_does_not_change_output(self):
        single = fastq_run("repro-t1", threads=1)
        multi = fastq_run("repro-t4", threads=4)
        for suffix in (".R1.fastq", ".R2.fastq"):
            self.assertEqual(
                single.file(suffix).read_bytes(), multi.file(suffix).read_bytes())
        self.assertEqual(truth(single), truth(multi))

    def test_master_seed_varies_reads_but_not_the_sample(self):
        base = fastq_run("repro-t1", threads=1)
        other = fastq_run("repro-seed", seed="9001")
        self.assertNotEqual(
            base.file(".R1.fastq").read_bytes(), other.file(".R1.fastq").read_bytes())
        self.assertEqual(truth(base), truth(other))

    def test_sequencing_settings_do_not_change_the_sample(self):
        base = fastq_run("repro-t1", threads=1)
        shorter = fastq_run("repro-length", read_length=80)
        self.assertEqual(truth(base), truth(shorter))

    def test_fastq_records_pair_and_name_fragments(self):
        run = fastq_run("repro-t1", threads=1)
        read1 = read_fastq(run.file(".R1.fastq"))
        read2 = read_fastq(run.file(".R2.fastq"))
        self.assertEqual(len(read1), 3_000)
        self.assertEqual(len(read2), 3_000)
        for first, second in zip(read1, read2, strict=True):
            self.assertEqual(first[0][:-2], second[0][:-2])
            self.assertEqual((first[0][-2:], second[0][-2:]), ("/1", "/2"))
            for record in (first, second):
                self.assertEqual(record[2], "+")
                self.assertEqual(len(record[1]), 100)
                self.assertEqual(len(record[3]), 100)
            contig, start, end = support.envelope(first[0][1:-2])
            self.assertIn(contig, support.GENOME)
            self.assertTrue(150 <= end - start <= 500 + 4)

    def test_sam_holds_the_bam_records(self):
        options = ("-n", 4_000, "--mutation-rate", 0.003, "--indel-fraction", 0.5,
                   "--fragment-realization", *SAMPLE_SEEDS)
        bam = support.simulate("repro-bam", "wgbs", *options, "--format", "bam").bam
        sam_run = support.simulate("repro-sam", "wgbs", *options, "--format", "sam")
        text = sam_run.file(".sam").read_text().splitlines()
        self.assertEqual(
            [tuple(line.split("\t")[1:3]) for line in text if line.startswith("@SQ")],
            [(f"SN:{name}", f"LN:{length}") for name, length in bam.references])
        lines = [line.split("\t") for line in text if not line.startswith("@")]
        self.assertEqual(len(lines), len(bam.records))
        for fields, record in zip(lines, bam.records, strict=True):
            self.assertEqual(
                (fields[0], int(fields[1]), int(fields[3]) - 1, int(fields[4]), fields[5],
                 int(fields[7]) - 1, int(fields[8]), fields[9],
                 tuple(ord(q) - 33 for q in fields[10])),
                (record.query_name, record.flag, record.position, record.mapq,
                 "".join(f"{length}{op}" for length, op in record.cigar),
                 record.next_position, record.template_length, record.sequence,
                 record.quality))
            tags = dict(sam_tag(field) for field in fields[11:])
            self.assertEqual({tag: value for tag, value in tags.items() if tag != "RG"},
                             {tag: value for tag, value in record.tags.items() if tag != "RG"})


if __name__ == "__main__":
    unittest.main()
