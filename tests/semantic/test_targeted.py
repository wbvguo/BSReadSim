"""Targeted assays place fragments around strand-aware capture targets."""

from collections import Counter, defaultdict
import statistics
import unittest

from tests.semantic import support

CENTER_SD = 50
INSERT = ("--insert-mean", 250, "--insert-sd", 20,
          "--insert-min", 150, "--insert-max", 400)


def make_targets():
    targets = []
    for k, start in enumerate(range(10_000, 230_000, 13_000)):
        targets.append(("chrA", start, start + 200, "+-."[k % 3], (1, 3, 0)[k // 3 % 3]))
    for start in range(8_000, 118_000, 13_000):
        if 55_000 < start < 64_000:
            continue
        index = len(targets)
        targets.append(("chrB", start, start + 200, "+-."[index % 3], (1, 3, 0)[index // 3 % 3]))
    return targets


TARGETS = make_targets()
TARGET_BED = support.write_text("targets.bed", "".join(
    "{}\t{}\t{}\tt{}\t{}\t{}\n".format(contig, start, end, k, score, strand)
    for k, (contig, start, end, strand, score) in enumerate(TARGETS)
))


def nearest_target(name):
    contig, start, end = support.envelope(name)
    centre = (start + end) / 2
    index = min(
        (k for k, t in enumerate(TARGETS) if t[0] == contig),
        key=lambda k: abs((TARGETS[k][1] + TARGETS[k][2]) / 2 - centre),
    )
    target = TARGETS[index]
    return index, centre - (target[1] + target[2]) / 2


def targeted_run(name, assay, reads, *options):
    return support.simulate(
        name, assay, "-n", reads, "-l", 100, *INSERT, "--threads", 4,
        "--targets", TARGET_BED, "--center-sd", CENTER_SD,
        "--mutation-rate", 0, "-e", 0, "--format", "bam", *options,
    )


class TbsTests(unittest.TestCase):
    def test_uniform_sampling_treats_targets_equally(self):
        run = targeted_run("tbs-uniform", "tbs", 24_000, "--conversion-rate", 1)
        counts = Counter(nearest_target(name)[0] for name in run.fragments)
        expected = len(run.fragments) / len(TARGETS)
        for index in range(len(TARGETS)):
            self.assertAlmostEqual(counts[index], expected, delta=0.2 * expected)

    def test_fragment_centres_scatter_around_targets(self):
        run = targeted_run("tbs-uniform", "tbs", 24_000, "--conversion-rate", 1)
        offsets = [nearest_target(name)[1] for name in run.fragments]
        self.assertAlmostEqual(statistics.fmean(offsets), 0.0, delta=2.0)
        self.assertAlmostEqual(statistics.pstdev(offsets), CENTER_SD, delta=3.0)

    def test_capture_strand_selects_library_strand(self):
        run = targeted_run("tbs-uniform", "tbs", 24_000, "--conversion-rate", 1)
        strands = defaultdict(Counter)
        for name, records in run.fragments.items():
            index, _offset = nearest_target(name)
            strands[TARGETS[index][3]][records[0].tags["YS"]] += 1
        self.assertEqual(set(strands["+"]), {"OT"})
        self.assertEqual(set(strands["-"]), {"OB"})
        unknown = strands["."]
        self.assertAlmostEqual(
            unknown["OT"] / (unknown["OT"] + unknown["OB"]), 0.5, delta=0.05)

    def test_reads_match_converted_reference(self):
        run = targeted_run("tbs-uniform", "tbs", 24_000, "--conversion-rate", 1)
        audit = support.audit_bisulfite_reads(run, None)
        self.assertEqual(audit.base_mismatches, 0, audit.examples)

    def test_score_sampling_follows_target_scores(self):
        run = targeted_run(
            "tbs-score", "tbs", 24_000, "--conversion-rate", 1, "--sampling", "score")
        counts = Counter(nearest_target(name)[0] for name in run.fragments)
        by_score = {0: 0, 1: 0, 3: 0}
        sizes = Counter(t[4] for t in TARGETS)
        for index, count in counts.items():
            by_score[TARGETS[index][4]] += count
        self.assertEqual(by_score[0], 0)
        ratio = (by_score[3] / sizes[3]) / (by_score[1] / sizes[1])
        self.assertAlmostEqual(ratio, 3.0, delta=0.3)


class VariantTargetedTests(unittest.TestCase):
    def test_fragments_follow_their_haplotype_around_targets(self):
        run = support.simulate(
            "tbs-variants", "tbs", "-n", 24_000, "-l", 100, *INSERT, "--threads", 4,
            "--targets", TARGET_BED, "--center-sd", CENTER_SD,
            "--mutation-rate", 0.006, "--indel-fraction", 0.4,
            "--conversion-rate", 1, "-e", 0, "--format", "bam", "--save-vcf",
        )
        problems, operations = support.reconstruction_problems(run, run.haplotypes)
        self.assertEqual(problems[:10], [])
        self.assertGreater(operations["I"], 100)
        self.assertGreater(operations["D"], 100)
        offsets = [nearest_target(name)[1] for name in run.fragments]
        self.assertAlmostEqual(statistics.fmean(offsets), 0.0, delta=3.0)
        self.assertAlmostEqual(statistics.pstdev(offsets), CENTER_SD, delta=4.0)


class NonBisulfiteTargetedTests(unittest.TestCase):
    def test_wes_and_ts_use_targets_without_conversion(self):
        for assay in ("wes", "ts"):
            run = targeted_run(assay, assay, 4_000)
            self.assertEqual(len(run.bam.records), 4_000)
            for name, records in run.fragments.items():
                _index, offset = nearest_target(name)
                self.assertLess(abs(offset), 6 * CENTER_SD + 10)
                for record in records:
                    self.assertNotIn("XG", record.tags)
                    sequence = support.GENOME[run.contig_of(record)]
                    self.assertEqual(
                        record.sequence,
                        sequence[record.position:record.reference_end],
                    )


if __name__ == "__main__":
    unittest.main()
