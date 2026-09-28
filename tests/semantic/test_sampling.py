"""Whole-genome GC-profile sampling and depth-derived read counts."""

from collections import Counter
import unittest

from tests.semantic import support

GC_PROFILE = (0.0, 0.5, 0.2, 0.3, 0.0)
FRAGMENT_LENGTH = 200


class GcSamplingTests(unittest.TestCase):
    def test_fragment_gc_matches_target_profile(self):
        profile = support.write_text(
            "gc-profile.txt", "".join("{:g}\n".format(p) for p in GC_PROFILE))
        run = support.simulate(
            "wgs-gc", "wgs", "-n", 20_000, "-l", 50, "--threads", 4,
            "--sampling", "gc", "--gc-profile", profile,
            "--insert-mean", FRAGMENT_LENGTH, "--insert-sd", 0,
            "--insert-min", 100, "--insert-max", 300,
            "-e", 0, "--format", "bam",
        )
        bins = Counter()
        bin_count = len(GC_PROFILE)
        for name in run.fragments:
            contig, start, end = support.envelope(name)
            self.assertEqual(end - start, FRAGMENT_LENGTH)
            gc = sum(b in "CG" for b in support.GENOME[contig][start:end])
            bins[int(gc * (bin_count - 1) / FRAGMENT_LENGTH + 0.5)] += 1
        total = sum(bins.values())
        for index, p in enumerate(GC_PROFILE):
            self.assertAlmostEqual(bins[index] / total, p, delta=0.025, msg=str(index))

    def test_template_gc_matches_target_profile_on_haplotypes_with_indels(self):
        """With variants, the template (not its reference envelope) has the
        insert length and follows the GC profile."""
        profile = support.write_text(
            "gc-profile.txt", "".join("{:g}\n".format(p) for p in GC_PROFILE))
        run = support.simulate(
            "wgs-gc-variants", "wgs", "-n", 20_000, "-l", 50, "--threads", 4,
            "--sampling", "gc", "--gc-profile", profile,
            "--insert-mean", FRAGMENT_LENGTH, "--insert-sd", 0,
            "--insert-min", 100, "--insert-max", 300,
            "--mutation-rate", 0.01, "--indel-fraction", 0.5,
            "-e", 0, "--format", "bam", "--save-vcf",
        )
        haplotypes = {}
        bins = Counter()
        bin_count = len(GC_PROFILE)
        moved = 0
        for name, records in run.fragments.items():
            contig, start, end = support.envelope(name)
            key = (contig, support.haplotype_of(records[0]))
            if key not in haplotypes:
                sequence, origins = support.haplotype_sequence(contig, run.variants, key[1])
                offsets = {o: k for k, o in enumerate(origins) if not o[1]}
                haplotypes[key] = sequence, origins, offsets
            sequence, origins, offsets = haplotypes[key]
            first = offsets[(start, False)]
            template = sequence[first:first + FRAGMENT_LENGTH]
            self.assertEqual(origins[first + FRAGMENT_LENGTH - 1][0] + 1, end, name)
            moved += end - start != FRAGMENT_LENGTH
            gc = sum(b in "CG" for b in template)
            bins[int(gc * (bin_count - 1) / FRAGMENT_LENGTH + 0.5)] += 1
        self.assertGreater(moved, 500)
        total = sum(bins.values())
        for index, p in enumerate(GC_PROFILE):
            self.assertAlmostEqual(bins[index] / total, p, delta=0.025, msg=str(index))


class DepthTests(unittest.TestCase):
    def test_depth_sets_read_count_from_genome_size(self):
        run = support.simulate(
            "wgs-depth", "wgs", "-d", 3, "-l", 100, "--mutation-rate", 0,
            "--format", "bam",
        )
        genome_bases = sum(len(s) for s in support.GENOME.values())
        self.assertEqual(len(run.bam.records), genome_bases * 3 // 100)


if __name__ == "__main__":
    unittest.main()
