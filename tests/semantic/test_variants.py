"""Reads are drawn from a phased diploid genome with SNVs and short indels."""

from collections import Counter
import unittest

from tests.semantic import support

READS = 30_000
MUTATION_RATE = 0.006
INDEL_FRACTION = 0.4
EXTENSION = 0.3


def wgs_run():
    return support.simulate(
        "wgs-variants", "wgs",
        "-n", READS, "-l", 100, "--threads", 4,
        "--insert-mean", 300, "--insert-sd", 30,
        "--insert-min", 150, "--insert-max", 500,
        "--mutation-rate", MUTATION_RATE, "--indel-fraction", INDEL_FRACTION,
        "--indel-extension-probability", EXTENSION,
        "-e", 0, "-q", 35, "--format", "bam", "--fragment-summary", "--save-vcf",
    )


def reconstruction_problems(run):
    return support.reconstruction_problems(run, run.haplotypes)


class HaplotypeReadTests(unittest.TestCase):
    def test_reads_reconstruct_their_haplotype(self):
        problems, operations = reconstruction_problems(wgs_run())
        self.assertEqual(problems[:10], [])
        self.assertGreater(operations["I"], 200)
        self.assertGreater(operations["D"], 200)

    def test_mates_come_from_one_haplotype(self):
        for name, records in wgs_run().fragments.items():
            self.assertEqual(len(records), 2, name)
            self.assertEqual(len({support.haplotype_of(r) for r in records}), 1, name)

    def test_snv_alleles_follow_genotypes(self):
        run = wgs_run()
        snvs = {
            (v.contig, v.position): v for v in run.variants if v.kind == "snv"
        }
        counts = {1: [0, 0], 2: [0, 0]}  # zygosity -> [alt reads, reads]
        for record in run.bam.records:
            contig = run.contig_of(record)
            for query, position, operation in record.aligned_pairs():
                variant = snvs.get((contig, position)) if operation == "M" else None
                if variant is None:
                    continue
                entry = counts[len(variant.haplotypes)]
                entry[0] += record.sequence[query] == variant.alternate
                entry[1] += 1
        self.assertGreater(counts[1][1], 5_000)
        self.assertAlmostEqual(counts[1][0] / counts[1][1], 0.5, delta=0.03)
        self.assertEqual(counts[2][0], counts[2][1])

    def test_variant_bases_are_flagged(self):
        run = wgs_run()
        unflagged = 0
        for record in run.bam.records:
            key = (run.contig_of(record), support.haplotype_of(record))
            snvs = run.haplotypes.snvs.get(key, {})
            states = record.tags["zt"]
            for query, position, operation in record.aligned_pairs():
                if operation == "I" or (operation == "M" and position in snvs):
                    unflagged += not support.ZT_VALUES[states[query]] & 16
        self.assertEqual(unflagged, 0)

    def test_non_bisulfite_reads_carry_no_methylation_state(self):
        for record in wgs_run().bam.records:
            self.assertNotIn("XG", record.tags)
            self.assertNotIn("YS", record.tags)
            self.assertTrue(
                all(support.ZT_VALUES[c] & 15 == 0 for c in record.tags["zt"])
            )


class GeneratedVariantTests(unittest.TestCase):
    def test_indel_spectrum(self):
        variants = wgs_run().variants
        callable_bases = sum(
            len(s) - s.count("N") for s in support.GENOME.values()
        )
        expected = MUTATION_RATE * callable_bases
        self.assertAlmostEqual(len(variants), expected, delta=0.12 * expected)
        indels = [v for v in variants if v.kind != "snv"]
        self.assertAlmostEqual(len(indels) / len(variants), INDEL_FRACTION, delta=0.05)
        insertions = sum(v.kind == "insertion" for v in indels)
        self.assertAlmostEqual(insertions / len(indels), 0.5, delta=0.08)
        lengths = Counter(v.length for v in indels)
        self.assertLessEqual(max(lengths), 4)
        self.assertAlmostEqual(lengths[1] / len(indels), 1 - EXTENSION, delta=0.07)
        homozygous = sum(len(v.haplotypes) == 2 for v in variants)
        self.assertAlmostEqual(homozygous / len(variants), 1 / 3, delta=0.05)

    def test_variants_match_reference_and_do_not_overlap(self):
        variants = wgs_run().variants
        last_end = {}
        for variant in variants:
            sequence = support.GENOME[variant.contig]
            if variant.kind != "insertion":
                end = variant.position + len(variant.reference)
                self.assertEqual(sequence[variant.position:end], variant.reference)
                self.assertNotIn("N", variant.reference)
            else:
                end = variant.position
            self.assertGreaterEqual(variant.position, last_end.get(variant.contig, 0))
            last_end[variant.contig] = end


class VcfInputTests(unittest.TestCase):
    def test_vcf_variants_are_used_with_their_phase(self):
        source = wgs_run().directory / "truth" / "sim.variants.vcf.gz"
        run = support.simulate(
            "wgs-from-vcf", "wgs", "-n", 6_000, "-l", 100, "--vcf", source,
            "--insert-mean", 300, "--insert-sd", 30,
            "--insert-min", 150, "--insert-max", 500,
            "-e", 0, "--format", "bam", "--save-vcf",
        )
        self.assertEqual(run.variants, wgs_run().variants)
        problems, _operations = reconstruction_problems(run)
        self.assertEqual(problems[:10], [])


if __name__ == "__main__":
    unittest.main()
