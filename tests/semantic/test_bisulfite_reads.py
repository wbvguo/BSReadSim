"""WGBS reads realize the diploid genome, methylation profile, and library."""

from collections import Counter, defaultdict
import statistics
import unittest

from tests.semantic import support
from tests.semantic.support import CG, CHG, CHH

READS = 40_000
CONVERSION_RATE = 0.95
ERROR_RATE = 0.01
MUTATION_RATE = 0.004
INSERT_MEAN, INSERT_SD, INSERT_MIN, INSERT_MAX = 300, 30, 150, 500

# Bismark strand conventions: (XG, R1 XR, R2 XR, R1 reverse-mapped)
STRANDS = {
    "OT": ("CT", "CT", "GA", False),
    "OB": ("GA", "CT", "GA", True),
    "CTOT": ("CT", "GA", "CT", True),
    "CTOB": ("GA", "GA", "CT", False),
}


def wgbs_run():
    return support.simulate(
        "wgbs", "wgbs",
        "-n", READS, "-l", 100, "--threads", 4,
        "--insert-mean", INSERT_MEAN, "--insert-sd", INSERT_SD,
        "--insert-min", INSERT_MIN, "--insert-max", INSERT_MAX,
        "--mutation-rate", MUTATION_RATE, "--indel-fraction", 0,
        "--conversion-rate", CONVERSION_RATE, "--undirectional",
        "-e", ERROR_RATE, "-q", 30,
        "--format", "bam", "--fragment-realization", "--save-truth",
    )


_AUDIT = None


def wgbs_audit():
    global _AUDIT
    if _AUDIT is None:
        run = wgbs_run()
        _AUDIT = support.audit_bisulfite_reads(run, run.haplotypes)
    return _AUDIT


class ReadLayoutTests(unittest.TestCase):
    def test_read_count_and_mate_pairing(self):
        run = wgbs_run()
        self.assertEqual(len(run.bam.records), READS)
        self.assertEqual(len(run.fragments), READS // 2)
        for name, records in run.fragments.items():
            self.assertEqual(sorted(r.mate_number for r in records), [1, 2], name)
            self.assertEqual(len({support.haplotype_of(r) for r in records}), 1)
            self.assertEqual(len({r.tags["YS"] for r in records}), 1)
            self.assertEqual(len({r.tags["zx"] for r in records}), 1)

    def test_strand_tags_follow_bismark_conventions(self):
        for record in wgbs_run().bam.records:
            xg, r1_xr, r2_xr, r1_reverse = STRANDS[record.tags["YS"]]
            self.assertEqual(record.tags["XG"], xg)
            if record.mate_number == 1:
                self.assertEqual(record.tags["XR"], r1_xr)
                self.assertEqual(record.is_reverse, r1_reverse)
            else:
                self.assertEqual(record.tags["XR"], r2_xr)
                self.assertEqual(record.is_reverse, not r1_reverse)

    def test_undirectional_library_samples_four_strands_equally(self):
        counts = Counter(
            records[0].tags["YS"] for records in wgbs_run().fragments.values()
        )
        total = sum(counts.values())
        for strand in STRANDS:
            self.assertAlmostEqual(counts[strand] / total, 0.25, delta=0.02)

    def test_mates_span_the_fragment_envelope(self):
        lengths = []
        for name, records in wgbs_run().fragments.items():
            _contig, start, end = support.envelope(name)
            self.assertEqual(min(r.position for r in records), start, name)
            self.assertEqual(max(r.reference_end for r in records), end, name)
            lengths.append(end - start)
        self.assertGreaterEqual(min(lengths), INSERT_MIN)
        self.assertLessEqual(max(lengths), INSERT_MAX)
        self.assertAlmostEqual(statistics.fmean(lengths), INSERT_MEAN, delta=1.5)
        self.assertAlmostEqual(statistics.pstdev(lengths), INSERT_SD, delta=1.5)

    def test_reads_respect_ambiguous_base_limit(self):
        for record in wgbs_run().bam.records:
            self.assertLessEqual(record.sequence.count("N"), 5, record.query_name)

    def test_fragments_are_spread_uniformly(self):
        run = wgbs_run()
        per_contig = Counter()
        bins = Counter()
        for name in run.fragments:
            contig, start, _end = support.envelope(name)
            per_contig[contig] += 1
            if contig == "chrA":
                bins[start // 20_000] += 1
        total = sum(per_contig.values())
        expected_a = len(support.GENOME["chrA"]) / (
            len(support.GENOME["chrA"]) + len(support.GENOME["chrB"]) - 2_000
        )
        self.assertAlmostEqual(per_contig["chrA"] / total, expected_a, delta=0.02)
        expected = per_contig["chrA"] / 12
        for index in range(12):
            self.assertAlmostEqual(bins[index], expected, delta=0.15 * expected)

    def test_haplotypes_are_sampled_equally(self):
        counts = Counter(
            support.haplotype_of(records[0])
            for records in wgbs_run().fragments.values()
        )
        self.assertAlmostEqual(counts[0] / (counts[0] + counts[1]), 0.5, delta=0.02)


class ReadContentTests(unittest.TestCase):
    def test_every_base_matches_converted_haplotype(self):
        audit = wgbs_audit()
        self.assertGreater(audit.bases, READS * 99)
        self.assertEqual(audit.base_mismatches, 0, audit.examples)

    def test_zt_context_matches_haplotype_sequence(self):
        audit = wgbs_audit()
        self.assertGreater(audit.context_checked, 1_000_000)
        self.assertEqual(audit.context_mismatches, 0, audit.examples)
        self.assertEqual(audit.flags_on_non_cytosine, 0)

    def test_opposite_strand_cytosines_are_not_converted(self):
        self.assertEqual(wgbs_audit().opposite_strand_converted, 0)

    def test_variant_flag_marks_haplotype_snvs(self):
        self.assertEqual(wgbs_audit().variant_flag_mismatches, 0)

    def test_methylated_cytosines_are_never_converted(self):
        self.assertEqual(wgbs_audit().methylated_and_converted, 0)

    def test_unmethylated_cytosines_convert_at_configured_rate(self):
        audit = wgbs_audit()
        rate = audit.converted_unmethylated_sites / audit.unmethylated_sites
        self.assertAlmostEqual(rate, CONVERSION_RATE, delta=0.005)

    def test_mates_share_one_molecule_state(self):
        self.assertEqual(wgbs_audit().conflicting_mate_states, 0)

    def test_sequencing_errors_follow_uniform_rate(self):
        audit = wgbs_audit()
        self.assertAlmostEqual(
            audit.sequencing_errors / audit.bases, ERROR_RATE, delta=0.001
        )
        self.assertEqual(audit.errors_equal_to_truth, 0)
        for record in wgbs_run().bam.records:
            self.assertEqual(set(record.quality), {30})


class MethylationProfileTests(unittest.TestCase):
    def test_generated_profile_follows_context_beta_means(self):
        by_context = defaultdict(list)
        for site in wgbs_run().methdb_sites:
            by_context[site.context.split("-")[0]].append(site.probability)
        self.assertEqual(set(by_context), {"CG", "CHG", "CHH"})
        self.assertAlmostEqual(statistics.fmean(by_context["CG"]), 0.5, delta=0.02)
        self.assertAlmostEqual(statistics.pstdev(by_context["CG"]), 0.3536, delta=0.02)
        for context in ("CHG", "CHH"):
            self.assertAlmostEqual(
                statistics.fmean(by_context[context]), 1 / 6, delta=0.02
            )
        self.assertTrue(all(site.source == "beta" for site in wgbs_run().methdb_sites))

    def test_profile_covers_every_reference_cytosine(self):
        run = wgbs_run()
        shared = {
            (site.contig, site.position, site.strand): site.context
            for site in run.methdb_sites
            if site.set == "shared"
        }
        variant_positions = {
            (variant.contig, variant.position) for variant in run.variants
        }
        names = {CG: "CG", CHG: "CHG", CHH: "CHH"}
        missing = 0
        checked = 0
        for contig, sequence in support.GENOME.items():
            for position in range(2, len(sequence) - 2):
                if any(
                    (contig, near) in variant_positions
                    for near in range(position - 2, position + 3)
                ):
                    continue
                for strand, base in (("+", "C"), ("-", "G")):
                    if sequence[position] != base:
                        continue
                    context = support.context_at(
                        sequence.__getitem__, position, len(sequence), strand == "-"
                    )
                    if context is None:
                        continue
                    checked += 1
                    expected = "{}-{}".format(names[context], "C" if strand == "+" else "G")
                    missing += shared.get((contig, position, strand)) != expected
        self.assertGreater(checked, 150_000)
        self.assertEqual(missing, 0)

    def test_read_states_realize_site_probabilities(self):
        run = wgbs_run()
        probability = {
            (site.contig, site.position, site.strand): site.probability
            for site in run.methdb_sites
            if site.set == "shared" and site.origin == "reference"
        }
        bins = defaultdict(lambda: [0.0, 0, 0])  # sum p, methylated, observed
        for (contig, position, strand, _hap), (methylated, observed) in (
            wgbs_audit().site_states.items()
        ):
            p = probability.get((contig, position, strand))
            if p is None:
                continue
            entry = bins[min(int(p * 10), 9)]
            entry[0] += p * observed
            entry[1] += methylated
            entry[2] += observed
        checked = 0
        for index, (sum_p, methylated, observed) in sorted(bins.items()):
            if observed < 2_000:
                continue
            checked += 1
            self.assertAlmostEqual(
                methylated / observed, sum_p / observed, delta=0.03, msg=str(index)
            )
        self.assertGreaterEqual(checked, 8)


class GeneratedVariantTests(unittest.TestCase):
    def test_de_novo_snvs_follow_rate_and_genotype_model(self):
        variants = wgbs_run().variants
        self.assertTrue(all(variant.kind == "snv" for variant in variants))
        callable_bases = sum(
            len(sequence) - sequence.count("N") for sequence in support.GENOME.values()
        )
        expected = MUTATION_RATE * callable_bases
        self.assertAlmostEqual(len(variants), expected, delta=0.15 * expected)
        homozygous = sum(len(variant.haplotypes) == 2 for variant in variants)
        self.assertAlmostEqual(homozygous / len(variants), 1 / 3, delta=0.05)
        for variant in variants:
            self.assertNotEqual(variant.reference, variant.alternate)
            self.assertEqual(
                support.GENOME[variant.contig][variant.position], variant.reference
            )


class SingleEndTests(unittest.TestCase):
    def test_single_end_emits_one_read_per_fragment(self):
        run = support.simulate(
            "wgbs-single", "wgbs", "-n", 2_000, "-l", 100, "--single-end",
            "--mutation-rate", 0, "--format", "bam",
        )
        self.assertEqual(len(run.bam.records), 2_000)
        self.assertEqual(len(run.fragments), 2_000)
        strands = Counter(record.tags["YS"] for record in run.bam.records)
        self.assertEqual(set(strands), {"OT", "OB"})
        for record in run.bam.records:
            _contig, start, end = support.envelope(record.query_name)
            if record.is_reverse:
                self.assertEqual(record.reference_end, end)
            else:
                self.assertEqual(record.position, start)
        audit = support.audit_bisulfite_reads(run, None)
        self.assertEqual(audit.base_mismatches, 0, audit.examples)


if __name__ == "__main__":
    unittest.main()
