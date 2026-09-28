"""RRBS fragments are restriction fragments selected by size and score."""

from collections import Counter
import unittest

from tests.semantic import support

READ_LENGTH = 36
WINDOW = ("--insert-min", 40, "--insert-mean", 150, "--insert-max", 300)
MOTIF = "CCGG"  # default cut site C|CGG


def cut_positions(sequence):
    positions = []
    start = sequence.find(MOTIF)
    while start >= 0:
        positions.append(start + 1)
        start = sequence.find(MOTIF, start + 1)
    return positions


CUTS = {name: cut_positions(sequence) for name, sequence in support.GENOME.items()}
_CANDIDATES = {}


def candidates():
    if "rows" not in _CANDIDATES:
        path = support.WORK / "rrbs-candidates.bed"
        support.bsreadsim(
            "build", "rrbs", "-r", support.REFERENCE, "-o", path,
            "-l", READ_LENGTH, *WINDOW,
        )
        rows = []
        for line in path.read_text().splitlines():
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            rows.append(fields)
        _CANDIDATES["rows"] = rows
    return _CANDIDATES["rows"]


def eligible(sequence, start, end):
    """Each mate may contain at most 5% ambiguous bases."""
    limit = int(0.05 * READ_LENGTH)
    return (
        sequence[start:start + READ_LENGTH].count("N") <= limit
        and sequence[end - READ_LENGTH:end].count("N") <= limit
    )


def rrbs_run(name, *options):
    return support.simulate(
        name, "rrbs", "-n", 20_000, "-l", READ_LENGTH, *WINDOW, "--threads", 4,
        "--conversion-rate", 1, "-e", 0, "--format", "bam", *options,
    )


class CandidateTests(unittest.TestCase):
    def test_candidates_are_size_selected_restriction_fragments(self):
        rows = candidates()
        self.assertGreater(len(rows), 1_000)
        for fields in rows:
            contig, start, end = fields[0], int(fields[1]), int(fields[2])
            cuts = CUTS[contig]
            sequence = support.GENOME[contig]
            self.assertIn(start, cuts)
            self.assertIn(end, cuts)
            self.assertTrue(40 <= end - start <= 300)
            self.assertEqual(fields[4], "1")
            self.assertEqual(fields[6], "3")
            self.assertEqual(int(fields[7]), end - start)
            self.assertEqual(
                int(fields[8]), sum(b in "CG" for b in sequence[start:end]))
            self.assertEqual(
                int(fields[9]), sum(start <= c <= end for c in cuts))

    def test_candidate_domain_is_complete(self):
        expected = set()
        for contig, cuts in CUTS.items():
            sequence = support.GENOME[contig]
            for i, start in enumerate(cuts):
                for end in cuts[i + 1:]:
                    if end - start > 300:
                        break
                    if end - start >= 40 and eligible(sequence, start, end):
                        expected.add((contig, start, end))
        observed = {(f[0], int(f[1]), int(f[2])) for f in candidates()}
        self.assertEqual(observed, expected)


class SamplingTests(unittest.TestCase):
    def test_uniform_sampling_draws_candidates_evenly(self):
        run = rrbs_run("rrbs-uniform", "--mutation-rate", 0)
        domain = {(f[0], int(f[1]), int(f[2])) for f in candidates()}
        counts = Counter(support.envelope(name) for name in run.fragments)
        self.assertTrue(set(counts) <= domain)
        mean = sum(counts.values()) / len(domain)
        chi_square = sum((counts[key] - mean) ** 2 / mean for key in domain)
        self.assertAlmostEqual(chi_square / (len(domain) - 1), 1.0, delta=0.15)
        audit = support.audit_bisulfite_reads(run, None)
        self.assertEqual(audit.base_mismatches, 0, audit.examples)
        self.assertEqual({r.tags["YS"] for r in run.bam.records}, {"OT", "OB"})

    def test_score_sampling_follows_candidate_scores(self):
        weights = {}
        lines = []
        for index, fields in enumerate(candidates()):
            score = (0, 1, 3)[index % 3]
            weights[(fields[0], int(fields[1]), int(fields[2]))] = score
            lines.append("\t".join(fields[:4] + [str(score)] + fields[5:]))
        path = support.write_text("rrbs-scored.bed", "\n".join(lines) + "\n")
        run = rrbs_run(
            "rrbs-score", "--sampling", "score", "--rrbs-candidates", path,
        )
        counts = Counter(support.envelope(name) for name in run.fragments)
        by_score = {0: 0, 1: 0, 3: 0}
        for key, count in counts.items():
            by_score[weights[key]] += count
        sizes = Counter(weights.values())
        self.assertEqual(by_score[0], 0)
        ratio = (by_score[3] / sizes[3]) / (by_score[1] / sizes[1])
        self.assertAlmostEqual(ratio, 3.0, delta=0.3)


# Variants that create, remove, and move cut sites; candidates then belong to one haplotype.
VARIANTS = ("--mutation-rate", 0.02, "--indel-fraction", 0.5,
            "--indel-extension-probability", 0.3, "--seed-mut", 11, "--seed-phase", 12)
_VARIANT_RUN = {}


def variant_vcf():
    if "vcf" not in _VARIANT_RUN:
        path = support.WORK / "rrbs-variants.vcf.gz"
        support.bsreadsim("build", "variants", "-r", support.REFERENCE, "-o", path, *VARIANTS)
        _VARIANT_RUN["vcf"] = path
    return _VARIANT_RUN["vcf"]


def haplotype_candidates(contig, sequence, origins, mask):
    """Candidate BED rows (without IDs and scores) of one haplotype."""
    cuts = cut_positions(sequence)
    rows = []
    for i, start in enumerate(cuts):
        for end in cuts[i + 1:]:
            if end - start > 300:
                break
            if end - start < 40 or not eligible(sequence, start, end):
                continue
            span = origins[start:end]
            envelope_start = next(position for position, inserted in span if not inserted)
            rows.append((
                contig, envelope_start, span[-1][0] + 1, mask, end - start,
                sum(b in "CG" for b in sequence[start:end]),
                sum(start <= c <= end for c in cuts),
            ))
    return rows


class VariantTests(unittest.TestCase):
    def test_candidates_are_the_restriction_fragments_of_each_haplotype(self):
        variants = support.read_vcf(variant_vcf())
        self.assertGreater(sum(v.kind != "snv" for v in variants), 1_000)
        path = support.WORK / "rrbs-variant-candidates.bed"
        support.bsreadsim("build", "rrbs", "-r", support.REFERENCE, "-o", path,
                          "-l", READ_LENGTH, *WINDOW, *VARIANTS)
        observed = Counter()
        for line in path.read_text().splitlines():
            if not line.startswith("#"):
                f = line.split("\t")
                observed[(f[0], int(f[1]), int(f[2]), *map(int, f[6:10]))] += 1
        expected = Counter()
        for contig in support.GENOME:
            for haplotype in (0, 1):
                sequence, origins = support.haplotype_sequence(contig, variants, haplotype)
                expected.update(haplotype_candidates(contig, sequence, origins, haplotype + 1))
        changed = sum(n for row, n in expected.items() if row[4] != row[2] - row[1])
        self.assertGreater(changed, 1_000)
        self.assertEqual(observed, expected)

    def test_reads_come_from_haplotype_restriction_fragments(self):
        vcf = variant_vcf()
        run = support.simulate(
            "rrbs-variants", "rrbs", "-n", 20_000, "-l", READ_LENGTH, *WINDOW, "--threads", 4,
            "--vcf", vcf, "--conversion-rate", 0, "-e", 0, "--format", "bam",
        )
        problems, operations = support.reconstruction_problems(
            run, support.Haplotypes(support.read_vcf(vcf)))
        self.assertEqual(problems[:10], [])
        self.assertGreater(operations["I"], 100)
        self.assertGreater(operations["D"], 100)
        for name, records in run.fragments.items():
            records = sorted(records, key=lambda r: r.position)
            self.assertTrue(records[0].sequence.startswith("CGG"), name)
            self.assertTrue(records[-1].sequence.endswith("C"), name)


if __name__ == "__main__":
    unittest.main()
