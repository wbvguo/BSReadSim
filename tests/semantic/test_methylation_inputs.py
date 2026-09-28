"""Methylation profiles: text inputs, fallback, pooling, CpG-only, MethDB reuse."""

from collections import Counter
import math
import unittest

from tests.semantic import support

CG_VALUES = (0.0, 0.2, 0.5, 0.8, 1.0)
CHH_VALUES = (0.1, 0.3)


def u16(probability: float) -> int:
    """UNORM16 storage: round half away from zero, as documented."""
    return math.floor(probability * 65535 + 0.5)


def profile_sites():
    """Listed '+' strand cytosines on chrA: (position, context, probability)."""
    sequence = support.GENOME["chrA"]
    cg, chh = [], []
    for position in range(1_000, 230_000):
        if sequence[position] != "C":
            continue
        if sequence[position + 1] == "G":
            cg.append(position)
        elif sequence[position + 2] != "G":
            chh.append(position)
    sites = [(p, "CG", CG_VALUES[k % 5]) for k, p in enumerate(cg[::5][:2_500])]
    sites += [(p, "CHH", CHH_VALUES[k % 2]) for k, p in enumerate(chh[::40][:1_000])]
    return sorted(sites)


SITES = profile_sites()
LISTED = {position: probability for position, _context, probability in SITES}


def write_profiles():
    sequence = support.GENOME["chrA"]
    cgmap = ["#CHR\tNUC\tPOS\tCONTEXT\tDINUC\tMETH\tMC\tNC"]
    bedmethyl, methbg, methbed = [], [], []
    for position, context, p in SITES:
        cgmap.append("chrA\tC\t{}\t{}\t{}\t{:g}\t{}\t10".format(
            position + 1, context, sequence[position:position + 2], p, round(p * 10)))
        bedmethyl.append("chrA\t{0}\t{1}\tm\t{2}\t+\t{0}\t{1}\t0\t10\t{3:g}".format(
            position, position + 1, round(p * 1000), p * 100))
        methbg.append("chrA\t{}\t{}\t{:g}".format(position, position + 1, p))
        methbed.append("chrA\t{}\t{}\t.\t{}\t+".format(
            position, position + 1, round(p * 1000)))
    return {
        "cgmap": support.write_text("profile.CGmap", "\n".join(cgmap) + "\n"),
        "bedmethyl": support.write_text("profile.bedmethyl", "\n".join(bedmethyl) + "\n"),
        "methbg": support.write_text("profile.methbg", "\n".join(methbg) + "\n"),
        "methbed": support.write_text("profile.methbed", "\n".join(methbed) + "\n"),
    }


PROFILES = write_profiles()
_BUILDS = {}


def build_methdb(name, *options):
    if name not in _BUILDS:
        output = support.WORK / (name + ".methdb")
        support.bsreadsim(
            "build", "methdb", "-r", support.REFERENCE, "-o", output,
            "--mutation-rate", 0, "--seed-meth", 7, *options,
        )
        _BUILDS[name] = support.export_methdb(output)
    return _BUILDS[name]


def site_key(site):
    return (site.contig, site.position, site.strand, site.set, site.origin,
            site.context, site.probability_u16)


class TextProfileTests(unittest.TestCase):
    def test_all_text_formats_describe_the_same_profile(self):
        exports = {
            name: [site_key(s) for s in build_methdb(name, "--" + name, path)]
            for name, path in PROFILES.items()
        }
        reference = exports.pop("cgmap")
        for name, sites in exports.items():
            self.assertEqual(sites, reference, name)

    def test_listed_sites_keep_values_and_unlisted_sites_fall_back(self):
        listed = 0
        for site in build_methdb("cgmap", "--cgmap", PROFILES["cgmap"]):
            if site.contig == "chrA" and site.strand == "+" and site.position in LISTED:
                listed += 1
                self.assertEqual(site.probability_u16, u16(LISTED[site.position]))
                self.assertNotEqual(site.source, "beta")
            else:
                self.assertEqual(site.source, "beta")
        self.assertEqual(listed, len(LISTED))

    def test_reads_realize_listed_values(self):
        run = support.simulate(
            "wgbs-cgmap", "wgbs", "-n", 40_000, "-l", 100, "--threads", 4,
            "--cgmap", PROFILES["cgmap"], "--mutation-rate", 0,
            "--conversion-rate", 1, "-e", 0, "--format", "bam",
        )
        self.assertEqual(
            {r.tags["YS"] for r in run.bam.records}, {"OT", "OB"}
        )
        audit = support.audit_bisulfite_reads(run, None)
        self.assertEqual(audit.base_mismatches, 0, audit.examples)
        totals = {p: [0, 0] for p in CG_VALUES + CHH_VALUES}
        for (contig, position, strand, _h), (methylated, observed) in audit.site_states.items():
            if contig == "chrA" and strand == "+" and position in LISTED:
                entry = totals[LISTED[position]]
                entry[0] += methylated
                entry[1] += observed
        self.assertEqual(totals[0.0][0], 0)
        self.assertEqual(totals[1.0][0], totals[1.0][1])
        for p in (0.2, 0.5, 0.8, 0.1, 0.3):
            methylated, observed = totals[p]
            self.assertGreater(observed, 800)
            self.assertAlmostEqual(methylated / observed, p, delta=0.06, msg=str(p))


class PoolingTests(unittest.TestCase):
    def test_pooling_resamples_values_within_contig_and_context(self):
        sites = build_methdb("pooled", "--cgmap", PROFILES["cgmap"], "--pool-meth")
        cg_pool = Counter()
        for site in sites:
            context = site.context.split("-")[0]
            if site.contig != "chrA" or context == "CHG":
                self.assertEqual(site.source, "beta", site)
            elif context == "CG":
                self.assertIn(site.probability_u16, {u16(p) for p in CG_VALUES})
                cg_pool[site.probability_u16] += 1
            else:
                self.assertIn(site.probability_u16, {u16(p) for p in CHH_VALUES})
        total = sum(cg_pool.values())
        self.assertGreater(total, 10_000)
        for p in CG_VALUES:
            self.assertAlmostEqual(cg_pool[u16(p)] / total, 0.2, delta=0.02)


class CpgOnlyTests(unittest.TestCase):
    def test_cpg_only_profile_keeps_exactly_the_cpg_sites(self):
        full = build_methdb("full")
        cpg = build_methdb("cpg-only", "--cpg-only")
        self.assertTrue(all(s.context in ("CG-C", "CG-G") for s in cpg))
        self.assertEqual(
            [(s.contig, s.position, s.strand) for s in cpg],
            [(s.contig, s.position, s.strand) for s in full
             if s.context in ("CG-C", "CG-G")],
        )


def methdb_runs():
    first = support.simulate(
        "wgbs-save-methdb", "wgbs", "-n", 2_000, "-l", 100,
        "--mutation-rate", 0.003, "--indel-fraction", 0.2,
        "-e", 0, "--format", "bam", "--save-truth",
    )
    second = support.simulate(
        "wgbs-load-methdb", "wgbs", "-n", 2_000, "-l", 100,
        "--methdb", first.directory / "truth" / "sim.methdb",
        "-e", 0, "--format", "bam", "--save-truth", seed="77",
    )
    return first, second


class MethdbReuseTests(unittest.TestCase):
    def test_saved_methdb_reproduces_profile(self):
        first, second = methdb_runs()
        self.assertEqual(
            [site_key(s) for s in second.methdb_sites],
            [site_key(s) for s in first.methdb_sites],
        )

    def test_reads_use_the_embedded_variants(self):
        first, second = methdb_runs()
        self.assertGreater(len(first.variants), 500)
        problems, _operations = support.reconstruction_problems(
            second, first.haplotypes)
        self.assertEqual(problems[:10], [])

    def test_saved_vcf_reports_the_embedded_variants(self):
        first, second = methdb_runs()
        self.assertEqual(second.variants, first.variants)


if __name__ == "__main__":
    unittest.main()
