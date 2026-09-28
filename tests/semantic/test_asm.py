"""Allele-specific methylation follows the allele of the linked SNV."""

import unittest

from tests.semantic import support


def asm_records():
    """CpG targets on chrA, each linked to an A<->T SNV six bases downstream.

    A<->T substitutions cannot change any cytosine context, and bisulfite
    conversion leaves both alleles readable.
    """
    sequence = support.GENOME["chrA"]
    records = []
    last = -100
    for target in range(2_000, 230_000):
        if target - last < 40:
            continue
        snv = target + 6
        if sequence[target:target + 2] == "CG" and sequence[snv] in "AT":
            reference = sequence[snv]
            alternate = "T" if reference == "A" else "A"
            ref_meth, alt_meth = (0.0, 1.0) if len(records) % 2 == 0 else (1.0, 0.0)
            records.append((target, snv, reference, alternate, ref_meth, alt_meth))
            last = target
            if len(records) == 300:
                break
    return records


RECORDS = asm_records()


def write_asm_bed():
    lines = [
        "chrA\t{0}\t{1}\tasm{2}\t1000\t+\t{3}\t{4}\t{5}\t{6}\t{7:g}\t{8:g}".format(
            target, target + 1, k, snv, snv + 1, ref, alt, ref_meth, alt_meth)
        for k, (target, snv, ref, alt, ref_meth, alt_meth) in enumerate(RECORDS)
    ]
    return support.write_text("asm.bed", "\n".join(lines) + "\n")


def write_asm_ass():
    lines = [
        "Chr\tSNP_Pos\tRef\tAllele1\tAllele2\tC_Pos\tAllele1_linked_C\t"
        "Allele2_linked_C\tAllele1_linked_C_met\tAllele2_linked_C_met\t"
        "pvalue\tfdr\tASM"
    ]
    for k, (target, snv, ref, alt, ref_meth, alt_meth) in enumerate(RECORDS):
        alleles = ((ref, ref_meth), (alt, alt_meth))
        if k % 3 == 0:  # CGmapTools does not promise REF as Allele1
            alleles = alleles[::-1]
        (allele1, met1), (allele2, met2) = alleles
        lines.append("chrA\t{}\t{}\t{}\t{}\t{}\t8-2\t2-8\t{:g}\t{:g}\t0.001\t0.005\tTRUE".format(
            snv + 1, ref, allele1, allele2, target + 1, met1, met2))
    return support.write_text("asm.ass", "\n".join(lines) + "\n")


ASM_BED = write_asm_bed()
ASM_ASS = write_asm_ass()


def asm_run():
    return support.simulate(
        "wgbs-asm", "wgbs", "-n", 60_000, "-l", 100, "--threads", 4,
        "--asm-bed", ASM_BED, "--mutation-rate", 0,
        "--conversion-rate", 1, "-e", 0, "--format", "bam", "--save-truth",
    )


class AsmTests(unittest.TestCase):
    def test_linked_snvs_become_heterozygous_variants(self):
        variants = asm_run().variants
        self.assertEqual(
            [(v.position, v.reference, v.alternate) for v in variants],
            [(snv, ref, alt) for _t, snv, ref, alt, _r, _a in RECORDS],
        )
        self.assertTrue(all(len(v.haplotypes) == 1 for v in variants))
        self.assertEqual({v.contig for v in variants}, {"chrA"})

    def test_methylation_follows_the_linked_allele(self):
        run = asm_run()
        alternate_haplotype = {v.position: next(iter(v.haplotypes)) for v in run.variants}
        audit = support.audit_bisulfite_reads(run, run.haplotypes)
        self.assertEqual(audit.base_mismatches, 0, audit.examples)
        observed = violations = 0
        for target, snv, _ref, _alt, ref_meth, alt_meth in RECORDS:
            for haplotype in (0, 1):
                methylated, total = audit.site_states.get(
                    ("chrA", target, "+", haplotype), (0, 0))
                p = alt_meth if haplotype == alternate_haplotype[snv] else ref_meth
                observed += total
                violations += methylated if p == 0 else total - methylated
        self.assertGreater(observed, 1_000)
        self.assertEqual(violations, 0)

    def test_asm_formats_are_equivalent(self):
        exports = []
        for option, path in (("--asm-bed", ASM_BED), ("--asm", ASM_ASS)):
            output = support.WORK / "asm{}.methdb".format(option)
            support.bsreadsim(
                "build", "methdb", "-r", support.REFERENCE, "-o", output,
                option, path, "--seed-meth", 7, "--seed-phase", 8,
            )
            exports.append([
                (s.contig, s.position, s.strand, s.set, s.context, s.probability_u16)
                for s in support.export_methdb(output)
            ])
        self.assertEqual(exports[0], exports[1])


class AsmCheckTests(unittest.TestCase):
    def test_runs_check_asm_on_contigs_without_fragments(self):
        """Every contig is prepared as `build methdb` would prepare it, even one that
        gets no fragment: an ASM row on chrB, outside the only target, still fails."""
        sequence = support.GENOME["chrB"]
        target = next(t for t in range(1_000, len(sequence) - 10)
                      if sequence[t:t + 2] == "CG" and sequence[t + 6] in "AT")
        snv, ref = target + 6, sequence[target + 6]
        alt = "T" if ref == "A" else "A"
        asm = support.write_text("asm-chrB.bed", "chrB\t{}\t{}\tbad\t1000\t+\t{}\t{}\t{}\t{}\t0\t1\n".format(
            target, target + 1, snv, snv + 1, ref, alt))
        vcf = support.write_text("asm-chrB.vcf", (
            "##fileformat=VCFv4.2\n"
            '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n"
            "chrB\t{}\t.\t{}\t{}\t.\tPASS\t.\tGT\t1|1\n").format(snv + 1, ref, alt))
        targets = support.write_text("asm-chrA-target.bed", "chrA\t50000\t50200\tt\t1\t.\n")
        message = "must retain one reference haplotype"
        run = support.bsreadsim(
            "run", "tbs", "-r", support.REFERENCE, "-o", support.WORK / "asm-chrB-run",
            "--targets", targets, "--vcf", vcf, "--asm-bed", asm, "-n", 200, "-l", 100,
            "--format", "bam", check=False)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn(message, run.stderr)
        build = support.bsreadsim(
            "build", "methdb", "-r", support.REFERENCE, "-o", support.WORK / "asm-chrB.methdb",
            "--vcf", vcf, "--asm-bed", asm, check=False)
        self.assertIn(message, build.stderr)


if __name__ == "__main__":
    unittest.main()
