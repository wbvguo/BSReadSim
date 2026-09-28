"""Reference letters other than A, C, G, T, and N (such as IUPAC codes) are simulated as N."""

import hashlib
import json
import random
import unittest

from tests.semantic import support


def _sequence():
    rng = random.Random(7)
    bases = [rng.choice("ACGT") for _ in range(20_000)]
    other = {}
    for k, letter in enumerate("MRWSYKVHDBmrwsykvhdbXx"):
        position = 1_000 + 800 * k
        bases[position] = letter
        other[position] = letter
    return "".join(bases), other


SEQUENCE, OTHER = _sequence()
FASTA = support.write_text("other-letters.fa", ">chrI\n" + "".join(
    SEQUENCE[i:i + 60] + "\n" for i in range(0, len(SEQUENCE), 60)))


class OtherLetterTests(unittest.TestCase):
    def test_other_letters_are_simulated_as_n(self):
        directory = support.WORK / "other-letters"
        support.bsreadsim(
            "run", "wgs", "-r", FASTA, "-o", directory, "--seed", support.SEED,
            "-n", 4_000, "-l", 100, "--mutation-rate", 0, "-e", 0, "--format", "bam",
        )
        run = support.Run(directory)
        covered = 0
        for record in run.bam.records:
            for query, position, operation in record.aligned_pairs():
                if operation != "M":
                    continue
                expected = "N" if position in OTHER else SEQUENCE[position].upper()
                self.assertEqual(record.sequence[query], expected, (record.query_name, position))
                covered += position in OTHER
        self.assertGreater(covered, 100)

        # The contig digest is that of the sequence as written, like the SAM M5 tag.
        manifest = json.loads(run.file(".manifest.json").read_text())
        self.assertEqual(manifest["details"]["contigs"][0]["md5"],
                         hashlib.md5(SEQUENCE.upper().encode()).hexdigest())


if __name__ == "__main__":
    unittest.main()
