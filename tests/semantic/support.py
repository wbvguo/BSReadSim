"""Shared fixtures for the semantic test suite.

The semantic suite states what BSReadSim simulates, not how. It drives only
public interfaces -- the ``bsreadsim`` CLI, BAM records with their documented
truth tags, and exported truth artifacts (VCF, MethDB export, RRBS candidate
BED) -- so it stays valid when htsim or the Python processing stack is
rewritten. Every simulation uses fixed seeds, and statistical tolerances are
at least four standard deviations wide, so any correct implementation passes
regardless of its random-number streams.

Set ``BSREADSIM_HTSIM`` to the htsim executable under test. Set
``BSREADSIM_KEEP_SEMANTIC=1`` to keep the generated files for inspection and
``BSREADSIM_SEMANTIC_SEED`` to rerun every simulation with another seed.
"""

from __future__ import annotations

import atexit
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from functools import cached_property
import gzip
import itertools
import os
from pathlib import Path
import random
import shutil
import string
import subprocess
import sys
import tempfile

from tests.helpers.bam_reader import BamFile, read_bam

HTSIM = Path(os.environ["BSREADSIM_HTSIM"]).resolve()
# The CLI option that points bsreadsim at a specific htsim executable.
HTSIM_OPTION = "--core"
# Any seed must pass; override to check that tolerances do not rely on luck.
SEED = os.environ.get("BSREADSIM_SEMANTIC_SEED", "20261003")

WORK = Path(tempfile.mkdtemp(prefix="bsreadsim-semantic-"))
if os.environ.get("BSREADSIM_KEEP_SEMANTIC"):
    print("semantic outputs kept in", WORK, file=sys.stderr)
else:
    atexit.register(shutil.rmtree, WORK, ignore_errors=True)

COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
ZT_VALUES = {
    character: value
    for value, character in enumerate(
        string.ascii_uppercase + string.ascii_lowercase + string.digits + "-_"
    )
}

# zt context codes
NO_CONTEXT, CG, CHG, CHH = 0, 1, 2, 3
CONTEXT_CODES = {"CG": CG, "CHG": CHG, "CHH": CHH}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


# ---------------------------------------------------------------- reference

def _random_block(rng: random.Random, length: int, gc: float) -> str:
    at = (1.0 - gc) / 2.0
    return "".join(rng.choices("ACGT", weights=(at, gc / 2, gc / 2, at), k=length))


def _make_genome() -> dict[str, str]:
    """Two contigs; chrA alternates 8 kb blocks of 35/50/65% GC."""
    rng = random.Random(20261003)
    gc_cycle = itertools.cycle((0.35, 0.50, 0.65))
    chr_a = "".join(_random_block(rng, 8_000, next(gc_cycle)) for _ in range(30))
    chr_b = (
        _random_block(rng, 60_000, 0.45)
        + "N" * 2_000
        + _random_block(rng, 58_000, 0.45)
    )
    return {"chrA": chr_a, "chrB": chr_b}


GENOME = _make_genome()
REFERENCE = WORK / "reference.fa"
with REFERENCE.open("w") as handle:
    for name, sequence in GENOME.items():
        handle.write(">{}\n".format(name))
        for offset in range(0, len(sequence), 60):
            handle.write(sequence[offset:offset + 60] + "\n")


def context_at(base_at, position: int, length: int, reverse: bool):
    """Classify the cytosine at ``position`` on one strand.

    ``reverse`` selects the bottom strand, where the cytosine is a reference G.
    Returns None when a required flanking base is missing or ambiguous.
    """
    step = -1 if reverse else 1
    guanine = "C" if reverse else "G"
    first = position + step
    if not 0 <= first < length:
        return None
    first_base = base_at(first)
    if first_base == "N":
        return None
    if first_base == guanine:
        return CG
    second = position + 2 * step
    if not 0 <= second < length:
        return None
    second_base = base_at(second)
    if second_base == "N":
        return None
    return CHG if second_base == guanine else CHH


# ---------------------------------------------------------------- CLI runs

def bsreadsim(*arguments, check: bool = True) -> subprocess.CompletedProcess:
    command = [
        sys.executable, "-m", "bsreadsim",
        *(str(argument) for argument in arguments),
        HTSIM_OPTION, str(HTSIM),
    ]
    completed = subprocess.run(
        command, capture_output=True, text=True, cwd=str(WORK), check=False
    )
    if check and completed.returncode != 0:
        raise AssertionError(
            "bsreadsim {} failed ({}):\n{}".format(
                " ".join(str(a) for a in arguments[:2]),
                completed.returncode,
                completed.stderr,
            )
        )
    return completed


@dataclass(frozen=True)
class Variant:
    contig: str
    kind: str  # "snv", "insertion", or "deletion"
    position: int  # SNV/deletion: first affected base; insertion: base after the gap
    reference: str
    alternate: str
    haplotypes: frozenset

    @property
    def length(self) -> int:
        return max(len(self.reference), len(self.alternate))


def read_vcf(path: Path) -> list[Variant]:
    data = path.read_bytes()
    text = (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode()
    variants = []
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        position = int(fields[1]) - 1
        reference, alternate = fields[3], fields[4]
        prefix = 0
        while (
            prefix < min(len(reference), len(alternate))
            and reference[prefix] == alternate[prefix]
        ):
            prefix += 1
        reference, alternate = reference[prefix:], alternate[prefix:]
        while reference and alternate and reference[-1] == alternate[-1]:
            reference, alternate = reference[:-1], alternate[:-1]
        position += prefix
        if len(reference) == 1 and len(alternate) == 1:
            kind = "snv"
        elif not reference:
            kind = "insertion"
        elif not alternate:
            kind = "deletion"
        else:
            raise AssertionError("unexpected complex variant: " + line)
        genotype = fields[9].split(":")[0].replace("/", "|").split("|")
        haplotypes = frozenset(i for i, allele in enumerate(genotype) if allele == "1")
        variants.append(
            Variant(fields[0], kind, position, reference, alternate, haplotypes)
        )
    return variants


class Haplotypes:
    """Diploid genome = reference + phased variants (haplotype 0 and 1)."""

    def __init__(self, variants: list[Variant]):
        self.snvs = defaultdict(dict)  # (contig, h) -> {position: base}
        self.deleted = defaultdict(set)  # (contig, h) -> {position}
        self.insertions = defaultdict(dict)  # (contig, h) -> {position: bases}
        for variant in variants:
            for haplotype in variant.haplotypes:
                key = (variant.contig, haplotype)
                if variant.kind == "snv":
                    self.snvs[key][variant.position] = variant.alternate
                elif variant.kind == "deletion":
                    self.deleted[key].update(
                        range(variant.position, variant.position + len(variant.reference))
                    )
                else:
                    self.insertions[key][variant.position] = variant.alternate

    def base_getter(self, contig: str, haplotype: int):
        sequence = GENOME[contig]
        snvs = self.snvs.get((contig, haplotype), {})

        def base_at(position: int) -> str:
            return snvs.get(position, sequence[position])

        return base_at


@dataclass(frozen=True)
class MethdbSite:
    contig: str
    position: int
    strand: str
    set: str  # shared, haplotype-1, haplotype-2
    origin: str  # reference or insertion
    context: str  # e.g. CG-C
    source: str
    probability_u16: int

    @property
    def probability(self) -> float:
        return self.probability_u16 / 65535.0


def export_methdb(path: Path) -> list[MethdbSite]:
    destination = path.with_suffix(".export.bed")
    if not destination.exists():
        bsreadsim(
            "export", "methdb", "-i", path, "-o", destination, "--no-compression"
        )
    sites = []
    for line in destination.read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        sites.append(
            MethdbSite(
                contig=fields[0],
                position=int(fields[1]),
                strand=fields[5],
                set=fields[6],
                origin=fields[8],
                context=fields[11],
                source=fields[12],
                probability_u16=int(fields[14]),
            )
        )
    return sites


@dataclass
class Run:
    directory: Path
    prefix: str = "sim"

    @cached_property
    def bam(self) -> BamFile:
        return read_bam(self.directory / (self.prefix + ".bam"))

    @cached_property
    def variants(self) -> list[Variant]:
        return read_vcf(self.directory / "truth" / (self.prefix + ".variants.vcf.gz"))

    @cached_property
    def haplotypes(self) -> Haplotypes:
        return Haplotypes(self.variants)

    @cached_property
    def methdb_sites(self) -> list[MethdbSite]:
        return export_methdb(self.directory / "truth" / (self.prefix + ".methdb"))

    def file(self, suffix: str) -> Path:
        return self.directory / (self.prefix + suffix)

    @cached_property
    def fragments(self) -> dict[str, list]:
        grouped = defaultdict(list)
        for record in self.bam.records:
            grouped[record.query_name].append(record)
        return grouped

    def contig_of(self, record) -> str:
        return self.bam.references[record.reference_id][0]


_RUNS: dict[str, Run] = {}


def simulate(name: str, assay: str, *options, seed: str = SEED) -> Run:
    """Run ``bsreadsim run ASSAY`` once per test process and cache the result."""
    if name not in _RUNS:
        directory = WORK / name
        bsreadsim(
            "run", assay, "-r", REFERENCE, "-o", directory, "--seed", seed, *options
        )
        _RUNS[name] = Run(directory)
    return _RUNS[name]


def envelope(query_name: str) -> tuple[str, int, int]:
    """Return (contig, start, end) of the zero-based, half-open fragment envelope."""
    locus, _ordinal = query_name.rsplit(":", 1)
    contig, interval = locus.rsplit(":", 1)
    start, end = interval.split("-")
    return contig, int(start) - 1, int(end)


def haplotype_of(record) -> int:
    return record.tags["zr"][0] & 3


def haplotype_sequence(contig, variants, haplotype):
    """Bases of one haplotype, with the reference position of each (an inserted
    base takes the position of the base it follows) and whether it is inserted."""
    sequence = GENOME[contig]
    snvs, deleted, inserted = {}, set(), {}
    for v in variants:
        if v.contig != contig or haplotype not in v.haplotypes:
            continue
        if v.kind == "snv":
            snvs[v.position] = v.alternate
        elif v.kind == "deletion":
            deleted.update(range(v.position, v.position + len(v.reference)))
        else:
            inserted[v.position] = v.alternate
    bases, origins = [], []
    for position, base in enumerate(sequence):
        for inserted_base in inserted.get(position, ""):
            bases.append(inserted_base)
            origins.append((position - 1, True))
        if position not in deleted:
            bases.append(snvs.get(position, base))
            origins.append((position, False))
    return "".join(bases), origins


# ---------------------------------------------------------------- base audit

def reconstruction_problems(run: Run, haplotypes: Haplotypes):
    """Compare every read of an error-free run with its haplotype.

    M bases must equal the haplotype base, D operations must cover deleted
    reference bases, and I operations must carry the inserted alleles. Every
    insertion strictly inside a read's aligned span must appear as an I.
    Bisulfite reads may additionally show their XG conversion (C->T or G->A).
    Returns (problems, Counter of CIGAR operations).
    """
    problems = []
    operations = Counter()
    for record in run.bam.records:
        contig = run.contig_of(record)
        key = (contig, haplotype_of(record))
        sequence = GENOME[contig]
        snvs = haplotypes.snvs.get(key, {})
        deleted = haplotypes.deleted.get(key, set())
        insertions = haplotypes.insertions.get(key, {})
        query = 0
        reference = record.position
        inserted_at = set()
        last = len(record.cigar) - 1
        for index, (length, operation) in enumerate(record.cigar):
            operations[operation] += 1
            if operation == "M":
                for offset in range(length):
                    position = reference + offset
                    expected = snvs.get(position, sequence[position])
                    observed = record.sequence[query + offset]
                    if position in deleted or not _same_base(observed, expected, record):
                        problems.append((record.query_name, "M", position))
                query += length
                reference += length
            elif operation == "I":
                observed = record.sequence[query:query + length]
                expected = insertions.get(reference, "")
                # A read may begin or end inside an insertion.
                if index == 0:
                    expected = expected[len(expected) - length:]
                elif index == last:
                    expected = expected[:length]
                if len(expected) != length or not all(
                    _same_base(o, e, record) for o, e in zip(observed, expected, strict=True)
                ):
                    problems.append((record.query_name, "I", reference))
                inserted_at.add(reference)
                query += length
            elif operation == "D":
                for offset in range(length):
                    if reference + offset not in deleted:
                        problems.append((record.query_name, "D", reference + offset))
                reference += length
            else:
                problems.append((record.query_name, operation, reference))
        for position in insertions:
            if record.position < position < record.reference_end and position not in inserted_at:
                problems.append((record.query_name, "missing I", position))
    return problems, operations


def _same_base(observed: str, expected: str, record) -> bool:
    """Equality, allowing bisulfite conversion in bisulfite reads."""
    if observed == expected:
        return True
    conversion = record.tags.get("XG")
    return (conversion == "CT" and (expected, observed) == ("C", "T")) or (
        conversion == "GA" and (expected, observed) == ("G", "A"))


@dataclass
class BisulfiteAudit:
    """Per-base comparison of bisulfite reads against the diploid truth."""

    bases: int = 0
    base_mismatches: int = 0
    context_mismatches: int = 0
    context_checked: int = 0
    variant_flag_mismatches: int = 0
    flags_on_non_cytosine: int = 0
    opposite_strand_converted: int = 0
    methylated_and_converted: int = 0
    unmethylated_sites: int = 0
    converted_unmethylated_sites: int = 0
    sequencing_errors: int = 0
    errors_equal_to_truth: int = 0
    examples: list = field(default_factory=list)
    # (contig, position, strand, haplotype) -> [methylated, observed]
    site_states: dict = field(default_factory=lambda: defaultdict(lambda: [0, 0]))
    # query_name -> {(position, strand): (methylated, converted)}
    fragment_states: dict = field(default_factory=lambda: defaultdict(dict))
    conflicting_mate_states: int = 0

    def note(self, message: str) -> None:
        if len(self.examples) < 10:
            self.examples.append(message)


def audit_bisulfite_reads(run: Run, haplotypes: Haplotypes | None) -> BisulfiteAudit:
    """Check every aligned base of an SNV-only bisulfite run.

    Each base must equal the read's haplotype base, after C->T (XG=CT) or
    G->A (XG=GA) conversion where zt reports conversion, unless zt reports a
    sequencing error. zt context must equal the context of the cytosine on
    the read's haplotype, and zt variant flags must mark that haplotype's SNVs.
    zt also reports the molecule state of cytosines on the opposite strand;
    those carry their own context and are never converted in this read.
    """
    audit = BisulfiteAudit()
    getters = {}
    for record in run.bam.records:
        contig = run.contig_of(record)
        haplotype = haplotype_of(record)
        key = (contig, haplotype)
        if key not in getters:
            getters[key] = (
                haplotypes.base_getter(contig, haplotype)
                if haplotypes is not None
                else GENOME[contig].__getitem__
            )
        base_at = getters[key]
        snvs = haplotypes.snvs.get(key, {}) if haplotypes is not None else {}
        length = len(GENOME[contig])
        reverse_strand = record.tags["XG"] == "GA"
        cytosine = "G" if reverse_strand else "C"
        opposite_cytosine = "C" if reverse_strand else "G"
        converted_base = "A" if reverse_strand else "T"
        strand = "-" if reverse_strand else "+"
        opposite_strand = "+" if reverse_strand else "-"
        states = record.tags["zt"]
        fragment = audit.fragment_states[record.query_name]
        for query_index, position, operation in record.aligned_pairs():
            if operation != "M":
                audit.note("unexpected CIGAR operation {} in {}".format(
                    operation, record.query_name))
                audit.base_mismatches += 1
                continue
            audit.bases += 1
            state = ZT_VALUES[states[query_index]]
            context = state & 3
            methylated = bool(state & 4)
            converted = bool(state & 8)
            if bool(state & 16) != (position in snvs):
                audit.variant_flag_mismatches += 1
            truth = base_at(position)
            if truth == cytosine:
                expected_context = context_at(base_at, position, length, reverse_strand)
                if expected_context is not None:
                    audit.context_checked += 1
                    if context != expected_context:
                        audit.context_mismatches += 1
                        audit.note("context {}!={} at {}:{}{} in {}".format(
                            context, expected_context, contig, position, strand,
                            record.query_name))
                if methylated and converted:
                    audit.methylated_and_converted += 1
                if context:
                    if not methylated:
                        audit.unmethylated_sites += 1
                        audit.converted_unmethylated_sites += converted
                    counts = audit.site_states[(contig, position, strand, haplotype)]
                    counts[0] += methylated
                    counts[1] += 1
                    previous = fragment.setdefault(
                        (position, strand), (methylated, converted))
                    if previous != (methylated, converted):
                        audit.conflicting_mate_states += 1
                expected = converted_base if converted else truth
            elif truth == opposite_cytosine:
                expected_context = context_at(
                    base_at, position, length, not reverse_strand)
                if expected_context is not None:
                    audit.context_checked += 1
                    if context != expected_context:
                        audit.context_mismatches += 1
                        audit.note("context {}!={} at {}:{}{} in {}".format(
                            context, expected_context, contig, position,
                            opposite_strand, record.query_name))
                if converted:
                    audit.opposite_strand_converted += 1
                if context:
                    previous = fragment.setdefault(
                        (position, opposite_strand), (methylated, converted))
                    if previous != (methylated, converted):
                        audit.conflicting_mate_states += 1
                expected = truth
            else:
                if context or methylated or converted:
                    audit.flags_on_non_cytosine += 1
                expected = truth
            observed = record.sequence[query_index]
            if state & 32:
                audit.sequencing_errors += 1
                audit.errors_equal_to_truth += observed == expected
            elif observed != expected:
                audit.base_mismatches += 1
                audit.note("base {}!={} at {}:{} in {}".format(
                    observed, expected, contig, position, record.query_name))
    return audit


def write_text(name: str, text: str) -> Path:
    path = WORK / name
    path.write_text(text)
    return path
