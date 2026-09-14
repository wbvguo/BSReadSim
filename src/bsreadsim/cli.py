"""The ``bsreadsim`` command line: its parser, the ``Settings`` a parsed command
resolves to, and the commands.

Options default to None, so that resolution knows what was given; the field
defaults of ``Settings`` are the defaults shown in the help. The run imports
NumPy only when it starts, so ``--help`` and ``--version`` stay light.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Mapping, Sequence
import json
import os
from pathlib import Path
import secrets
import sys

from . import __version__
from .errors import BSReadSimError
from .htsim import HtsimCore
from .input import INPUT_FORMATS, validate_inputs
from .settings import (
    BISULFITE, DEFAULT_CUT_SITES, DEFAULTS, WHOLE_GENOME, Settings, derive_stage_seed,
)

UINT64_MAX = (1 << 64) - 1
UINT32_MAX = (1 << 32) - 1
TECHNOLOGIES = {
    "wgbs": "simulate whole-genome bisulfite sequencing",
    "rrbs": "simulate reduced-representation bisulfite sequencing",
    "tbs": "simulate targeted bisulfite sequencing",
    "wgs": "simulate whole-genome sequencing",
    "wes": "simulate whole-exome capture sequencing",
    "ts": "simulate targeted capture sequencing",
}
CUT_SITE_HELP = ("RRBS motif with | marking the cut position (default: C|CGG); "
                 "separate multiple motifs with commas")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the BSReadSim command-line interface."""
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = build_parser()
    arguments = parser.parse_args(argv)
    if arguments.command is None:
        parser.print_help()
        return 2
    commands = {"run": _run, "validate": _validate, "build": _build, "export": _export}
    try:
        return commands[arguments.command](arguments, argv)
    except (BSReadSimError, OSError) as error:
        print(f"bsreadsim: error: {error}", file=sys.stderr)
        return 1


# ---------------------------------------------------------------- commands

def _run(arguments: argparse.Namespace, argv: Sequence[str]) -> int:
    from .pipeline import run

    result = run(from_arguments(arguments), ["bsreadsim", *argv], arguments.core)
    print(result.manifest_path)
    return 0


def _validate(arguments: argparse.Namespace, argv: Sequence[str]) -> int:
    summary = validate_inputs(from_arguments(arguments), HtsimCore.find(arguments.core))
    print(json.dumps(summary, indent=2, sort_keys=True) if arguments.json
          else _validation_text(summary))
    skipped = summary["vcf"]["skipped"]["total"] if summary.get("vcf") else 0
    if arguments.strict and skipped:
        print(f"bsreadsim: error: strict validation failed: {skipped} unsupported VCF records "
              "would be skipped", file=sys.stderr)
        return 1
    return 0


def _build(arguments: argparse.Namespace, argv: Sequence[str]) -> int:
    if arguments.build_target == "variants":
        _require_suffix(arguments.output, ".vcf.gz", "variant VCF output")
    print(HtsimCore.find(arguments.core).build(
        from_arguments(arguments), arguments.build_target, arguments.output))
    return 0


def _export(arguments: argparse.Namespace, argv: Sequence[str]) -> int:
    _require_suffix(arguments.output, ".bed" if arguments.no_compression else ".bed.gz",
                    "MethDB BED output")
    print(HtsimCore.find(arguments.core).export_methdb(arguments.input, arguments.output))
    return 0


def _require_suffix(path: Path, suffix: str, what: str) -> None:
    if not str(path).lower().endswith(suffix):
        raise BSReadSimError(f"{what} path must end in {suffix}")


# ---------------------------------------------------------------- settings

def from_arguments(arguments: argparse.Namespace,
                   entropy: Callable[[int], int] = secrets.randbits) -> Settings:
    """The settings of a parsed ``run``, ``build``, or ``validate`` command: defaults filled
    in, paths absolute, and omitted seeds drawn (the master seed) or derived from it."""
    given = {name: value for name, value in vars(arguments).items() if value is not None}
    run = arguments.command == "run"
    target = getattr(arguments, "build_target", None)
    technology = (arguments.run_technology.upper() if run
                  else "RRBS" if target == "rrbs" else "WGBS")
    fields: dict = {"technology": technology}

    for name in INPUT_FORMATS:
        if name in given:
            fields[name] = _absolute(given[name])
    if run:
        fields["output"] = _absolute(given["output_directory"])
    for name in ("max_ambiguous_fraction", "indel_fraction",
                 "indel_extension_probability", "homozygous_only", "cpg_only", "pool_meth",
                 "meth_model", "sampling", "center_sd", "conversion_rate", "threads",
                 "prefix", "format", "gzip_level", "depth", "beta_cg", "beta_chg", "beta_chh",
                 "seed_mut", "seed_phase", "seed_meth"):
        if name in given:
            fields[name] = tuple(given[name]) if name.startswith("beta") else given[name]
    if "read_length" in given:  # two equal lengths are one
        lengths = tuple(given["read_length"])
        fields["read_length"] = lengths[:1] if len(set(lengths)) == 1 else lengths
    fields["paired_end"] = not given.get("single_end", False)
    fields["directional"] = not given.get("undirectional", False)
    if technology not in BISULFITE:
        fields["conversion_rate"] = 0.0
    if technology == "RRBS":
        fields["cut_sites"] = _cut_sites(given.get("cut_sites"))

    # A fixed insert (SD 0) needs no bounds beyond its length.
    mean = fields["insert_mean"] = given.get("insert_mean", DEFAULTS.insert_mean)
    fields["insert_sd"] = given.get("insert_sd", DEFAULTS.insert_sd)
    collapse = run and technology != "RRBS" and fields["insert_sd"] == 0
    fields["insert_min"] = given.get("insert_min", mean if collapse else DEFAULTS.insert_min)
    fields["insert_max"] = given.get("insert_max", mean if collapse else DEFAULTS.insert_max)

    fields["mutation_rate"] = given.get("mutation_rate", _default_mutation_rate(given, run, target))

    if "depth" in given:
        fields["reads"] = None
    elif run:
        reads = given.get("reads", DEFAULTS.reads)
        mates = 2 if fields["paired_end"] else 1
        if reads < 1:
            raise BSReadSimError("--reads must be a positive integer")
        if reads % mates:
            raise BSReadSimError("--reads must be even for paired-end simulation")
        if reads // mates > UINT32_MAX:
            raise BSReadSimError(f"--reads exceeds {UINT32_MAX * mates}")
        fields["reads"] = reads

    # a model replaces the uniform value
    fields["phred"] = None if "quality_model" in given else given.get("phred", DEFAULTS.phred)
    fields["error_rate"] = (None if "error_model" in given
                            else given.get("error_rate", DEFAULTS.error_rate))

    realization = given.get("fragment_realization", False)
    fields["fragment_realization"] = realization
    fields["fragment_summary"] = given.get("fragment_summary", False) or realization
    truth = given.get("save_truth", False)
    fields["save_methdb"] = technology in BISULFITE and (given.get("save_methdb", False) or truth)
    fields["save_vcf"] = given.get("save_vcf", False) or truth

    if run:
        seed = given["seed"] if "seed" in given else entropy(64)
        fields["seed"] = seed
        for name, domain in (("seed_mut", "mutation"), ("seed_phase", "phasing"),
                             ("seed_meth", "methylation")):
            if name not in given:
                fields[name] = derive_stage_seed(seed, domain)
    return Settings(**fields)


def _default_mutation_rate(given: dict, run: bool, target: str | None) -> float:
    variant_source = any(name in given for name in ("vcf", "asm", "asm_bed", "methdb"))
    if run:
        fixed_domain = "rrbs_candidates" in given or "gc_profile" in given
        return 0.0 if variant_source or fixed_domain else 0.001
    if target in ("variants", "methdb"):
        return 0.0 if variant_source else 0.001
    return 0.0


def _cut_sites(values: Sequence[str] | None) -> tuple[str, ...]:
    if not values:
        return DEFAULT_CUT_SITES
    if len(values) != 1:
        raise BSReadSimError("--cut-site may be supplied once; separate multiple sites with commas")
    sites = [site.upper() for site in values[0].split(",")]
    if not all(sites):
        raise BSReadSimError("--cut-site contains an empty value")
    if len(set(sites)) != len(sites):
        raise BSReadSimError("--cut-site repeats a site")
    return tuple(sites)


def _absolute(path: os.PathLike | str) -> Path:
    return Path(path).expanduser().resolve()


def _validation_text(summary: Mapping[str, object]) -> str:
    """The text form of a ``validate`` summary."""
    contigs = summary["reference"]["contigs"]
    lines = ["status: valid",
             f"reference: {contigs} contigs, {summary['reference']['bases']} bases"]
    vcf = summary.get("vcf")
    if vcf:
        skipped = vcf["skipped"]
        lines.append(f"VCF: {vcf['rows']} rows on {vcf['contigs']}/{contigs} contigs; "
                     f"{vcf['retained']} retained, {vcf['reference_genotypes']} "
                     f"reference-genotype, {skipped['total']} unsupported skipped")
        if skipped["total"]:
            lines.append(f"  skipped: {skipped['mnp']} MNP, {skipped['complex_replacement']} "
                         f"complex replacement, {skipped['long_indel']} >4 bp indel")
    methylation = summary.get("methylation")
    if methylation:
        lines.append(f"{methylation['format']}: {methylation['rows']} rows on "
                     f"{methylation['contigs']}/{contigs} contigs; "
                     f"{methylation['defined_probabilities']} defined probabilities")
    asm = summary.get("asm")
    if asm:
        lines.append(f"{asm['format']}: {asm['rows']} rows on {asm['contigs']}/{contigs} contigs")
    return "\n".join(lines)


# ---------------------------------------------------------------- parser

def build_parser() -> argparse.ArgumentParser:
    """The public command-line parser."""
    parser = argparse.ArgumentParser(
        prog="bsreadsim", allow_abbrev=False,
        description="Generate biological fragments with htsim and apply BSReadSim "
                    "optional chemistry and sequencing effects.")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    commands = parser.add_subparsers(dest="command")

    run = _subparser(commands, "run", help="simulate directly from command-line arguments",
                     description="Simulate reads from CLI arguments.")
    assays = run.add_subparsers(dest="run_technology", required=True)
    for name, help_text in TECHNOLOGIES.items():
        _add_run_arguments(_subparser(assays, name, help=help_text), name)

    validate = _subparser(
        commands, "validate",
        help="validate reference-coordinate inputs without generating output",
        description="Validate a reference and optional VCF, methylation profile, and ASM "
                    "inputs through the same native boundaries used for generation.")
    _add_validate_arguments(validate)

    build = _subparser(commands, "build", help="build a reusable truth or sampling artifact")
    targets = build.add_subparsers(dest="build_target", required=True)
    _add_build_rrbs_arguments(_subparser(targets, "rrbs", help="export exact RRBS candidates",
                                         description="Build the htsim RRBS candidate BED."))
    _add_build_variants_arguments(_subparser(
        targets, "variants", help="build a normalized, phased VCF.gz",
        description="Build the exact deterministic de novo variant set, or normalize and "
                    "phase --vcf input, as a gzip-compressed one-sample VCF."))
    _add_build_methdb_arguments(_subparser(
        targets, "methdb", help="build a reusable prepared methylation profile as MethDB"))

    export = _subparser(commands, "export", help="decode a BSReadSim artifact")
    export_targets = export.add_subparsers(dest="export_target", required=True)
    _add_export_methdb_arguments(_subparser(
        export_targets, "methdb",
        help="decode a MethDB snapshot as a human-readable extended BED"))
    return parser


def _subparser(commands, name: str, **keywords) -> argparse.ArgumentParser:
    return commands.add_parser(name, allow_abbrev=False, **keywords)


def _add_run_arguments(parser: argparse.ArgumentParser, name: str) -> None:
    technology = name.upper()
    bisulfite = technology in BISULFITE
    required = parser.add_argument_group("required run inputs")
    _add_reference(required)
    required.add_argument("-o", "--output", dest="output_directory", type=Path, required=True,
                          help="new output directory")
    count = parser.add_argument_group("read count").add_mutually_exclusive_group()
    count.add_argument("-n", "--reads", type=int, metavar="N",
                       help="total number of read records; paired-end values must be even "
                            f"(default: {DEFAULTS.reads})")
    count.add_argument("-d", "--depth", type=float,
                       help="requested mean depth over the technology target region")

    if bisulfite:
        _add_methylation_inputs(parser)
    _add_seed_arguments(parser, run=True, methylation=bisulfite)

    sampling = parser.add_argument_group("sampling")
    if technology in WHOLE_GENOME:
        sampling.add_argument("--sampling", choices=("uniform", "gc"),
                              help=f"fragment sampling (default: {DEFAULTS.sampling})")
        sampling.add_argument("--gc-profile", type=Path,
                              help="target fragment-GC distribution for --sampling gc")
    elif technology == "RRBS":
        sampling.add_argument("--sampling", choices=("uniform", "score"),
                              help=f"fragment sampling (default: {DEFAULTS.sampling})")
        sampling.add_argument("--cut-site", action="append", dest="cut_sites",
                              metavar="CUT_SITE", help=CUT_SITE_HELP)
        sampling.add_argument("--rrbs-candidates", type=Path,
                              help="candidate BED produced by 'bsreadsim build rrbs'")
    else:
        sampling.add_argument("--sampling", choices=("uniform", "score"),
                              help=f"fragment sampling (default: {DEFAULTS.sampling})")
        sampling.add_argument("--targets", type=Path, required=True, help="BED6 capture targets")
        sampling.add_argument("--center-sd", type=float,
                              help="fragment-center displacement SD "
                                   f"(default: {DEFAULTS.center_sd:g})")

    _add_fragment_arguments(parser)
    _add_variant_arguments(parser)
    if bisulfite:
        _add_methylation_arguments(parser)

    sequencing = parser.add_argument_group("sequencing")
    if bisulfite:
        sequencing.add_argument("--conversion-rate", type=float,
                                help=f"bisulfite conversion rate "
                                     f"(default: {DEFAULTS.conversion_rate:g})")
        sequencing.add_argument("--undirectional", action="store_true",
                                help="sample OT, OB, CTOT, and CTOB instead of directional OT/OB")
    quality = sequencing.add_mutually_exclusive_group()
    quality.add_argument("-q", "--phred", type=int,
                         help=f"uniform Phred score (default: {DEFAULTS.phred})")
    quality.add_argument("--quality-model", type=Path, help="quality Markov JSON")
    errors = sequencing.add_mutually_exclusive_group()
    errors.add_argument("-e", "--error-rate", type=float,
                        help=f"uniform substitution rate (default: {DEFAULTS.error_rate:g})")
    errors.add_argument("--error-model", type=Path, help="quality-confusion JSON")

    execution = parser.add_argument_group("execution")
    execution.add_argument("-t", "--threads", type=int, metavar="N",
                           help=f"number of threads to use (default: {DEFAULTS.threads})")

    output = parser.add_argument_group("output")
    output.add_argument("-p", "--prefix", help=f"output file prefix (default: {DEFAULTS.prefix})")
    output.add_argument("-f", "--format", choices=("fastq", "fastq.gz", "bam", "sam"),
                        help=f"read output format (default: {DEFAULTS.format})")
    output.add_argument("--gzip-level", type=int,
                        help=f"gzip and BAM compression level (default: {DEFAULTS.gzip_level})")
    output.add_argument("--fragment-summary", action="store_true",
                        help="add optional full-fragment zf summaries to BAM/SAM records; "
                             "zt and zr are always present")
    if bisulfite:
        output.add_argument("--fragment-realization", action="store_true",
                            help="emit complete-fragment methylation/conversion state in "
                                 "BAM/SAM zx")
        output.add_argument("--save-methdb", action="store_true",
                            help="save the run's prepared methylation profile and variants as "
                                 "a MethDB truth artifact")
    output.add_argument("--save-vcf", action="store_true",
                        help="save the prepared variant set as a phased VCF truth artifact")
    output.add_argument("--save-truth", action="store_true",
                        help="save the prepared variant set and methylation profile as "
                             "simulation truth artifacts" if bisulfite
                        else "save the prepared variant set as a simulation truth artifact")
    _add_core_option(parser)


def _add_validate_arguments(parser: argparse.ArgumentParser) -> None:
    _add_reference(parser)
    parser.add_argument_group("variants").add_argument(
        "--vcf", type=Path, help="one-sample diploid VCF input")
    _add_methylation_inputs(parser, methdb=False)
    checks = parser.add_argument_group("methylation validation")
    checks.add_argument("--cpg-only", action="store_true",
                        help="validate ASM targets under a CG-only methylation domain")
    checks.add_argument("--pool-meth", action="store_true",
                        help="require a defined value in the text methylation profile")
    parser.add_argument("--seed-phase", type=_parse_seed,
                        help="seed for phasing unphased VCF/ASM heterozygotes (default: 0)")
    report = parser.add_argument_group("output")
    report.add_argument("--json", action="store_true", help="write the validation summary as JSON")
    report.add_argument("--strict", action="store_true",
                        help="fail when unsupported VCF records would be skipped")
    _add_core_option(parser)


def _add_build_rrbs_arguments(parser: argparse.ArgumentParser) -> None:
    _add_reference(parser)
    parser.add_argument("-o", "--output", type=Path, required=True, help="new candidate BED path")
    parser.add_argument_group("RRBS").add_argument(
        "--cut-site", action="append", dest="cut_sites", metavar="CUT_SITE", help=CUT_SITE_HELP)
    _add_fragment_arguments(parser)
    _add_variant_arguments(parser)
    _add_seed_arguments(parser, run=False, methylation=False)
    _add_core_option(parser)


def _add_build_variants_arguments(parser: argparse.ArgumentParser) -> None:
    _add_reference(parser)
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="new gzip-compressed .vcf.gz path")
    _add_variant_arguments(parser)
    _add_seed_arguments(parser, run=False, methylation=False)
    _add_core_option(parser)


def _add_build_methdb_arguments(parser: argparse.ArgumentParser) -> None:
    _add_reference(parser)
    parser.add_argument("-o", "--output", type=Path, required=True, help="new .methdb path")
    _add_methylation_inputs(parser, methdb=False)
    _add_variant_arguments(parser)
    _add_methylation_arguments(parser)
    _add_seed_arguments(parser, run=False)
    _add_core_option(parser)


def _add_export_methdb_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-i", "--input", type=Path, required=True, help="input MethDB path")
    parser.add_argument("-o", "--output", type=Path, required=True,
                        help="new .bed.gz path, or .bed with --no-compression")
    parser.add_argument("--no-compression", action="store_true",
                        help="write plain .bed instead of BGZF")
    _add_core_option(parser)


# ---------------------------------------------------------------- option groups

def _add_core_option(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--core", type=Path, help="override the bundled htsim executable")


def _add_reference(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("-r", "--reference", type=Path, required=True, help="reference FASTA")


def _add_fragment_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("fragments")
    group.add_argument("--single-end", action="store_true",
                       help="emit one read per fragment instead of paired-end reads")
    group.add_argument("-l", "--read-length", type=_read_lengths, metavar="LENGTH[,LENGTH]",
                       help="read length; READ1,READ2 gives paired mates different lengths "
                            f"(default: {DEFAULTS.read_length[0]})")
    group.add_argument("--insert-min", type=int,
                       help=f"lower insert bound (default: {DEFAULTS.insert_min})")
    group.add_argument("--insert-mean", type=int,
                       help="fixed insert length or variable-distribution mean "
                            f"(default: {DEFAULTS.insert_mean})")
    group.add_argument("--insert-max", type=int,
                       help=f"upper insert bound (default: {DEFAULTS.insert_max})")
    group.add_argument("--insert-sd", type=float,
                       help="insert-length SD; 0 selects a fixed insert for non-RRBS "
                            f"assays (default: {DEFAULTS.insert_sd:g})")
    group.add_argument("--max-ambiguous-fraction", type=float,
                       help="maximum N fraction in each emitted mate "
                            f"(default: {DEFAULTS.max_ambiguous_fraction:g})")


def _add_variant_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("variants")
    source = group.add_mutually_exclusive_group()
    source.add_argument("--vcf", type=Path,
                        help="one-sample diploid VCF input instead of de novo mutations")
    source.add_argument("--mutation-rate", type=float,
                        help="total de novo mutation-event rate (run/build variants default: "
                             "0.001; build RRBS default: 0)")
    group.add_argument("--indel-fraction", type=float,
                       help=f"share of indels among mutations "
                            f"(default: {DEFAULTS.indel_fraction:g})")
    group.add_argument("--indel-extension-probability", type=float,
                       help=f"chance an indel grows by one more base "
                            f"(default: {DEFAULTS.indel_extension_probability:g})")
    group.add_argument("--homozygous-only", action="store_true")


def _add_seed_arguments(parser: argparse.ArgumentParser, *, run: bool,
                        methylation: bool = True) -> None:
    group = parser.add_argument_group("random seeds")
    if run:
        group.add_argument("--seed", type=_parse_seed,
                           help="master unsigned 64-bit seed; omit to generate and record one")
    origin = "omit to derive from --seed" if run else "default: 0"
    group.add_argument("--seed-mut", type=_parse_seed,
                       help=f"de novo mutation seed ({origin})")
    group.add_argument("--seed-phase", type=_parse_seed,
                       help=f"seed for phasing unphased VCF heterozygotes ({origin})")
    if methylation:
        group.add_argument("--seed-meth", type=_parse_seed,
                           help=f"methylation-probability seed ({origin})")


def _add_methylation_inputs(parser: argparse.ArgumentParser, *, methdb: bool = True) -> None:
    group = parser.add_argument_group("methylation inputs")
    profile = group.add_mutually_exclusive_group()
    profile.add_argument("--cgmap", type=Path, help="CGmap methylation input")
    profile.add_argument("--bedmethyl", dest="bed_methyl", type=Path,
                         help="UCSC/ENCODE bedMethyl (BED9+2 or BED9+9) input")
    profile.add_argument("--methbg", type=Path, help="four-column MethBG methylation input")
    profile.add_argument("--methbed", type=Path, help="MethBED methylation input")
    if methdb:
        profile.add_argument("--methdb", type=Path,
                             help="prepared methylation profile and embedded variants (MethDB)")
    asm = group.add_mutually_exclusive_group()
    asm.add_argument("--asm", type=Path,
                     help="CGmapTools 'asm -m ass' allele-specific methylation output")
    asm.add_argument("--asm-bed", type=Path,
                     help="BSReadSim ASM BED6+6 or BED6+10 allele-specific methylation input")


def _add_methylation_arguments(parser: argparse.ArgumentParser) -> None:
    group = parser.add_argument_group("methylation")
    for context in ("cg", "chg", "chh"):
        alpha, beta = getattr(DEFAULTS, "beta_" + context)
        group.add_argument("--beta-" + context, type=_beta_pair, metavar="ALPHA,BETA",
                           help=f"Beta shape of {context.upper()} levels "
                                f"(default: {alpha:g},{beta:g})")
    group.add_argument("--cpg-only", action="store_true",
                       help="omit CHG and CHH sites from the methylation profile")
    group.add_argument("--pool-meth", action="store_true",
                       help="pool text-profile values by contig and context")
    group.add_argument("--meth-model", choices=("bernoulli", "bilstm"),
                       help="fragment-level methylation state model "
                            f"(default: {DEFAULTS.meth_model})")


# ---------------------------------------------------------------- value types

def _parse_seed(text: str) -> int:
    """An unsigned 64-bit decimal seed."""
    if not text.isdecimal() or (len(text) > 1 and text.startswith("0")) \
            or int(text) > UINT64_MAX:
        raise argparse.ArgumentTypeError("must be an unsigned 64-bit decimal integer")
    return int(text)


def _beta_pair(value: str) -> tuple[float, float]:
    parts = value.split(",")
    try:
        if len(parts) == 2 and all(parts):
            return float(parts[0]), float(parts[1])
    except ValueError:
        pass
    raise argparse.ArgumentTypeError("must be ALPHA,BETA")


def _read_lengths(value: str) -> tuple[int, ...]:
    """LENGTH for every read, or READ1,READ2 for paired mates of different lengths."""
    parts = value.split(",")
    if 1 <= len(parts) <= 2 and all(part.isdecimal() for part in parts):
        return tuple(int(part) for part in parts)
    raise argparse.ArgumentTypeError("must be LENGTH or READ1,READ2")
