"""Command-line entry point for BSReadSim."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys
from collections.abc import Mapping, Sequence

from . import __version__


_BISULFITE_TECHNOLOGIES = frozenset(("WGBS", "RRBS", "TBS"))
_WHOLE_GENOME_TECHNOLOGIES = frozenset(("WGBS", "WGS"))
_TARGETED_TECHNOLOGIES = frozenset(("TBS", "WES", "TS"))
_UINT32_MAX = (1 << 32) - 1


class CommandLineError(ValueError):
    """Direct simulation arguments cannot form one valid run."""


def _parse_beta_pair(value: str) -> list[float]:
    parts = value.split(",")
    if len(parts) != 2 or not all(parts):
        raise argparse.ArgumentTypeError("must be ALPHA,BETA")
    try:
        return [float(parts[0]), float(parts[1])]
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be ALPHA,BETA") from error


def _add_runtime_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--core",
        type=Path,
        help="override the bundled htsim-core executable",
    )


def _add_fragment_domain_arguments(
    parser: argparse.ArgumentParser, *, targeted: bool = False
) -> None:
    """Add fragment-shape arguments shared by simulation and build commands."""
    fragments = parser.add_argument_group("fragments")
    fragments.add_argument(
        "--single-end",
        action="store_true",
        help="emit one read per fragment instead of paired-end reads",
    )
    fragments.add_argument(
        "-l", "--read-length", type=int, default=100, help="read length (default: 100)"
    )
    fragments.add_argument(
        "--insert-min", type=int, help="lower insert bound (default: 100)"
    )
    fragments.add_argument(
        "--insert-mean",
        type=int,
        help="fixed insert length or variable-distribution mean (default: 400)",
    )
    fragments.add_argument(
        "--insert-max", type=int, help="upper insert bound (default: 1000)"
    )
    fragments.add_argument(
        "--insert-sd",
        type=float,
        default=0.0 if targeted else None,
        help=(
            "insert-length SD; 0 selects a fixed whole-genome/targeted "
            "insert (default: {} for this command)".format(
                0 if targeted else 25
            )
        ),
    )
    fragments.add_argument(
        "--max-ambiguous-fraction",
        type=float,
        default=0.05,
        help="maximum N fraction in each emitted mate (default: 0.05)",
    )


def _add_variant_domain_arguments(parser: argparse.ArgumentParser) -> None:
    """Add de novo mutation arguments shared by run and build commands."""
    mutation = parser.add_argument_group("variants")
    source = mutation.add_mutually_exclusive_group()
    source.add_argument(
        "--vcf",
        type=Path,
        help="one-sample diploid VCF input instead of de novo mutations",
    )
    source.add_argument(
        "--mutation-rate",
        type=float,
        help=(
            "total de novo mutation-event rate (run/build variants default: "
            "0.001; build RRBS default: 0)"
        ),
    )
    mutation.add_argument("--indel-fraction", type=float, default=0.15)
    mutation.add_argument(
        "--indel-extension-probability", type=float, default=0.15
    )
    mutation.add_argument("--homozygous-only", action="store_true")


def _add_seed_arguments(
    parser: argparse.ArgumentParser,
    *,
    include_master: bool,
    include_methylation: bool = True,
) -> None:
    seeds = parser.add_argument_group("random seeds")
    if include_master:
        seeds.add_argument(
            "-s",
            "--seed",
            help="read-simulation unsigned 64-bit seed; omit to generate one",
        )
    seeds.add_argument(
        "--seed-mut",
        default="0",
        help="de novo mutation seed (default: 0)",
    )
    seeds.add_argument(
        "--seed-phase",
        default="0",
        help="seed for phasing unphased VCF heterozygotes (default: 0)",
    )
    if include_methylation:
        seeds.add_argument(
            "--seed-meth",
            default="0",
            help="methylation-probability seed (default: 0)",
        )


def _add_methylation_input_arguments(
    parser: argparse.ArgumentParser, *, include_methdb: bool = True
) -> None:
    biology = parser.add_argument_group("methylation inputs")
    methylation_input = biology.add_mutually_exclusive_group()
    methylation_input.add_argument(
        "--cgmap", type=Path, help="CGmap methylation input"
    )
    methylation_input.add_argument(
        "--bed-methyl",
        type=Path,
        help="UCSC/ENCODE bedMethyl (BED9+2 or BED9+9) input",
    )
    if include_methdb:
        methylation_input.add_argument(
            "--methdb",
            type=Path,
            help=(
                "prepared methylation profile and embedded variants in "
                "MethDB v2 format"
            ),
        )
    asm_input = biology.add_mutually_exclusive_group()
    asm_input.add_argument(
        "--asm", type=Path, help="htsim allele-specific methylation input"
    )
    asm_input.add_argument(
        "--asm-bed",
        type=Path,
        help="htsim ASM BED6+6 allele-specific methylation input",
    )


def _add_methylation_arguments(parser: argparse.ArgumentParser) -> None:
    methylation = parser.add_argument_group("methylation")
    methylation.add_argument(
        "--beta-cg", type=_parse_beta_pair, metavar="ALPHA,BETA"
    )
    methylation.add_argument(
        "--beta-chg", type=_parse_beta_pair, metavar="ALPHA,BETA"
    )
    methylation.add_argument(
        "--beta-chh", type=_parse_beta_pair, metavar="ALPHA,BETA"
    )
    methylation.add_argument(
        "--cpg-only",
        action="store_true",
        help="omit CHG and CHH sites from the generated MethDB",
    )
    methylation.add_argument("--cgmap-pool", action="store_true")
    methylation.add_argument(
        "--methylation-model",
        choices=("bernoulli", "bilstm"),
        default="bernoulli",
        help="methylation state model (BiLSTM currently warns and falls back)",
    )
    methylation.add_argument(
        "--no-update-variant-boundaries",
        action="store_true",
        help="do not recompute methylation contexts affected by variants",
    )


def _add_direct_run_arguments(
    parser: argparse.ArgumentParser, technology: str
) -> None:
    bisulfite = technology in _BISULFITE_TECHNOLOGIES
    required = parser.add_argument_group("required run inputs")
    required.add_argument(
        "-r", "--reference", type=Path, required=True, help="reference FASTA"
    )
    required.add_argument(
        "-o",
        "--output",
        dest="output_directory",
        type=Path,
        required=True,
        help="new output directory",
    )
    count = required.add_mutually_exclusive_group(required=True)
    count.add_argument(
        "-n",
        "--reads",
        type=int,
        metavar="N",
        help="total number of read records; paired-end values must be even",
    )
    count.add_argument(
        "-d",
        "--depth",
        type=float,
        help="requested mean depth over the technology target region",
    )

    if bisulfite:
        _add_methylation_input_arguments(parser)
    _add_seed_arguments(
        parser,
        include_master=True,
        include_methylation=bisulfite,
    )

    sampling = parser.add_argument_group("sampling")
    if technology in _WHOLE_GENOME_TECHNOLOGIES:
        sampling.add_argument(
            "--gc-profile",
            type=Path,
            help="target fragment-GC distribution; omit for uniform sampling",
        )
    elif technology == "RRBS":
        sampling.add_argument(
            "--cut-site",
            action="append",
            dest="cut_sites",
            required=True,
            help="RRBS cut motif, for example C|CGG; repeat for multiple enzymes",
        )
        sampling.add_argument(
            "--rrbs-candidates",
            type=Path,
            help="candidate BED produced by 'bsreadsim build rrbs'",
        )
        sampling.add_argument(
            "--sampling", choices=("uniform", "score"), default="uniform"
        )
    elif technology in _TARGETED_TECHNOLOGIES:
        sampling.add_argument(
            "--targets", type=Path, required=True, help="BED6 capture targets"
        )
        sampling.add_argument(
            "--sampling", choices=("uniform", "score"), default="uniform"
        )
        sampling.add_argument(
            "--fragment-center-stddev",
            type=float,
            default=50.0,
            help="capture fragment-center standard deviation (default: 50)",
        )
    else:
        raise CommandLineError("unsupported run technology")

    _add_fragment_domain_arguments(
        parser, targeted=technology in _TARGETED_TECHNOLOGIES
    )
    _add_variant_domain_arguments(parser)

    if bisulfite:
        _add_methylation_arguments(parser)

    sequencing = parser.add_argument_group("sequencing")
    if bisulfite:
        sequencing.add_argument(
            "--conversion-rate", type=float, default=0.998
        )
        sequencing.add_argument(
            "--undirectional",
            action="store_true",
            help="sample OT, OB, CTOT, and CTOB instead of directional OT/OB",
        )
    sequencing.add_argument(
        "-q", "--phred", type=int, help="uniform Phred score (default: 40)"
    )
    sequencing.add_argument(
        "--quality-model",
        type=Path,
        help="quality Markov JSON; identity is hashed automatically",
    )
    sequencing.add_argument(
        "-e", "--error-rate",
        type=float,
        help="uniform substitution rate (default: 0.005)",
    )
    sequencing.add_argument(
        "--error-model",
        type=Path,
        help="quality-confusion JSON; identity is hashed automatically",
    )

    execution = parser.add_argument_group("execution")
    execution.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        metavar="N",
        help="total CPU thread budget for the simulation pipeline (default: 1)",
    )

    output = parser.add_argument_group("output")
    output.add_argument("-p", "--prefix", default="sim")
    output.add_argument(
        "-f",
        "--format",
        choices=("fastq", "fastq.gz", "bam"),
        default="fastq.gz",
        help="read output format (default: fastq.gz)",
    )
    output.add_argument("--gzip-level", type=int, default=6)
    output.add_argument(
        "--fragment-summary",
        action="store_true",
        help=(
            "add optional full-fragment zf summaries to BAM records; "
            "zt and zr are always present"
        ),
    )
    if bisulfite:
        output.add_argument(
            "--fragment-realization",
            action="store_true",
            help="emit complete-fragment methylation/conversion state in BAM zx",
        )
        output.add_argument(
            "--save-methdb",
            action="store_true",
            help=(
                "save the run's prepared methylation/variant world as a "
                "MethDB truth artifact"
            ),
        )
    output.add_argument(
        "--save-vcf",
        action="store_true",
        help="save the prepared variant set as a phased VCF truth artifact",
    )
    output.add_argument(
        "--save-truth",
        action="store_true",
        help=(
            "save the prepared variant set and methylation profile as "
            "simulation truth artifacts"
            if bisulfite
            else "save the prepared variant set as a simulation truth artifact"
        ),
    )


def _add_rrbs_build_arguments(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("required build inputs")
    required.add_argument(
        "-r", "--reference", type=Path, required=True, help="reference FASTA"
    )
    required.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="new candidate BED path",
    )
    required.add_argument(
        "--cut-site",
        action="append",
        dest="cut_sites",
        required=True,
        help="RRBS cut motif, for example C|CGG; repeat for multiple enzymes",
    )

    _add_fragment_domain_arguments(parser)
    _add_variant_domain_arguments(parser)
    _add_seed_arguments(
        parser, include_master=False, include_methylation=False
    )
    _add_runtime_options(parser)


def _add_variant_build_arguments(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("required build inputs")
    required.add_argument(
        "-r", "--reference", type=Path, required=True, help="reference FASTA"
    )
    required.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="new gzip-compressed .vcf.gz path",
    )
    _add_variant_domain_arguments(parser)
    _add_seed_arguments(
        parser, include_master=False, include_methylation=False
    )
    _add_runtime_options(parser)


def _add_methdb_build_arguments(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("required build inputs")
    required.add_argument(
        "-r", "--reference", type=Path, required=True, help="reference FASTA"
    )
    required.add_argument(
        "-o", "--output", type=Path, required=True, help="new .methdb path"
    )
    _add_methylation_input_arguments(parser, include_methdb=False)
    _add_variant_domain_arguments(parser)
    _add_methylation_arguments(parser)
    _add_seed_arguments(parser, include_master=False)
    _add_runtime_options(parser)


def build_parser() -> argparse.ArgumentParser:
    """Build the public command-line parser without importing heavy runtime modules."""
    parser = argparse.ArgumentParser(
        prog="bsreadsim",
        description=(
            "Generate biological fragments with htsim and apply BSReadSim "
            "optional chemistry and sequencing effects."
        ),
    )
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
    )
    commands = parser.add_subparsers(dest="command")

    run_parser = commands.add_parser(
        "run",
        help="simulate directly from command-line arguments",
        description=(
            "Simulate reads from CLI arguments; JSON run-configuration "
            "files are not accepted."
        ),
    )
    run_commands = run_parser.add_subparsers(dest="run_technology", required=True)
    for name, technology, help_text in (
        ("wgbs", "WGBS", "simulate whole-genome bisulfite sequencing"),
        ("rrbs", "RRBS", "simulate reduced-representation bisulfite sequencing"),
        ("tbs", "TBS", "simulate targeted bisulfite sequencing"),
        ("wgs", "WGS", "simulate whole-genome sequencing"),
        ("wes", "WES", "simulate whole-exome capture sequencing"),
        ("ts", "TS", "simulate targeted capture sequencing"),
    ):
        technology_parser = run_commands.add_parser(name, help=help_text)
        _add_direct_run_arguments(technology_parser, technology)
        _add_runtime_options(technology_parser)

    build_parser = commands.add_parser(
        "build", help="build a reusable truth or sampling artifact"
    )
    build_commands = build_parser.add_subparsers(dest="build_target", required=True)
    rrbs_build_parser = build_commands.add_parser(
        "rrbs",
        help="export exact RRBS candidates directly from CLI arguments",
        description=(
            "Build the htsim RRBS candidate BED directly; no JSON "
            "configuration file is read or written."
        ),
    )
    _add_rrbs_build_arguments(rrbs_build_parser)
    variant_build_parser = build_commands.add_parser(
        "variants",
        help="build a normalized, phased VCF.gz",
        description=(
            "Build the exact deterministic de novo variant set, or normalize and "
            "phase --vcf input, as a gzip-compressed one-sample VCF."
        ),
    )
    _add_variant_build_arguments(variant_build_parser)
    methdb_build_parser = build_commands.add_parser(
        "methdb", help="build a reusable methylation profile as MethDB"
    )
    _add_methdb_build_arguments(methdb_build_parser)

    export_parser = commands.add_parser(
        "export",
        help="export bundled data or decode a BSReadSim artifact",
    )
    export_commands = export_parser.add_subparsers(
        dest="export_target", required=True
    )
    test_fasta_export = export_commands.add_parser(
        "test-fasta",
        help="export the bundled synthetic test FASTA",
    )
    test_fasta_export.add_argument(
        "-o", "--output", type=Path, required=True, help="new FASTA path"
    )
    methdb_export = export_commands.add_parser(
        "methdb",
        help="decode a MethDB methylation profile as a human-readable BED",
    )
    methdb_export.add_argument(
        "-i", "--input", type=Path, required=True, help="input MethDB path"
    )
    methdb_export.add_argument(
        "-o",
        "--output",
        type=Path,
        required=True,
        help="new .bed.gz path, or .bed with --no-compression",
    )
    methdb_export.add_argument(
        "--no-compression",
        action="store_true",
        help="write plain .bed instead of deterministic BGZF",
    )
    methdb_export.add_argument(
        "--core", type=Path, help="override the bundled htsim-core executable"
    )
    return parser


def _resolved_input(path: Path, base_directory: Path) -> Path:
    expanded = path.expanduser()
    if not expanded.is_absolute():
        expanded = base_directory / expanded
    return expanded.resolve(strict=False)


def _artifact_sha256(path: Path, base_directory: Path) -> str:
    resolved = _resolved_input(path, base_directory)
    digest = hashlib.sha256()
    try:
        with resolved.open("rb") as input_file:
            while True:
                block = input_file.read(1024 * 1024)
                if not block:
                    break
                digest.update(block)
    except OSError as error:
        raise CommandLineError(
            "cannot read artifact {}: {}".format(resolved, error)
        ) from error
    return digest.hexdigest()


def _artifact(
    path: Path,
    base_directory: Path,
) -> dict[str, str]:
    return {
        "path": str(path.expanduser()),
        "sha256": _artifact_sha256(path, base_directory),
    }


def _insert_parameters(
    arguments: argparse.Namespace,
    *,
    collapse_fixed_bounds: bool = False,
) -> Mapping[str, object]:
    insert_mean = 400 if arguments.insert_mean is None else arguments.insert_mean
    insert_sd = 25.0 if arguments.insert_sd is None else arguments.insert_sd
    fixed_insert = collapse_fixed_bounds and insert_sd == 0
    insert_min = 100 if arguments.insert_min is None else arguments.insert_min
    insert_max = 1000 if arguments.insert_max is None else arguments.insert_max
    if fixed_insert:
        if arguments.insert_min is None:
            insert_min = insert_mean
        if arguments.insert_max is None:
            insert_max = insert_mean
    return {
        "insert_min": insert_min,
        "insert_mean": insert_mean,
        "insert_max": insert_max,
        "insert_sd": insert_sd,
    }


def _resolve_mutation_rate(
    arguments: argparse.Namespace, *, default: float
) -> float:
    if arguments.vcf is not None:
        if arguments.mutation_rate is not None:
            raise CommandLineError(
                "--vcf and --mutation-rate are mutually exclusive"
            )
        return 0
    return default if arguments.mutation_rate is None else arguments.mutation_rate


def _read_quantity(arguments: argparse.Namespace) -> dict[str, object]:
    """Project the public read quantity without leaking fragment semantics."""
    if arguments.reads is None:
        return {"depth": arguments.depth}

    read_count = arguments.reads
    if read_count < 1:
        raise CommandLineError("--reads must be a positive integer")
    reads_per_fragment = 1 if arguments.single_end else 2
    if read_count % reads_per_fragment != 0:
        raise CommandLineError(
            "--reads must be even for paired-end simulation"
        )
    maximum = _UINT32_MAX * reads_per_fragment
    if read_count > maximum:
        raise CommandLineError(
            "--reads exceeds the maximum for {} simulation ({})".format(
                "single-end" if arguments.single_end else "paired-end",
                maximum,
            )
        )
    return {"count": read_count}


def build_rrbs_document(
    arguments: argparse.Namespace,
    base_directory: Path,
) -> dict[str, object]:
    """Project ``build rrbs`` arguments into an in-memory core config."""
    if (
        arguments.command != "build"
        or arguments.build_target != "rrbs"
    ):
        raise CommandLineError("build rrbs arguments are required")

    base = base_directory.expanduser().resolve(strict=False)
    mutation_rate = _resolve_mutation_rate(arguments, default=0)
    fragments: dict[str, object] = {
        "paired_end": not arguments.single_end,
        "read_length_1": arguments.read_length,
        **_insert_parameters(arguments),
        "max_ambiguous_fraction": arguments.max_ambiguous_fraction,
    }
    if not arguments.single_end:
        fragments["read_length_2"] = arguments.read_length

    document: dict[str, object] = {
        "reference": str(arguments.reference.expanduser()),
        "inputs": {},
        "technology": "RRBS",
        # Artifact construction emits no reads. This minimum complete quantity
        # exists only because build and run share one normalized core shape.
        "reads": {"count": 1 if arguments.single_end else 2},
        "rrbs": {"cut_sites": list(arguments.cut_sites)},
        "mutation": {
            "rate": mutation_rate,
            "indel_fraction": arguments.indel_fraction,
            "indel_extension_probability": (
                arguments.indel_extension_probability
            ),
            "homozygous_only": arguments.homozygous_only,
        },
        "seeds": {
            "mutation": arguments.seed_mut,
            "phasing": arguments.seed_phase,
            "methylation": getattr(arguments, "seed_meth", "0"),
        },
        "fragments": fragments,
        "methylation": {
            "state_model": "bernoulli",
            "beta": {
                "CG": [0.5, 0.5],
                "CHG": [0.01, 0.05],
                "CHH": [0.01, 0.05],
            }
        },
        "coverage": {"kind": "uniform"},
        "sequencing": {
            "conversion_rate": 1,
            "directional": True,
            "quality": {"kind": "uniform", "phred": 40},
            "error": {"kind": "uniform", "rate": 0},
        },
        "execution": {},
        # Required by the shared schema, but catalog generation never creates
        # this directory or publishes simulation outputs.
        "output": {
            "directory": str(base),
            "prefix": "rrbs_build",
            "format": "fastq",
        },
    }
    if arguments.vcf is not None:
        inputs = document["inputs"]
        if not isinstance(inputs, dict):
            raise CommandLineError("internal RRBS input projection failed")
        inputs["vcf"] = str(arguments.vcf.expanduser())
    return document


def build_variant_document(
    arguments: argparse.Namespace,
    base_directory: Path,
) -> dict[str, object]:
    """Project ``build variants`` arguments into an in-memory core config."""
    if (
        arguments.command != "build"
        or arguments.build_target != "variants"
    ):
        raise CommandLineError("build variants arguments are required")

    base = base_directory.expanduser().resolve(strict=False)
    mutation_rate = _resolve_mutation_rate(arguments, default=0.001)
    document: dict[str, object] = {
        "reference": str(arguments.reference.expanduser()),
        "inputs": {},
        "technology": "WGBS",
        "reads": {"count": 1},
        "mutation": {
            "rate": mutation_rate,
            "indel_fraction": arguments.indel_fraction,
            "indel_extension_probability": (
                arguments.indel_extension_probability
            ),
            "homozygous_only": arguments.homozygous_only,
        },
        "seeds": {
            "mutation": arguments.seed_mut,
            "phasing": arguments.seed_phase,
            "methylation": getattr(arguments, "seed_meth", "0"),
        },
        # Catalog-only generation never samples fragments, but the shared
        # normalized contract deliberately retains one complete shape.
        "fragments": {
            "paired_end": False,
            "read_length_1": 1,
            "insert_min": 1,
            "insert_mean": 1,
            "insert_max": 1,
            "insert_sd": 0,
            "max_ambiguous_fraction": 0,
        },
        "methylation": {
            "state_model": "bernoulli",
            "beta": {
                "CG": [0.5, 0.5],
                "CHG": [0.01, 0.05],
                "CHH": [0.01, 0.05],
            },
        },
        "coverage": {"kind": "uniform"},
        "sequencing": {
            "conversion_rate": 1,
            "directional": True,
            "quality": {"kind": "uniform", "phred": 40},
            "error": {"kind": "uniform", "rate": 0},
        },
        "execution": {},
        "output": {
            "directory": str(base),
            "prefix": "variant_build",
            "format": "fastq",
        },
    }
    if arguments.vcf is not None:
        inputs = document["inputs"]
        if not isinstance(inputs, dict):
            raise CommandLineError("internal variant input projection failed")
        inputs["vcf"] = str(arguments.vcf.expanduser())
    return document


def build_methdb_document(
    arguments: argparse.Namespace,
    base_directory: Path,
) -> dict[str, object]:
    """Project ``build methdb`` arguments into an in-memory core config."""
    if arguments.command != "build" or arguments.build_target != "methdb":
        raise CommandLineError("build methdb arguments are required")
    if (arguments.asm is not None or arguments.asm_bed is not None) and (
        arguments.vcf is None
    ):
        raise CommandLineError("--asm/--asm-bed requires --vcf")

    mutation_rate = _resolve_mutation_rate(arguments, default=0.001)

    base = base_directory.expanduser().resolve(strict=False)
    document: dict[str, object] = {
        "reference": str(arguments.reference.expanduser()),
        "inputs": {},
        "technology": "WGBS",
        "reads": {"count": 1},
        "mutation": {
            "rate": mutation_rate,
            "indel_fraction": arguments.indel_fraction,
            "indel_extension_probability": arguments.indel_extension_probability,
            "homozygous_only": arguments.homozygous_only,
        },
        "seeds": {
            "mutation": arguments.seed_mut,
            "phasing": arguments.seed_phase,
            "methylation": arguments.seed_meth,
        },
        "fragments": {
            "paired_end": False,
            "read_length_1": 1,
            "insert_min": 1,
            "insert_mean": 1,
            "insert_max": 1,
            "insert_sd": 0,
            "max_ambiguous_fraction": 0,
        },
        "methylation": {
            "state_model": arguments.methylation_model,
            "collect_non_cpg": not arguments.cpg_only,
            "cgmap_pool": arguments.cgmap_pool,
            "update_variant_boundaries": not arguments.no_update_variant_boundaries,
            "beta": {
                "CG": list(arguments.beta_cg or (0.5, 0.5)),
                "CHG": list(arguments.beta_chg or (0.01, 0.05)),
                "CHH": list(arguments.beta_chh or (0.01, 0.05)),
            },
        },
        "coverage": {"kind": "uniform"},
        "sequencing": {
            "conversion_rate": 1,
            "directional": True,
            "quality": {"kind": "uniform", "phred": 40},
            "error": {"kind": "uniform", "rate": 0},
        },
        "execution": {},
        "output": {
            "directory": str(base),
            "prefix": "methdb_build",
            "format": "fastq",
        },
    }
    inputs = document["inputs"]
    if not isinstance(inputs, dict):
        raise CommandLineError("internal MethDB input projection failed")
    for name in ("vcf", "cgmap", "bed_methyl", "asm", "asm_bed"):
        value = getattr(arguments, name)
        if value is not None:
            inputs[name] = str(value.expanduser())
    return document


def build_run_document(
    arguments: argparse.Namespace,
    base_directory: Path,
) -> dict[str, object]:
    """Project direct CLI arguments into the sole normalized config contract."""
    if arguments.command != "run" or arguments.run_technology not in (
        "wgbs",
        "rrbs",
        "tbs",
        "wgs",
        "wes",
        "ts",
    ):
        raise CommandLineError("direct run arguments are required")
    base = base_directory.expanduser().resolve(strict=False)
    technology = arguments.run_technology.upper()
    bisulfite = technology in _BISULFITE_TECHNOLOGIES
    asm = getattr(arguments, "asm", None)
    asm_bed = getattr(arguments, "asm_bed", None)
    methdb = getattr(arguments, "methdb", None)
    cgmap_pool = getattr(arguments, "cgmap_pool", False)
    if (asm is not None or asm_bed is not None) and (
        arguments.vcf is None
    ):
        raise CommandLineError("--asm/--asm-bed requires --vcf")
    if methdb is not None and (
        asm is not None
        or asm_bed is not None
        or cgmap_pool
    ):
        raise CommandLineError(
            "--methdb cannot be combined with ASM or --cgmap-pool"
        )
    if methdb is not None and arguments.vcf is not None:
        raise CommandLineError(
            "--methdb embeds variants and cannot be combined with --vcf"
        )
    if (
        methdb is not None
        and arguments.mutation_rate is not None
        and arguments.mutation_rate != 0
    ):
        raise CommandLineError(
            "--methdb embeds variants and forbids de novo mutations"
        )
    default_mutation_rate = (
        0.0
        if methdb is not None
        or getattr(arguments, "rrbs_candidates", None) is not None
        or getattr(arguments, "gc_profile", None) is not None
        else 0.001
    )
    mutation_rate = _resolve_mutation_rate(
        arguments, default=default_mutation_rate
    )
    fragment_realization = getattr(arguments, "fragment_realization", False)
    if arguments.format != "bam" and (
        arguments.fragment_summary or fragment_realization
    ):
        raise CommandLineError(
            "--fragment-summary/--fragment-realization requires --format bam"
        )

    document: dict[str, object] = {
        "reference": str(arguments.reference.expanduser()),
        "inputs": {},
        "technology": technology,
        "reads": _read_quantity(arguments),
        "mutation": {
            "rate": mutation_rate,
            "indel_fraction": arguments.indel_fraction,
            "indel_extension_probability": arguments.indel_extension_probability,
            "homozygous_only": arguments.homozygous_only,
        },
        "seeds": {
            "mutation": arguments.seed_mut,
            "phasing": arguments.seed_phase,
            "methylation": getattr(arguments, "seed_meth", "0"),
        },
        "fragments": {
            "paired_end": not arguments.single_end,
            "read_length_1": arguments.read_length,
            **_insert_parameters(
                arguments,
                collapse_fixed_bounds=(
                    technology in _WHOLE_GENOME_TECHNOLOGIES
                    or technology in _TARGETED_TECHNOLOGIES
                ),
            ),
            "max_ambiguous_fraction": arguments.max_ambiguous_fraction,
        },
        "methylation": {
            "state_model": getattr(arguments, "methylation_model", "bernoulli"),
            "collect_non_cpg": not getattr(arguments, "cpg_only", False),
            "cgmap_pool": cgmap_pool,
            "update_variant_boundaries": not getattr(
                arguments, "no_update_variant_boundaries", False
            ),
            "beta": {
                "CG": list(getattr(arguments, "beta_cg", None) or (0.5, 0.5)),
                "CHG": list(getattr(arguments, "beta_chg", None) or (0.01, 0.05)),
                "CHH": list(getattr(arguments, "beta_chh", None) or (0.01, 0.05)),
            },
        },
        "coverage": {"kind": "uniform"},
        "sequencing": {
            "conversion_rate": (
                arguments.conversion_rate if bisulfite else 0.0
            ),
            "directional": (
                not arguments.undirectional if bisulfite else True
            ),
        },
        "execution": {},
        "output": {
            "directory": str(arguments.output_directory.expanduser()),
            "prefix": arguments.prefix,
            "format": arguments.format,
            "gzip_level": arguments.gzip_level,
            "save_methdb": (
                getattr(arguments, "save_methdb", False) or arguments.save_truth
                if bisulfite
                else False
            ),
            "save_vcf": arguments.save_vcf or arguments.save_truth,
            "fragment_summary": (
                arguments.fragment_summary or fragment_realization
            ),
            "fragment_realization": fragment_realization,
        },
    }
    if arguments.seed is not None:
        document["seed"] = arguments.seed

    inputs = document["inputs"]
    if not isinstance(inputs, dict):
        raise CommandLineError("internal CLI input projection failed")
    for name in ("vcf", "cgmap", "bed_methyl", "methdb", "asm", "asm_bed"):
        value = getattr(arguments, name, None)
        if value is not None:
            inputs[name] = str(value.expanduser())

    fragments = document["fragments"]
    if not isinstance(fragments, dict):
        raise CommandLineError("internal CLI fragment projection failed")
    if not arguments.single_end:
        fragments["read_length_2"] = arguments.read_length
    if technology == "RRBS":
        cut_sites = arguments.cut_sites
        if not cut_sites:
            raise CommandLineError("RRBS requires at least one --cut-site")
        rrbs = {"cut_sites": cut_sites}
        if arguments.rrbs_candidates is not None:
            rrbs["candidate_bed"] = str(arguments.rrbs_candidates.expanduser())
        document["rrbs"] = rrbs

    if technology in _TARGETED_TECHNOLOGIES:
        document["tbs"] = {
            "bed": str(arguments.targets.expanduser()),
            "fragment_center_stddev": arguments.fragment_center_stddev,
        }

    gc_profile = getattr(arguments, "gc_profile", None)
    if gc_profile is not None:
        variable_insert = fragments["insert_sd"] != 0
        if variable_insert and (arguments.vcf is not None or mutation_rate != 0):
            raise CommandLineError(
                "variable-insert --gc-profile does not yet support variants"
            )
        document["coverage"] = {
            "kind": "profile",
            "artifact": _artifact(gc_profile, base),
        }
    elif technology in _TARGETED_TECHNOLOGIES and arguments.sampling == "score":
        document["coverage"] = {"kind": "target-score"}
    elif technology == "RRBS" and arguments.sampling == "score":
        if arguments.rrbs_candidates is None:
            raise CommandLineError(
                "--sampling score requires --rrbs-candidates"
            )
        document["coverage"] = {"kind": "profile"}

    sequencing = document["sequencing"]
    if not isinstance(sequencing, dict):
        raise CommandLineError("internal CLI sequencing projection failed")
    if arguments.quality_model is not None:
        if arguments.phred is not None:
            raise CommandLineError("--quality-model cannot be combined with --phred")
        sequencing["quality"] = {
            "kind": "markov",
            "artifact": _artifact(arguments.quality_model, base),
        }
    else:
        sequencing["quality"] = {
            "kind": "uniform",
            "phred": 40 if arguments.phred is None else arguments.phred,
        }
    if arguments.error_model is not None:
        if arguments.error_rate is not None:
            raise CommandLineError("--error-model cannot be combined with --error-rate")
        sequencing["error"] = {
            "kind": "quality-confusion",
            "artifact": _artifact(arguments.error_model, base),
        }
    else:
        sequencing["error"] = {
            "kind": "uniform",
            "rate": 0.005 if arguments.error_rate is None else arguments.error_rate,
        }

    execution = document["execution"]
    if not isinstance(execution, dict):
        raise CommandLineError("internal CLI execution projection failed")
    execution["threads"] = arguments.threads
    return document


def main(argv: Sequence[str] | None = None) -> int:
    """Run the BSReadSim command-line interface."""
    parser = build_parser()
    command_arguments = list(sys.argv[1:] if argv is None else argv)
    arguments = parser.parse_args(command_arguments)
    if arguments.command is None:
        parser.print_help()
        return 2

    if arguments.command == "export" and arguments.export_target == "test-fasta":
        from .resources import ResourceError, copy_resource

        try:
            print(copy_resource("test-fasta", arguments.output))
        except ResourceError as error:
            print("bsreadsim: error: {}".format(error), file=sys.stderr)
            return 1
        return 0

    # Keep artifact building independent from optional array/runtime modules.
    from .run.config import ConfigError
    from .run.catalog import (
        CatalogError,
        export_methdb_bed,
        build_methdb_snapshot,
        export_rrbs_catalog,
        export_variant_catalog,
    )
    from .htsim.launch import CoreArgvError
    from .run.prepare import PreparationError

    if arguments.command == "build":
        try:
            if arguments.build_target == "rrbs":
                output_path = export_rrbs_catalog(
                    build_rrbs_document(arguments, Path.cwd()),
                    arguments.output,
                    base_directory=Path.cwd(),
                    core_executable=arguments.core,
                )
            elif arguments.build_target == "variants":
                output_path = export_variant_catalog(
                    build_variant_document(arguments, Path.cwd()),
                    arguments.output,
                    base_directory=Path.cwd(),
                    core_executable=arguments.core,
                )
            else:
                output_path = build_methdb_snapshot(
                    build_methdb_document(arguments, Path.cwd()),
                    arguments.output,
                    base_directory=Path.cwd(),
                    core_executable=arguments.core,
                )
            print(output_path)
            return 0
        except (
            CatalogError,
            CommandLineError,
            ConfigError,
            CoreArgvError,
            PreparationError,
        ) as error:
            print("bsreadsim: error: {}".format(error), file=sys.stderr)
            return 1

    if arguments.command == "export":
        try:
            print(
                export_methdb_bed(
                    arguments.input,
                    arguments.output,
                    compressed=not arguments.no_compression,
                    core_executable=arguments.core,
                )
            )
            return 0
        except CatalogError as error:
            print("bsreadsim: error: {}".format(error), file=sys.stderr)
            return 1

    from .output import BamError, OutputError
    from .htsim.subprocess import CoreProcessError
    from .run.manifest import ManifestError
    from .run.execute import PipelineError, run_document
    from .process import ProcessError

    try:
        document = build_run_document(arguments, Path.cwd())
        result = run_document(
            document,
            base_directory=Path.cwd(),
            core_executable=arguments.core,
            invocation_argv=("bsreadsim", *command_arguments),
        )
    except (
        CatalogError,
        CommandLineError,
        ConfigError,
        CoreArgvError,
        CoreProcessError,
        ManifestError,
        OutputError,
        PipelineError,
        ProcessError,
        PreparationError,
        BamError,
    ) as error:
        print("bsreadsim: error: {}".format(error), file=sys.stderr)
        return 1

    print(result.manifest_path)
    return 0


__all__ = [
    "CommandLineError",
    "build_methdb_document",
    "build_parser",
    "build_rrbs_document",
    "build_run_document",
    "build_variant_document",
    "main",
]
