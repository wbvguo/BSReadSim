"""Command-line entry point for BSReadSim."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path
import sys
from collections.abc import Mapping, Sequence

from . import __version__


class CommandLineError(ValueError):
    """Direct simulation arguments cannot form one valid run."""


def _add_runtime_options(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--core",
        type=Path,
        help="override the bundled htsim-core executable",
    )


def _add_fragment_domain_arguments(parser: argparse.ArgumentParser) -> None:
    """Add fragment-shape arguments shared by simulation and catalog export."""
    fragments = parser.add_argument_group("fragments")
    fragments.add_argument(
        "--single-end",
        action="store_true",
        help="emit one read per fragment instead of paired-end reads",
    )
    fragments.add_argument(
        "--read-length", type=int, default=100, help="read length (default: 100)"
    )
    fragments.add_argument(
        "--insert-size",
        type=int,
        help=(
            "fixed insert size; cannot be combined with insert range options"
        ),
    )
    fragments.add_argument("--insert-min", type=int)
    fragments.add_argument("--insert-mean", type=int)
    fragments.add_argument("--insert-max", type=int)
    fragments.add_argument("--insert-stddev", type=float)
    fragments.add_argument(
        "--max-ambiguous-fraction",
        type=float,
        default=0.05,
        help="maximum N fraction in each emitted mate (default: 0.05)",
    )


def _add_variant_domain_arguments(parser: argparse.ArgumentParser) -> None:
    """Add de novo mutation arguments shared by run and RRBS export."""
    mutation = parser.add_argument_group("variants")
    mutation.add_argument(
        "--mutation-rate",
        type=float,
        help="de novo mutation rate (run default: 0.001; catalog default: 0)",
    )
    mutation.add_argument("--indel-fraction", type=float, default=0.15)
    mutation.add_argument(
        "--indel-extension-probability", type=float, default=0.15
    )
    mutation.add_argument("--homozygous-only", action="store_true")


def _add_direct_run_arguments(parser: argparse.ArgumentParser) -> None:
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
        "-n", "--read-pairs", type=int, help="number of fragments/read pairs"
    )
    count.add_argument(
        "-d",
        "--depth",
        type=float,
        help="requested mean depth over the technology target region",
    )

    biology = parser.add_argument_group("biological inputs and technology")
    biology.add_argument(
        "--technology",
        choices=("WGBS", "RRBS", "TBS"),
        default="WGBS",
    )
    biology.add_argument(
        "--seed",
        help="simulation-only unsigned 64-bit seed; omit to generate one",
    )
    biology.add_argument(
        "--methdb-seed",
        default="0",
        help="catalog-only unsigned 64-bit seed (default: 0)",
    )
    biology.add_argument("--vcf", type=Path, help="phased VCF input")
    methylation_input = biology.add_mutually_exclusive_group()
    methylation_input.add_argument(
        "--cgmap", type=Path, help="CGmap methylation input"
    )
    methylation_input.add_argument(
        "--bed-methyl",
        type=Path,
        help="UCSC/ENCODE bedMethyl (BED9+2 or BED9+9) input",
    )
    methylation_input.add_argument(
        "--methdb",
        type=Path,
        help="fixed normalized BSReadSim MethDB snapshot",
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
    biology.add_argument(
        "--cut-site",
        action="append",
        dest="cut_sites",
        help="RRBS cut motif, for example C|CGG; repeat for multiple enzymes",
    )
    biology.add_argument(
        "--rrbs-candidates",
        type=Path,
        help=(
            "RRBS candidate BED exported by 'bsreadsim catalog rrbs'; "
            "scores are ignored unless --rrbs-score is selected"
        ),
    )
    biology.add_argument("--targets", type=Path, help="TBS BED6 targets")
    biology.add_argument(
        "--fragment-center-stddev",
        type=float,
        default=50.0,
        help="TBS fragment-center standard deviation (default: 50)",
    )

    coverage = parser.add_argument_group("coverage")
    coverage_mode = coverage.add_mutually_exclusive_group()
    coverage_mode.add_argument(
        "--coverage-profile",
        type=Path,
        help=(
            "WGBS target GC distribution with one probability per line; "
            "identity is hashed automatically"
        ),
    )
    coverage_mode.add_argument(
        "--target-score",
        action="store_true",
        help="sample TBS targets by exact BED6 output weights",
    )
    coverage_mode.add_argument(
        "--rrbs-score",
        action="store_true",
        help="sample RRBS fragments by candidate BED score weights",
    )

    _add_fragment_domain_arguments(parser)
    _add_variant_domain_arguments(parser)

    methylation = parser.add_argument_group("methylation")
    methylation.add_argument(
        "--beta-cg", type=float, nargs=2, metavar=("ALPHA", "BETA")
    )
    methylation.add_argument(
        "--beta-chg", type=float, nargs=2, metavar=("ALPHA", "BETA")
    )
    methylation.add_argument(
        "--beta-chh", type=float, nargs=2, metavar=("ALPHA", "BETA")
    )
    methylation.add_argument(
        "--cpg-only",
        action="store_true",
        help="omit CHG and CHH sites from the generated methylation catalog",
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

    sequencing = parser.add_argument_group("sequencing")
    sequencing.add_argument(
        "--conversion-rate", type=float, default=0.998
    )
    sequencing.add_argument(
        "--undirectional",
        action="store_true",
        help="sample both directional bisulfite library orientations",
    )
    sequencing.add_argument(
        "--phred", type=int, help="uniform Phred score (default: 40)"
    )
    sequencing.add_argument(
        "--quality-model",
        type=Path,
        help="quality Markov JSON; identity is hashed automatically",
    )
    sequencing.add_argument(
        "--error-rate",
        type=float,
        help="uniform substitution rate (default: 0.005)",
    )
    sequencing.add_argument(
        "--error-model",
        type=Path,
        help="quality-confusion JSON; identity is hashed automatically",
    )

    execution = parser.add_argument_group("execution")
    execution.add_argument("--workers", type=int)
    execution.add_argument("--core-workers", type=int)
    execution.add_argument("--chunk-size", type=int)
    execution.add_argument("--max-in-flight-fragments", type=int)

    output = parser.add_argument_group("output")
    output.add_argument("-p", "--prefix", default="sim")
    output.add_argument(
        "--compression", choices=("none", "gzip"), default="gzip"
    )
    output.add_argument("--gzip-level", type=int, default=6)
    output.add_argument(
        "--bam",
        action="store_true",
        help=(
            "also emit an unsorted annotated BAM; Full Details projection "
            "is retained internally"
        ),
    )
    output.add_argument(
        "--fragment-summary",
        action="store_true",
        help=(
            "add optional full-fragment zf summaries to BAM records; "
            "zt and zr are always present"
        ),
    )
    output.add_argument(
        "--fragment-realization",
        action="store_true",
        help="emit complete-fragment methylation/conversion state in BAM zx",
    )
    output.add_argument(
        "--save-methdb",
        type=Path,
        help="save the fixed normalized MethDB catalog used by the run",
    )


def _add_rrbs_catalog_arguments(parser: argparse.ArgumentParser) -> None:
    required = parser.add_argument_group("required catalog inputs")
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

    biology = parser.add_argument_group("biological inputs")
    biology.add_argument(
        "--methdb-seed",
        dest="catalog_seed",
        default="0",
        help=(
            "catalog-only unsigned 64-bit decimal seed (default: 0)"
        ),
    )
    biology.add_argument("--vcf", type=Path, help="phased VCF input")

    _add_fragment_domain_arguments(parser)
    _add_variant_domain_arguments(parser)
    parser.add_argument(
        "--core", type=Path, help="override the bundled htsim-core executable"
    )


def build_parser() -> argparse.ArgumentParser:
    """Build the public command-line parser without importing heavy runtime modules."""
    parser = argparse.ArgumentParser(
        prog="bsreadsim",
        description=(
            "Generate biological fragments with htsim and apply BSReadSim "
            "bisulfite conversion and sequencing effects."
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
    _add_direct_run_arguments(run_parser)
    _add_runtime_options(run_parser)

    catalog_parser = commands.add_parser(
        "catalog",
        help="export an htsim fragment-candidate catalog",
    )
    catalog_commands = catalog_parser.add_subparsers(
        dest="catalog_technology", required=True
    )
    rrbs_catalog_parser = catalog_commands.add_parser(
        "rrbs",
        help="export exact RRBS candidates directly from CLI arguments",
        description=(
            "Export the htsim RRBS candidate BED directly; no JSON "
            "configuration file is read or written."
        ),
    )
    _add_rrbs_catalog_arguments(rrbs_catalog_parser)
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
    *,
    format_name: str,
    version: str,
) -> dict[str, str]:
    return {
        "path": str(path.expanduser()),
        "format": format_name,
        "version": version,
        "sha256": _artifact_sha256(path, base_directory),
    }


def _insert_parameters(arguments: argparse.Namespace) -> Mapping[str, object]:
    ranged = (
        arguments.insert_min,
        arguments.insert_mean,
        arguments.insert_max,
        arguments.insert_stddev,
    )
    if arguments.insert_size is not None:
        if any(value is not None for value in ranged):
            raise CommandLineError(
                "--insert-size cannot be combined with insert range options"
            )
        return {
            "insert_min": arguments.insert_size,
            "insert_mean": arguments.insert_size,
            "insert_max": arguments.insert_size,
            "insert_stddev": 0,
        }
    return {
        "insert_min": 100 if arguments.insert_min is None else arguments.insert_min,
        "insert_mean": 400 if arguments.insert_mean is None else arguments.insert_mean,
        "insert_max": 1000 if arguments.insert_max is None else arguments.insert_max,
        "insert_stddev": (
            25 if arguments.insert_stddev is None else arguments.insert_stddev
        ),
    }


def build_rrbs_catalog_document(
    arguments: argparse.Namespace,
    base_directory: Path,
) -> dict[str, object]:
    """Project direct RRBS catalog arguments into an in-memory core config."""
    if (
        arguments.command != "catalog"
        or arguments.catalog_technology != "rrbs"
    ):
        raise CommandLineError("direct RRBS catalog arguments are required")

    base = base_directory.expanduser().resolve(strict=False)
    mutation_rate = (
        0 if arguments.mutation_rate is None else arguments.mutation_rate
    )
    if arguments.vcf is not None and mutation_rate != 0:
        raise CommandLineError("--vcf requires --mutation-rate 0")
    fragments: dict[str, object] = {
        "paired_end": not arguments.single_end,
        "read_length_1": arguments.read_length,
        **_insert_parameters(arguments),
        # The htsim config contract requires a count even though catalog-only
        # generation never samples or emits protocol fragments.
        "read_pairs": 1,
        "max_ambiguous_fraction": arguments.max_ambiguous_fraction,
    }
    if not arguments.single_end:
        fragments["read_length_2"] = arguments.read_length

    document: dict[str, object] = {
        "schema_version": "1.1",
        "reference": str(arguments.reference.expanduser()),
        "inputs": {},
        "technology": "RRBS",
        "rrbs": {"cut_sites": list(arguments.cut_sites)},
        "mutation": {
            "rate": mutation_rate,
            "indel_fraction": arguments.indel_fraction,
            "indel_extension_probability": (
                arguments.indel_extension_probability
            ),
            "homozygous_only": arguments.homozygous_only,
        },
        "fragments": fragments,
        "methylation": {
            "catalog_seed": arguments.catalog_seed,
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
            "prefix": "rrbs_catalog",
            "compression": "none",
        },
    }
    if arguments.vcf is not None:
        inputs = document["inputs"]
        if not isinstance(inputs, dict):
            raise CommandLineError("internal catalog input projection failed")
        inputs["vcf"] = str(arguments.vcf.expanduser())
    return document


def build_run_document(
    arguments: argparse.Namespace,
    base_directory: Path,
) -> dict[str, object]:
    """Project direct CLI arguments into the sole normalized config contract."""
    if arguments.command != "run":
        raise CommandLineError("direct run arguments are required")
    base = base_directory.expanduser().resolve(strict=False)
    technology = arguments.technology
    if (arguments.asm is not None or arguments.asm_bed is not None) and (
        arguments.vcf is None
    ):
        raise CommandLineError("--asm/--asm-bed requires --vcf")
    if arguments.methdb is not None and (
        arguments.asm is not None
        or arguments.asm_bed is not None
        or arguments.cgmap_pool
    ):
        raise CommandLineError(
            "--methdb cannot be combined with ASM or --cgmap-pool"
        )
    if arguments.methdb is not None and arguments.save_methdb is not None:
        raise CommandLineError("--methdb and --save-methdb are mutually exclusive")

    mutation_rate = arguments.mutation_rate
    if mutation_rate is None:
        mutation_rate = (
            0
            if arguments.vcf is not None
            or arguments.rrbs_candidates is not None
            or arguments.coverage_profile is not None
            else 0.001
        )
    if arguments.vcf is not None and mutation_rate != 0:
        raise CommandLineError("--vcf requires --mutation-rate 0")

    document: dict[str, object] = {
        "schema_version": "1.1",
        "reference": str(arguments.reference.expanduser()),
        "inputs": {},
        "technology": technology,
        "mutation": {
            "rate": mutation_rate,
            "indel_fraction": arguments.indel_fraction,
            "indel_extension_probability": arguments.indel_extension_probability,
            "homozygous_only": arguments.homozygous_only,
        },
        "fragments": {
            "paired_end": not arguments.single_end,
            "read_length_1": arguments.read_length,
            **_insert_parameters(arguments),
            "max_ambiguous_fraction": arguments.max_ambiguous_fraction,
        },
        "methylation": {
            "catalog_seed": arguments.methdb_seed,
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
            "conversion_rate": arguments.conversion_rate,
            "directional": not arguments.undirectional,
        },
        "execution": {},
        "output": {
            "directory": str(arguments.output_directory.expanduser()),
            "prefix": arguments.prefix,
            "compression": arguments.compression,
            "gzip_level": arguments.gzip_level,
            "bam": arguments.bam or arguments.fragment_realization,
            "fragment_summary": (
                arguments.fragment_summary or arguments.fragment_realization
            ),
            "fragment_realization": arguments.fragment_realization,
        },
    }
    if arguments.seed is not None:
        document["seed"] = arguments.seed

    inputs = document["inputs"]
    if not isinstance(inputs, dict):
        raise CommandLineError("internal CLI input projection failed")
    for name in ("vcf", "cgmap", "bed_methyl", "methdb", "asm", "asm_bed"):
        value = getattr(arguments, name)
        if value is not None:
            inputs[name] = str(value.expanduser())

    fragments = document["fragments"]
    if not isinstance(fragments, dict):
        raise CommandLineError("internal CLI fragment projection failed")
    if not arguments.single_end:
        fragments["read_length_2"] = arguments.read_length
    if arguments.read_pairs is not None:
        fragments["read_pairs"] = arguments.read_pairs
    else:
        fragments["depth"] = arguments.depth

    cut_sites = arguments.cut_sites or []
    if technology == "RRBS":
        if not cut_sites:
            raise CommandLineError("RRBS requires at least one --cut-site")
        if arguments.targets is not None:
            raise CommandLineError("RRBS forbids --targets")
        rrbs = {"cut_sites": cut_sites}
        if arguments.rrbs_candidates is not None:
            rrbs["candidate_bed"] = str(arguments.rrbs_candidates.expanduser())
        document["rrbs"] = rrbs
    elif cut_sites:
        raise CommandLineError("--cut-site requires --technology RRBS")
    elif arguments.rrbs_candidates is not None:
        raise CommandLineError("--rrbs-candidates requires --technology RRBS")

    if technology == "TBS":
        if arguments.targets is None:
            raise CommandLineError("TBS requires --targets")
        document["tbs"] = {
            "bed": str(arguments.targets.expanduser()),
            "fragment_center_stddev": arguments.fragment_center_stddev,
        }
    elif arguments.targets is not None:
        raise CommandLineError("--targets requires --technology TBS")

    if arguments.coverage_profile is not None:
        if technology != "WGBS":
            raise CommandLineError("--coverage-profile supports WGBS only")
        variable_insert = not (
            fragments["insert_min"]
            == fragments["insert_mean"]
            == fragments["insert_max"]
            and fragments["insert_stddev"] == 0
        )
        if variable_insert and (arguments.vcf is not None or mutation_rate != 0):
            raise CommandLineError(
                "variable-insert --coverage-profile does not yet support variants"
            )
        from .run.config import WGBS_GC_PROFILE_FORMAT, WGBS_GC_PROFILE_VERSION

        document["coverage"] = {
            "kind": "profile",
            "artifact": _artifact(
                arguments.coverage_profile,
                base,
                format_name=WGBS_GC_PROFILE_FORMAT,
                version=WGBS_GC_PROFILE_VERSION,
            ),
        }
    elif arguments.target_score:
        if technology != "TBS":
            raise CommandLineError("--target-score supports TBS only")
        document["coverage"] = {"kind": "target-score"}
    elif arguments.rrbs_score:
        if technology != "RRBS":
            raise CommandLineError("--rrbs-score supports RRBS only")
        if arguments.rrbs_candidates is None:
            raise CommandLineError(
                "--rrbs-score requires --rrbs-candidates"
            )
        document["coverage"] = {"kind": "profile"}

    sequencing = document["sequencing"]
    if not isinstance(sequencing, dict):
        raise CommandLineError("internal CLI sequencing projection failed")
    if arguments.quality_model is not None:
        if arguments.phred is not None:
            raise CommandLineError("--quality-model cannot be combined with --phred")
        from .process.sequencing import QUALITY_MARKOV_FORMAT, QUALITY_MARKOV_VERSION

        sequencing["quality"] = {
            "kind": "markov",
            "artifact": _artifact(
                arguments.quality_model,
                base,
                format_name=QUALITY_MARKOV_FORMAT,
                version=QUALITY_MARKOV_VERSION,
            ),
        }
    else:
        sequencing["quality"] = {
            "kind": "uniform",
            "phred": 40 if arguments.phred is None else arguments.phred,
        }
    if arguments.error_model is not None:
        if arguments.error_rate is not None:
            raise CommandLineError("--error-model cannot be combined with --error-rate")
        from .process.sequencing import (
            QUALITY_CONFUSION_FORMAT,
            QUALITY_CONFUSION_VERSION,
        )

        sequencing["error"] = {
            "kind": "quality-confusion",
            "artifact": _artifact(
                arguments.error_model,
                base,
                format_name=QUALITY_CONFUSION_FORMAT,
                version=QUALITY_CONFUSION_VERSION,
            ),
        }
    else:
        sequencing["error"] = {
            "kind": "uniform",
            "rate": 0.005 if arguments.error_rate is None else arguments.error_rate,
        }

    execution = document["execution"]
    if not isinstance(execution, dict):
        raise CommandLineError("internal CLI execution projection failed")
    for name in (
        "workers",
        "core_workers",
        "chunk_size",
        "max_in_flight_fragments",
    ):
        value = getattr(arguments, name)
        if value is not None:
            execution[name] = value
    return document


def main(argv: Sequence[str] | None = None) -> int:
    """Run the BSReadSim command-line interface."""
    parser = build_parser()
    arguments = parser.parse_args(argv)
    if arguments.command is None:
        parser.print_help()
        return 2

    # Keep metadata-only commands independent from optional runtime
    # dependencies and the simulation modules.
    from .run.config import ConfigError
    from .output import BamError
    from .run.catalog import (
        CatalogError,
        export_methdb_catalog,
        export_rrbs_catalog,
    )
    from .htsim.launch import CoreArgvError
    from .htsim.subprocess import CoreProcessError
    from .run.manifest import ManifestError
    from .output import OutputError
    from .run.execute import PipelineError, run_document
    from .process import ProcessError
    from .run.prepare import PreparationError

    try:
        if arguments.command == "run":
            document = build_run_document(arguments, Path.cwd())
            if arguments.save_methdb is not None:
                methdb_path = export_methdb_catalog(
                    document,
                    arguments.save_methdb,
                    base_directory=Path.cwd(),
                    core_executable=arguments.core,
                )
                inputs = document["inputs"]
                methylation = document["methylation"]
                if not isinstance(inputs, dict) or not isinstance(methylation, dict):
                    raise CommandLineError("internal MethDB projection failed")
                for name in ("cgmap", "bed_methyl", "asm", "asm_bed"):
                    inputs.pop(name, None)
                inputs["methdb"] = str(methdb_path)
                methylation["cgmap_pool"] = False
            result = run_document(
                document,
                base_directory=Path.cwd(),
                core_executable=arguments.core,
            )
        else:
            output_path = export_rrbs_catalog(
                build_rrbs_catalog_document(arguments, Path.cwd()),
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
    "build_parser",
    "build_rrbs_catalog_document",
    "build_run_document",
    "main",
]
