"""Canonical CLI reconstruction for a normalized simulation configuration."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any


_BISULFITE_TECHNOLOGIES = frozenset(("WGBS", "RRBS", "TBS"))
_WHOLE_GENOME_TECHNOLOGIES = frozenset(("WGBS", "WGS"))
_TARGETED_TECHNOLOGIES = frozenset(("TBS", "WES", "TS"))
_TECHNOLOGIES = (
    _BISULFITE_TECHNOLOGIES
    | _WHOLE_GENOME_TECHNOLOGIES
    | _TARGETED_TECHNOLOGIES
)


class FullCommandError(ValueError):
    """An effective run configuration cannot be represented by the run CLI."""


def build_full_run_argv(
    config: Mapping[str, Any], invocation_argv: Sequence[str]
) -> list[str]:
    """Return a canonical run command containing every effective CLI setting."""
    if not isinstance(config, Mapping):
        raise FullCommandError("normalized config must be an object")
    if isinstance(invocation_argv, (str, bytes)):
        raise FullCommandError("invocation argv must be a sequence of strings")
    received = list(invocation_argv)
    if not received or any(not isinstance(token, str) for token in received):
        raise FullCommandError("invocation argv must contain strings")

    technology = _string(config, "technology")
    if technology not in _TECHNOLOGIES:
        raise FullCommandError("normalized technology is not a released run command")
    bisulfite = technology in _BISULFITE_TECHNOLOGIES
    reads = _mapping(config, "reads")
    fragments = _mapping(config, "fragments")
    inputs = _mapping(config, "inputs")
    mutation = _mapping(config, "mutation")
    seeds = _mapping(config, "seeds")
    coverage = _mapping(config, "coverage")
    sequencing = _mapping(config, "sequencing")
    execution = _mapping(config, "execution")
    output = _mapping(config, "output")

    argv = [received[0], "run", technology.lower()]
    _option(argv, "--reference", _string(config, "reference"))
    _option(argv, "--output", _string(output, "directory"))
    if "count" in reads:
        _option(argv, "--reads", reads["count"])
    elif "depth" in reads:
        _option(argv, "--depth", reads["depth"])
    else:
        raise FullCommandError("normalized reads require count or depth")

    _option(argv, "--seed", _string(config, "seed"))
    _option(argv, "--seed-mut", _string(seeds, "mutation"))
    _option(argv, "--seed-phase", _string(seeds, "phasing"))
    if bisulfite:
        _option(argv, "--seed-meth", _string(seeds, "methylation"))
        _methylation_input_options(argv, inputs)
    elif _string(seeds, "methylation") != "0":
        raise FullCommandError(
            "standard sequencing cannot represent a nonzero methylation seed"
        )

    _sampling_options(argv, config, coverage, technology)
    _fragment_options(argv, fragments)
    _variant_options(argv, inputs, mutation)
    if bisulfite:
        _methylation_options(argv, _mapping(config, "methylation"))
        _option(argv, "--conversion-rate", sequencing["conversion_rate"])
        if not _boolean(sequencing, "directional"):
            argv.append("--undirectional")
    _sequencing_options(argv, sequencing)
    _execution_options(argv, execution)
    _output_options(argv, output, bisulfite=bisulfite)

    core = _received_option(received, "--core")
    if core is not None:
        _option(argv, "--core", core)
    return argv


def _methylation_input_options(
    argv: list[str], inputs: Mapping[str, Any]
) -> None:
    for name, option in (
        ("cgmap", "--cgmap"),
        ("bed_methyl", "--bedmethyl"),
        ("methbg", "--methbg"),
        ("methbed", "--methbed"),
        ("methdb", "--methdb"),
        ("asm", "--asm"),
        ("asm_bed", "--asm-bed"),
    ):
        if name in inputs:
            _option(argv, option, inputs[name])


def _sampling_options(
    argv: list[str],
    config: Mapping[str, Any],
    coverage: Mapping[str, Any],
    technology: str,
) -> None:
    kind = _string(coverage, "kind")
    if technology in _WHOLE_GENOME_TECHNOLOGIES:
        if kind == "uniform":
            _option(argv, "--sampling", "uniform")
            return
        if kind != "profile" or "artifact" not in coverage:
            raise FullCommandError(
                "whole-genome coverage cannot be represented by the run CLI"
            )
        artifact = _mapping(coverage, "artifact")
        _option(argv, "--sampling", "gc")
        _option(argv, "--gc-profile", _string(artifact, "path"))
        return

    if technology == "RRBS":
        if kind == "uniform":
            sampling = "uniform"
        elif kind == "profile" and "artifact" not in coverage:
            sampling = "score"
        else:
            raise FullCommandError(
                "normalized RRBS coverage cannot be represented by the run CLI"
            )
        _option(argv, "--sampling", sampling)
        rrbs = _mapping(config, "rrbs")
        cut_sites = rrbs.get("cut_sites")
        if (
            not isinstance(cut_sites, list)
            or not cut_sites
            or any(not isinstance(site, str) for site in cut_sites)
        ):
            raise FullCommandError("normalized RRBS cut sites are invalid")
        _option(argv, "--cut-site", ",".join(cut_sites))
        if "candidate_bed" in rrbs:
            _option(argv, "--rrbs-candidates", rrbs["candidate_bed"])
        return

    if technology in _TARGETED_TECHNOLOGIES:
        if kind == "uniform":
            sampling = "uniform"
        elif kind == "target-score":
            sampling = "score"
        else:
            raise FullCommandError(
                "normalized target coverage cannot be represented by the run CLI"
            )
        _option(argv, "--sampling", sampling)
        tbs = _mapping(config, "tbs")
        _option(argv, "--targets", _string(tbs, "bed"))
        _option(
            argv,
            "--center-sd",
            tbs["center_sd"],
        )
        return

    raise FullCommandError("normalized technology has no sampling projection")


def _fragment_options(argv: list[str], fragments: Mapping[str, Any]) -> None:
    paired_end = _boolean(fragments, "paired_end")
    read_length = fragments["read_length_1"]
    if paired_end:
        if fragments.get("read_length_2") != read_length:
            raise FullCommandError(
                "the run CLI cannot represent different R1 and R2 lengths"
            )
    else:
        argv.append("--single-end")
    _option(argv, "--read-length", read_length)
    _option(argv, "--insert-min", fragments["insert_min"])
    _option(argv, "--insert-mean", fragments["insert_mean"])
    _option(argv, "--insert-max", fragments["insert_max"])
    _option(argv, "--insert-sd", fragments["insert_sd"])
    _option(
        argv,
        "--max-ambiguous-fraction",
        fragments["max_ambiguous_fraction"],
    )


def _variant_options(
    argv: list[str],
    inputs: Mapping[str, Any],
    mutation: Mapping[str, Any],
) -> None:
    if "vcf" in inputs:
        _option(argv, "--vcf", inputs["vcf"])
    else:
        _option(argv, "--mutation-rate", mutation["rate"])
    _option(argv, "--indel-fraction", mutation["indel_fraction"])
    _option(
        argv,
        "--indel-extension-probability",
        mutation["indel_extension_probability"],
    )
    if _boolean(mutation, "homozygous_only"):
        argv.append("--homozygous-only")


def _methylation_options(
    argv: list[str], methylation: Mapping[str, Any]
) -> None:
    beta = _mapping(methylation, "beta")
    for context, option in (
        ("CG", "--beta-cg"),
        ("CHG", "--beta-chg"),
        ("CHH", "--beta-chh"),
    ):
        pair = beta.get(context)
        if not isinstance(pair, list) or len(pair) != 2:
            raise FullCommandError("normalized beta parameters are invalid")
        _option(argv, option, "{},{}".format(*map(_value, pair)))
    if not _boolean(methylation, "collect_non_cpg"):
        argv.append("--cpg-only")
    if _boolean(methylation, "cgmap_pool"):
        argv.append("--pool-meth")
    _option(argv, "--meth-model", methylation["state_model"])
    if not _boolean(methylation, "update_variant_boundaries"):
        argv.append("--no-update-variant-boundaries")


def _sequencing_options(
    argv: list[str], sequencing: Mapping[str, Any]
) -> None:
    quality = _mapping(sequencing, "quality")
    quality_kind = _string(quality, "kind")
    if quality_kind == "uniform":
        _option(argv, "--phred", quality["phred"])
    elif quality_kind == "markov":
        artifact = _mapping(quality, "artifact")
        _option(argv, "--quality-model", _string(artifact, "path"))
    else:
        raise FullCommandError("normalized quality model has no CLI projection")

    error = _mapping(sequencing, "error")
    error_kind = _string(error, "kind")
    if error_kind == "uniform":
        _option(argv, "--error-rate", error["rate"])
    elif error_kind == "quality-confusion":
        artifact = _mapping(error, "artifact")
        _option(argv, "--error-model", _string(artifact, "path"))
    else:
        raise FullCommandError("normalized error model has no CLI projection")


def _execution_options(argv: list[str], execution: Mapping[str, Any]) -> None:
    _option(argv, "--threads", execution["threads"])


def _output_options(
    argv: list[str], output: Mapping[str, Any], *, bisulfite: bool
) -> None:
    _option(argv, "--prefix", output["prefix"])
    _option(argv, "--format", output["format"])
    _option(argv, "--gzip-level", output["gzip_level"])
    if _boolean(output, "fragment_summary"):
        argv.append("--fragment-summary")
    if _boolean(output, "fragment_realization"):
        if not bisulfite:
            raise FullCommandError(
                "standard sequencing cannot represent fragment realization"
            )
        argv.append("--fragment-realization")
    if _boolean(output, "save_methdb"):
        if not bisulfite:
            raise FullCommandError(
                "standard sequencing cannot represent MethDB truth output"
            )
        argv.append("--save-methdb")
    if _boolean(output, "save_vcf"):
        argv.append("--save-vcf")


def _received_option(argv: Sequence[str], option: str) -> str | None:
    value = None
    index = 1
    while index < len(argv):
        token = argv[index]
        if token == "--":
            break
        if token == option:
            if index + 1 >= len(argv):
                raise FullCommandError("{} is missing its value".format(option))
            value = argv[index + 1]
            index += 2
            continue
        prefix = option + "="
        if token.startswith(prefix):
            value = token[len(prefix) :]
        index += 1
    return value


def _option(argv: list[str], name: str, value: object) -> None:
    argv.extend((name, _value(value)))


def _value(value: object) -> str:
    if isinstance(value, bool) or not isinstance(value, (str, int, float)):
        raise FullCommandError("normalized CLI value has an invalid type")
    return str(value)


def _mapping(value: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    result = value.get(name)
    if not isinstance(result, Mapping):
        raise FullCommandError("normalized {} section is invalid".format(name))
    return result


def _string(value: Mapping[str, Any], name: str) -> str:
    result = value.get(name)
    if not isinstance(result, str) or not result:
        raise FullCommandError("normalized {} value is invalid".format(name))
    return result


def _boolean(value: Mapping[str, Any], name: str) -> bool:
    result = value.get(name)
    if not isinstance(result, bool):
        raise FullCommandError("normalized {} value is invalid".format(name))
    return result


__all__ = ["FullCommandError", "build_full_run_argv"]
