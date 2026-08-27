#include "core.h"

#include <cstring>

#include <array>
#include <charconv>
#include <cmath>
#include <limits>
#include <string_view>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <algorithm>
#include <cfenv>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "wgbs.h"
#include "fragment.h"
#include "methdb.h"
#include "variant.h"
#include "types.h"
#include "protocol.h"
#include "reference.h"
#include "utilities.h"
#include "rrbs.h"
#include "tbs.h"

// ---- config --------------------------------------------------------

namespace htsim::core {
namespace {

using Values = std::unordered_map<std::string, std::vector<std::string>>;

constexpr std::array<std::string_view, 47> known_options = {{
    "--emit-details",
    "--protocol-batch-fragments",
    "--run-id",
    "--config-sha256",
    "--seed",
    "--seed-mut",
    "--seed-phase",
    "--seed-meth",
    "--reference",
    "--vcf",
    "--cgmap",
    "--bed-methyl",
    "--methdb",
    "--asm",
    "--asm-bed",
    "--technology",
    "--directional",
    "--paired-end",
    "--read-length-1",
    "--read-length-2",
    "--insert-min",
    "--insert-mean",
    "--insert-max",
    "--insert-sd",
    "--depth",
    "--fragments",
    "--max-ambiguous-fraction",
    "--chunk-size",
    "--core-workers",
    "--mutation-rate",
    "--indel-fraction",
    "--indel-extension-probability",
    "--homozygous-only",
    "--collect-non-cpg",
    "--cgmap-pool",
    "--update-variant-boundaries",
    "--beta-cg",
    "--beta-chg",
    "--beta-chh",
    "--coverage",
    "--coverage-profile",
    "--rrbs-cut-site",
    "--rrbs-candidate-bed",
    "--tbs-bed",
    "--tbs-center-stddev",
}};

bool is_known_option(std::string_view option)
{
    for (const std::string_view known : known_options) {
        if (option == known) {return true;}
    }
    return false;
}

Values collect_values(const std::vector<std::string> &arguments)
{
    Values values;
    for (std::size_t index = 0; index < arguments.size();) {
        const std::string &option = arguments[index];
        if (!is_known_option(option)) {
            throw CoreConfigError("unknown core option: " + option);
        }
        if (index + 1 >= arguments.size()) {
            throw CoreConfigError("missing value for core option: " + option);
        }
        const std::string &value = arguments[index + 1];
        if (option != "--rrbs-cut-site" && values.find(option) != values.end()) {
            throw CoreConfigError("duplicate core option: " + option);
        }
        values[option].push_back(value);
        index += 2;
    }
    return values;
}

const std::string &required(const Values &values, const std::string &name)
{
    const auto found = values.find(name);
    if (found == values.end()) {
        throw CoreConfigError("missing required core option: " + name);
    }
    return found->second.front();
}

std::optional<std::string> optional_text(
    const Values &values,
    const std::string &name)
{
    const auto found = values.find(name);
    if (found == values.end()) {return std::nullopt;}
    if (found->second.front().empty()) {
        throw CoreConfigError(name + " must not be empty");
    }
    return found->second.front();
}

template <typename Integer>
Integer parse_unsigned(std::string_view text, std::string_view name)
{
    static_assert(std::is_unsigned<Integer>::value, "unsigned type required");
    if (text.empty()) {throw CoreConfigError(std::string(name) + " is empty");}
    Integer value = 0;
    const auto result = std::from_chars(
        text.data(), text.data() + text.size(), value, 10);
    if (result.ec != std::errc() || result.ptr != text.data() + text.size()) {
        throw CoreConfigError(std::string(name) + " must be an unsigned decimal integer");
    }
    return value;
}

double parse_number(std::string_view text, std::string_view name)
{
    if (text.empty() || text.find('\0') != std::string_view::npos) {
        throw CoreConfigError(std::string(name) + " must be a finite number");
    }
    std::size_t index = text.front() == '-' ? 1 : 0;
    if (index == text.size()) {
        throw CoreConfigError(std::string(name) + " must be a finite number");
    }
    if (text[index] == '0') {
        ++index;
    } else if (text[index] >= '1' && text[index] <= '9') {
        while (index < text.size() && text[index] >= '0' && text[index] <= '9') {
            ++index;
        }
    } else {
        throw CoreConfigError(std::string(name) + " must be a finite number");
    }
    if (index < text.size() && text[index] == '.') {
        ++index;
        const std::size_t fraction_begin = index;
        while (index < text.size() && text[index] >= '0' && text[index] <= '9') {
            ++index;
        }
        if (index == fraction_begin) {
            throw CoreConfigError(std::string(name) + " must be a finite number");
        }
    }
    if (index < text.size() && (text[index] == 'e' || text[index] == 'E')) {
        ++index;
        if (index < text.size() && (text[index] == '+' || text[index] == '-')) {
            ++index;
        }
        const std::size_t exponent_begin = index;
        while (index < text.size() && text[index] >= '0' && text[index] <= '9') {
            ++index;
        }
        if (index == exponent_begin) {
            throw CoreConfigError(std::string(name) + " must be a finite number");
        }
    }
    if (index != text.size()) {
        throw CoreConfigError(std::string(name) + " must be a finite number");
    }
    double value = 0.0;
    const auto result = std::from_chars(
        text.data(), text.data() + text.size(), value, std::chars_format::general);
    if (result.ec != std::errc() || result.ptr != text.data() + text.size()
        || !std::isfinite(value)) {
        throw CoreConfigError(std::string(name) + " must be a finite number");
    }
    return value;
}

bool parse_boolean(std::string_view text, std::string_view name)
{
    if (text == "true") {return true;}
    if (text == "false") {return false;}
    throw CoreConfigError(std::string(name) + " must be true or false");
}

crypto::Sha256Digest parse_digest(std::string_view text, std::string_view name)
{
    if (text.size() != 64) {
        throw CoreConfigError(std::string(name) +
                              " must contain 64 lowercase hex digits");
    }
    crypto::Sha256Digest digest = {};
    for (std::size_t index = 0; index < digest.size(); ++index) {
        const auto nibble = [](char character) -> int {
            if (character >= '0' && character <= '9') {return character - '0';}
            if (character >= 'a' && character <= 'f') {
                return character - 'a' + 10;
            }
            return -1;
        };
        const int high = nibble(text[index * 2]);
        const int low = nibble(text[index * 2 + 1]);
        if (high < 0 || low < 0) {
            throw CoreConfigError(std::string(name) +
                                  " must contain 64 lowercase hex digits");
        }
        digest[index] = static_cast<std::uint8_t>((high << 4) | low);
    }
    return digest;
}

bool canonical_uuid(std::string_view value)
{
    if (value.size() != 36) {return false;}
    for (std::size_t index = 0; index < value.size(); ++index) {
        if (index == 8 || index == 13 || index == 18 || index == 23) {
            if (value[index] != '-') {return false;}
            continue;
        }
        const char character = value[index];
        if (!((character >= '0' && character <= '9')
              || (character >= 'a' && character <= 'f'))) {
            return false;
        }
    }
    return true;
}

bool valid_rrbs_cut_site(std::string_view value)
{
    if (value.size() > 1025) {return false;}
    const std::size_t separator = value.find('|');
    if (separator == std::string_view::npos || separator + 1 == value.size()
        || value.find('|', separator + 1) != std::string_view::npos) {
        return false;
    }
    for (std::size_t index = 0; index < value.size(); ++index) {
        if (index == separator) {continue;}
        const char base = value[index];
        if (base != 'A' && base != 'C' && base != 'G'
            && base != 'T' && base != 'N') {
            return false;
        }
    }
    return true;
}

BetaShape parse_beta_shape(std::string_view text, std::string_view name)
{
    const std::size_t separator = text.find(',');
    if (separator == std::string_view::npos
        || text.find(',', separator + 1) != std::string_view::npos) {
        throw CoreConfigError(std::string(name) + " must be ALPHA,BETA");
    }
    const double alpha = parse_number(text.substr(0, separator), name);
    const double beta = parse_number(text.substr(separator + 1), name);
    return BetaShape{alpha, beta};
}

Technology parse_technology(std::string_view value)
{
    if (value == "WGBS") {return Technology::wgbs;}
    if (value == "RRBS") {return Technology::rrbs;}
    if (value == "TBS") {return Technology::tbs;}
    if (value == "WGS") {return Technology::wgs;}
    if (value == "WES") {return Technology::wes;}
    if (value == "TS") {return Technology::ts;}
    throw CoreConfigError(
        "--technology must be WGBS, RRBS, TBS, WGS, WES, or TS");
}

bool whole_genome_technology(Technology technology) noexcept
{
    return technology == Technology::wgbs || technology == Technology::wgs;
}

bool targeted_technology(Technology technology) noexcept
{
    return technology == Technology::tbs || technology == Technology::wes
        || technology == Technology::ts;
}

bool bisulfite_technology(Technology technology) noexcept
{
    return technology == Technology::wgbs || technology == Technology::rrbs
        || technology == Technology::tbs;
}

CoverageMode parse_coverage(std::string_view value)
{
    if (value == "uniform") {return CoverageMode::uniform;}
    if (value == "profile") {return CoverageMode::profile;}
    if (value == "target-score") {return CoverageMode::target_score;}
    throw CoreConfigError(
        "--coverage must be uniform, profile, or target-score");
}

void require_nonempty(
    const std::optional<std::string> &value,
    const char *name)
{
    if (value && value->empty()) {
        throw CoreConfigError(std::string(name) + " must not be empty");
    }
}

bool valid_probability(double value) noexcept
{
    return std::isfinite(value) && value >= 0.0 && value <= 1.0;
}

bool valid_shape(const BetaShape &shape) noexcept
{
    return std::isfinite(shape.alpha) && std::isfinite(shape.beta)
        && shape.alpha > 0.0 && shape.beta > 0.0;
}

} // namespace

void validate_core_config(const CoreConfig &config)
{
    if (!canonical_uuid(config.run_id)) {
        throw CoreConfigError("run_id must be canonical lowercase UUID text");
    }
    if (config.reference_path.empty()) {
        throw CoreConfigError("reference path must not be empty");
    }
    require_nonempty(config.vcf_path, "VCF path");
    require_nonempty(config.cgmap_path, "CGmap path");
    require_nonempty(config.bed_methyl_path, "bedMethyl path");
    require_nonempty(config.methdb_path, "MethDB path");
    require_nonempty(config.asm_path, "ASM path");
    require_nonempty(config.asm_bed_path, "ASM BED path");
    require_nonempty(config.coverage_profile_path, "coverage-profile path");
    require_nonempty(config.rrbs_candidate_bed_path, "RRBS candidate BED path");
    require_nonempty(config.tbs_bed_path, "TBS BED path");

    const bool has_vcf = config.vcf_path.has_value();
    const bool has_cgmap = config.cgmap_path.has_value();
    const bool has_bed_methyl = config.bed_methyl_path.has_value();
    const bool has_methdb = config.methdb_path.has_value();
    const bool has_asm = config.asm_path.has_value();
    const bool has_asm_bed = config.asm_bed_path.has_value();
    if (has_cgmap && has_bed_methyl) {
        throw CoreConfigError(
            "CGmap and bedMethyl inputs are mutually exclusive");
    }
    if (has_methdb
        && (has_cgmap || has_bed_methyl || has_asm || has_asm_bed)) {
        throw CoreConfigError(
            "MethDB and methylation overlays are mutually exclusive");
    }
    if (has_asm && has_asm_bed) {
        throw CoreConfigError(
            "ASM and ASM BED inputs are mutually exclusive");
    }
    if ((has_asm || has_asm_bed) && !has_vcf) {
        throw CoreConfigError("ASM generation requires a VCF input");
    }
    if (config.cgmap_pool && !has_cgmap && !has_bed_methyl) {
        throw CoreConfigError(
            "cgmap_pool=true requires a CGmap or bedMethyl input");
    }
    if (config.cgmap_pool && has_methdb) {
        throw CoreConfigError("MethDB forbids cgmap_pool=true");
    }

    if (config.read_length_1 == 0 || config.insert_min == 0
        || config.insert_mean == 0 || config.insert_max == 0
        || config.chunk_size == 0 || config.core_workers == 0
        || config.protocol_batch_fragments == 0) {
        throw CoreConfigError("read, insert, and chunk sizes must be positive");
    }
    if (config.core_workers > 64U) {
        throw CoreConfigError("core_workers must be in [1, 64]");
    }
    if (config.protocol_batch_fragments > 64U) {
        throw CoreConfigError("protocol_batch_fragments must be in [1, 64]");
    }
    if (!(config.insert_min <= config.insert_mean
          && config.insert_mean <= config.insert_max)) {
        throw CoreConfigError("insert_min <= insert_mean <= insert_max must hold");
    }
    if (!std::isfinite(config.insert_sd) || config.insert_sd < 0.0) {
        throw CoreConfigError("insert_sd must be finite and non-negative");
    }
    const bool fixed_mean_insert =
        (whole_genome_technology(config.technology)
         || targeted_technology(config.technology))
        && config.insert_sd == 0.0;
    const std::uint32_t read_boundary =
        fixed_mean_insert ? config.insert_mean : config.insert_min;
    if (config.read_length_1 > read_boundary) {
        throw CoreConfigError(
            "read length exceeds the fixed mean or minimum variable insert");
    }
    if (config.paired_end) {
        if (!config.read_length_2) {
            throw CoreConfigError("paired-end mode requires read_length_2");
        }
        if (*config.read_length_2 == 0
            || *config.read_length_2 > read_boundary) {
            throw CoreConfigError("read length 2 is outside the insert boundary");
        }
        if (*config.read_length_2 != config.read_length_1) {
            throw CoreConfigError(
                "the current generator requires equal R1 and R2 lengths");
        }
    } else if (config.read_length_2) {
        throw CoreConfigError("single-end mode forbids read_length_2");
    }
    if (config.depth.has_value() == config.fragment_count.has_value()) {
        throw CoreConfigError("exactly one of depth and fragment_count is required");
    }
    if (config.depth
        && (!std::isfinite(*config.depth) || *config.depth <= 0.0)) {
        throw CoreConfigError("depth must be finite and positive");
    }
    if (config.fragment_count && *config.fragment_count == 0U) {
        throw CoreConfigError("fragment_count must be positive");
    }

    if (!valid_probability(config.max_ambiguous_fraction)) {
        throw CoreConfigError("max_ambiguous_fraction must be in [0, 1]");
    }
    if (!valid_probability(config.mutation_rate)
        || !valid_probability(config.indel_fraction)
        || !valid_probability(config.indel_extension_probability)) {
        throw CoreConfigError(
            "mutation probabilities must be finite and in [0, 1]");
    }
    if (!valid_shape(config.beta_cg) || !valid_shape(config.beta_chg)
        || !valid_shape(config.beta_chh)) {
        throw CoreConfigError("Beta shapes must be finite and positive");
    }
    if (has_vcf && config.mutation_rate > 0.0) {
        throw CoreConfigError(
            "VCF and de novo mutation generation are mutually exclusive");
    }
    if (bisulfite_technology(config.technology)
        && (has_vcf || config.mutation_rate > 0.0)
        && !config.update_variant_boundaries) {
        throw CoreConfigError(
            "variant generation requires update_variant_boundaries=true");
    }

    std::unordered_set<std::string> cut_sites;
    for (const std::string &site : config.rrbs_cut_sites) {
        if (!valid_rrbs_cut_site(site)) {
            throw CoreConfigError(
                "RRBS cut site must match ^[ACGTN]*\\|[ACGTN]+$");
        }
        if (!cut_sites.insert(site).second) {
            throw CoreConfigError("duplicate RRBS cut site: " + site);
        }
    }
    if (config.technology == Technology::rrbs) {
        if (config.rrbs_cut_sites.empty()) {
            throw CoreConfigError("RRBS requires at least one cut site");
        }
        if (config.tbs_bed_path || config.tbs_center_stddev) {
            throw CoreConfigError("RRBS forbids TBS inputs");
        }
    } else if (targeted_technology(config.technology)) {
        if (!config.rrbs_cut_sites.empty()) {
            throw CoreConfigError("targeted technologies forbid RRBS cut sites");
        }
        if (config.rrbs_candidate_bed_path) {
            throw CoreConfigError(
                "targeted technologies forbid an RRBS candidate BED");
        }
        if (!config.tbs_bed_path || !config.tbs_center_stddev) {
            throw CoreConfigError(
                "targeted technologies require a BED path and center standard deviation");
        }
        if (!std::isfinite(*config.tbs_center_stddev)
            || *config.tbs_center_stddev < 0.0) {
            throw CoreConfigError(
                "TBS center standard deviation must be finite and non-negative");
        }
        if (config.insert_sd != 0.0) {
            throw CoreConfigError(
                "targeted technologies require --insert-sd 0");
        }
    } else if (whole_genome_technology(config.technology)) {
        if (!config.rrbs_cut_sites.empty()) {
            throw CoreConfigError(
                "whole-genome technologies forbid RRBS cut sites");
        }
        if (config.rrbs_candidate_bed_path) {
            throw CoreConfigError(
                "whole-genome technologies forbid an RRBS candidate BED");
        }
        if (config.tbs_bed_path || config.tbs_center_stddev) {
            throw CoreConfigError(
                "whole-genome technologies forbid targeted inputs");
        }
    } else {
        throw CoreConfigError("technology is outside the core contract");
    }
    if (!bisulfite_technology(config.technology)
        && (has_cgmap || has_bed_methyl || has_methdb || has_asm
            || has_asm_bed || config.cgmap_pool)) {
        throw CoreConfigError(
            "standard sequencing forbids methylation inputs and pooling");
    }
    if (!bisulfite_technology(config.technology) && !config.directional) {
        throw CoreConfigError(
            "standard sequencing requires directional=true as an inert value");
    }
    const bool has_coverage_artifact = config.coverage_profile_path.has_value();
    if (config.coverage == CoverageMode::profile) {
        if (whole_genome_technology(config.technology)) {
            if (!config.coverage_profile_path) {
                throw CoreConfigError(
                    "whole-genome profile coverage requires a profile path");
            }
            if (config.insert_sd != 0.0
                && (has_vcf || config.mutation_rate != 0.0)) {
                throw CoreConfigError(
                    "variable-insert target GC does not yet support variants");
            }
        } else if (config.technology == Technology::rrbs) {
            if (!config.rrbs_candidate_bed_path) {
                throw CoreConfigError(
                    "RRBS profile coverage requires a candidate BED");
            }
            if (has_coverage_artifact) {
                throw CoreConfigError(
                    "RRBS profile coverage forbids whole-genome profile inputs");
            }
        } else {
            throw CoreConfigError(
                "profile coverage supports WGBS, WGS, or RRBS only");
        }
    } else if (config.coverage == CoverageMode::target_score) {
        if (!targeted_technology(config.technology)) {
            throw CoreConfigError("target-score coverage requires TBS, WES, or TS");
        }
        if (has_coverage_artifact) {
            throw CoreConfigError(
                "target-score coverage forbids coverage-profile inputs");
        }
    } else if (config.coverage == CoverageMode::uniform) {
        if (has_coverage_artifact) {
            throw CoreConfigError(
                "uniform coverage forbids coverage-profile inputs");
        }
    } else {
        throw CoreConfigError("coverage mode is outside the core contract");
    }
}

CoreConfig parse_core_config(const std::vector<std::string> &arguments)
{
    const Values values = collect_values(arguments);
    CoreConfig config;
    if (values.find("--emit-details") != values.end()) {
        config.emit_details = parse_boolean(
            required(values, "--emit-details"), "--emit-details");
    }
    if (values.find("--protocol-batch-fragments") != values.end()) {
        config.protocol_batch_fragments = parse_unsigned<std::uint32_t>(
            required(values, "--protocol-batch-fragments"),
            "--protocol-batch-fragments");
    }
    config.run_id = required(values, "--run-id");
    if (!canonical_uuid(config.run_id)) {
        throw CoreConfigError("--run-id must be canonical lowercase UUID text");
    }
    config.normalized_config_sha256 = parse_digest(
        required(values, "--config-sha256"), "--config-sha256");
    config.master_seed = parse_unsigned<std::uint64_t>(
        required(values, "--seed"), "--seed");
    config.mutation_seed = parse_unsigned<std::uint64_t>(
        required(values, "--seed-mut"), "--seed-mut");
    config.phasing_seed = parse_unsigned<std::uint64_t>(
        required(values, "--seed-phase"), "--seed-phase");
    config.methylation_seed = parse_unsigned<std::uint64_t>(
        required(values, "--seed-meth"), "--seed-meth");
    config.reference_path = required(values, "--reference");
    config.vcf_path = optional_text(values, "--vcf");
    config.cgmap_path = optional_text(values, "--cgmap");
    config.bed_methyl_path = optional_text(values, "--bed-methyl");
    config.methdb_path = optional_text(values, "--methdb");
    config.asm_path = optional_text(values, "--asm");
    config.asm_bed_path = optional_text(values, "--asm-bed");
    config.technology = parse_technology(required(values, "--technology"));
    config.directional = parse_boolean(
        required(values, "--directional"), "--directional");

    config.paired_end = parse_boolean(
        required(values, "--paired-end"), "--paired-end");
    config.read_length_1 = parse_unsigned<std::uint32_t>(
        required(values, "--read-length-1"), "--read-length-1");
    if (values.find("--read-length-2") != values.end()) {
        config.read_length_2 = parse_unsigned<std::uint32_t>(
            required(values, "--read-length-2"), "--read-length-2");
    }
    config.insert_min = parse_unsigned<std::uint32_t>(
        required(values, "--insert-min"), "--insert-min");
    config.insert_mean = parse_unsigned<std::uint32_t>(
        required(values, "--insert-mean"), "--insert-mean");
    config.insert_max = parse_unsigned<std::uint32_t>(
        required(values, "--insert-max"), "--insert-max");
    config.insert_sd = parse_number(
        required(values, "--insert-sd"), "--insert-sd");
    if (values.find("--depth") != values.end()) {
        config.depth = parse_number(required(values, "--depth"), "--depth");
    }
    if (values.find("--fragments") != values.end()) {
        config.fragment_count = parse_unsigned<std::uint32_t>(
            required(values, "--fragments"), "--fragments");
    }
    config.max_ambiguous_fraction = parse_number(
        required(values, "--max-ambiguous-fraction"),
        "--max-ambiguous-fraction");
    config.chunk_size = parse_unsigned<std::uint32_t>(
        required(values, "--chunk-size"), "--chunk-size");
    if (values.find("--core-workers") != values.end()) {
        config.core_workers = parse_unsigned<std::uint32_t>(
            required(values, "--core-workers"), "--core-workers");
    }

    config.mutation_rate = parse_number(
        required(values, "--mutation-rate"), "--mutation-rate");
    config.indel_fraction = parse_number(
        required(values, "--indel-fraction"), "--indel-fraction");
    config.indel_extension_probability = parse_number(
        required(values, "--indel-extension-probability"),
        "--indel-extension-probability");
    config.homozygous_only = parse_boolean(
        required(values, "--homozygous-only"), "--homozygous-only");
    config.collect_non_cpg = parse_boolean(
        required(values, "--collect-non-cpg"), "--collect-non-cpg");
    config.cgmap_pool = parse_boolean(
        required(values, "--cgmap-pool"), "--cgmap-pool");
    config.update_variant_boundaries = parse_boolean(
        required(values, "--update-variant-boundaries"),
        "--update-variant-boundaries");
    config.beta_cg = parse_beta_shape(
        required(values, "--beta-cg"), "--beta-cg");
    config.beta_chg = parse_beta_shape(
        required(values, "--beta-chg"), "--beta-chg");
    config.beta_chh = parse_beta_shape(
        required(values, "--beta-chh"), "--beta-chh");

    config.coverage = parse_coverage(required(values, "--coverage"));
    config.coverage_profile_path = optional_text(values, "--coverage-profile");
    const auto rrbs = values.find("--rrbs-cut-site");
    if (rrbs != values.end()) {
        for (const std::string &site : rrbs->second) {
            config.rrbs_cut_sites.push_back(site);
        }
    }
    config.rrbs_candidate_bed_path = optional_text(
        values, "--rrbs-candidate-bed");
    config.tbs_bed_path = optional_text(values, "--tbs-bed");
    if (values.find("--tbs-center-stddev") != values.end()) {
        config.tbs_center_stddev = parse_number(
            required(values, "--tbs-center-stddev"),
            "--tbs-center-stddev");
    }

    validate_core_config(config);
    return config;
}

CoreConfig parse_core_config(int argc, char *argv[])
{
    if (argc < 1 || argv == nullptr) {
        throw CoreConfigError("invalid process argv");
    }
    std::vector<std::string> arguments;
    arguments.reserve(static_cast<std::size_t>(argc - 1));
    for (int index = 1; index < argc; ++index) {
        if (argv[index] == nullptr) {
            throw CoreConfigError("argv contains a null argument");
        }
        arguments.emplace_back(argv[index]);
    }
    return parse_core_config(arguments);
}

} // namespace htsim::core

// ---- generator --------------------------------------------------------

namespace htsim::core {
namespace {

bool uses_variable_wgbs_insert(const CoreConfig &config) noexcept
{
    return whole_genome_technology(config.technology)
        && config.insert_sd != 0.0;
}

tbs::SamplingMode tbs_sampling_mode(const CoreConfig &config) noexcept
{
    return config.coverage == CoverageMode::target_score
        ? tbs::SamplingMode::output_weight
        : tbs::SamplingMode::uniform;
}

struct LibraryOrientation {
    model::CaptureStrand informative_strand;
    bool reverse_molecule;
};

LibraryOrientation sample_library_orientation(
    const CoreConfig &config,
    std::uint64_t key,
    std::uint64_t fragment_ordinal,
    model::CaptureStrand constraint)
{
    if (!bisulfite_technology(config.technology)) {
        return {constraint, false};
    }
    if (constraint != model::CaptureStrand::unknown
        && constraint != model::CaptureStrand::forward
        && constraint != model::CaptureStrand::reverse) {
        throw CoreGeneratorError("fragment capture strand is invalid");
    }

    const bool informative_reverse =
        constraint == model::CaptureStrand::reverse
        || (constraint == model::CaptureStrand::unknown
            && rng::bernoulli(key, fragment_ordinal, 0U, 0.5, 0U));
    const bool complementary = !config.directional
        && rng::bernoulli(key, fragment_ordinal, 0U, 0.5, 1U);
    return {
        informative_reverse
            ? model::CaptureStrand::reverse
            : model::CaptureStrand::forward,
        informative_reverse != complementary,
    };
}

void require_generation_environment(const CoreConfig &config)
{
    validate_core_config(config);
    if (std::fegetround() != FE_TONEAREST) {
        throw CoreGeneratorError("minimal core requires round-to-nearest floating point");
    }
    fragment_builder::require_payload_fits_protocol({
        config.technology == Technology::rrbs
            || uses_variable_wgbs_insert(config)
            ? config.insert_max
            : config.insert_mean,
        config.read_length_1,
        config.paired_end});
}

std::uint64_t checked_add(
    std::uint64_t left,
    std::uint64_t right,
    const char *field)
{
    if (right > std::numeric_limits<std::uint64_t>::max() - left) {
        throw CoreGeneratorError(std::string(field) + " size exceeds uint64");
    }
    return left + right;
}

crypto::Sha256Digest methdb_binding(
    const CoreConfig &config,
    const crypto::Sha256Digest &reference_sha256,
    const std::optional<crypto::Sha256Digest> &vcf_sha256)
{
    crypto::Sha256 hash;
    const auto bytes = [&](const void *data, std::size_t size) {
        hash.update(static_cast<const std::uint8_t *>(data), size);
    };
    const auto u64 = [&](std::uint64_t value) {
        std::uint8_t encoded[8];
        for (unsigned index = 0U; index < 8U; ++index) {
            encoded[index] = static_cast<std::uint8_t>(value >> (index * 8U));
        }
        bytes(encoded, sizeof(encoded));
    };
    const auto f64 = [&](double value) {
        std::uint64_t encoded = 0U;
        static_assert(sizeof(encoded) == sizeof(value));
        std::memcpy(&encoded, &value, sizeof(encoded));
        u64(encoded);
    };
    static constexpr char contract[] = "methdb-binding";
    bytes(contract, sizeof(contract) - 1U);
    bytes(reference_sha256.data(), reference_sha256.size());
    if (vcf_sha256) {
        u64(config.methylation_seed);
        u64(config.phasing_seed);
        const std::uint8_t marker = 1U;
        bytes(&marker, 1U);
        bytes(vcf_sha256->data(), vcf_sha256->size());
    } else {
        u64(config.methylation_seed);
        const std::uint8_t marker = 0U;
        bytes(&marker, 1U);
        if (config.mutation_rate > 0.0) {u64(config.mutation_seed);}
        f64(config.mutation_rate);
        f64(config.indel_fraction);
        f64(config.indel_extension_probability);
        const std::uint8_t homozygous = config.homozygous_only ? 1U : 0U;
        bytes(&homozygous, 1U);
    }
    return hash.digest();
}

template <typename Range, typename Start, typename End>
std::uint64_t interval_union_bases(
    const Range &values,
    Start start_of,
    End end_of)
{
    std::vector<std::pair<std::uint32_t, std::uint32_t>> intervals;
    intervals.reserve(values.size());
    for (const auto &value : values) {
        const std::uint32_t begin = start_of(value);
        const std::uint32_t end = end_of(value);
        if (end > begin) {intervals.emplace_back(begin, end);}
    }
    std::sort(intervals.begin(), intervals.end());
    std::uint64_t total = 0U;
    std::uint32_t begin = 0U;
    std::uint32_t end = 0U;
    bool active = false;
    for (const auto &interval : intervals) {
        if (!active) {
            begin = interval.first;
            end = interval.second;
            active = true;
        } else if (interval.first <= end) {
            end = std::max(end, interval.second);
        } else {
            total = checked_add(total, end - begin, "target interval union");
            begin = interval.first;
            end = interval.second;
        }
    }
    return active
        ? checked_add(total, end - begin, "target interval union")
        : 0U;
}

std::uint64_t variant_payload_extra(
    const variant::Variant &event)
{
    std::uint64_t bytes = UINT64_C(32);
    bytes = checked_add(bytes, event.ref_bases.size(), "variant payload");
    bytes = checked_add(bytes, event.alt_bases.size(), "variant payload");
    if (event.kind == model::VariantKind::insertion) {
        bytes = checked_add(
            bytes,
            UINT64_C(37) * event.alt_bases.size(),
            "variant payload");
    }
    return bytes;
}

std::uint64_t maximum_variant_payload_extra(
    const variant::ContigVariants &variants,
    std::uint32_t reference_span,
    std::uint8_t haplotype)
{
    if (haplotype > 1U || reference_span == 0U) {
        throw CoreGeneratorError("variant payload preflight input is invalid");
    }
    if (reference_span > variants.reference_length()) {return 0U;}

    struct WeightedEvent {
        std::uint32_t position;
        std::uint64_t bytes;
    };
    std::vector<WeightedEvent> weighted;
    weighted.reserve(variants.variants().size());
    for (const variant::Variant &event : variants.variants()) {
        if (!model::mask_contains(event.alt_haplotypes, haplotype)) {
            continue;
        }
        std::uint32_t position = event.reference_start;
        if (event.kind == model::VariantKind::insertion
            && position == variants.reference_length()) {
            if (position == 0U) {
                throw CoreGeneratorError(
                    "terminal insertion cannot address an empty contig");
            }
            --position;
        }
        weighted.push_back({position, variant_payload_extra(event)});
    }

    std::size_t first = 0U;
    std::uint64_t current = 0U;
    std::uint64_t maximum = 0U;
    for (std::size_t last = 0U; last < weighted.size(); ++last) {
        current = checked_add(current, weighted[last].bytes, "variant payload");
        while (weighted[last].position - weighted[first].position
               >= reference_span) {
            current -= weighted[first].bytes;
            ++first;
        }
        maximum = std::max(maximum, current);
    }
    return maximum;
}

void require_variant_payload_fits_protocol(
    const variant::ContigVariants &variants,
    const fragment_builder::ReadLayout &layout)
{
    const std::uint64_t baseline =
        fragment_builder::maximum_payload_bytes(layout);
    const std::uint64_t extra = std::max(
        maximum_variant_payload_extra(variants, layout.insert_length, 0U),
        maximum_variant_payload_extra(variants, layout.insert_length, 1U));
    if (checked_add(baseline, extra, "variant payload")
        > protocol::maximum_frame_payload) {
        throw CoreGeneratorError(
            "a VCF-projected fragment can exceed the protocol frame limit");
    }
}

void require_haplotype_payload_fits_protocol(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    const fragment_builder::ReadLayout &layout)
{
    std::array<std::uint64_t, 2> extra = {0U, 0U};
    for (std::uint8_t haplotype = 0U; haplotype < 2U; ++haplotype) {
        const haplotype::HaplotypeLayout physical_layout(
            contig, variants, haplotype, false);
        extra[haplotype] =
            physical_layout.maximum_variant_payload_bytes(
                variants, layout.insert_length);
    }
    const std::uint64_t baseline =
        fragment_builder::maximum_payload_bytes(layout);
    if (checked_add(
            baseline,
            std::max(extra[0], extra[1]),
            "haplotype variant payload")
        > protocol::maximum_frame_payload) {
        throw CoreGeneratorError(
            "a haplotype-fragmented record can exceed the protocol frame limit");
    }
}

void require_header_fits_protocol(
    const CoreConfig &config,
    const std::vector<reference::ContigMetadata> &catalog)
{
    // Three strings, seed, config digest, and contig count.
    std::uint64_t bytes = UINT64_C(4) + config.run_id.size();
    bytes = checked_add(bytes, UINT64_C(4) + sizeof(core_version) - 1U, "header");
    bytes = checked_add(
        bytes, UINT64_C(4) + protocol::rng_contract.size(), "header");
    bytes = checked_add(bytes, UINT64_C(8) + 32U, "header");
    // technology/details/mates/encoding/ambiguity/reserved, two read lengths,
    // and the contig count.
    bytes = checked_add(bytes, UINT64_C(8) + 8U + 4U, "header");
    for (const reference::ContigMetadata &contig : catalog) {
        bytes = checked_add(bytes, UINT64_C(4) + contig.name.size(), "header");
        bytes = checked_add(bytes, UINT64_C(4) + 32U, "header");
    }
    bytes = checked_add(bytes, (4U - bytes % 4U) % 4U, "header padding");
    if (bytes > protocol::maximum_frame_payload) {
        throw CoreGeneratorError("protocol header exceeds the frame limit");
    }
}

wgbs::FixedFragmentShape sampling_shape(const CoreConfig &config)
{
    return {
        config.insert_mean,
        config.read_length_1,
        config.paired_end,
        config.max_ambiguous_fraction,
    };
}

methdb::ContextShapes methylation_shapes(const CoreConfig &config)
{
    return {
        {config.beta_cg.alpha, config.beta_cg.beta},
        {config.beta_chg.alpha, config.beta_chg.beta},
        {config.beta_chh.alpha, config.beta_chh.beta},
    };
}

std::uint8_t choose_haplotype(
    model::HaplotypeMask mask,
    std::uint64_t key,
    std::uint64_t fragment_ordinal)
{
    if (!model::is_haplotype_mask(
            static_cast<std::uint8_t>(mask))) {
        throw CoreGeneratorError("fragment candidate has an invalid haplotype mask");
    }
    if (mask == model::HaplotypeMask::haplotype_1) {return 0;}
    if (mask == model::HaplotypeMask::haplotype_2) {return 1;}
    return static_cast<std::uint8_t>(
        rng::bernoulli(key, fragment_ordinal, UINT64_C(0), 0.5) ? 1 : 0);
}

protocol::Header make_header(
    const CoreConfig &config,
    const std::vector<reference::ContigMetadata> &catalog)
{
    protocol::Header header;
    header.run_id = config.run_id;
    header.core_version = core_version;
    header.rng_contract = std::string(protocol::rng_contract);
    header.master_seed = config.master_seed;
    header.normalized_config_sha256 = config.normalized_config_sha256;
    switch (config.technology) {
    case Technology::wgbs:
        header.technology = protocol::Technology::wgbs;
        break;
    case Technology::rrbs:
        header.technology = protocol::Technology::rrbs;
        break;
    case Technology::tbs:
        header.technology = protocol::Technology::tbs;
        break;
    case Technology::wgs:
        header.technology = protocol::Technology::wgs;
        break;
    case Technology::wes:
        header.technology = protocol::Technology::wes;
        break;
    case Technology::ts:
        header.technology = protocol::Technology::ts;
        break;
    }
    header.has_details = config.emit_details;
    header.mates_per_fragment = config.paired_end ? 2U : 1U;
    header.base_encoding = protocol::BaseEncoding::acgtn_u8;
    header.ambiguity_policy = protocol::AmbiguityPolicy::preserve_n;
    header.read_length_r1 = config.read_length_1;
    header.read_length_r2 = config.paired_end ? *config.read_length_2 : 0U;
    header.contigs.reserve(catalog.size());
    for (const reference::ContigMetadata &contig : catalog) {
        header.contigs.push_back({
            contig.name,
            static_cast<std::uint32_t>(contig.length),
            contig.reference_sha256,
        });
    }
    return header;
}

} // namespace

void generate_rrbs_candidate_bed(
    const CoreConfig &config,
    std::ostream &sink)
{
    try {
        require_generation_environment(config);
        if (config.technology != Technology::rrbs) {
            throw CoreGeneratorError(
                "RRBS candidate export requires technology RRBS");
        }
        if (config.coverage != CoverageMode::uniform
            || config.rrbs_candidate_bed_path) {
            throw CoreGeneratorError(
                "RRBS candidate export requires uniform coverage and forbids a candidate BED input");
        }
        reference::ReferenceSnapshot snapshot(config.reference_path);
        const auto &catalog = snapshot.catalog();
        for (const reference::ContigMetadata &contig : catalog) {
            if (contig.length > std::numeric_limits<std::uint32_t>::max()) {
                throw CoreGeneratorError(
                    "RRBS candidate export requires contig length <= UINT32_MAX");
            }
        }
        std::unique_ptr<variant::VariantFile> variant_file;
        if (config.vcf_path) {
            variant_file = std::make_unique<variant::VariantFile>(
                *config.vcf_path,
                catalog,
                config.phasing_seed);
        }
        const bool generate_mutations = config.mutation_rate > 0.0;
        const variant::MutationParameters mutation_parameters{
            config.mutation_rate,
            config.indel_fraction,
            config.indel_extension_probability,
            config.homozygous_only,
        };
        const std::vector<rrbs::CutSite> cut_sites =
            rrbs::parse_cut_sites(config.rrbs_cut_sites);

        rrbs::write_candidate_bed_header(sink);
        snapshot.visit_contigs([&](const reference::Contig &contig) {
            std::vector<variant::Variant> generated_variants;
            const std::vector<variant::Variant> *variant_records = nullptr;
            if (variant_file) {
                variant_records = &variant_file->variants(contig.index);
            } else if (generate_mutations) {
                generated_variants = variant::generate_de_novo_events(
                    contig, config.mutation_seed, mutation_parameters);
                variant_records = &generated_variants;
            }
            std::unique_ptr<variant::ContigVariants> contig_variants;
            if (variant_records != nullptr) {
                contig_variants = std::make_unique<variant::ContigVariants>(
                    contig.bases, *variant_records, contig.index);
            }
            if (contig_variants && !contig_variants->variants().empty()) {
                const rrbs::DiploidCandidateCatalog candidate_catalog(
                    contig,
                    *contig_variants,
                    cut_sites,
                    config.insert_min,
                    config.insert_max,
                    config.read_length_1,
                    config.paired_end,
                    config.max_ambiguous_fraction);
                rrbs::write_candidate_bed_contig(
                    sink, contig.name, candidate_catalog.candidates());
            } else {
                const rrbs::CandidateCatalog candidate_catalog(
                    contig.bases,
                    cut_sites,
                    config.insert_min,
                    config.insert_max,
                    config.read_length_1,
                    config.paired_end,
                    config.max_ambiguous_fraction);
                rrbs::write_candidate_bed_contig(
                    sink, contig.name, candidate_catalog.candidates());
            }
        });
        if (!sink) {
            throw CoreGeneratorError(
                "failed while flushing the RRBS candidate BED");
        }
    } catch (const CoreGeneratorError &) {
        throw;
    } catch (const std::exception &error) {
        throw CoreGeneratorError(error.what());
    }
}

void generate_methdb_catalog(
    const CoreConfig &config,
    std::ostream &sink)
{
    try {
        require_generation_environment(config);
        if (config.methdb_path) {
            throw CoreGeneratorError(
                "MethDB export cannot load another MethDB snapshot");
        }
        reference::ReferenceSnapshot snapshot(config.reference_path);
        const auto &catalog = snapshot.catalog();
        if (catalog.size() > std::numeric_limits<std::uint32_t>::max()) {
            throw CoreGeneratorError("MethDB contig count exceeds uint32");
        }
        std::unique_ptr<variant::VariantFile> variant_file;
        if (config.vcf_path) {
            variant_file = std::make_unique<variant::VariantFile>(
                *config.vcf_path,
                catalog,
                config.phasing_seed);
        }
        std::optional<crypto::Sha256Digest> vcf_sha256;
        if (variant_file) {vcf_sha256 = variant_file->file_sha256();}
        const bool generate_mutations = config.mutation_rate > 0.0;
        const variant::MutationParameters mutation_parameters{
            config.mutation_rate,
            config.indel_fraction,
            config.indel_extension_probability,
            config.homozygous_only,
        };
        std::unique_ptr<methdb::CgmapProfile> cgmap_profile;
        if (config.cgmap_path || config.bed_methyl_path) {
            const bool bed_methyl = config.bed_methyl_path.has_value();
            cgmap_profile = std::make_unique<methdb::CgmapProfile>(
                bed_methyl ? *config.bed_methyl_path : *config.cgmap_path,
                catalog,
                bed_methyl
                    ? methdb::MethylationProfileFormat::bed_methyl
                    : methdb::MethylationProfileFormat::cgmap);
            if (config.cgmap_pool
                && cgmap_profile->defined_probability_count() == 0U) {
                throw CoreGeneratorError(
                    "methylation-profile pooling requires a defined probability");
            }
        }
        std::unique_ptr<methdb::AsmProfile> asm_profile;
        if (config.asm_path || config.asm_bed_path) {
            const bool asm_bed = config.asm_bed_path.has_value();
            asm_profile = std::make_unique<methdb::AsmProfile>(
                asm_bed ? *config.asm_bed_path : *config.asm_path,
                catalog,
                asm_bed
                    ? methdb::AsmProfileFormat::bed
                    : methdb::AsmProfileFormat::htsim);
        }
        const methdb::ContextShapes shapes = methylation_shapes(config);
        methdb::SnapshotWriter writer(
            sink,
            methdb_binding(config, snapshot.file_sha256(), vcf_sha256),
            static_cast<std::uint32_t>(catalog.size()));
        snapshot.visit_contigs([&](const reference::Contig &contig) {
            std::vector<methdb::CgmapRecord> cgmap_records;
            if (cgmap_profile) {
                cgmap_profile->validate_contig(contig);
                cgmap_records = cgmap_profile->records(contig);
            }
            const auto *cgmap_or_null = cgmap_profile ? &cgmap_records : nullptr;
            std::vector<methdb::AsmRecord> asm_records;
            if (asm_profile) {
                asm_profile->validate_contig(contig);
                asm_records = asm_profile->records(contig);
            }
            const auto *asm_or_null = asm_profile ? &asm_records : nullptr;
            std::vector<variant::Variant> generated_variants;
            const std::vector<variant::Variant> *variant_records = nullptr;
            if (variant_file) {
                variant_records = &variant_file->variants(contig.index);
            } else if (generate_mutations) {
                generated_variants = variant::generate_de_novo_events(
                    contig, config.mutation_seed, mutation_parameters);
                variant_records = &generated_variants;
            }
            if (variant_records) {
                const variant::ContigVariants contig_variants(
                    contig.bases, *variant_records, contig.index);
                const methdb::DiploidMethylationCatalog methylation(
                    contig,
                    contig_variants,
                    config.methylation_seed,
                    config.collect_non_cpg,
                    shapes,
                    cgmap_or_null,
                    asm_or_null,
                    config.cgmap_pool);
                writer.write_diploid(catalog.at(contig.index), methylation);
            } else {
                const methdb::MethylationCatalog methylation(
                    contig.bases,
                    contig.index,
                    config.methylation_seed,
                    config.collect_non_cpg,
                    shapes,
                    cgmap_or_null,
                    config.cgmap_pool);
                writer.write_reference(catalog.at(contig.index), methylation);
            }
        });
        writer.finish();
    } catch (const CoreGeneratorError &) {
        throw;
    } catch (const std::exception &error) {
        throw CoreGeneratorError(error.what());
    }
}

void generate_variant_catalog_vcf(
    const CoreConfig &config,
    std::ostream &sink)
{
    try {
        require_generation_environment(config);
        reference::ReferenceSnapshot snapshot(config.reference_path);
        const auto &catalog = snapshot.catalog();
        std::unique_ptr<variant::VariantFile> variant_file;
        if (config.vcf_path) {
            variant_file = std::make_unique<variant::VariantFile>(
                *config.vcf_path, catalog, config.phasing_seed);
        }
        const variant::MutationParameters mutation_parameters{
            config.mutation_rate,
            config.indel_fraction,
            config.indel_extension_probability,
            config.homozygous_only,
        };

        variant::write_vcf_header(sink);
        snapshot.visit_contigs([&](const reference::Contig &contig) {
            if (variant_file) {
                variant::write_vcf_contig(
                    sink, contig, variant_file->variants(contig.index));
            } else {
                const std::vector<variant::Variant> events =
                    variant::generate_de_novo_events(
                        contig, config.mutation_seed, mutation_parameters);
                variant::write_vcf_contig(sink, contig, events);
            }
        });
        if (!sink) {
            throw CoreGeneratorError("failed while flushing the variant VCF");
        }
    } catch (const CoreGeneratorError &) {
        throw;
    } catch (const std::exception &error) {
        throw CoreGeneratorError(error.what());
    }
}

protocol::Trailer generate_core_stream(
    const CoreConfig &config,
    std::ostream &sink)
{
    try {
        require_generation_environment(config);
        reference::ReferenceSnapshot snapshot(config.reference_path);
        const auto &catalog = snapshot.catalog();
        for (const reference::ContigMetadata &contig : catalog) {
            if (contig.length > std::numeric_limits<std::uint32_t>::max()) {
                throw CoreGeneratorError(
                    "minimal core requires every contig length <= UINT32_MAX");
            }
        }
        require_header_fits_protocol(config, catalog);

        std::unique_ptr<variant::VariantFile> variant_file;
        if (config.vcf_path) {
            variant_file = std::make_unique<variant::VariantFile>(
                *config.vcf_path,
                catalog,
                config.phasing_seed);
        }
        std::optional<crypto::Sha256Digest> vcf_sha256;
        if (variant_file) {vcf_sha256 = variant_file->file_sha256();}
        const bool generate_mutations = config.mutation_rate > 0.0;
        const variant::MutationParameters mutation_parameters{
            config.mutation_rate,
            config.indel_fraction,
            config.indel_extension_probability,
            config.homozygous_only,
        };
        std::unique_ptr<methdb::CgmapProfile> cgmap_profile;
        if (config.cgmap_path || config.bed_methyl_path) {
            const bool bed_methyl = config.bed_methyl_path.has_value();
            cgmap_profile = std::make_unique<methdb::CgmapProfile>(
                bed_methyl ? *config.bed_methyl_path : *config.cgmap_path,
                catalog,
                bed_methyl
                    ? methdb::MethylationProfileFormat::bed_methyl
                    : methdb::MethylationProfileFormat::cgmap);
            if (config.cgmap_pool
                && cgmap_profile->defined_probability_count() == 0U) {
                throw CoreGeneratorError(
                    "methylation-profile pooling requires at least one defined probability");
            }
        }
        std::unique_ptr<methdb::AsmProfile> asm_profile;
        if (config.asm_path || config.asm_bed_path) {
            const bool asm_bed = config.asm_bed_path.has_value();
            asm_profile = std::make_unique<methdb::AsmProfile>(
                asm_bed ? *config.asm_bed_path : *config.asm_path,
                catalog,
                asm_bed
                    ? methdb::AsmProfileFormat::bed
                    : methdb::AsmProfileFormat::htsim);
        }
        const methdb::ContextShapes context_shapes =
            methylation_shapes(config);
        const wgbs::FixedFragmentShape sampler_shape = sampling_shape(config);
        const bool variable_wgbs_insert = uses_variable_wgbs_insert(config);
        const wgbs::FixedFragmentShape target_calibration_shape{
            config.insert_mean,
            config.read_length_1,
            config.paired_end,
            config.max_ambiguous_fraction,
        };
        const insert_length::Parameters insert_parameters{
            config.insert_min,
            config.insert_mean,
            config.insert_max,
            config.insert_sd,
        };
        const std::vector<rrbs::CutSite> rrbs_cut_sites =
            config.technology == Technology::rrbs
            ? rrbs::parse_cut_sites(config.rrbs_cut_sites)
            : std::vector<rrbs::CutSite>{};
        const bool rrbs_profile = config.technology == Technology::rrbs
            && config.coverage == CoverageMode::profile;
        std::unique_ptr<rrbs::CandidateBed> rrbs_candidate_bed;
        if (config.rrbs_candidate_bed_path) {
            rrbs_candidate_bed = std::make_unique<rrbs::CandidateBed>(
                *config.rrbs_candidate_bed_path, catalog);
        }
        std::unique_ptr<tbs::TargetFile> tbs_targets;
        if (targeted_technology(config.technology)) {
            tbs_targets = std::make_unique<tbs::TargetFile>(
                *config.tbs_bed_path, catalog);
        }
        std::unique_ptr<wgbs::WgbsGcProfile> coverage_profile;
        if (whole_genome_technology(config.technology)
            && config.coverage == CoverageMode::profile) {
            coverage_profile = std::make_unique<wgbs::WgbsGcProfile>(
                *config.coverage_profile_path);
        }
        const bool haplotype_gc_profile = coverage_profile
            && (variant_file || generate_mutations);
        std::unique_ptr<methdb::Snapshot> fixed_methdb;
        if (config.methdb_path) {
            fixed_methdb = std::make_unique<methdb::Snapshot>(
                *config.methdb_path,
                methdb_binding(
                    config, snapshot.file_sha256(), vcf_sha256),
                catalog);
        }
        std::vector<std::uint32_t> candidate_weights(catalog.size(), 0);
        std::vector<double> rrbs_profile_weights(catalog.size(), 0.0);
        std::vector<std::vector<std::uint32_t>> target_bin_counts(
            catalog.size());
        std::vector<std::uint64_t> target_reference_bases(catalog.size(), 0U);
        snapshot.visit_contigs([&](const reference::Contig &contig) {
            if (contig.index >= candidate_weights.size()) {
                throw CoreGeneratorError("reference contig index exceeds its catalog");
            }
            if (cgmap_profile) {
                cgmap_profile->validate_contig(contig);
            }
            if (asm_profile) {
                asm_profile->validate_contig(contig);
            }
            std::vector<variant::Variant> generated_events;
            const std::vector<variant::Variant> *variants = nullptr;
            if (variant_file) {
                variants = &variant_file->variants(contig.index);
            } else if (generate_mutations) {
                generated_events = variant::generate_de_novo_events(
                    contig, config.mutation_seed, mutation_parameters);
                variants = &generated_events;
            }
            std::unique_ptr<variant::ContigVariants> planned_variants;
            if (variants != nullptr) {
                planned_variants = std::make_unique<variant::ContigVariants>(
                    contig.bases,
                    *variants,
                    contig.index);
                const fragment_builder::ReadLayout payload_layout{
                    variable_wgbs_insert
                        ? config.insert_max
                        : (config.technology == Technology::rrbs
                            ? config.insert_max
                            : config.insert_mean),
                    config.read_length_1,
                    config.paired_end,
                    whole_genome_technology(config.technology)
                        ? fragment_builder::ReadLayout::InsertCoordinate::reference
                        : fragment_builder::ReadLayout::InsertCoordinate::haplotype,
                };
                if (whole_genome_technology(config.technology)) {
                    require_variant_payload_fits_protocol(
                        *planned_variants, payload_layout);
                } else {
                    require_haplotype_payload_fits_protocol(
                        contig, *planned_variants, payload_layout);
                }
            }
            if (asm_profile) {
                const std::vector<methdb::AsmRecord> asm_records =
                    asm_profile->records(contig);
                if (!asm_records.empty()) {
                    if (!planned_variants) {
                        throw CoreGeneratorError(
                            "ASM planning lost its required VCF variant set");
                    }
                    std::vector<methdb::CgmapRecord> cgmap_records;
                    if (cgmap_profile) {
                        cgmap_records = cgmap_profile->records(contig);
                    }
                    const auto *cgmap_records_or_null = cgmap_profile
                        ? &cgmap_records
                        : nullptr;
                    (void)methdb::DiploidMethylationCatalog(
                        contig,
                        *planned_variants,
                        config.methylation_seed,
                        config.collect_non_cpg,
                        context_shapes,
                        cgmap_records_or_null,
                        &asm_records,
                        config.cgmap_pool);
                }
            }
            if (whole_genome_technology(config.technology)) {
                if (coverage_profile) {
                    if (variable_wgbs_insert) {
                        if (planned_variants) {
                            throw CoreGeneratorError(
                                "variable target GC escaped its variant gate");
                        }
                        const wgbs::WgbsGcSampler target_domain(
                            contig.bases,
                            target_calibration_shape,
                            *coverage_profile);
                        candidate_weights[contig.index] =
                            wgbs::VariableWgbsSampler(
                                contig.bases,
                                contig.index,
                                config.master_seed,
                                insert_parameters,
                                config.read_length_1,
                                config.paired_end,
                                config.max_ambiguous_fraction)
                                .allocation_weight();
                        target_bin_counts[contig.index].assign(
                            coverage_profile->bin_count(), 0U);
                        if (candidate_weights[contig.index] > 0U) {
                            target_bin_counts[contig.index] =
                                target_domain.bin_opportunity_counts();
                        }
                    } else if (planned_variants) {
                        const wgbs::HaplotypeGcSampler target_domain(
                            contig,
                            *planned_variants,
                            sampler_shape,
                            *coverage_profile);
                        candidate_weights[contig.index] =
                            target_domain.physical_candidate_count();
                        target_bin_counts[contig.index] =
                            target_domain.category_opportunity_counts();
                    } else {
                        const wgbs::WgbsGcSampler target_domain(
                            contig.bases,
                            target_calibration_shape,
                            *coverage_profile);
                        candidate_weights[contig.index] =
                            target_domain.valid_start_count();
                        target_bin_counts[contig.index] =
                            target_domain.bin_opportunity_counts();
                    }
                } else if (variable_wgbs_insert) {
                    if (planned_variants) {
                        candidate_weights[contig.index] =
                            wgbs::VariableHaplotypeSampler(
                                contig,
                                *planned_variants,
                                config.master_seed,
                                insert_parameters,
                                config.read_length_1,
                                config.paired_end,
                                config.max_ambiguous_fraction)
                                .allocation_weight();
                    } else {
                        candidate_weights[contig.index] =
                            wgbs::VariableWgbsSampler(
                                contig.bases,
                                contig.index,
                                config.master_seed,
                                insert_parameters,
                                config.read_length_1,
                                config.paired_end,
                                config.max_ambiguous_fraction)
                                .allocation_weight();
                    }
                } else if (planned_variants) {
                    candidate_weights[contig.index] =
                        wgbs::HaplotypeStartIndex(
                            contig,
                            *planned_variants,
                            sampler_shape).valid_start_count();
                } else {
                    candidate_weights[contig.index] =
                        wgbs::count_valid_starts(
                            contig.bases, sampler_shape);
                }
            } else if (config.technology == Technology::rrbs) {
                const bool diploid = planned_variants
                    && !planned_variants->variants().empty();
                std::unique_ptr<rrbs::CandidateCatalog> reference_catalog;
                std::unique_ptr<rrbs::DiploidCandidateCatalog> diploid_catalog;
                if (diploid) {
                    diploid_catalog =
                        std::make_unique<rrbs::DiploidCandidateCatalog>(
                            contig,
                            *planned_variants,
                            rrbs_cut_sites,
                            config.insert_min,
                            config.insert_max,
                            config.read_length_1,
                            config.paired_end,
                            config.max_ambiguous_fraction);
                } else {
                    reference_catalog =
                        std::make_unique<rrbs::CandidateCatalog>(
                            contig.bases,
                            rrbs_cut_sites,
                            config.insert_min,
                            config.insert_max,
                            config.read_length_1,
                            config.paired_end,
                            config.max_ambiguous_fraction);
                }
                const auto &candidates = diploid_catalog
                    ? diploid_catalog->candidates()
                    : reference_catalog->candidates();
                target_reference_bases[contig.index] = interval_union_bases(
                    candidates,
                    [](const auto &candidate) {
                        return candidate.reference_start;
                    },
                    [](const auto &candidate) {
                        return candidate.reference_end;
                    });
                std::vector<double> scores;
                if (rrbs_candidate_bed) {
                    scores = rrbs_candidate_bed->match_scores(
                        contig.index, contig.name, candidates, rrbs_profile);
                }
                if (rrbs_profile) {
                    rrbs_profile_weights[contig.index] =
                        rrbs::ProfileSampler(candidates, scores)
                            .allocation_weight();
                } else {
                    candidate_weights[contig.index] = diploid_catalog
                        ? diploid_catalog->allocation_weight()
                        : reference_catalog->allocation_weight();
                }
            } else {
                target_reference_bases[contig.index] = interval_union_bases(
                    tbs_targets->targets(contig.index),
                    [](const auto &target) {return target.interval_start;},
                    [](const auto &target) {return target.interval_end;});
                candidate_weights[contig.index] =
                    planned_variants && !planned_variants->variants().empty()
                    ? tbs::DiploidCandidateCatalog(
                          contig,
                          *planned_variants,
                          tbs_targets->targets(contig.index),
                          *config.tbs_center_stddev,
                          config.insert_mean,
                          config.read_length_1,
                          config.paired_end,
                          config.max_ambiguous_fraction,
                          tbs_sampling_mode(config)).allocation_weight()
                    : tbs::CandidateCatalog(
                          contig.bases,
                          tbs_targets->targets(contig.index),
                          contig.index,
                          *config.tbs_center_stddev,
                          config.insert_mean,
                          config.read_length_1,
                          config.paired_end,
                          config.max_ambiguous_fraction,
                          tbs_sampling_mode(config)).allocation_weight();
            }
        });
        std::optional<wgbs::WgbsGcTargetCalibration> target_calibration;
        if (coverage_profile) {
            target_calibration = haplotype_gc_profile
                ? wgbs::calibrate_haplotype_gc_target(
                      *coverage_profile, target_bin_counts)
                : wgbs::calibrate_gc_target(
                      *coverage_profile,
                      target_bin_counts,
                      variable_wgbs_insert
                          ? wgbs::UnreachableTargetPolicy::drop_and_renormalize
                          : wgbs::UnreachableTargetPolicy::reject);
        }

        std::uint32_t requested_fragment_count = 0U;
        if (config.fragment_count) {
            requested_fragment_count = *config.fragment_count;
        } else {
            std::uint64_t effective_reference_bases = 0U;
            for (std::size_t index = 0U; index < catalog.size(); ++index) {
                std::uint64_t contribution = target_reference_bases[index];
                if (whole_genome_technology(config.technology)) {
                    const bool eligible = target_calibration
                        ? target_calibration->contig_allocation_weights[index] > 0.0
                        : candidate_weights[index] > 0U;
                    contribution = eligible ? catalog[index].length : 0U;
                }
                effective_reference_bases = checked_add(
                    effective_reference_bases,
                    contribution,
                    "target reference");
            }
            requested_fragment_count = depth_count::fragments(
                *config.depth,
                effective_reference_bases,
                config.read_length_1,
                config.paired_end);
        }
        const std::vector<std::uint32_t> fragment_counts = target_calibration
            ? allocation::largest_remainder_real(
                  target_calibration->contig_allocation_weights,
                  requested_fragment_count)
            : (rrbs_profile
                ? allocation::largest_remainder_real(
                      rrbs_profile_weights, requested_fragment_count)
                : allocation::largest_remainder(
                      candidate_weights, requested_fragment_count));

        const protocol::Header header = make_header(config, catalog);
        protocol::Writer writer(sink);
        writer.write_header(header);
        protocol::BatchEmitter emitter(
            writer,
            header,
            config.core_workers,
            config.protocol_batch_fragments);

        std::uint64_t fragment_ordinal = 0;
        std::uint64_t skipped_fragment_count = 0;
        const fragment_builder::FragmentDetail fragment_detail =
            !config.emit_details
            ? fragment_builder::FragmentDetail::common_columns
            : fragment_builder::FragmentDetail::full;
        snapshot.visit_contigs([&](const reference::Contig &contig) {
            const std::uint32_t requested = fragment_counts.at(contig.index);
            if (requested == 0) {return;}
            std::vector<methdb::CgmapRecord> cgmap_records;
            if (cgmap_profile) {
                cgmap_records = cgmap_profile->records(contig);
            }
            const auto *cgmap_records_or_null = cgmap_profile
                ? &cgmap_records
                : nullptr;
            std::vector<methdb::AsmRecord> asm_records;
            if (asm_profile) {
                asm_records = asm_profile->records(contig);
            }
            const auto *asm_records_or_null = asm_profile
                ? &asm_records
                : nullptr;
            std::vector<variant::Variant> generated_events;
            const std::vector<variant::Variant> *variants = nullptr;
            if (variant_file) {
                variants = &variant_file->variants(contig.index);
            } else if (generate_mutations) {
                generated_events = variant::generate_de_novo_events(
                    contig, config.mutation_seed, mutation_parameters);
                variants = &generated_events;
            }
            std::unique_ptr<variant::ContigVariants> contig_variants;
            std::unique_ptr<methdb::DiploidMethylationCatalog>
                diploid_methylation_catalog;
            std::unique_ptr<methdb::MethylationCatalog>
                reference_methylation_catalog;
            if (variants != nullptr) {
                contig_variants = std::make_unique<variant::ContigVariants>(
                    contig.bases,
                    *variants,
                    contig.index);
                if (!bisulfite_technology(config.technology)) {
                    diploid_methylation_catalog =
                        std::make_unique<methdb::DiploidMethylationCatalog>(
                            contig.index,
                            static_cast<std::uint32_t>(contig.length),
                            std::vector<methdb::DiploidSite>{},
                            std::array<std::vector<methdb::DiploidSite>, 2>{});
                } else if (fixed_methdb) {
                    const auto &saved = fixed_methdb->contig(contig.index);
                    if (!saved.diploid) {
                        throw CoreGeneratorError(
                            "MethDB contig is not diploid for the prepared "
                            "variant set");
                    }
                    diploid_methylation_catalog =
                        std::make_unique<methdb::DiploidMethylationCatalog>(
                            contig.index,
                            static_cast<std::uint32_t>(contig.length),
                            saved.shared_sites,
                            saved.haplotype_sites);
                } else {
                    diploid_methylation_catalog =
                        std::make_unique<methdb::DiploidMethylationCatalog>(
                            contig,
                            *contig_variants,
                            config.methylation_seed,
                            config.collect_non_cpg,
                            context_shapes,
                            cgmap_records_or_null,
                            asm_records_or_null,
                            config.cgmap_pool);
                }
            } else {
                if (!bisulfite_technology(config.technology)) {
                    reference_methylation_catalog =
                        std::make_unique<methdb::MethylationCatalog>(
                            static_cast<std::uint32_t>(contig.length),
                            std::vector<methdb::CatalogSite>{});
                } else if (fixed_methdb) {
                    const auto &saved = fixed_methdb->contig(contig.index);
                    if (saved.diploid) {
                        throw CoreGeneratorError(
                            "MethDB contig is diploid without a prepared "
                            "variant set");
                    }
                    reference_methylation_catalog =
                        std::make_unique<methdb::MethylationCatalog>(
                            static_cast<std::uint32_t>(contig.length),
                            saved.reference_sites);
                } else {
                    reference_methylation_catalog =
                        std::make_unique<methdb::MethylationCatalog>(
                            contig.bases,
                            contig.index,
                            config.methylation_seed,
                            config.collect_non_cpg,
                            context_shapes,
                            cgmap_records_or_null,
                            config.cgmap_pool);
                }
            }
            const std::uint64_t haplotype_key = rng::derive_key(
                config.master_seed, rng::Stage::haplotype, contig.index);
            const std::uint64_t library_orientation_key = rng::derive_key(
                config.master_seed,
                rng::Stage::library_orientation,
                contig.index);
            const auto build_fragment = [&](
                std::uint32_t start,
                std::uint8_t haplotype,
                model::CaptureStrand capture_strand,
                const fragment_builder::ReadLayout &layout) {
                const LibraryOrientation orientation = sample_library_orientation(
                    config,
                    library_orientation_key,
                    fragment_ordinal,
                    capture_strand);
                if (!contig_variants) {
                    return fragment_builder::build_fragment(
                        contig,
                        *reference_methylation_catalog,
                        fragment_ordinal,
                        start,
                        haplotype,
                        orientation.informative_strand,
                        layout,
                        fragment_detail,
                        orientation.reverse_molecule);
                }
                if (layout.insert_length
                    > std::numeric_limits<std::uint32_t>::max() - start) {
                    throw CoreGeneratorError(
                        "variant fragment reference interval exceeds uint32");
                }
                auto projection = haplotype::project_interval(
                    contig,
                    *contig_variants,
                    haplotype,
                    start,
                    start + layout.insert_length);
                return fragment_builder::build_fragment(
                    std::move(projection),
                    *diploid_methylation_catalog,
                    fragment_ordinal,
                    orientation.informative_strand,
                    layout,
                    fragment_detail,
                    orientation.reverse_molecule);
            };
            const auto build_haplotype_fragment = [&](
                std::uint32_t reference_start,
                std::uint32_t reference_end,
                bool include_start_anchor_insertion,
                bool include_end_anchor_insertion,
                std::uint8_t haplotype,
                model::CaptureStrand capture_strand,
                const fragment_builder::ReadLayout &layout) {
                if (!contig_variants || !diploid_methylation_catalog) {
                    throw CoreGeneratorError(
                        "haplotype fragment lost its prepared variant set");
                }
                auto projection = haplotype::project_interval(
                    contig,
                    *contig_variants,
                    haplotype,
                    reference_start,
                    reference_end,
                    haplotype::ProjectionBoundaryPolicy{
                        include_start_anchor_insertion,
                        include_end_anchor_insertion,
                    });
                const LibraryOrientation orientation = sample_library_orientation(
                    config,
                    library_orientation_key,
                    fragment_ordinal,
                    capture_strand);
                return fragment_builder::build_fragment(
                    std::move(projection),
                    *diploid_methylation_catalog,
                    fragment_ordinal,
                    orientation.informative_strand,
                    layout,
                    fragment_detail,
                    orientation.reverse_molecule);
            };

            std::uint32_t emitted_for_contig = 0;
            std::uint64_t candidate_ordinal = 0;
            if (whole_genome_technology(config.technology)) {
                std::unique_ptr<wgbs::VariableHaplotypeSampler>
                    variable_haplotype_starts;
                std::unique_ptr<wgbs::VariableWgbsSampler>
                    variable_starts;
                std::unique_ptr<wgbs::VariableWgbsGcSampler>
                    profiled_variable_starts;
                std::unique_ptr<wgbs::ValidStartIndex> uniform_starts;
                std::unique_ptr<wgbs::HaplotypeStartIndex> variant_starts;
                std::unique_ptr<wgbs::WgbsGcSampler> profiled_starts;
                std::unique_ptr<wgbs::HaplotypeGcSampler>
                    profiled_haplotype_starts;
                if (coverage_profile) {
                    if (!target_calibration) {
                        throw CoreGeneratorError(
                            "target GC profile escaped its generation gate");
                    }
                    if (variable_wgbs_insert) {
                        if (contig_variants) {
                            throw CoreGeneratorError(
                                "variable target GC escaped its variant gate");
                        }
                        profiled_variable_starts = std::make_unique<
                            wgbs::VariableWgbsGcSampler>(
                            contig.bases,
                            contig.index,
                            config.master_seed,
                            insert_parameters,
                            config.read_length_1,
                            config.paired_end,
                            config.max_ambiguous_fraction,
                            *coverage_profile);
                        if (profiled_variable_starts->allocation_weight()
                            != candidate_weights[contig.index]
                            || candidate_weights[contig.index] == 0U) {
                            throw CoreGeneratorError(
                                "profiled variable WGBS domain changed between planning and generation");
                        }
                    } else if (contig_variants) {
                        profiled_haplotype_starts = std::make_unique<
                            wgbs::HaplotypeGcSampler>(
                            contig,
                            *contig_variants,
                            sampler_shape,
                            *coverage_profile);
                        if (profiled_haplotype_starts->physical_candidate_count()
                                != candidate_weights[contig.index]
                            || profiled_haplotype_starts
                                   ->category_opportunity_counts()
                                != target_bin_counts[contig.index]) {
                            throw CoreGeneratorError(
                                "haplotype target GC domain changed between planning and generation");
                        }
                    } else {
                        profiled_starts =
                            std::make_unique<wgbs::WgbsGcSampler>(
                                contig.bases,
                                sampler_shape,
                                *coverage_profile);
                        if (profiled_starts->valid_start_count()
                                != candidate_weights[contig.index]
                            || profiled_starts->bin_opportunity_counts()
                                != target_bin_counts[contig.index]) {
                            throw CoreGeneratorError(
                                "target GC opportunity domain changed between planning and generation");
                        }
                    }
                } else if (variable_wgbs_insert) {
                    if (contig_variants) {
                        variable_haplotype_starts = std::make_unique<
                            wgbs::VariableHaplotypeSampler>(
                            contig,
                            *contig_variants,
                            config.master_seed,
                            insert_parameters,
                            config.read_length_1,
                            config.paired_end,
                            config.max_ambiguous_fraction);
                    } else {
                        variable_starts =
                            std::make_unique<wgbs::VariableWgbsSampler>(
                            contig.bases,
                            contig.index,
                            config.master_seed,
                            insert_parameters,
                            config.read_length_1,
                            config.paired_end,
                            config.max_ambiguous_fraction);
                    }
                    std::uint32_t observed_weight = 0U;
                    if (variable_haplotype_starts) {
                        observed_weight =
                            variable_haplotype_starts->allocation_weight();
                    } else {
                        observed_weight = variable_starts->allocation_weight();
                    }
                    if (observed_weight
                        != candidate_weights[contig.index]) {
                        throw CoreGeneratorError(
                            "variable WGBS domain changed between planning and generation");
                    }
                } else if (contig_variants) {
                    variant_starts =
                        std::make_unique<wgbs::HaplotypeStartIndex>(
                            contig, *contig_variants, sampler_shape);
                    if (variant_starts->valid_start_count()
                        != candidate_weights[contig.index]) {
                        throw CoreGeneratorError(
                            "variant start domain changed between planning and generation");
                    }
                } else {
                    uniform_starts =
                        std::make_unique<wgbs::ValidStartIndex>(
                            contig.bases, sampler_shape);
                    if (uniform_starts->valid_start_count()
                        != candidate_weights[contig.index]) {
                        throw CoreGeneratorError(
                            "valid-start count changed between planning and generation");
                    }
                }
                const fragment_builder::ReadLayout fixed_read_layout{
                    config.insert_mean,
                    config.read_length_1,
                    config.paired_end};
                while (emitted_for_contig < requested) {
                    const std::uint32_t remaining =
                        requested - emitted_for_contig;
                    const std::uint32_t chunk =
                        std::min(remaining, config.chunk_size);
                    if (variable_haplotype_starts) {
                        const auto batch = variable_haplotype_starts->sample(
                            candidate_ordinal, chunk);
                        if (batch.skipped_count
                            > std::numeric_limits<std::uint64_t>::max()
                                - skipped_fragment_count) {
                            throw CoreGeneratorError(
                                "variable haplotype skipped fragment count exceeds uint64");
                        }
                        skipped_fragment_count += batch.skipped_count;
                        for (const wgbs::VariableHaplotypeCandidate &candidate
                             : batch.candidates) {
                            const fragment_builder::ReadLayout read_layout{
                                candidate.reference_span,
                                config.read_length_1,
                                config.paired_end};
                            emitter.write(build_fragment(
                                candidate.reference_start,
                                choose_haplotype(
                                    candidate.eligible_haplotypes,
                                    haplotype_key,
                                    fragment_ordinal),
                                model::CaptureStrand::unknown,
                                read_layout));
                            ++fragment_ordinal;
                        }
                        candidate_ordinal = batch.next_candidate_ordinal;
                    } else if (profiled_variable_starts || variable_starts) {
                        const auto batch = profiled_variable_starts
                            ? profiled_variable_starts->sample(
                                  candidate_ordinal,
                                  chunk,
                                  target_calibration->acceptance_probabilities)
                            : variable_starts->sample(candidate_ordinal, chunk);
                        if (batch.skipped_count
                            > std::numeric_limits<std::uint64_t>::max()
                                - skipped_fragment_count) {
                            throw CoreGeneratorError(
                                "variable WGBS skipped fragment count exceeds uint64");
                        }
                        skipped_fragment_count += batch.skipped_count;
                        for (const wgbs::VariableWgbsCandidate &candidate
                             : batch.candidates) {
                            const fragment_builder::ReadLayout read_layout{
                                candidate.insert_length,
                                config.read_length_1,
                                config.paired_end};
                            emitter.write(build_fragment(
                                candidate.reference_start,
                                choose_haplotype(
                                    model::HaplotypeMask::both,
                                    haplotype_key,
                                    fragment_ordinal),
                                model::CaptureStrand::unknown,
                                read_layout));
                            ++fragment_ordinal;
                        }
                        candidate_ordinal = batch.next_candidate_ordinal;
                    } else if (profiled_haplotype_starts) {
                        const auto batch = profiled_haplotype_starts->sample(
                            contig.index,
                            config.master_seed,
                            candidate_ordinal,
                            chunk,
                            target_calibration->acceptance_probabilities);
                        if (batch.skipped_count
                            > std::numeric_limits<std::uint64_t>::max()
                                - skipped_fragment_count) {
                            throw CoreGeneratorError(
                                "haplotype coverage skipped fragment count exceeds uint64");
                        }
                        skipped_fragment_count += batch.skipped_count;
                        for (const wgbs::HaplotypeCandidate &candidate
                             : batch.candidates) {
                            emitter.write(build_fragment(
                                candidate.reference_start,
                                candidate.haplotype,
                                model::CaptureStrand::unknown,
                                fixed_read_layout));
                            ++fragment_ordinal;
                        }
                        candidate_ordinal += chunk;
                    } else {
                        std::vector<std::uint32_t> starts;
                        if (profiled_starts) {
                            auto batch = profiled_starts->sample(
                                contig.index,
                                config.master_seed,
                                candidate_ordinal,
                                chunk,
                                target_calibration->acceptance_probabilities);
                            if (batch.skipped_count
                                > std::numeric_limits<std::uint64_t>::max()
                                    - skipped_fragment_count) {
                                throw CoreGeneratorError(
                                    "coverage skipped fragment count exceeds uint64");
                            }
                            skipped_fragment_count += batch.skipped_count;
                            starts = std::move(batch.starts);
                        } else if (variant_starts) {
                            starts = variant_starts->sample(
                                contig.index,
                                config.master_seed,
                                candidate_ordinal,
                                chunk);
                        } else {
                            starts = uniform_starts->sample(
                                contig.index,
                                config.master_seed,
                                candidate_ordinal,
                                chunk);
                        }
                        for (const std::uint32_t start : starts) {
                            emitter.write(build_fragment(
                                start,
                                choose_haplotype(
                                    variant_starts
                                        ? variant_starts->haplotype_mask(start)
                                        : model::HaplotypeMask::both,
                                    haplotype_key,
                                    fragment_ordinal),
                                model::CaptureStrand::unknown,
                                fixed_read_layout));
                            ++fragment_ordinal;
                        }
                        candidate_ordinal += chunk;
                    }
                    emitted_for_contig += chunk;
                }
            } else if (config.technology == Technology::rrbs) {
                const bool diploid_rrbs = contig_variants
                    && !contig_variants->variants().empty();
                std::unique_ptr<rrbs::CandidateCatalog> reference_catalog;
                std::unique_ptr<rrbs::DiploidCandidateCatalog> diploid_catalog;
                if (diploid_rrbs) {
                        diploid_catalog =
                            std::make_unique<rrbs::DiploidCandidateCatalog>(
                            contig,
                            *contig_variants,
                            rrbs_cut_sites,
                            config.insert_min,
                            config.insert_max,
                            config.read_length_1,
                            config.paired_end,
                            config.max_ambiguous_fraction);
                } else {
                    reference_catalog = std::make_unique<rrbs::CandidateCatalog>(
                        contig.bases,
                        rrbs_cut_sites,
                        config.insert_min,
                        config.insert_max,
                        config.read_length_1,
                        config.paired_end,
                        config.max_ambiguous_fraction);
                }
                const auto &candidates = diploid_catalog
                    ? diploid_catalog->candidates()
                    : reference_catalog->candidates();
                std::vector<double> scores;
                if (rrbs_candidate_bed) {
                    scores = rrbs_candidate_bed->match_scores(
                        contig.index, contig.name, candidates, rrbs_profile);
                }
                std::unique_ptr<rrbs::ProfileSampler> profile_sampler;
                if (rrbs_profile) {
                    profile_sampler = std::make_unique<rrbs::ProfileSampler>(
                        candidates, scores);
                    if (profile_sampler->allocation_weight()
                        != rrbs_profile_weights[contig.index]) {
                        throw CoreGeneratorError(
                            "RRBS profile weight changed between planning and generation");
                    }
                } else {
                    const std::uint32_t observed_allocation_weight =
                        diploid_catalog
                        ? diploid_catalog->allocation_weight()
                        : reference_catalog->allocation_weight();
                    if (observed_allocation_weight
                        != candidate_weights[contig.index]) {
                        throw CoreGeneratorError(
                            "RRBS candidate set changed between planning "
                            "and generation");
                    }
                }
                while (emitted_for_contig < requested) {
                    const std::uint32_t remaining =
                        requested - emitted_for_contig;
                    const std::uint32_t chunk =
                        std::min(remaining, config.chunk_size);
                    std::vector<std::uint32_t> indices;
                    if (profile_sampler) {
                        indices = profile_sampler->sample_indices(
                            contig.index,
                            config.master_seed,
                            candidate_ordinal,
                            chunk);
                    } else if (diploid_catalog) {
                        indices = diploid_catalog->sample_indices(
                            contig.index,
                            config.master_seed,
                            candidate_ordinal,
                            chunk);
                    } else {
                        indices = reference_catalog->sample_indices(
                            contig.index,
                            config.master_seed,
                            candidate_ordinal,
                            chunk);
                    }
                    for (const std::uint32_t index : indices) {
                        const rrbs::Candidate &candidate = diploid_catalog
                            ? diploid_catalog->candidate(index)
                            : reference_catalog->candidate(index);
                        const fragment_builder::ReadLayout read_layout{
                            candidate.template_length,
                            config.read_length_1,
                            config.paired_end,
                            diploid_catalog
                                ? fragment_builder::ReadLayout::InsertCoordinate::haplotype
                                : fragment_builder::ReadLayout::InsertCoordinate::reference,
                        };
                        const std::uint8_t selected_haplotype = choose_haplotype(
                            candidate.haplotype_mask,
                            haplotype_key,
                            fragment_ordinal);
                        emitter.write(diploid_catalog
                            ? build_haplotype_fragment(
                                  candidate.reference_start,
                                  candidate.reference_end,
                                  candidate.include_start_anchor_insertion,
                                  candidate.include_end_anchor_insertion,
                                  selected_haplotype,
                                  model::CaptureStrand::unknown,
                                  read_layout)
                            : build_fragment(
                                  candidate.reference_start,
                                  selected_haplotype,
                                  model::CaptureStrand::unknown,
                                  read_layout));
                        ++fragment_ordinal;
                    }
                    emitted_for_contig += chunk;
                    candidate_ordinal += chunk;
                }
            } else {
                const bool diploid_tbs = contig_variants
                    && !contig_variants->variants().empty();
                std::unique_ptr<tbs::CandidateCatalog> reference_catalog;
                std::unique_ptr<tbs::DiploidCandidateCatalog> diploid_catalog;
                if (diploid_tbs) {
                        diploid_catalog =
                            std::make_unique<tbs::DiploidCandidateCatalog>(
                            contig,
                            *contig_variants,
                            tbs_targets->targets(contig.index),
                            *config.tbs_center_stddev,
                            config.insert_mean,
                            config.read_length_1,
                            config.paired_end,
                            config.max_ambiguous_fraction,
                            tbs_sampling_mode(config));
                } else {
                    reference_catalog = std::make_unique<tbs::CandidateCatalog>(
                        contig.bases,
                        tbs_targets->targets(contig.index),
                        contig.index,
                        *config.tbs_center_stddev,
                        config.insert_mean,
                        config.read_length_1,
                        config.paired_end,
                        config.max_ambiguous_fraction,
                        tbs_sampling_mode(config));
                }
                const std::uint32_t observed_allocation_weight = diploid_catalog
                    ? diploid_catalog->allocation_weight()
                    : reference_catalog->allocation_weight();
                if (observed_allocation_weight
                    != candidate_weights[contig.index]) {
                    throw CoreGeneratorError(
                        "TBS target set changed between planning and "
                        "generation");
                }
                while (emitted_for_contig < requested) {
                    const std::uint32_t remaining =
                        requested - emitted_for_contig;
                    const std::uint32_t chunk =
                        std::min(remaining, config.chunk_size);
                    const auto batch = diploid_catalog
                        ? diploid_catalog->sample(
                              contig.index,
                              config.master_seed,
                              candidate_ordinal,
                              fragment_ordinal,
                              chunk)
                        : reference_catalog->sample(
                              contig.index,
                              config.master_seed,
                              candidate_ordinal,
                              chunk);
                    if (batch.skipped_count
                        > std::numeric_limits<std::uint64_t>::max()
                            - skipped_fragment_count) {
                        throw CoreGeneratorError(
                            "TBS skipped fragment count exceeds uint64");
                    }
                    skipped_fragment_count += batch.skipped_count;
                    for (const tbs::Candidate &candidate : batch.candidates) {
                        const fragment_builder::ReadLayout read_layout{
                            candidate.template_length,
                            config.read_length_1,
                            config.paired_end,
                            diploid_catalog
                                ? fragment_builder::ReadLayout::InsertCoordinate::haplotype
                                : fragment_builder::ReadLayout::InsertCoordinate::reference,
                        };
                        const std::uint8_t selected_haplotype = choose_haplotype(
                            candidate.haplotype_mask,
                            haplotype_key,
                            fragment_ordinal);
                        emitter.write(diploid_catalog
                            ? build_haplotype_fragment(
                                  candidate.reference_start,
                                  candidate.reference_end,
                                  candidate.include_start_anchor_insertion,
                                  candidate.include_end_anchor_insertion,
                                  selected_haplotype,
                                  candidate.capture_strand,
                                  read_layout)
                            : build_fragment(
                                  candidate.reference_start,
                                  selected_haplotype,
                                  candidate.capture_strand,
                                  read_layout));
                        ++fragment_ordinal;
                    }
                    emitted_for_contig += chunk;
                    candidate_ordinal += chunk;
                }
            }
        });
        if (fragment_ordinal
            != static_cast<std::uint64_t>(requested_fragment_count)) {
            throw CoreGeneratorError("generated fragment count disagrees with plan");
        }
        emitter.finish();
        return writer.finish(skipped_fragment_count);
    } catch (const CoreGeneratorError &) {
        throw;
    } catch (const std::exception &error) {
        throw CoreGeneratorError(error.what());
    }
}

} // namespace htsim::core
