#include "core.h"

#include <cstdlib>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

using htsim::core::CoreConfigError;
using htsim::core::CoverageMode;
using htsim::core::Technology;
using htsim::core::parse_core_config;
using htsim::core::validate_core_config;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Callable>
void require_error(Callable &&callable, const std::string &message)
{
    try {
        callable();
    } catch (const CoreConfigError &) {
        return;
    }
    throw std::runtime_error(message);
}

std::vector<std::string> base_arguments()
{
    return {
        "--run-id", "12345678-1234-4234-8234-123456789abc",
        "--config-sha256", std::string(64, '1'),
        "--seed", "0",
        "--seed-mut", "1",
        "--seed-phase", "2",
        "--seed-meth", "3",
        "--reference", "/data/reference.fa",
        "--technology", "WGBS",
        "--paired-end", "true",
        "--read-length-1", "100",
        "--read-length-2", "100",
        "--insert-min", "100",
        "--insert-mean", "400",
        "--insert-max", "1000",
        "--insert-sd", "25",
        "--fragments", "20",
        "--max-ambiguous-fraction", "0.05",
        "--chunk-size", "10000",
        "--mutation-rate", "0.001",
        "--indel-fraction", "0.15",
        "--indel-extension-probability", "0.15",
        "--homozygous-only", "false",
        "--collect-non-cpg", "true",
        "--cgmap-pool", "false",
        "--update-variant-boundaries", "true",
        "--beta-cg", "0.5,0.5",
        "--beta-chg", "0.05,0.10",
        "--beta-chh", "0.01,0.20",
        "--coverage", "uniform",
    };
}

void replace_value(
    std::vector<std::string> &arguments,
    const std::string &option,
    const std::string &value)
{
    for (std::size_t index = 0; index + 1 < arguments.size(); index += 2) {
        if (arguments[index] == option) {
            arguments[index + 1] = value;
            return;
        }
    }
    throw std::runtime_error("test option not found: " + option);
}

void remove_option(
    std::vector<std::string> &arguments,
    const std::string &option)
{
    for (auto iterator = arguments.begin(); iterator != arguments.end(); iterator += 2) {
        if (*iterator == option) {
            arguments.erase(iterator, iterator + 2);
            return;
        }
    }
    throw std::runtime_error("test option not found: " + option);
}

void test_valid_wgbs_projection()
{
    const auto config = parse_core_config(base_arguments());
    require(!config.emit_details,
            "direct invocation did not default to no-Details");
    require(config.technology == Technology::wgbs, "WGBS technology changed");
    require(config.master_seed == 0, "seed zero was not preserved");
    require(config.paired_end && config.read_length_2 == 100,
            "paired-end projection was not parsed");
    require(config.fragment_count == 20 && !config.depth.has_value(),
            "read-count mode was not parsed");
    require(config.coverage == CoverageMode::uniform,
            "uniform coverage was not parsed");
    require(config.core_workers == 1,
            "core worker default changed");
    require(config.beta_chh.alpha == 0.01 && config.beta_chh.beta == 0.20,
            "beta shape projection changed");
    require(config.normalized_config_sha256.front() == 0x11,
            "config digest was not decoded");
}

void test_output_controls()
{
    auto arguments = base_arguments();
    arguments.insert(
        arguments.begin(),
        {"--emit-details", "true", "--protocol-batch-fragments", "7"});
    const auto config = parse_core_config(arguments);
    require(config.emit_details,
            "Full Details selection was not parsed");
    require(config.protocol_batch_fragments == 7,
            "protocol batch bound was not parsed");

    auto retired_selection = arguments;
    retired_selection.insert(
        retired_selection.begin(), {"--protocol-major", "1"});
    require_error(
        [&] {parse_core_config(retired_selection);},
        "retired protocol selector was accepted");

    for (const std::string value : {"", "TRUE", "details"}) {
        auto invalid = arguments;
        replace_value(invalid, "--emit-details", value);
        require_error(
            [&] {parse_core_config(invalid);},
            "invalid Details-column mode was accepted: " + value);
    }
    for (const std::string value : {"0", "65", "-1"}) {
        auto invalid = arguments;
        replace_value(invalid, "--protocol-batch-fragments", value);
        require_error(
            [&] {parse_core_config(invalid);},
            "invalid protocol batch bound was accepted: " + value);
    }

}

void test_valid_optional_input_projection()
{
    auto arguments = base_arguments();
    replace_value(arguments, "--mutation-rate", "0");
    arguments.insert(
        arguments.end(),
        {"--vcf", "/data/sample.vcf.gz",
         "--cgmap", "/data/sample.cgmap.gz",
         "--asm", "/data/sample.asm.gz"});
    const auto config = parse_core_config(arguments);
    require(config.vcf_path == "/data/sample.vcf.gz",
            "VCF path projection was lost");
    require(config.cgmap_path == "/data/sample.cgmap.gz",
            "CGmap path projection was lost");
    require(config.asm_path == "/data/sample.asm.gz",
            "ASM path projection was lost");

    auto bed_arguments = base_arguments();
    replace_value(bed_arguments, "--mutation-rate", "0");
    bed_arguments.insert(
        bed_arguments.end(),
        {"--vcf", "/data/sample.vcf.gz",
         "--bed-methyl", "/data/sample.bedmethyl.gz",
         "--asm-bed", "/data/sample.asm.bed.gz"});
    const auto bed_config = parse_core_config(bed_arguments);
    require(
        bed_config.bed_methyl_path == "/data/sample.bedmethyl.gz"
            && bed_config.asm_bed_path == "/data/sample.asm.bed.gz",
        "bedMethyl/ASM BED path projection was lost");
}

void test_cgmap_pool_requires_its_input()
{
    auto missing = base_arguments();
    replace_value(missing, "--cgmap-pool", "true");
    require_error(
        [&] {parse_core_config(missing);},
        "CGmap pooling without a CGmap input was accepted");

    missing.insert(
        missing.end(),
        {"--cgmap", "/data/sample.cgmap.gz"});
    const auto pooled = parse_core_config(missing);
    require(pooled.cgmap_pool && pooled.cgmap_path.has_value(),
            "valid CGmap pooling projection was lost");

    auto bed = base_arguments();
    replace_value(bed, "--cgmap-pool", "true");
    bed.insert(
        bed.end(),
        {"--bed-methyl", "/data/sample.bedmethyl"});
    const auto bed_pooled = parse_core_config(bed);
    require(bed_pooled.cgmap_pool && bed_pooled.bed_methyl_path.has_value(),
            "valid bedMethyl pooling projection was lost");
}

void test_valid_rrbs_tbs_and_profile_projection()
{
    auto rrbs_arguments = base_arguments();
    replace_value(rrbs_arguments, "--technology", "RRBS");
    rrbs_arguments.insert(
        rrbs_arguments.end(),
        {"--rrbs-cut-site", "C|CGG", "--rrbs-cut-site", "CCTN|AGG"});
    const auto rrbs = parse_core_config(rrbs_arguments);
    require(rrbs.rrbs_cut_sites.size() == 2, "repeated RRBS cut sites were lost");
    rrbs_arguments.insert(
        rrbs_arguments.end(),
        {"--rrbs-candidate-bed", "/data/rrbs-candidates.bed"});
    replace_value(rrbs_arguments, "--coverage", "profile");
    const auto rrbs_profile = parse_core_config(rrbs_arguments);
    require(rrbs_profile.coverage == CoverageMode::profile
                && rrbs_profile.rrbs_candidate_bed_path
                    == "/data/rrbs-candidates.bed",
            "RRBS candidate profile projection was lost");

    auto tbs_arguments = base_arguments();
    replace_value(tbs_arguments, "--technology", "TBS");
    replace_value(tbs_arguments, "--insert-mean", "300");
    replace_value(tbs_arguments, "--insert-sd", "0");
    tbs_arguments.insert(
        tbs_arguments.end(),
        {"--tbs-bed", "/data/targets.bed",
         "--tbs-center-stddev", "50"});
    const auto tbs = parse_core_config(tbs_arguments);
    require(tbs.tbs_bed_path == "/data/targets.bed"
                && tbs.tbs_center_stddev == 50.0
                && tbs.insert_min == 100U
                && tbs.insert_mean == 300U
                && tbs.insert_max == 1000U
                && tbs.insert_sd == 0.0,
            "TBS projection was not parsed");
    replace_value(tbs_arguments, "--coverage", "target-score");
    const auto target_score = parse_core_config(tbs_arguments);
    require(target_score.coverage == CoverageMode::target_score,
            "TBS target-score coverage projection was lost");

    auto profile_arguments = base_arguments();
    replace_value(profile_arguments, "--coverage", "profile");
    replace_value(profile_arguments, "--mutation-rate", "0");
    profile_arguments.insert(
        profile_arguments.end(),
        {"--coverage-profile", "/data/coverage.tsv"});
    const auto profile = parse_core_config(profile_arguments);
    require(profile.coverage_profile_path == "/data/coverage.tsv",
            "coverage profile path was lost");
    require(profile.insert_min == 100U
                && profile.insert_mean == 400U
                && profile.insert_max == 1000U
                && profile.insert_sd == 25.0,
            "profile coverage did not retain variable insert parameters");
}

void test_standard_technology_projection()
{
    auto wgs_arguments = base_arguments();
    replace_value(wgs_arguments, "--technology", "WGS");
    const auto wgs = parse_core_config(wgs_arguments);
    require(wgs.technology == Technology::wgs,
            "WGS technology projection was lost");

    for (const auto &[name, expected] : {
             std::pair<std::string, Technology>{"WES", Technology::wes},
             std::pair<std::string, Technology>{"TS", Technology::ts}}) {
        auto arguments = base_arguments();
        replace_value(arguments, "--technology", name);
        replace_value(arguments, "--insert-sd", "0");
        arguments.insert(
            arguments.end(),
            {"--tbs-bed", "/data/targets.bed",
             "--tbs-center-stddev", "50"});
        require(parse_core_config(arguments).technology == expected,
                name + " technology projection was lost");
    }

    auto methylation_input = wgs_arguments;
    methylation_input.insert(
        methylation_input.end(), {"--cgmap", "/data/sample.cgmap"});
    require_error(
        [&] {parse_core_config(methylation_input);},
        "standard sequencing accepted a methylation input");
}

void test_unknown_duplicate_and_missing_options_fail()
{
    auto unknown = base_arguments();
    unknown.insert(unknown.end(), {"--not-real", "value"});
    require_error([&] {parse_core_config(unknown);}, "unknown option was accepted");

    auto duplicate = base_arguments();
    duplicate.insert(duplicate.end(), {"--seed", "1"});
    require_error([&] {parse_core_config(duplicate);}, "duplicate option was accepted");

    auto missing = base_arguments();
    remove_option(missing, "--reference");
    require_error([&] {parse_core_config(missing);}, "missing option was accepted");

    auto no_value = base_arguments();
    no_value.push_back("--vcf");
    require_error([&] {parse_core_config(no_value);}, "missing option value was accepted");
}

void test_core_worker_projection_and_boundaries()
{
    auto configured = base_arguments();
    configured.insert(configured.end(), {"--core-workers", "4"});
    require(parse_core_config(configured).core_workers == 4,
            "core worker projection was lost");

    for (const std::string value : {"0", "65"}) {
        auto invalid = base_arguments();
        invalid.insert(invalid.end(), {"--core-workers", value});
        require_error(
            [&] {parse_core_config(invalid);},
            "invalid core worker count was accepted");
    }
}

void test_number_boolean_and_identity_boundaries()
{
    const std::vector<std::pair<std::string, std::string>> invalid = {
        {"--seed", "-1"},
        {"--seed", "18446744073709551616"},
        {"--paired-end", "1"},
        {"--mutation-rate", "nan"},
        {"--mutation-rate", "+0.1"},
        {"--mutation-rate", "0x1p-1"},
        {"--mutation-rate", "01e-2"},
        {"--indel-fraction", "1.1"},
        {"--insert-sd", "-0.1"},
        {"--chunk-size", "0"},
        {"--config-sha256", std::string(64, 'A')},
        {"--run-id", "not-a-uuid"},
    };
    for (const auto &[option, value] : invalid) {
        auto arguments = base_arguments();
        replace_value(arguments, option, value);
        require_error(
            [&] {parse_core_config(arguments);},
            "invalid boundary was accepted for " + option);
    }

    auto nil_uuid = base_arguments();
    replace_value(nil_uuid, "--run-id", "00000000-0000-0000-0000-000000000000");
    require(
        parse_core_config(nil_uuid).run_id ==
            "00000000-0000-0000-0000-000000000000",
        "canonical UUID was rejected");
}

void test_cross_field_boundaries()
{
    auto read_longer_than_minimum_insert = base_arguments();
    replace_value(read_longer_than_minimum_insert, "--insert-min", "99");
    require_error(
        [&] {parse_core_config(read_longer_than_minimum_insert);},
        "read length greater than the minimum insert was accepted");

    auto fixed_insert = base_arguments();
    replace_value(fixed_insert, "--read-length-1", "150");
    replace_value(fixed_insert, "--read-length-2", "150");
    replace_value(fixed_insert, "--insert-mean", "300");
    replace_value(fixed_insert, "--insert-sd", "0");
    const auto fixed = parse_core_config(fixed_insert);
    require(fixed.read_length_1 == 150U && fixed.insert_mean == 300U,
            "fixed mean insert incorrectly used insert_min as its read boundary");

    auto asymmetric = base_arguments();
    replace_value(asymmetric, "--read-length-2", "101");
    require_error(
        [&] {parse_core_config(asymmetric);},
        "unsupported asymmetric paired reads were accepted");

    auto both_counts = base_arguments();
    both_counts.insert(both_counts.end(), {"--depth", "10"});
    require_error(
        [&] {parse_core_config(both_counts);},
        "depth and fragments were accepted together");

    auto asm_without_vcf = base_arguments();
    asm_without_vcf.insert(
        asm_without_vcf.end(),
        {"--asm", "/data/a.asm"});
    require_error(
        [&] {parse_core_config(asm_without_vcf);},
        "ASM without VCF was accepted");

    auto asm_bed_without_vcf = base_arguments();
    asm_bed_without_vcf.insert(
        asm_bed_without_vcf.end(),
        {"--asm-bed", "/data/a.asm.bed"});
    require_error(
        [&] {parse_core_config(asm_bed_without_vcf);},
        "ASM BED without VCF was accepted");

    auto duplicate_methylation_formats = base_arguments();
    duplicate_methylation_formats.insert(
        duplicate_methylation_formats.end(),
        {"--cgmap", "/data/a.cgmap",
         "--bed-methyl", "/data/a.bedmethyl"});
    require_error(
        [&] {parse_core_config(duplicate_methylation_formats);},
        "CGmap and bedMethyl were accepted together");

    auto duplicate_asm_formats = base_arguments();
    replace_value(duplicate_asm_formats, "--mutation-rate", "0");
    duplicate_asm_formats.insert(
        duplicate_asm_formats.end(),
        {"--vcf", "/data/a.vcf",
         "--asm", "/data/a.asm",
         "--asm-bed", "/data/a.asm.bed"});
    require_error(
        [&] {parse_core_config(duplicate_asm_formats);},
        "ASM and ASM BED were accepted together");

    auto rrbs_without_sites = base_arguments();
    replace_value(rrbs_without_sites, "--technology", "RRBS");
    require_error(
        [&] {parse_core_config(rrbs_without_sites);},
        "RRBS without cut sites was accepted");

    auto profile_without_artifact = base_arguments();
    replace_value(profile_without_artifact, "--coverage", "profile");
    require_error(
        [&] {parse_core_config(profile_without_artifact);},
        "profile coverage without an artifact was accepted");

    auto rrbs_profile_without_candidates = base_arguments();
    replace_value(rrbs_profile_without_candidates, "--technology", "RRBS");
    replace_value(rrbs_profile_without_candidates, "--coverage", "profile");
    rrbs_profile_without_candidates.insert(
        rrbs_profile_without_candidates.end(),
        {"--rrbs-cut-site", "C|CGG"});
    require_error(
        [&] {parse_core_config(rrbs_profile_without_candidates);},
        "RRBS profile without a candidate BED was accepted");

    auto rrbs_hash_option = rrbs_profile_without_candidates;
    rrbs_hash_option.insert(
        rrbs_hash_option.end(),
        {"--rrbs-candidate-bed", "/data/candidates.bed",
         "--rrbs-candidate-bed-sha256", std::string(64, '9')});
    require_error(
        [&] {parse_core_config(rrbs_hash_option);},
        "an RRBS candidate BED hash option unexpectedly entered the contract");

    auto target_score_wgbs = base_arguments();
    replace_value(target_score_wgbs, "--coverage", "target-score");
    require_error(
        [&] {parse_core_config(target_score_wgbs);},
        "TBS target-score coverage was accepted for WGBS");

    auto target_score_with_artifact = base_arguments();
    replace_value(target_score_with_artifact, "--technology", "TBS");
    replace_value(target_score_with_artifact, "--coverage", "target-score");
    target_score_with_artifact.insert(
        target_score_with_artifact.end(),
        {"--tbs-bed", "/data/targets.bed",
         "--tbs-center-stddev", "50",
         "--coverage-profile", "/data/coverage.tsv"});
    require_error(
        [&] {parse_core_config(target_score_with_artifact);},
        "TBS target-score coverage accepted a second profile artifact");

    for (const std::string retired_digest_option : {
             "--reference-sha256", "--vcf-sha256", "--cgmap-sha256",
             "--bed-methyl-sha256", "--methdb-sha256", "--asm-sha256",
             "--asm-bed-sha256", "--coverage-profile-sha256",
             "--tbs-bed-sha256"}) {
        auto retired = base_arguments();
        retired.insert(
            retired.end(), {retired_digest_option, std::string(64, '6')});
        require_error(
            [&] {parse_core_config(retired);},
            "retired input digest option was accepted: "
                + retired_digest_option);
    }

    auto profile_without_path = base_arguments();
    replace_value(profile_without_path, "--coverage", "profile");
    require_error(
        [&] {parse_core_config(profile_without_path);},
        "coverage profile without a path was accepted");

    auto tbs_without_center = base_arguments();
    replace_value(tbs_without_center, "--technology", "TBS");
    tbs_without_center.insert(
        tbs_without_center.end(), {"--tbs-bed", "/data/targets.bed"});
    require_error(
        [&] {parse_core_config(tbs_without_center);},
        "TBS BED path without center standard deviation was accepted");

    auto invalid_rrbs_site = base_arguments();
    replace_value(invalid_rrbs_site, "--technology", "RRBS");
    invalid_rrbs_site.insert(
        invalid_rrbs_site.end(), {"--rrbs-cut-site", "ccgg"});
    require_error(
        [&] {parse_core_config(invalid_rrbs_site);},
        "invalid RRBS cut-site grammar was accepted");

    auto duplicate_rrbs_site = base_arguments();
    replace_value(duplicate_rrbs_site, "--technology", "RRBS");
    duplicate_rrbs_site.insert(
        duplicate_rrbs_site.end(),
        {"--rrbs-cut-site", "C|CGG", "--rrbs-cut-site", "C|CGG"});
    require_error(
        [&] {parse_core_config(duplicate_rrbs_site);},
        "duplicate RRBS cut site was accepted");

    auto oversized_rrbs_site = base_arguments();
    replace_value(oversized_rrbs_site, "--technology", "RRBS");
    oversized_rrbs_site.insert(
        oversized_rrbs_site.end(),
        {"--rrbs-cut-site", std::string(1025, 'A') + "|A"});
    require_error(
        [&] {parse_core_config(oversized_rrbs_site);},
        "oversized RRBS cut site was accepted");
}

void test_programmatic_config_uses_the_same_validator()
{
    auto config = parse_core_config(base_arguments());
    validate_core_config(config);

    config.core_workers = 0U;
    require_error(
        [&] {validate_core_config(config);},
        "programmatic configuration bypassed worker validation");

    config = parse_core_config(base_arguments());
    config.mutation_rate = 0.0;
    config.vcf_path = "/data/sample.vcf";
    validate_core_config(config);
    config.vcf_path.reset();
    config.asm_path = "/data/sample.asm";
    require_error(
        [&] {validate_core_config(config);},
        "programmatic configuration bypassed ASM/VCF validation");
}

} // namespace

int main()
{
    try {
        test_valid_wgbs_projection();
        test_output_controls();
        test_valid_optional_input_projection();
        test_cgmap_pool_requires_its_input();
        test_valid_rrbs_tbs_and_profile_projection();
        test_standard_technology_projection();
        test_unknown_duplicate_and_missing_options_fail();
        test_core_worker_projection_and_boundaries();
        test_number_boolean_and_identity_boundaries();
        test_cross_field_boundaries();
        test_programmatic_config_uses_the_same_validator();
    } catch (const std::exception &error) {
        std::cerr << "core_config_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
