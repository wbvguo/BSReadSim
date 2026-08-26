#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include "core.h"
#include "protocol.h"

namespace {

using htsim::core::CoreConfig;
using htsim::core::CoreGeneratorError;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-core-generator-XXXXXX";
        const int descriptor = mkstemp(pattern);
        if (descriptor < 0) {throw std::runtime_error("mkstemp failed");}
        path_ = pattern;
        if (close(descriptor) != 0) {
            throw std::runtime_error("temporary descriptor close failed");
        }
    }

    ~TempFile() {if (!path_.empty()) {(void)unlink(path_.c_str());}}
    const std::string &path() const noexcept {return path_;}

private:
    std::string path_;
};

std::vector<std::uint8_t> bytes_of(const std::string &text)
{
    return {text.begin(), text.end()};
}

void write_text(const std::string &path, const std::string &text)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output.write(text.data(), static_cast<std::streamsize>(text.size()));
    if (!output) {throw std::runtime_error("temporary FASTA write failed");}
}

CoreConfig baseline_config(const TempFile &reference)
{
    CoreConfig config;
    config.emit_details = true;
    config.run_id = "00000000-0000-0000-0000-000000000001";
    config.normalized_config_sha256 = htsim::crypto::sha256(bytes_of("config"));
    config.master_seed = UINT64_C(0x123456789abcdef0);
    config.reference_path = reference.path();
    config.technology = htsim::core::Technology::wgbs;
    config.paired_end = true;
    config.read_length_1 = 3;
    config.read_length_2 = 3;
    config.insert_min = 5;
    config.insert_mean = 5;
    config.insert_max = 5;
    config.insert_sd = 0.0;
    config.fragment_count = 7;
    config.max_ambiguous_fraction = 0.0;
    config.chunk_size = 2;
    config.mutation_rate = 0.0;
    config.collect_non_cpg = true;
    config.cgmap_pool = false;
    config.beta_cg = {2.0, 5.0};
    config.beta_chg = {3.0, 4.0};
    config.beta_chh = {5.0, 2.0};
    config.coverage = htsim::core::CoverageMode::uniform;
    return config;
}

std::string vcf_header()
{
    return
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
}

template <typename Operation>
void require_empty_failure(Operation operation, const std::string &message)
{
    std::ostringstream sink(std::ios::binary);
    try {
        operation(sink);
    } catch (const CoreGeneratorError &) {
        require(sink.str().empty(), message + ": failure wrote output bytes");
        return;
    }
    throw std::runtime_error(message + ": failure was not reported");
}

void test_valid_generation_and_chunk_independence()
{
    const std::string fasta = ">chr1\nACGTCGTAACGT\n>chr2\nAACGTACGTTAA\n";
    TempFile reference;
    write_text(reference.path(), fasta);
    CoreConfig config = baseline_config(reference);

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(first_trailer.fragment_count == 7
                && first_trailer.mate_count == 14
                && first_trailer.skipped_fragment_count == 0,
            "valid generation returned wrong trailer counts");
    require(first_trailer.per_contig_fragment_counts.size() == 2
                && first_trailer.per_contig_fragment_counts[0]
                    + first_trailer.per_contig_fragment_counts[1] == 7,
            "per-contig allocation did not preserve the total");

    std::ostringstream repeated(std::ios::binary);
    const auto repeated_trailer =
        htsim::core::generate_core_stream(config, repeated);
    require(repeated.str() == first.str()
                && repeated_trailer.stream_sha256 == first_trailer.stream_sha256,
            "same configuration did not produce identical protocol bytes");

    CoreConfig bounded_fixed = config;
    bounded_fixed.insert_min = 3U;
    bounded_fixed.insert_max = 8U;
    std::ostringstream bounded(std::ios::binary);
    const auto bounded_trailer =
        htsim::core::generate_core_stream(bounded_fixed, bounded);
    require(bounded.str() == first.str()
                && bounded_trailer.stream_sha256 == first_trailer.stream_sha256,
            "fixed WGBS output depended on insert bounds instead of insert_mean");

    config.chunk_size = 5;
    std::ostringstream rechunked(std::ios::binary);
    const auto rechunked_trailer =
        htsim::core::generate_core_stream(config, rechunked);
    require(rechunked.str() == first.str()
                && rechunked_trailer.stream_sha256 == first_trailer.stream_sha256,
            "chunk size changed protocol output");
}

void test_standard_technologies_emit_no_methylation_sites()
{
    TempFile reference;
    write_text(reference.path(), ">chr1\nACGTCGTAACGT\n");

    CoreConfig wgs = baseline_config(reference);
    wgs.technology = htsim::core::Technology::wgs;
    std::ostringstream wgs_stream(std::ios::binary);
    const auto wgs_trailer =
        htsim::core::generate_core_stream(wgs, wgs_stream);
    require(wgs_trailer.fragment_count == 7U
                && wgs_trailer.methylation_site_count == 0U,
            "WGS generated methylation rows");

    TempFile targets;
    write_text(targets.path(), "chr1\t4\t5\ttarget\t1\t+\n");
    for (const auto technology : {
             htsim::core::Technology::wes,
             htsim::core::Technology::ts}) {
        CoreConfig targeted = baseline_config(reference);
        targeted.technology = technology;
        targeted.tbs_bed_path = targets.path();
        targeted.tbs_center_stddev = 0.0;
        std::ostringstream stream(std::ios::binary);
        const auto trailer =
            htsim::core::generate_core_stream(targeted, stream);
        require(trailer.fragment_count == 7U
                    && trailer.methylation_site_count == 0U,
                "standard targeted sequencing generated methylation rows");
    }
}

void test_variable_wgbs_generation_and_capability_gate()
{
    std::string sequence;
    sequence.reserve(80U);
    constexpr char bases[] = "ACGT";
    for (std::size_t index = 0U; index < 80U; ++index) {
        sequence.push_back(index % 7U == 0U ? 'N' : bases[index % 4U]);
    }
    const std::string fasta = ">chrVariable\n" + sequence + "\n";
    TempFile reference;
    write_text(reference.path(), fasta);
    CoreConfig config = baseline_config(reference);
    config.master_seed = 91U;
    config.read_length_1 = 4U;
    config.read_length_2 = 4U;
    config.insert_min = 5U;
    config.insert_mean = 9U;
    config.insert_max = 16U;
    config.insert_sd = 4.0;
    config.fragment_count = 40U;
    config.chunk_size = 13U;

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(
        first_trailer.fragment_count == 40U
            && first_trailer.mate_count == 80U
            && first_trailer.skipped_fragment_count > 0U
            && first_trailer.per_contig_fragment_counts
                == std::vector<std::uint64_t>({40U}),
        "variable WGBS generation returned wrong trailer counts");

    config.chunk_size = 27U;
    std::ostringstream rechunked(std::ios::binary);
    const auto rechunked_trailer =
        htsim::core::generate_core_stream(config, rechunked);
    require(
        rechunked.str() == first.str()
            && rechunked_trailer.skipped_fragment_count
                == first_trailer.skipped_fragment_count,
        "chunk size changed variable WGBS protocol output");

    const std::string profile =
        "0\n0\n1\n0\n0\n";
    TempFile profile_file;
    write_text(profile_file.path(), profile);
    CoreConfig profiled = config;
    profiled.coverage = htsim::core::CoverageMode::profile;
    profiled.coverage_profile_path = profile_file.path();
    std::ostringstream profiled_stream(std::ios::binary);
    const auto profiled_trailer =
        htsim::core::generate_core_stream(profiled, profiled_stream);
    require(profiled_trailer.fragment_count == 40U
                && profiled_trailer.mate_count == 80U
                && profiled_trailer.skipped_fragment_count > 0U
                && profiled_stream.str() != first.str(),
            "variable-insert target GC generation lost counts or rejection");
    profiled.chunk_size = 27U;
    std::ostringstream profiled_rechunked(std::ios::binary);
    const auto profiled_rechunked_trailer =
        htsim::core::generate_core_stream(profiled, profiled_rechunked);
    require(profiled_rechunked.str() == profiled_stream.str()
                && profiled_rechunked_trailer.skipped_fragment_count
                    == profiled_trailer.skipped_fragment_count,
            "chunk size changed variable-insert target GC output");

    const std::string vcf = vcf_header()
        + "chrVariable\t2\t.\tC\tT\t.\tPASS\t.\tGT\t1|0\n"
          "chrVariable\t5\t.\tA\tACC\t.\tPASS\t.\tGT\t1|1\n"
          "chrVariable\t10\t.\tCG\tC\t.\tPASS\t.\tGT\t0|1\n"
          "chrVariable\t42\t.\tC\tT\t.\tPASS\t.\tGT\t1|1\n";
    TempFile variants;
    write_text(variants.path(), vcf);
    CoreConfig with_vcf = config;
    with_vcf.chunk_size = 13U;
    with_vcf.vcf_path = variants.path();
    std::ostringstream vcf_stream(std::ios::binary);
    const auto vcf_trailer =
        htsim::core::generate_core_stream(with_vcf, vcf_stream);
    require(
        vcf_trailer.fragment_count == 40U
            && vcf_trailer.mate_count == 80U
            && vcf_stream.str() != first.str(),
        "variable WGBS VCF generation lost counts or typed events");
    with_vcf.chunk_size = 27U;
    std::ostringstream vcf_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(with_vcf, vcf_rechunked);
    require(vcf_rechunked.str() == vcf_stream.str(),
            "chunk size changed variable WGBS VCF output");

    const std::string variable_asm =
        "chrVariable\tC\t18\tCG\tCG\t0.5\t2\tC\tT\t0.2\t0.8\t4\t0.01\tfixture\n"
        "chrVariable\tC\t26\tCG\tCG\t0.5\t2\tC\tT\t0.2\t0.8\t4\t0.01\tfixture\n"
        "chrVariable\tC\t34\tCG\tCG\t0.5\t2\tC\tT\t0.2\t0.8\t4\t0.01\tfixture\n"
        "chrVariable\tC\t58\tCG\tCG\t0.5\t2\tC\tT\t0.2\t0.8\t4\t0.01\tfixture\n"
        "chrVariable\tC\t66\tCG\tCG\t0.5\t2\tC\tT\t0.2\t0.8\t4\t0.01\tfixture\n"
        "chrVariable\tC\t74\tCG\tCG\t0.5\t2\tC\tT\t0.2\t0.8\t4\t0.01\tfixture\n";
    TempFile asm_file;
    write_text(asm_file.path(), variable_asm);
    CoreConfig with_asm = with_vcf;
    with_asm.fragment_count = 128U;
    with_asm.chunk_size = 11U;
    with_asm.asm_path = asm_file.path();
    std::ostringstream asm_stream(std::ios::binary);
    const auto asm_trailer =
        htsim::core::generate_core_stream(with_asm, asm_stream);
    require(
        asm_trailer.fragment_count == 128U
            && asm_trailer.mate_count == 256U
            && asm_stream.str() != vcf_stream.str(),
        "variable WGBS ASM generation lost counts or its overlay");
    with_asm.chunk_size = 29U;
    std::ostringstream asm_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(with_asm, asm_rechunked);
    require(asm_rechunked.str() == asm_stream.str(),
            "chunk size changed variable WGBS ASM output");

    CoreConfig with_mutations = config;
    with_mutations.chunk_size = 13U;
    with_mutations.mutation_rate = 0.8;
    with_mutations.indel_fraction = 0.65;
    with_mutations.indel_extension_probability = 0.6;
    std::ostringstream mutation_stream(std::ios::binary);
    const auto mutation_trailer =
        htsim::core::generate_core_stream(
            with_mutations, mutation_stream);
    require(
        mutation_trailer.fragment_count == 40U
            && mutation_trailer.mate_count == 80U
            && mutation_stream.str() != first.str(),
        "variable WGBS de novo generation lost counts or typed events");
    with_mutations.chunk_size = 27U;
    std::ostringstream mutation_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(
        with_mutations, mutation_rechunked);
    require(mutation_rechunked.str() == mutation_stream.str(),
            "chunk size changed variable WGBS de novo output");

    const std::string short_fasta = ">chrShort\nACGTACGTACGT\n";
    TempFile short_reference;
    write_text(short_reference.path(), short_fasta);
    CoreConfig unsupported = baseline_config(short_reference);
    unsupported.read_length_1 = 4U;
    unsupported.read_length_2 = 4U;
    unsupported.insert_min = 5U;
    unsupported.insert_mean = 9U;
    unsupported.insert_max = 16U;
    unsupported.insert_sd = 4.0;
    unsupported.fragment_count = 1U;
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(unsupported, sink);
        },
        "variable WGBS without a maximum-span allocation domain");
}

void test_wgbs_depth_conversion_and_preflight_rejections()
{
    const std::string fasta =
        ">chr1\nACGTCGTAACGT\n>ineligible\nNNNNNNNNNNNN\n";
    TempFile reference;
    write_text(reference.path(), fasta);
    CoreConfig config = baseline_config(reference);
    config.fragment_count.reset();
    config.depth = 2.0;

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(
        first_trailer.fragment_count == 4U
            && first_trailer.mate_count == 8U
            && first_trailer.per_contig_fragment_counts
                == std::vector<std::uint64_t>({4U, 0U}),
        "WGBS depth conversion used the wrong effective reference length");

    config.chunk_size = 5U;
    std::ostringstream rechunked(std::ios::binary);
    const auto rechunked_trailer =
        htsim::core::generate_core_stream(config, rechunked);
    require(
        rechunked.str() == first.str()
            && rechunked_trailer.stream_sha256
                == first_trailer.stream_sha256,
        "chunk size changed depth-derived WGBS protocol output");

    config.depth = 0.01;
    std::ostringstream minimum_depth(std::ios::binary);
    const auto minimum_depth_trailer =
        htsim::core::generate_core_stream(config, minimum_depth);
    require(minimum_depth_trailer.fragment_count == 1U,
            "positive WGBS depth did not ceil to one fragment");

    config.depth = 2.0;
    config.technology = htsim::core::Technology::rrbs;
    config.rrbs_cut_sites = {"|C"};
    std::ostringstream rrbs_depth(std::ios::binary);
    const auto rrbs_depth_trailer =
        htsim::core::generate_core_stream(config, rrbs_depth);
    require(rrbs_depth_trailer.fragment_count == 2U
                && rrbs_depth_trailer.mate_count == 4U,
            "RRBS depth did not use the eligible restriction-region denominator");

    const std::string target_text = "chr1\t1\t7\tprobe\t1\t+\n";
    TempFile targets;
    write_text(targets.path(), target_text);
    config.technology = htsim::core::Technology::tbs;
    config.rrbs_cut_sites.clear();
    config.tbs_bed_path = targets.path();
    config.tbs_center_stddev = 0.0;
    std::ostringstream tbs_depth(std::ios::binary);
    const auto tbs_depth_trailer =
        htsim::core::generate_core_stream(config, tbs_depth);
    require(tbs_depth_trailer.fragment_count == 2U
                && tbs_depth_trailer.mate_count == 4U,
            "TBS depth did not use the target-region union denominator");
}

void test_valid_wgbs_profile_generation_and_chunk_independence()
{
    const std::string fasta =
        ">chr1\nACGTCGTAACGT\n>chr2\nAACGTACGTTAA\n";
    const std::string profile_text = "0\n0\n1\n";
    TempFile reference;
    TempFile profile;
    write_text(reference.path(), fasta);
    write_text(profile.path(), profile_text);
    CoreConfig config = baseline_config(reference);
    config.coverage = htsim::core::CoverageMode::profile;
    config.coverage_profile_path = profile.path();

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(first_trailer.fragment_count == 7
                && first_trailer.mate_count == 14
                && first_trailer.skipped_fragment_count > 0
                && first_trailer.per_contig_fragment_counts
                    == std::vector<std::uint64_t>({7, 0}),
            "WGBS GC profile generation returned wrong trailer counts");

    config.chunk_size = 5;
    std::ostringstream rechunked(std::ios::binary);
    const auto rechunked_trailer =
        htsim::core::generate_core_stream(config, rechunked);
    require(rechunked.str() == first.str()
                && rechunked_trailer.skipped_fragment_count
                    == first_trailer.skipped_fragment_count,
            "chunk size changed profiled WGBS protocol output");

    config.technology = htsim::core::Technology::rrbs;
    config.rrbs_cut_sites = {"|C"};
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "RRBS with a WGBS GC profile");
}

void test_valid_wgbs_vcf_generation_and_chunk_independence()
{
    const std::string fasta = ">chr1\nAACGTAACGTT\n";
    const std::string vcf = vcf_header()
        + "chr1\t4\t.\tG\tA\t.\tPASS\t.\tGT\t1|0\n"
          "chr1\t6\t.\tA\tACG\t.\tPASS\t.\tGT\t1|1\n";
    TempFile reference;
    TempFile variants;
    write_text(reference.path(), fasta);
    write_text(variants.path(), vcf);
    CoreConfig config = baseline_config(reference);
    config.vcf_path = variants.path();

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(first_trailer.fragment_count == 7U
                && first_trailer.mate_count == 14U
                && first_trailer.skipped_fragment_count == 0U,
            "WGBS VCF generation returned wrong trailer counts");

    config.chunk_size = 5U;
    std::ostringstream rechunked(std::ios::binary);
    const auto rechunked_trailer =
        htsim::core::generate_core_stream(config, rechunked);
    require(rechunked.str() == first.str()
                && rechunked_trailer.stream_sha256
                    == first_trailer.stream_sha256,
            "chunk size changed the WGBS VCF protocol stream");

    CoreConfig baseline = baseline_config(reference);
    std::ostringstream reference_only(std::ios::binary);
    (void)htsim::core::generate_core_stream(baseline, reference_only);
    require(first.str() != reference_only.str(),
            "VCF events did not affect the emitted protocol stream");

    const std::string empty_vcf = vcf_header();
    write_text(variants.path(), empty_vcf);
    baseline.vcf_path = variants.path();
    std::ostringstream no_events(std::ios::binary);
    (void)htsim::core::generate_core_stream(baseline, no_events);
    require(no_events.str() == reference_only.str(),
            "an event-free VCF changed reference-only protocol bytes");
}

void test_wgbs_vcf_deletions_and_preflight_rejections()
{
    const std::string fasta = ">chr1\nAACGTAACGTT\n";
    const std::string deletion_vcf = vcf_header()
        + "chr1\t7\t.\tACG\tA\t.\tPASS\t.\tGT\t1|0\n";
    TempFile reference;
    TempFile variants;
    write_text(reference.path(), fasta);
    write_text(variants.path(), deletion_vcf);
    CoreConfig config = baseline_config(reference);
    config.vcf_path = variants.path();
    std::ostringstream deletion_stream(std::ios::binary);
    const auto deletion_trailer =
        htsim::core::generate_core_stream(config, deletion_stream);
    require(deletion_trailer.fragment_count == 7U
                && deletion_trailer.mate_count == 14U
                && deletion_trailer.skipped_fragment_count == 0U,
            "WGBS VCF deletion generation returned wrong counts");
    config.chunk_size = 5U;
    std::ostringstream deletion_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(
        config, deletion_rechunked);
    require(deletion_rechunked.str() == deletion_stream.str(),
            "chunk size changed the WGBS deletion protocol stream");

    const std::string all_one_profile = "0.5\n0.5\n";
    TempFile coverage_profile;
    write_text(coverage_profile.path(), all_one_profile);
    config.coverage = htsim::core::CoverageMode::profile;
    config.coverage_profile_path = coverage_profile.path();
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "target GC profile with a VCF");

    config = baseline_config(reference);
    config.vcf_path = variants.path();
    config.update_variant_boundaries = false;
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "VCF with disabled variant-boundary updates");
}

void test_de_novo_mutation_generation_and_preflight_rejections()
{
    const std::string fasta =
        ">chr1\nACGTACGTACGTACGTACGTACGTACGTACGT\n";
    TempFile reference;
    write_text(reference.path(), fasta);
    CoreConfig config = baseline_config(reference);
    config.mutation_rate = 0.8;
    config.indel_fraction = 0.65;
    config.indel_extension_probability = 0.6;

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(
        first_trailer.fragment_count == 7U
            && first_trailer.mate_count == 14U
            && first_trailer.skipped_fragment_count == 0U,
        "de novo mutation generation returned wrong counts");
    config.chunk_size = 5U;
    std::ostringstream rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(config, rechunked);
    require(rechunked.str() == first.str(),
            "chunk size changed the de novo mutation protocol stream");

    CoreConfig reference_only = config;
    reference_only.mutation_rate = 0.0;
    std::ostringstream reference_stream(std::ios::binary);
    (void)htsim::core::generate_core_stream(
        reference_only, reference_stream);
    require(reference_stream.str() != first.str(),
            "de novo mutations did not affect the protocol stream");

    const std::string all_one_profile = "0.5\n0.5\n";
    TempFile coverage_profile;
    write_text(coverage_profile.path(), all_one_profile);
    CoreConfig profiled = config;
    profiled.coverage = htsim::core::CoverageMode::profile;
    profiled.coverage_profile_path = coverage_profile.path();
    std::ostringstream profiled_stream(std::ios::binary);
    const auto profiled_trailer =
        htsim::core::generate_core_stream(profiled, profiled_stream);
    require(profiled_trailer.fragment_count == 7U
                && profiled_trailer.mate_count == 14U,
            "target GC profile lost de novo mutation fragment counts");
    profiled.chunk_size = 2U;
    std::ostringstream profiled_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(
        profiled, profiled_rechunked);
    require(profiled_rechunked.str() == profiled_stream.str(),
            "chunk size changed target GC de novo mutation output");

    const std::string empty_vcf = vcf_header();
    TempFile variants;
    write_text(variants.path(), empty_vcf);
    config.vcf_path = variants.path();
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "VCF with de novo mutations");

    config.vcf_path.reset();
    config.update_variant_boundaries = false;
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "de novo mutations with disabled variant-boundary updates");

    config.update_variant_boundaries = true;
    config.technology = htsim::core::Technology::rrbs;
    config.rrbs_cut_sites = {"|C"};
    config.insert_min = 4U;
    config.insert_mean = 5U;
    config.insert_max = 8U;
    config.insert_sd = 0.0;
    std::ostringstream rrbs_mutations(std::ios::binary);
    const auto rrbs_mutation_trailer =
        htsim::core::generate_core_stream(config, rrbs_mutations);
    require(rrbs_mutation_trailer.fragment_count == 7U
                && rrbs_mutation_trailer.mate_count == 14U,
            "RRBS de novo haplotypes did not reach fragmentation");
    config.chunk_size = 2U;
    std::ostringstream rrbs_mutations_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(
        config, rrbs_mutations_rechunked);
    require(rrbs_mutations_rechunked.str() == rrbs_mutations.str(),
            "chunk size changed RRBS de novo haplotype fragmentation");

    const std::string tbs_bed =
        "chr1\t8\t9\ta\t1\t+\n"
        "chr1\t16\t17\tb\t1\t-\n"
        "chr1\t24\t25\tc\t1\t.\n";
    TempFile tbs_targets;
    write_text(tbs_targets.path(), tbs_bed);
    CoreConfig tbs_mutations = baseline_config(reference);
    tbs_mutations.technology = htsim::core::Technology::tbs;
    tbs_mutations.mutation_rate = 0.2;
    tbs_mutations.indel_fraction = 0.65;
    tbs_mutations.indel_extension_probability = 0.6;
    tbs_mutations.tbs_bed_path = tbs_targets.path();
    tbs_mutations.tbs_center_stddev = 0.0;
    std::ostringstream tbs_mutation_stream(std::ios::binary);
    const auto tbs_mutation_trailer =
        htsim::core::generate_core_stream(
            tbs_mutations, tbs_mutation_stream);
    require(tbs_mutation_trailer.fragment_count == 7U
                && tbs_mutation_trailer.mate_count == 14U,
            "TBS de novo haplotypes did not reach fragmentation");
    tbs_mutations.chunk_size = 5U;
    std::ostringstream tbs_mutation_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(
        tbs_mutations, tbs_mutation_rechunked);
    require(tbs_mutation_rechunked.str() == tbs_mutation_stream.str(),
            "chunk size changed TBS de novo haplotype fragmentation");

    config = baseline_config(reference);
    config.mutation_rate = std::numeric_limits<double>::quiet_NaN();
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "non-finite de novo mutation probability");
}

void test_valid_cgmap_generation_and_preflight_rejections()
{
    const std::string fasta = ">chr1\nAACGTAACGTT\n";
    const std::string cgmap =
        "chr1\tC\t3\tCG\tCG\t1\t8\t8\n"
        "chr1\tG\t4\tCG\tCG\tna\t0\t0\n"
        "chr1\tC\t8\tCG\tCG\t0\t0\t8\n"
        "chr1\tG\t9\tCG\tCG\t0.25\t2\t8\n";
    TempFile reference;
    TempFile profile;
    write_text(reference.path(), fasta);
    write_text(profile.path(), cgmap);
    CoreConfig config = baseline_config(reference);
    config.cgmap_path = profile.path();

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(
        first_trailer.fragment_count == 7U
            && first_trailer.methylation_site_count > 0U,
        "CGmap generation returned wrong counts");
    config.chunk_size = 5U;
    std::ostringstream rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(config, rechunked);
    require(
        rechunked.str() == first.str(),
        "chunk size changed the CGmap protocol stream");

    CoreConfig beta = baseline_config(reference);
    std::ostringstream beta_stream(std::ios::binary);
    (void)htsim::core::generate_core_stream(beta, beta_stream);
    require(
        first.str() != beta_stream.str(),
        "CGmap values did not affect the protocol stream");

    const std::string bed_methyl =
        "chr1\t2\t3\tm\t8\t+\t2\t3\t255,0,0\t8\t100\n"
        "chr1\t7\t8\tm\t8\t+\t7\t8\t255,0,0\t8\t0\n"
        "chr1\t8\t9\tm\t8\t-\t8\t9\t255,0,0\t8\t25\n";
    TempFile bed_profile;
    write_text(bed_profile.path(), bed_methyl);
    CoreConfig bed_config = config;
    bed_config.chunk_size = 5U;
    bed_config.cgmap_path.reset();
    bed_config.bed_methyl_path = bed_profile.path();
    std::ostringstream bed_stream(std::ios::binary);
    (void)htsim::core::generate_core_stream(bed_config, bed_stream);
    require(
        bed_stream.str() == first.str(),
        "bedMethyl did not preserve the equivalent CGmap MethDB overlay");

    config.asm_path = profile.path();
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "ASM input without VCF");

    config = baseline_config(reference);
    config.cgmap_path = profile.path();
    const std::string mismatch =
        "chr1\tC\t3\tCHH\tCA\t1\t1\t1\n";
    write_text(profile.path(), mismatch);
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "CGmap context/reference mismatch");
}

void test_valid_asm_generation_and_preflight_rejections()
{
    const std::string fasta = ">chr1\nACGTA\n";
    const std::string vcf = vcf_header()
        + "chr1\t5\t.\tA\tC\t.\tPASS\t.\tGT\t1|0\n";
    const std::string cgmap =
        "chr1\tC\t2\tCG\tCG\t0.5\t4\t8\n";
    const std::string asm_profile =
        "chr1\tC\t2\tCG\tCG\t0.5\t5\tA\tC\t0.2\t0.8\t4\t0.01\tfixture\n";
    TempFile reference;
    TempFile variants;
    TempFile cgmap_file;
    TempFile asm_file;
    write_text(reference.path(), fasta);
    write_text(variants.path(), vcf);
    write_text(cgmap_file.path(), cgmap);
    write_text(asm_file.path(), asm_profile);

    CoreConfig config = baseline_config(reference);
    config.vcf_path = variants.path();
    config.cgmap_path = cgmap_file.path();
    config.asm_path = asm_file.path();

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(
        first_trailer.fragment_count == 7U
            && first_trailer.methylation_site_count >= 7U,
        "ASM generation returned wrong counts");
    config.chunk_size = 5U;
    std::ostringstream rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(config, rechunked);
    require(
        rechunked.str() == first.str(),
        "chunk size changed the ASM protocol stream");

    CoreConfig cgmap_only = config;
    cgmap_only.asm_path.reset();
    std::ostringstream cgmap_stream(std::ios::binary);
    (void)htsim::core::generate_core_stream(cgmap_only, cgmap_stream);
    require(
        cgmap_stream.str() != first.str(),
        "ASM did not override the CGmap probability");

    const std::string asm_bed =
        "chr1\t1\t2\tfixture\t0\t+\t4\t5\tA\tC\t0.2\t0.8\n";
    TempFile asm_bed_file;
    write_text(asm_bed_file.path(), asm_bed);
    CoreConfig bed_config = config;
    bed_config.asm_path.reset();
    bed_config.asm_bed_path = asm_bed_file.path();
    std::ostringstream bed_stream(std::ios::binary);
    (void)htsim::core::generate_core_stream(bed_config, bed_stream);
    require(
        bed_stream.str() == first.str(),
        "ASM BED did not preserve the equivalent typed ASM overlay");

    const std::string unresolved_asm =
        "chr1\tC\t2\tCG\tCG\t0.5\t4\tT\tC\t0.2\t0.8\t4\t0.01\tfixture\n";
    write_text(asm_file.path(), unresolved_asm);
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "ASM row without its exact linked VCF SNV");

    config = baseline_config(reference);
    config.asm_path = asm_file.path();
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "ASM without VCF");
}

void test_valid_rrbs_generation_and_chunk_independence()
{
    const std::string fasta = ">chrR\nCAACAACAAC\n";
    TempFile reference;
    write_text(reference.path(), fasta);
    CoreConfig config = baseline_config(reference);
    config.technology = htsim::core::Technology::rrbs;
    config.read_length_1 = 2;
    config.read_length_2 = 2;
    config.insert_min = 3;
    config.insert_mean = 5;
    config.insert_max = 9;
    config.insert_sd = 25.0;
    config.fragment_count = 5;
    config.rrbs_cut_sites = {"|C"};

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(first_trailer.fragment_count == 5
                && first_trailer.mate_count == 10
                && first_trailer.per_contig_fragment_counts
                    == std::vector<std::uint64_t>({5}),
            "RRBS generation returned wrong trailer counts");

    config.chunk_size = 4;
    std::ostringstream rechunked(std::ios::binary);
    const auto rechunked_trailer =
        htsim::core::generate_core_stream(config, rechunked);
    require(rechunked.str() == first.str()
                && rechunked_trailer.stream_sha256
                    == first_trailer.stream_sha256,
            "chunk size changed RRBS protocol output");

    std::ostringstream candidate_bed;
    htsim::core::generate_rrbs_candidate_bed(config, candidate_bed);
    require(candidate_bed.str().find(
                "#chrom\tstart\tend\tcandidate_id\tscore") == 0U
                && candidate_bed.str().find("chrR:") != std::string::npos,
            "RRBS catalog-only generation did not emit candidate BED rows");

    CoreConfig recursive_export = config;
    recursive_export.rrbs_candidate_bed_path = "/data/candidates.bed";
    require_empty_failure(
        [&](std::ostringstream &sink) {
            htsim::core::generate_rrbs_candidate_bed(
                recursive_export, sink);
        },
        "RRBS catalog export with a candidate BED input");
}

void test_rrbs_and_tbs_vcf_haplotype_fragmentation()
{
    const std::string fasta = ">chr1\nAACCGGAAACTGGAA\n";
    TempFile reference;
    write_text(reference.path(), fasta);

    const std::string rrbs_vcf = vcf_header()
        + "chr1\t11\t.\tT\tC\t.\tPASS\t.\tGT\t1|0\n";
    TempFile rrbs_variants;
    write_text(rrbs_variants.path(), rrbs_vcf);
    CoreConfig rrbs = baseline_config(reference);
    rrbs.technology = htsim::core::Technology::rrbs;
    rrbs.read_length_1 = 2U;
    rrbs.read_length_2 = 2U;
    rrbs.insert_min = 7U;
    rrbs.insert_mean = 7U;
    rrbs.insert_max = 7U;
    rrbs.insert_sd = 0.0;
    rrbs.rrbs_cut_sites = {"CCG|G"};
    rrbs.vcf_path = rrbs_variants.path();
    std::ostringstream rrbs_stream(std::ios::binary);
    const auto rrbs_trailer =
        htsim::core::generate_core_stream(rrbs, rrbs_stream);
    require(rrbs_trailer.fragment_count == 7U
                && rrbs_trailer.mate_count == 14U,
            "RRBS VCF haplotype motifs did not generate fragments");
    rrbs.chunk_size = 5U;
    std::ostringstream rrbs_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(rrbs, rrbs_rechunked);
    require(rrbs_rechunked.str() == rrbs_stream.str(),
            "chunk size changed RRBS VCF haplotype fragmentation");

    const std::string tbs_vcf = vcf_header()
        + "chr1\t7\t.\tA\tATT\t.\tPASS\t.\tGT\t1|0\n";
    const std::string bed = "chr1\t7\t8\tprobe\t1\t+\n";
    TempFile tbs_variants;
    TempFile targets;
    write_text(tbs_variants.path(), tbs_vcf);
    write_text(targets.path(), bed);
    CoreConfig tbs = baseline_config(reference);
    tbs.technology = htsim::core::Technology::tbs;
    tbs.tbs_bed_path = targets.path();
    tbs.tbs_center_stddev = 0.0;
    tbs.vcf_path = tbs_variants.path();
    std::ostringstream tbs_stream(std::ios::binary);
    const auto tbs_trailer =
        htsim::core::generate_core_stream(tbs, tbs_stream);
    require(tbs_trailer.fragment_count == 7U
                && tbs_trailer.mate_count == 14U,
            "TBS VCF haplotype centers did not generate fragments");
    tbs.chunk_size = 5U;
    std::ostringstream tbs_rechunked(std::ios::binary);
    (void)htsim::core::generate_core_stream(tbs, tbs_rechunked);
    require(tbs_rechunked.str() == tbs_stream.str(),
            "chunk size changed TBS VCF haplotype fragmentation");
}

void test_valid_tbs_generation_and_chunk_independence()
{
    const std::string fasta =
        ">chr1\nACGTCGTAACGT\n>chr2\nAACGTACGTTAA\n";
    const std::string bed =
        "chr1\t4\t5\tforward\t1\t+\n"
        "chr1\t8\t9\treverse\t2\t-\n"
        "chr2\t5\t6\tunknown\t0\t.\n";
    TempFile reference;
    TempFile targets;
    write_text(reference.path(), fasta);
    write_text(targets.path(), bed);
    CoreConfig config = baseline_config(reference);
    config.technology = htsim::core::Technology::tbs;
    config.tbs_bed_path = targets.path();
    config.tbs_center_stddev = 0.0;
    config.insert_min = 3U;
    config.insert_mean = 5U;
    config.insert_max = 8U;
    config.insert_sd = 0.0;

    std::ostringstream first(std::ios::binary);
    const auto first_trailer =
        htsim::core::generate_core_stream(config, first);
    require(first_trailer.fragment_count == 7
                && first_trailer.mate_count == 14
                && first_trailer.per_contig_fragment_counts.size() == 2,
            "TBS generation returned wrong trailer counts");

    config.chunk_size = 5;
    std::ostringstream rechunked(std::ios::binary);
    const auto rechunked_trailer =
        htsim::core::generate_core_stream(config, rechunked);
    require(rechunked.str() == first.str()
                && rechunked_trailer.stream_sha256
                    == first_trailer.stream_sha256,
            "chunk size changed TBS protocol output");

    config.coverage = htsim::core::CoverageMode::target_score;
    config.chunk_size = 2;
    std::ostringstream weighted;
    const auto weighted_trailer =
        htsim::core::generate_core_stream(config, weighted);
    require(weighted_trailer.per_contig_fragment_counts
                == std::vector<std::uint64_t>({7U, 0U})
                && weighted.str() != first.str(),
            "TBS target scores did not control exact contig allocation");
    config.chunk_size = 5;
    std::ostringstream weighted_rechunked;
    const auto weighted_rechunked_trailer =
        htsim::core::generate_core_stream(config, weighted_rechunked);
    require(weighted_rechunked.str() == weighted.str()
                && weighted_rechunked_trailer.stream_sha256
                    == weighted_trailer.stream_sha256,
            "chunk size changed target-score TBS protocol output");

    config.tbs_center_stddev = 4.0;
    config.chunk_size = 2;
    std::ostringstream displaced(std::ios::binary);
    const auto displaced_trailer =
        htsim::core::generate_core_stream(config, displaced);
    require(displaced_trailer.fragment_count == 7
                && displaced_trailer.skipped_fragment_count > 0,
            "TBS normal center sampling did not report accepted/skipped counts");
    config.chunk_size = 5;
    std::ostringstream displaced_rechunked(std::ios::binary);
    const auto displaced_rechunked_trailer =
        htsim::core::generate_core_stream(config, displaced_rechunked);
    require(displaced_rechunked.str() == displaced.str()
                && displaced_rechunked_trailer.skipped_fragment_count
                    == displaced_trailer.skipped_fragment_count,
            "chunk size changed displaced TBS protocol output");

    const std::string all_zero =
        "chr1\t4\t5\tforward\t0\t+\n"
        "chr2\t5\t6\tunknown\t0\t.\n";
    write_text(targets.path(), all_zero);
    config.tbs_center_stddev = 0.0;
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "all-zero TBS target-score profile");

    const std::string fractional =
        "chr1\t4\t5\tforward\t0.5\t+\n";
    write_text(targets.path(), fractional);
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "fractional TBS target-score weight");
}

void test_preflight_failures_write_nothing()
{
    const std::string fasta = ">chr1\nACGTACGT\n";
    TempFile reference;
    write_text(reference.path(), fasta);
    CoreConfig config = baseline_config(reference);

    config.coverage = htsim::core::CoverageMode::target_score;
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "TBS target-score coverage on WGBS");

    config = baseline_config(reference);

    config.technology = htsim::core::Technology::tbs;
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "TBS without its required input projection");

    config = baseline_config(reference);
    config.technology = htsim::core::Technology::rrbs;
    config.rrbs_cut_sites = {"|A"};
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "RRBS reference with fewer than two cut positions");

    config = baseline_config(reference);
    config.max_ambiguous_fraction = 0.0;
    write_text(reference.path(), ">chr1\nNNNNNNNN\n");
    require_empty_failure(
        [&](std::ostringstream &sink) {
            (void)htsim::core::generate_core_stream(config, sink);
        },
        "no eligible fragment starts");
}

} // namespace

int main()
{
    try {
        test_valid_generation_and_chunk_independence();
        test_standard_technologies_emit_no_methylation_sites();
        test_variable_wgbs_generation_and_capability_gate();
        test_wgbs_depth_conversion_and_preflight_rejections();
        test_valid_wgbs_profile_generation_and_chunk_independence();
        test_valid_wgbs_vcf_generation_and_chunk_independence();
        test_wgbs_vcf_deletions_and_preflight_rejections();
        test_de_novo_mutation_generation_and_preflight_rejections();
        test_valid_cgmap_generation_and_preflight_rejections();
        test_valid_asm_generation_and_preflight_rejections();
        test_valid_rrbs_generation_and_chunk_independence();
        test_rrbs_and_tbs_vcf_haplotype_fragmentation();
        test_valid_tbs_generation_and_chunk_independence();
        test_preflight_failures_write_nothing();
    } catch (const std::exception &error) {
        std::cerr << "core_generator_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
