#include "wgbs.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>
#include <zlib.h>

namespace {

using htsim::wgbs::CoverageProfileError;
using htsim::wgbs::WgbsGcProfile;
using htsim::wgbs::WgbsGcSampler;
using htsim::wgbs::VariableWgbsCandidate;
using htsim::wgbs::VariableWgbsGcSampler;
using htsim::wgbs::UnreachableTargetPolicy;
using htsim::wgbs::calibrate_gc_target;
using htsim::wgbs::FixedFragmentShape;

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-coverage-profile-XXXXXX";
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

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const CoverageProfileError &) {
        return;
    }
    throw std::runtime_error(message);
}

std::vector<std::uint8_t> bytes_of(const std::string &text)
{
    return {text.begin(), text.end()};
}

void write_bytes(const std::string &path, const std::vector<std::uint8_t> &bytes)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output.write(
        reinterpret_cast<const char *>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
    if (!output) {throw std::runtime_error("temporary profile write failed");}
}

std::vector<std::uint8_t> gzip_bytes(const std::string &text)
{
    z_stream stream = {};
    require(deflateInit2(
                &stream,
                Z_BEST_SPEED,
                Z_DEFLATED,
                MAX_WBITS + 16,
                8,
                Z_DEFAULT_STRATEGY) == Z_OK,
            "deflateInit2 failed");
    std::vector<std::uint8_t> result(text.size() + 128U);
    stream.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(text.data()));
    stream.avail_in = static_cast<uInt>(text.size());
    stream.next_out = result.data();
    stream.avail_out = static_cast<uInt>(result.size());
    const int status = deflate(&stream, Z_FINISH);
    require(status == Z_STREAM_END, "deflate failed");
    result.resize(result.size() - stream.avail_out);
    require(deflateEnd(&stream) == Z_OK, "deflateEnd failed");
    return result;
}

WgbsGcProfile load_profile(const std::string &text)
{
    TempFile file;
    const auto bytes = bytes_of(text);
    write_bytes(file.path(), bytes);
    return WgbsGcProfile(file.path(), htsim::crypto::sha256(bytes));
}

htsim::model::Bases encode(const std::string &text)
{
    htsim::model::Bases result;
    for (const char base : text) {
        switch (base) {
        case 'A': result.push_back(0); break;
        case 'C': result.push_back(1); break;
        case 'G': result.push_back(2); break;
        case 'T': result.push_back(3); break;
        case 'N': result.push_back(4); break;
        default: throw std::runtime_error("invalid test base");
        }
    }
    return result;
}

void test_profile_parse_and_exact_bin_mapping()
{
    const std::string text = "0\r\n0.25\r\n0.75";
    TempFile plain;
    const auto bytes = bytes_of(text);
    write_bytes(plain.path(), bytes);
    const WgbsGcProfile profile(
        plain.path(), htsim::crypto::sha256(bytes));
    require(profile.bin_count() == 3U,
            "profile did not retain every contiguous bin");
    require(profile.bin_for_counts(0, 4) == 0U
                && profile.bin_for_counts(1, 4) == 1U
                && profile.bin_for_counts(2, 4) == 1U
                && profile.bin_for_counts(3, 4) == 2U
                && profile.bin_for_counts(4, 4) == 2U,
            "exact half-up GC bin mapping changed");
    require(profile.target_probability_for_counts(1, 4) == 0.25
                && profile.target_probability_for_counts(4, 4) == 0.75
                && profile.file_sha256() == htsim::crypto::sha256(bytes),
            "target probability or raw digest changed");
    require_error(
        [&] {(void)profile.bin_for_counts(5, 4);},
        "GC count beyond the fragment was accepted");
    require_error(
        [&] {(void)profile.bin_for_counts(0, 0);},
        "zero fragment length was accepted");

    TempFile compressed_file;
    const auto compressed = gzip_bytes(text);
    write_bytes(compressed_file.path(), compressed);
    const WgbsGcProfile compressed_profile(
        compressed_file.path(), htsim::crypto::sha256(compressed));
    require(compressed_profile.bin_count() == profile.bin_count()
                && compressed_profile.target_probability_for_counts(3, 4)
                    == 0.75
                && compressed_profile.file_sha256()
                    == htsim::crypto::sha256(compressed),
            "gzip profile decoding or raw identity changed");
}

void test_profile_rejections()
{
    const std::vector<std::string> invalid = {
        "",
        "1\n",
        "0\t1\n1\t1\n",
        "0 1\n1\n",
        "1\textra\n1\n",
        "# comment\n1\n1\n",
        "1\n\n1\n",
        "nan\n1\n",
        "-0.1\n1\n",
        "1.1\n1\n",
        "0\n0\n",
        "0.5\n0.6\n",
        "0.5\n0.499999\n",
    };
    for (const std::string &text : invalid) {
        TempFile file;
        const auto bytes = bytes_of(text);
        write_bytes(file.path(), bytes);
        require_error(
            [&] {
                (void)WgbsGcProfile(
                    file.path(), htsim::crypto::sha256(bytes));
            },
            "invalid coverage profile was accepted: " + text);
    }

    TempFile file;
    const auto valid = bytes_of("0.5\n0.5\n");
    write_bytes(file.path(), valid);
    require_error(
        [&] {(void)WgbsGcProfile(file.path(), {});},
        "coverage profile digest mismatch was accepted");
}

void test_profiled_sampling_and_chunk_independence()
{
    const WgbsGcProfile high_gc = load_profile("0\n0\n1\n");
    const FixedFragmentShape shape{4, 4, false, 0.0};
    const WgbsGcSampler sampler(
        encode("AAAACCCCGGGGTTTT"), shape, high_gc);
    require(sampler.valid_start_count() == 13U
                && sampler.bin_opportunity_counts()
                    == std::vector<std::uint32_t>({2U, 4U, 7U}),
            "target sampler changed the GC opportunity domain");
    const auto calibration = calibrate_gc_target(
        high_gc, {sampler.bin_opportunity_counts()});
    require(calibration.acceptance_probabilities
                == std::vector<double>({0.0, 0.0, 1.0})
                && calibration.contig_allocation_weights
                    == std::vector<double>({7.0}),
            "high-GC target calibration changed");

    const auto all = sampler.sample(
        0U, 77, 10, 40, calibration.acceptance_probabilities);
    require(all.skipped_count == 67U
                && all.starts == std::vector<std::uint32_t>({
                    3, 8, 3, 8, 9, 6, 3, 4, 8, 4,
                    5, 6, 8, 3, 3, 5, 4, 7, 6, 9,
                    8, 4, 7, 3, 8, 3, 4, 5, 8, 4,
                    9, 4, 5, 8, 6, 7, 4, 8, 6, 9,
                }),
            "profiled WGBS golden vector changed");
    auto joined = sampler.sample(
        0U, 77, 10, 13, calibration.acceptance_probabilities);
    const auto tail = sampler.sample(
        0U, 77, 23, 27, calibration.acceptance_probabilities);
    joined.starts.insert(joined.starts.end(), tail.starts.begin(), tail.starts.end());
    joined.skipped_count += tail.skipped_count;
    require(all.starts == joined.starts
                && all.skipped_count == joined.skipped_count
                && all.skipped_count > 0U,
            "profiled WGBS sampling changed across chunk boundaries");
    require(std::all_of(
                all.starts.begin(),
                all.starts.end(),
                [](std::uint32_t start) {return start >= 3U && start <= 9U;}),
            "profiled WGBS accepted a zero-probability GC bin");

    const WgbsGcSampler unreachable(
        encode("AAAAAAAAAAAA"), shape, high_gc);
    require(unreachable.valid_start_count() == 9U,
            "unreachable target changed the valid-start domain");
    require_error(
        [&] {
            (void)calibrate_gc_target(
                high_gc, {unreachable.bin_opportunity_counts()});
        },
        "positive target mass without an opportunity was accepted");

    const WgbsGcProfile partly_unreachable = load_profile("0.5\n0\n0.5\n");
    const auto projected = calibrate_gc_target(
        partly_unreachable,
        {unreachable.bin_opportunity_counts()},
        UnreachableTargetPolicy::drop_and_renormalize);
    require(projected.acceptance_probabilities
                    == std::vector<double>({1.0, 0.0, 0.0})
                && projected.contig_allocation_weights
                    == std::vector<double>({9.0})
                && projected.dropped_target_probability == 0.5,
            "approximate calibration did not project unreachable target mass");
    require_error(
        [&] {
            (void)calibrate_gc_target(
                high_gc,
                {unreachable.bin_opportunity_counts()},
                UnreachableTargetPolicy::drop_and_renormalize);
        },
        "projection accepted a target with no reachable positive mass");

    const WgbsGcProfile certain = load_profile("0.5\n0.5\n");
    const WgbsGcSampler certain_sampler(
        encode("ACGT"), {1U, 1U, false, 0.0}, certain);
    const auto certain_calibration = calibrate_gc_target(
        certain, {certain_sampler.bin_opportunity_counts()});
    const auto certain_batch = certain_sampler.sample(
        1U,
        77,
        0,
        20,
        certain_calibration.acceptance_probabilities);
    require(certain_batch.starts.size() == 20U
                && certain_batch.skipped_count == 0U,
            "opportunity-matched target unexpectedly rejected candidates");
    require(certain_sampler.sample(
                1U,
                77,
                0,
                0,
                certain_calibration.acceptance_probabilities).starts.empty(),
            "zero-output request returned candidates");
    require_error(
        [&] {
            (void)certain_sampler.sample(
                1U,
                77,
                std::numeric_limits<std::uint64_t>::max(),
                2,
                certain_calibration.acceptance_probabilities);
        },
        "overflowing coverage candidate ordinals were accepted");

    const WgbsGcProfile mixed = load_profile("0.2\n0.3\n0.5\n");
    const auto mixed_calibration = calibrate_gc_target(
        mixed, {{2U, 4U, 0U}, {0U, 0U, 7U}});
    require(std::abs(mixed_calibration.acceptance_probabilities[0] - 1.0)
                    < 1e-15
                && std::abs(
                    mixed_calibration.acceptance_probabilities[1] - 0.75)
                    < 1e-15
                && std::abs(
                    mixed_calibration.acceptance_probabilities[2]
                        - 5.0 / 7.0) < 1e-15
                && std::abs(
                    mixed_calibration.contig_allocation_weights[0] - 5.0)
                    < 1e-15
                && std::abs(
                    mixed_calibration.contig_allocation_weights[1] - 5.0)
                    < 1e-15,
            "global target calibration or contig mixing changed");
}

void test_variable_insert_profile_rejection_and_chunk_independence()
{
    const auto bases = encode("AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT");
    const WgbsGcProfile high_gc = load_profile("0\n0\n1\n");
    const VariableWgbsGcSampler sampler(
        bases,
        0U,
        991U,
        {4U, 7U, 10U, 2.0},
        4U,
        false,
        0.0,
        high_gc);
    require(sampler.allocation_weight() > 0U,
            "variable GC sampler lost its proposal domain");

    const std::vector<double> high_gc_only{0.0, 0.0, 1.0};
    const auto whole = sampler.sample(0U, 80U, high_gc_only);
    const auto first = sampler.sample(0U, 23U, high_gc_only);
    const auto second = sampler.sample(
        first.next_candidate_ordinal, 57U, high_gc_only);
    std::vector<VariableWgbsCandidate> joined = first.candidates;
    joined.insert(
        joined.end(), second.candidates.begin(), second.candidates.end());
    require(whole.candidates.size() == 80U
                && whole.candidates.size() == joined.size()
                && whole.skipped_count > 0U
                && whole.skipped_count
                    == first.skipped_count + second.skipped_count
                && whole.next_candidate_ordinal
                    == second.next_candidate_ordinal,
            "profiled variable-insert accounting changed across chunks");

    bool saw_different_insert_lengths = false;
    const std::uint32_t first_insert = whole.candidates.front().insert_length;
    for (std::size_t index = 0U; index < whole.candidates.size(); ++index) {
        const auto &candidate = whole.candidates[index];
        require(candidate.reference_start == joined[index].reference_start
                    && candidate.insert_length == joined[index].insert_length,
                "profiled variable-insert candidate changed across chunks");
        saw_different_insert_lengths = saw_different_insert_lengths
            || candidate.insert_length != first_insert;
        const auto begin = bases.begin() + candidate.reference_start;
        const auto end = begin + candidate.insert_length;
        const std::uint32_t gc_count = static_cast<std::uint32_t>(std::count_if(
            begin,
            end,
            [](std::uint8_t base) {return base == 1U || base == 2U;}));
        require(high_gc.bin_for_counts(gc_count, candidate.insert_length) == 2U,
                "profiled variable-insert sampler accepted a zero-mass bin");
    }
    require(saw_different_insert_lengths,
            "GC rejection accidentally fixed the insert length");

    require(sampler.sample(0U, 0U, high_gc_only).candidates.empty(),
            "zero-output variable GC request returned candidates");
    require_error(
        [&] {(void)sampler.sample(0U, 1U, {1.0});},
        "variable GC sampler accepted the wrong bin count");
    require_error(
        [&] {(void)sampler.sample(0U, 1U, {0.0, -0.1, 1.1});},
        "variable GC sampler accepted invalid probabilities");
    require_error(
        [&] {(void)sampler.sample(0U, 1U, {0.0, 0.0, 0.0});},
        "variable GC sampler accepted an all-zero calibration");
}

} // namespace

int main()
{
    try {
        test_profile_parse_and_exact_bin_mapping();
        test_profile_rejections();
        test_profiled_sampling_and_chunk_independence();
        test_variable_insert_profile_rejection_and_chunk_independence();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "coverage_profile_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
