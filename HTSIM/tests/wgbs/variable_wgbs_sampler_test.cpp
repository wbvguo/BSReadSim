#include "wgbs.h"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using htsim::model::Bases;
using htsim::wgbs::VariableWgbsCandidate;
using htsim::wgbs::VariableWgbsOptions;
using htsim::wgbs::VariableWgbsSampler;
using htsim::wgbs::VariableWgbsSamplingError;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const VariableWgbsSamplingError &) {
        return;
    }
    throw std::runtime_error(message);
}

Bases repeating_bases(std::size_t length)
{
    Bases bases;
    bases.reserve(length);
    for (std::size_t index = 0U; index < length; ++index) {
        bases.push_back(static_cast<std::uint8_t>(index % 4U));
    }
    return bases;
}

void test_frozen_candidate_vector()
{
    const Bases bases = repeating_bases(60U);
    const VariableWgbsSampler sampler(
        bases,
        0U,
        UINT64_C(0x123456789abcdef0),
        {5U, 8U, 12U, 2.0},
        3U,
        true,
        0.0);
    require(sampler.allocation_weight() == 49U,
            "maximum-insert allocation weight changed");
    const auto first_candidate = sampler.candidate_at(0U);
    require(first_candidate
                && first_candidate->reference_start == 45U
                && first_candidate->insert_length == 9U,
            "single addressed variable proposal changed");
    const auto batch = sampler.sample(0U, 12U);
    const std::vector<VariableWgbsCandidate> expected = {
        {45U, 9U}, {27U, 5U}, {52U, 8U}, {20U, 8U},
        {16U, 8U}, {48U, 9U}, {16U, 7U}, {36U, 8U},
        {8U, 7U},  {50U, 8U}, {40U, 11U}, {23U, 8U},
    };
    require(batch.candidates.size() == expected.size(),
            "variable WGBS candidate count changed");
    for (std::size_t index = 0U; index < expected.size(); ++index) {
        require(
            batch.candidates[index].reference_start
                    == expected[index].reference_start
                && batch.candidates[index].insert_length
                    == expected[index].insert_length,
            "variable WGBS candidate vector changed");
    }
    require(batch.skipped_count == 0U
                && batch.next_candidate_ordinal == 12U,
            "candidate accounting changed without rejection");
}

void test_rejection_is_chunk_independent()
{
    Bases bases = repeating_bases(80U);
    for (std::size_t index = 0U; index < bases.size(); index += 7U) {
        bases[index] = 4U;
    }
    const VariableWgbsSampler sampler(
        bases, 0U, 91U, {5U, 9U, 16U, 4.0}, 4U, true, 0.0);
    const auto whole = sampler.sample(0U, 40U);
    const auto first = sampler.sample(0U, 13U);
    const auto second = sampler.sample(first.next_candidate_ordinal, 27U);
    std::vector<VariableWgbsCandidate> combined = first.candidates;
    combined.insert(
        combined.end(), second.candidates.begin(), second.candidates.end());
    require(whole.skipped_count > 0U,
            "N fixture did not exercise variable-insert rejection");
    require(whole.skipped_count
                    == first.skipped_count + second.skipped_count
                && whole.next_candidate_ordinal
                    == second.next_candidate_ordinal
                && whole.candidates.size() == combined.size(),
            "chunking changed variable-insert accounting");
    for (std::size_t index = 0U; index < combined.size(); ++index) {
        require(
            whole.candidates[index].reference_start
                    == combined[index].reference_start
                && whole.candidates[index].insert_length
                    == combined[index].insert_length,
            "chunking changed a variable-insert candidate");
    }
}

void test_single_end_eligibility_uses_only_the_emitted_mate()
{
    Bases bases(12U, 4U);
    bases[0] = 0U;
    bases[1] = 0U;
    bases[2] = 0U;
    const VariableWgbsSampler single_end(
        bases, 0U, 17U, {5U, 5U, 10U, 0.0}, 3U, false, 0.0);
    require(single_end.allocation_weight() == 1U,
            "single-end maximum-span eligibility changed");
    const auto batch = single_end.sample(0U, 3U);
    require(batch.candidates.size() == 3U,
            "single-end variable sampler returned the wrong count");
    for (const VariableWgbsCandidate &candidate : batch.candidates) {
        require(candidate.reference_start == 0U
                    && candidate.insert_length == 5U,
                "single-end variable sampler inspected an unsequenced mate");
    }

    const VariableWgbsSampler paired_end(
        bases, 0U, 17U, {5U, 5U, 10U, 0.0}, 3U, true, 0.0);
    require(paired_end.allocation_weight() == 0U,
            "paired-end eligibility ignored the reverse mate");
}

void test_attempt_cap_and_invalid_domains_fail_closed()
{
    Bases cap_bases(20U, 4U);
    for (const std::size_t index : {0U, 1U, 2U, 7U, 8U, 9U}) {
        cap_bases[index] = 0U;
    }
    const VariableWgbsSampler impossible(
        cap_bases,
        0U,
        3U,
        {10U, 10U, 10U, 0.0},
        3U,
        true,
        0.0,
        VariableWgbsOptions{1U});
    require(impossible.allocation_weight() > 0U,
            "attempt-cap fixture lost its maximum-length domain");
    require_error(
        [&] {(void)impossible.sample(0U, 1U);},
        "impossible insert distribution exhausted no attempt cap");

    const Bases all_n(20U, 4U);
    const VariableWgbsSampler no_starts(
        all_n, 0U, 0U, {5U, 5U, 5U, 0.0}, 3U, true, 0.0);
    require(no_starts.allocation_weight() == 0U,
            "all-N contig gained an allocation weight");
    require_error(
        [&] {(void)no_starts.sample(0U, 1U);},
        "zero-weight contig was sampled");

    require_error(
        [&] {
            (void)VariableWgbsSampler(
                cap_bases,
                0U,
                0U,
                {3U, 3U, 3U, 0.0},
                4U,
                true,
                0.0);
        },
        "read longer than minimum insert was accepted");
    require_error(
        [&] {
            (void)VariableWgbsSampler(
                cap_bases,
                0U,
                0U,
                {3U, 3U, 3U, 0.0},
                2U,
                false,
                0.0,
                VariableWgbsOptions{0U});
        },
        "zero attempt cap was accepted");
}

} // namespace

int main()
{
    try {
        test_frozen_candidate_vector();
        test_rejection_is_chunk_independent();
        test_single_end_eligibility_uses_only_the_emitted_mate();
        test_attempt_cap_and_invalid_domains_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "variable_wgbs_sampler_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
