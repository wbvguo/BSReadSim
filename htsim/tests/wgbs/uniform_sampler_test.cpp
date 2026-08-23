#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#include "protocol.h"
#include "wgbs.h"

namespace {

using htsim::wgbs::FixedFragmentShape;
using htsim::wgbs::UniformSamplingError;
using htsim::wgbs::ValidStartIndex;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const UniformSamplingError &) {
        return;
    }
    throw std::runtime_error(message);
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

void test_exact_validity_rules()
{
    using htsim::wgbs::count_valid_starts;
    const auto bases = encode("AANAANAAAA");
    require(count_valid_starts(bases, {4, 2, false, 0.0}) == 3,
            "single-end N window eligibility is wrong");
    require(count_valid_starts(bases, {4, 2, true, 0.0}) == 1,
            "paired-end per-mate N eligibility is wrong");
    require(count_valid_starts(bases, {4, 2, true, 0.5}) == 7,
            "ambiguous-fraction equality boundary is wrong");
    require(count_valid_starts(encode("AAAA"), {4, 2, true, 0.0}) == 1,
            "exact-length insert did not produce one start");
    require(count_valid_starts(encode("AAA"), {4, 2, true, 0.0}) == 0,
            "short contig produced a start");
}

void test_sampling_vectors_and_replacement()
{
    using htsim::wgbs::count_valid_starts;
    using htsim::wgbs::sample_valid_starts;
    const auto bases = encode("AACNAACGTAAC");
    const FixedFragmentShape shape{4, 2, true, 0.0};
    const std::uint32_t valid = count_valid_starts(bases, shape);
    require(valid == 5, "sampling fixture valid-start count changed");
    const auto starts = sample_valid_starts(
        bases, 4U, UINT64_C(0xfedcba9876543210), 20, 12,
        shape, valid);
    require(starts == std::vector<std::uint32_t>(
                {5, 4, 7, 8, 8, 4, 5, 7, 6, 7, 7, 5}),
            "uniform valid-rank golden vector changed");
    require(sample_valid_starts(
                bases, 4U, UINT64_C(0xfedcba9876543210), 20, 12,
                shape, valid) == starts,
            "uniform sampling is not deterministic");
}

void test_reusable_rank_index()
{
    const auto bases = encode("AANAANAAAA");
    const FixedFragmentShape shape{4, 2, true, 0.0};
    const ValidStartIndex index(bases, shape);
    require(index.possible_start_count() == 7
                && index.valid_start_count() == 1,
            "rank index counts changed");
    require(!index.is_valid_start(0) && !index.is_valid_start(5)
                && index.is_valid_start(6) && !index.is_valid_start(7),
            "valid-start membership index changed");
    require(index.start_for_rank(0) == 6,
            "rank index selected the wrong valid start");
    require(index.sample(5U, 7, 0, 4)
                == std::vector<std::uint32_t>({6, 6, 6, 6}),
            "rank index did not sample with replacement");
    require_error(
        [&] {(void)index.start_for_rank(1);},
        "out-of-range valid-start rank was accepted");

    const ValidStartIndex empty(encode("AAA"), {4, 2, true, 0.0});
    require(empty.possible_start_count() == 0
                && empty.valid_start_count() == 0,
            "short-contig rank index is not empty");
    require_error(
        [&] {(void)empty.sample(0U, 7, 0, 1);},
        "empty rank index was sampled");
}

void test_order_and_chunk_independence()
{
    using htsim::wgbs::count_valid_starts;
    using htsim::wgbs::sample_valid_starts;
    const auto bases = encode("ACGTACGTACGTACGT");
    const FixedFragmentShape shape{6, 4, true, 0.0};
    const std::uint32_t valid = count_valid_starts(bases, shape);
    const auto all = sample_valid_starts(bases, 0U, 19, 100, 10, shape, valid);
    auto first = sample_valid_starts(bases, 0U, 19, 100, 4, shape, valid);
    const auto second = sample_valid_starts(bases, 0U, 19, 104, 6, shape, valid);
    first.insert(first.end(), second.begin(), second.end());
    require(first == all, "sampling changed across chunk boundaries");

    std::vector<std::uint32_t> concurrent;
    std::thread worker([&] {
        concurrent = sample_valid_starts(
            bases, 0U, 19, 100, 10, shape, valid);
    });
    worker.join();
    require(concurrent == all, "sampling changed across threads");
    require(sample_valid_starts(bases, 1U, 19, 100, 10, shape, valid)
                != all,
            "contig RNG domains were not isolated");
}

void test_invalid_inputs_fail_closed()
{
    using htsim::wgbs::count_valid_starts;
    using htsim::wgbs::sample_valid_starts;
    const auto bases = encode("ACGT");
    for (const FixedFragmentShape shape : {
             FixedFragmentShape{0, 1, false, 0.0},
             {2, 0, false, 0.0},
             {2, 3, false, 0.0},
             {2, 1, false, -0.1},
             {2, 1, false, 1.1},
         }) {
        require_error(
            [&] {(void)count_valid_starts(bases, shape);},
            "invalid fragment shape was accepted");
    }
    require_error(
        [&] {(void)count_valid_starts({0, 5}, {1, 1, false, 0.0});},
        "invalid protocol base was accepted");
    require_error(
        [&] {
            (void)sample_valid_starts(
                bases, 0U, 1, 0, 1, {5, 1, false, 0.0}, 0);
        },
        "sampling from zero valid starts was accepted");
    require_error(
        [&] {
            (void)sample_valid_starts(
                bases, 0U, 1, 0, 1, {2, 1, false, 0.0}, 4);
        },
        "a stale valid-start count was accepted");
}

} // namespace

int main()
{
    try {
        test_exact_validity_rules();
        test_sampling_vectors_and_replacement();
        test_reusable_rank_index();
        test_order_and_chunk_independence();
        test_invalid_inputs_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "uniform_sampler_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
