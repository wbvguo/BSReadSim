#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

#include "fragment.h"

namespace {

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const htsim::allocation::AllocationError &) {
        return;
    }
    throw std::runtime_error(message);
}

void test_exact_proportions_and_remainders()
{
    using htsim::allocation::largest_remainder;
    require(largest_remainder({1, 2, 3}, 10) ==
                std::vector<std::uint32_t>({2, 3, 5}),
            "largest remainders were allocated incorrectly");
    require(largest_remainder({1, 1, 1}, 2) ==
                std::vector<std::uint32_t>({1, 1, 0}),
            "equal remainders did not use FASTA order");
    require(largest_remainder({0, 5, 5}, 1) ==
                std::vector<std::uint32_t>({0, 1, 0}),
            "zero-weight contig received a fragment or tie order changed");
}

void test_u32_products_and_u64_total_weight()
{
    using htsim::allocation::largest_remainder;
    constexpr std::uint32_t requested =
        std::numeric_limits<std::uint32_t>::max();
    const auto single = largest_remainder(
        {0, std::numeric_limits<std::uint32_t>::max()}, requested);
    require(single == std::vector<std::uint32_t>({0, requested}),
            "maximum u32 exact allocation overflowed");

    const auto equal_max = largest_remainder(
        {std::numeric_limits<std::uint32_t>::max(),
         std::numeric_limits<std::uint32_t>::max()},
        requested);
    require(equal_max == std::vector<std::uint32_t>(
                             {UINT32_C(2147483648), UINT32_C(2147483647)}),
            "u64 total weight or FASTA-order tie handling changed");

    const auto skewed = largest_remainder(
        {std::numeric_limits<std::uint32_t>::max(), 1}, requested);
    require(std::accumulate(skewed.begin(), skewed.end(), UINT64_C(0))
                == requested,
            "u32 allocation did not preserve its total");
    require(skewed == std::vector<std::uint32_t>({requested - 1, 1}),
            "maximum u32 product selected the wrong remainder seat");
}

void test_real_weights()
{
    using htsim::allocation::largest_remainder_real;
    require(largest_remainder_real({0.1, 0.2, 0.7}, 10U)
                == std::vector<std::uint32_t>({1U, 2U, 7U}),
            "real allocation changed exact quotas");
    require(largest_remainder_real({1.0, 1.0, 1.0}, 2U)
                == std::vector<std::uint32_t>({1U, 1U, 0U}),
            "real allocation tie order changed");
    require(largest_remainder_real({0.0, 0.25, 0.75}, 5U)
                == std::vector<std::uint32_t>({0U, 1U, 4U}),
            "real allocation used a zero-weight contig");
    require_error(
        [&] {largest_remainder_real({0.0, 0.0}, 1U);},
        "zero real weights were accepted");
    require_error(
        [&] {largest_remainder_real({1.0, -1.0}, 1U);},
        "negative real weight was accepted");
}

void test_invalid_inputs_fail_closed()
{
    using htsim::allocation::largest_remainder;
    require_error(
        [&] {largest_remainder({}, 1);},
        "empty weight vector was accepted");
    require_error(
        [&] {largest_remainder({1}, 0);},
        "zero requested total was accepted");
    require_error(
        [&] {largest_remainder({0, 0}, 1);},
        "all-zero weights were accepted");
}

} // namespace

int main()
{
    try {
        test_exact_proportions_and_remainders();
        test_u32_products_and_u64_total_weight();
        test_real_weights();
        test_invalid_inputs_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "allocation_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
