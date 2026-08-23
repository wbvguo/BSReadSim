#include "utilities.h"

#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {

using htsim::normal_sampler::SamplingError;
using htsim::normal_sampler::standard_normal;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const SamplingError &) {
        return;
    }
    throw std::runtime_error(message);
}

std::uint64_t double_bits(double value) noexcept
{
    std::uint64_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "binary64 size changed");
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

void test_identity_and_golden_vectors()
{
    require(htsim::normal_sampler::algorithm_id == "box-muller-normal-v1",
            "normal algorithm identifier changed");
    struct Golden {
        std::uint64_t key;
        std::uint64_t entity;
        std::uint64_t local;
        std::uint64_t expected_bits;
    };
    const std::array<Golden, 3> cases = {{
        {0, 0, 0, UINT64_C(0xbfd973628c131c06)},
        {UINT64_C(0x0123456789abcdef), 17, 2,
         UINT64_C(0x3fb5d995cc23b49b)},
        {UINT64_MAX, UINT64_MAX, UINT64_MAX,
         UINT64_C(0xbff7a7bdeaf3036d)},
    }};
    for (std::size_t index = 0; index < cases.size(); ++index) {
        const Golden &test = cases[index];
        const std::uint64_t observed = double_bits(
            standard_normal(test.key, test.entity, test.local));
        if (observed != test.expected_bits) {
            std::cerr << "normal_golden[" << index << "]=0x" << std::hex
                      << observed << std::dec << '\n';
            throw std::runtime_error("normal sampler golden vector changed");
        }
    }
}

void test_address_and_thread_independence()
{
    constexpr std::size_t count = 4096;
    constexpr std::size_t workers = 4;
    std::vector<double> baseline(count);
    for (std::size_t index = 0; index < count; ++index) {
        baseline[index] = standard_normal(77, index, 9);
    }
    std::vector<double> parallel(count);
    std::array<std::exception_ptr, workers> failures = {};
    std::array<std::thread, workers> threads;
    for (std::size_t worker = 0; worker < workers; ++worker) {
        threads[worker] = std::thread([&, worker] {
            try {
                for (std::size_t index = worker; index < count; index += workers) {
                    parallel[index] = standard_normal(77, index, 9);
                }
            } catch (...) {
                failures[worker] = std::current_exception();
            }
        });
    }
    for (auto &thread : threads) {thread.join();}
    for (const auto &failure : failures) {
        if (failure) {std::rethrow_exception(failure);}
    }
    for (std::size_t index = 0; index < count; ++index) {
        require(double_bits(parallel[index]) == double_bits(baseline[index]),
                "thread scheduling changed a normal draw");
    }
    require(double_bits(standard_normal(78, 1, 9))
                != double_bits(standard_normal(77, 1, 9)),
            "normal key did not isolate the stream");
    require(double_bits(standard_normal(77, 2, 9))
                != double_bits(standard_normal(77, 1, 9)),
            "normal entity did not isolate the stream");
    require(double_bits(standard_normal(77, 1, 10))
                != double_bits(standard_normal(77, 1, 9)),
            "normal local index did not isolate the stream");
}

void test_distribution_moments()
{
    constexpr std::uint64_t count = UINT64_C(100000);
    double mean = 0.0;
    double squared_deviation_sum = 0.0;
    for (std::uint64_t index = 0; index < count; ++index) {
        const double value = standard_normal(8675309, index, 12);
        const double delta = value - mean;
        mean += delta / static_cast<double>(index + 1U);
        squared_deviation_sum += delta * (value - mean);
    }
    const double variance =
        squared_deviation_sum / static_cast<double>(count - 1U);
    require(std::abs(mean) < 0.01, "normal mean exceeded tolerance");
    require(std::abs(variance - 1.0) < 0.02,
            "normal variance exceeded tolerance");
}

void test_rounding_mode_gate()
{
    const int original = std::fegetround();
    require(original != -1, "could not query floating-point rounding mode");
    require(std::fesetround(FE_DOWNWARD) == 0,
            "could not install non-contract rounding mode");
    try {
        require_error(
            [] {(void)standard_normal(1, 2, 3);},
            "non-nearest rounding mode was accepted");
    } catch (...) {
        (void)std::fesetround(original);
        throw;
    }
    require(std::fesetround(original) == 0,
            "could not restore floating-point rounding mode");
}

} // namespace

int main()
{
    try {
        test_identity_and_golden_vectors();
        test_address_and_thread_independence();
        test_distribution_moments();
        test_rounding_mode_gate();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "normal_sampler_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
