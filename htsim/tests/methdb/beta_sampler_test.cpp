#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "methdb.h"
#include "utilities.h"

namespace {

using htsim::beta_sampler::Options;
using htsim::beta_sampler::SamplingError;
using htsim::beta_sampler::sample_beta;
using htsim::beta_sampler::sample_beta_for_site;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_sampling_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const SamplingError &) {
        return;
    }
    throw std::runtime_error(message);
}

std::uint32_t float_bits(float value) noexcept
{
    std::uint32_t bits = 0;
    static_assert(sizeof(bits) == sizeof(value), "binary32 size changed");
    std::memcpy(&bits, &value, sizeof(bits));
    return bits;
}

template <typename Alpha, typename Beta, typename = void>
struct can_sample_beta : std::false_type {};

template <typename Alpha, typename Beta>
struct can_sample_beta<
    Alpha,
    Beta,
    std::void_t<decltype(sample_beta(
        UINT64_C(1),
        UINT32_C(0),
        UINT64_C(2),
        std::declval<Alpha>(),
        std::declval<Beta>()))>> : std::true_type {};

static_assert(can_sample_beta<double, double>::value,
              "numeric Beta shapes must be accepted");
static_assert(can_sample_beta<int, int>::value,
              "integer-valued numeric Beta shapes must convert to double");
static_assert(!can_sample_beta<bool, double>::value,
              "boolean alpha must be rejected at compile time");
static_assert(!can_sample_beta<double, bool>::value,
              "boolean beta must be rejected at compile time");

template <typename Alpha, typename Beta, typename = void>
struct can_sample_beta_for_site : std::false_type {};

template <typename Alpha, typename Beta>
struct can_sample_beta_for_site<
    Alpha,
    Beta,
    std::void_t<decltype(sample_beta_for_site(
        UINT64_C(1),
        UINT32_C(0),
        std::declval<htsim::methdb::SiteEntity>(),
        std::declval<Alpha>(),
        std::declval<Beta>()))>> : std::true_type {};

static_assert(can_sample_beta_for_site<double, double>::value,
              "numeric site-entity Beta shapes must be accepted");
static_assert(!can_sample_beta_for_site<bool, double>::value,
              "boolean site-entity alpha must be rejected at compile time");
static_assert(!can_sample_beta_for_site<double, bool>::value,
              "boolean site-entity beta must be rejected at compile time");

void test_algorithm_identity_and_branches()
{
    using htsim::beta_sampler::algorithm_id;
    using htsim::beta_sampler::site_entity_algorithm_id;
    require(algorithm_id == "marsaglia-tsang-box-muller-beta",
            "Beta sampler algorithm identifier changed");
    require(site_entity_algorithm_id
                == "marsaglia-tsang-box-muller-beta-site-entity",
            "site-entity Beta sampler algorithm identifier changed");

    const float below_one = sample_beta(11, 0U, 3, 0.25, 0.75);
    const float equal_one = sample_beta(11, 0U, 4, 1.0, 1.0);
    const float above_one = sample_beta(11, 0U, 5, 2.5, 7.25);
    for (const float value : {below_one, equal_one, above_one}) {
        require(std::isfinite(value) && value >= 0.0F && value <= 1.0F,
                "a Gamma shape branch returned a value outside [0, 1]");
    }
}

void test_site_entity_vectors_and_reference_identity()
{
    using htsim::model::HaplotypeMask;
    using namespace htsim::methdb;
    const SiteEntity unchanged = reference_site_entity(700U);
    require(float_bits(sample_beta_for_site(
                99, 1U, unchanged, 0.5, 3.0))
                == float_bits(sample_beta(99, 1U, 700U, 0.5, 3.0)),
            "tag-0 site entity changed a reference-only Beta draw");

    struct SiteGolden {
        SiteEntity entity;
        std::uint32_t expected_bits;
    };
    const std::array<SiteGolden, 4> cases = {{
        {variant_reference_site_entity(700U, HaplotypeMask::both, 0U),
         UINT32_C(0x3d6862eb)},
        {variant_reference_site_entity(
             700U, HaplotypeMask::haplotype_1, 0U),
         UINT32_C(0x3dc2ba9d)},
        {variant_reference_site_entity(
             700U, HaplotypeMask::haplotype_2, 1U),
         UINT32_C(0x3ca9258f)},
        {insertion_site_entity(
             7U, 2U, HaplotypeMask::haplotype_2, 1U),
         UINT32_C(0x3e568c7d)},
    }};
    bool all_matched = true;
    for (std::size_t index = 0; index < cases.size(); ++index) {
        const std::uint32_t observed = float_bits(sample_beta_for_site(
            99, 1U, cases[index].entity, 0.5, 3.0));
        if (observed != cases[index].expected_bits) {
            std::cerr << "site-entity-golden[" << index << "] observed=0x"
                      << std::hex << observed << std::dec << '\n';
            all_matched = false;
        }
    }
    require(all_matched,
            "exact site-entity Beta sampler golden vector changed");
}

struct GoldenCase {
    std::uint64_t seed;
    std::uint32_t contig_index;
    std::uint64_t reference_position;
    double alpha;
    double beta;
    std::uint32_t expected_bits;
};

void test_exact_golden_vectors()
{
    const std::array<GoldenCase, 3> cases = {{
        {0, 0U, 0, 0.25, 0.75, UINT32_C(0x3f1ac2f8)},
        {42, 3U, UINT64_C(123456789), 1.0, 1.0,
         UINT32_C(0x3f72308e)},
        {UINT64_C(12345678901234567890), 20U,
         std::numeric_limits<std::uint64_t>::max(), 2.5, 7.25,
         UINT32_C(0x3e8e2d92)},
    }};

    bool all_matched = true;
    for (std::size_t index = 0; index < cases.size(); ++index) {
        const GoldenCase &test = cases[index];
        const std::uint32_t observed = float_bits(sample_beta(
            test.seed,
            test.contig_index,
            test.reference_position,
            test.alpha,
            test.beta));
        if (observed != test.expected_bits) {
            std::cerr << "golden[" << index << "] observed=0x" << std::hex
                      << observed << std::dec << '\n';
            all_matched = false;
        }
    }
    require(all_matched, "exact Beta sampler golden vector changed");
}

void test_address_determinism_and_domain_separation()
{
    const float baseline = sample_beta(99, 1U, 700, 0.5, 3.0);
    require(float_bits(sample_beta(99, 1U, 700, 0.5, 3.0))
                == float_bits(baseline),
            "identical Beta addresses were not deterministic");
    require(float_bits(sample_beta(100, 1U, 700, 0.5, 3.0))
                != float_bits(baseline),
            "master seed did not isolate the Beta stream");
    require(float_bits(sample_beta(99, 2U, 700, 0.5, 3.0))
                != float_bits(baseline),
            "contig did not isolate the Beta stream");
    require(float_bits(sample_beta(99, 1U, 701, 0.5, 3.0))
                != float_bits(baseline),
            "reference position did not isolate the Beta entity");
    require(float_bits(sample_beta(99, 1U, 700, 1.5, 3.0))
                != float_bits(baseline),
            "shape selection did not affect the Beta result");

    for (std::uint64_t position = 0; position < UINT64_C(10000); ++position) {
        const float value = sample_beta(1234, 4U, position, 0.05, 0.20);
        require(std::isfinite(value) && value >= 0.0F && value <= 1.0F,
                "Beta sample left its finite output domain");
    }
}

void test_order_and_thread_independence()
{
    constexpr std::size_t sample_count = 4096;
    constexpr std::size_t thread_count = 4;
    std::vector<float> baseline(sample_count);
    for (std::size_t index = 0; index < sample_count; ++index) {
        baseline[index] = sample_beta(2026, 5U, index, 0.7, 4.2);
    }

    std::vector<float> reverse(sample_count);
    for (std::size_t offset = sample_count; offset != 0; --offset) {
        const std::size_t index = offset - 1;
        reverse[index] = sample_beta(2026, 5U, index, 0.7, 4.2);
    }
    for (std::size_t index = 0; index < sample_count; ++index) {
        require(float_bits(reverse[index]) == float_bits(baseline[index]),
                "reverse call order changed an addressed Beta draw");
    }

    std::vector<float> parallel(sample_count);
    std::array<std::exception_ptr, thread_count> failures = {};
    std::array<std::thread, thread_count> threads;
    for (std::size_t worker = 0; worker < thread_count; ++worker) {
        threads[worker] = std::thread([&, worker] {
            try {
                for (std::size_t index = worker;
                     index < sample_count;
                     index += thread_count) {
                    parallel[index] =
                        sample_beta(2026, 5U, index, 0.7, 4.2);
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
    for (std::size_t index = 0; index < sample_count; ++index) {
        require(float_bits(parallel[index]) == float_bits(baseline[index]),
                "thread scheduling changed an addressed Beta draw");
    }
}

void test_distribution_moments()
{
    struct ShapeCase {
        double alpha;
        double beta;
        std::uint32_t contig_index;
    };
    const std::array<ShapeCase, 3> shapes = {{
        {0.25, 0.75, 6U},
        {1.0, 1.0, 7U},
        {2.0, 5.0, 8U},
    }};
    constexpr std::uint64_t sample_count = UINT64_C(100000);
    for (const ShapeCase &shape : shapes) {
        double mean = 0.0;
        double squared_deviation_sum = 0.0;
        for (std::uint64_t index = 0; index < sample_count; ++index) {
            const double value = sample_beta(
                8675309, shape.contig_index, index, shape.alpha, shape.beta);
            const double delta = value - mean;
            mean += delta / static_cast<double>(index + 1);
            squared_deviation_sum += delta * (value - mean);
        }
        const double variance =
            squared_deviation_sum / static_cast<double>(sample_count - 1);
        const double expected_mean = shape.alpha / (shape.alpha + shape.beta);
        const double expected_variance =
            shape.alpha * shape.beta
            / ((shape.alpha + shape.beta) * (shape.alpha + shape.beta)
               * (shape.alpha + shape.beta + 1.0));
        require(std::abs(mean - expected_mean) < 0.004,
                "100k-sample Beta mean exceeded its statistical tolerance");
        require(std::abs(variance - expected_variance) < 0.002,
                "100k-sample Beta variance exceeded its statistical tolerance");
    }
}

void test_invalid_and_extreme_inputs()
{
    const double infinity = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const double invalid : {0.0, -1.0, infinity, -infinity, nan}) {
        require_sampling_error(
            [&] {sample_beta(1, 0U, 2, invalid, 1.0);},
            "invalid alpha shape was accepted");
        require_sampling_error(
            [&] {sample_beta(1, 0U, 2, 1.0, invalid);},
            "invalid beta shape was accepted");
    }

    require_sampling_error(
        [] {sample_beta(1, 0U, 2, 1.0, 1.0, Options{0});},
        "zero Gamma attempt cap was accepted");
    require_sampling_error(
        [] {
            sample_beta(
                1,
                0U,
                2,
                1.0,
                1.0,
                Options{htsim::beta_sampler::default_max_gamma_attempts + 1});
        },
        "an attempt cap above the algorithm hard maximum was accepted");

    const int original_rounding = std::fegetround();
    require(original_rounding != -1, "could not query floating-point rounding mode");
    require(std::fesetround(FE_DOWNWARD) == 0,
            "could not install a non-contract rounding mode");
    try {
        require_sampling_error(
            [] {sample_beta(1, 0U, 2, 1.0, 1.0);},
            "non-nearest floating-point rounding mode was accepted");
    } catch (...) {
        (void)std::fesetround(original_rounding);
        throw;
    }
    require(std::fesetround(original_rounding) == 0,
            "could not restore floating-point rounding mode");

    const double maximum = std::numeric_limits<double>::max();
    const float large = sample_beta(7, 9U, 8, maximum, maximum);
    require(std::isfinite(large) && large >= 0.0F && large <= 1.0F,
            "maximum finite shapes did not produce a finite Beta value");

    const double minimum = std::numeric_limits<double>::denorm_min();
    require_sampling_error(
        [&] {sample_beta(7, 9U, 8, minimum, minimum);},
        "indeterminate double-underflow Beta ratio did not fail closed");
}

void test_attempt_cap_exhaustion()
{
    bool found_exhaustion = false;
    for (std::uint64_t position = 0; position < UINT64_C(10000); ++position) {
        try {
            (void)sample_beta(314159, 10U, position, 1.0, 1.0, Options{1});
        } catch (const SamplingError &) {
            found_exhaustion = true;
            const float recovered =
                sample_beta(314159, 10U, position, 1.0, 1.0);
            require(std::isfinite(recovered),
                    "default hard cap did not recover from a first rejection");
            break;
        }
    }
    require(found_exhaustion,
            "could not exercise deterministic Gamma attempt-cap exhaustion");
}

} // namespace

int main()
{
    try {
        test_algorithm_identity_and_branches();
        test_site_entity_vectors_and_reference_identity();
        test_exact_golden_vectors();
        test_address_determinism_and_domain_separation();
        test_order_and_thread_independence();
        test_distribution_moments();
        test_invalid_and_extreme_inputs();
        test_attempt_cap_exhaustion();
    } catch (const std::exception &error) {
        std::cerr << "beta_sampler_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
