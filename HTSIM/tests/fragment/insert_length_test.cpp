#include "fragment.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <future>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using htsim::insert_length::InsertLengthError;
using htsim::insert_length::Parameters;
using htsim::insert_length::sample;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const InsertLengthError &) {
        return;
    }
    throw std::runtime_error(message);
}

void test_frozen_default_vector()
{
    constexpr Parameters parameters{100U, 400U, 1000U, 25.0};
    constexpr std::array<std::uint32_t, 12> expected = {{
        420U, 355U, 402U, 407U, 400U, 414U,
        387U, 400U, 387U, 403U, 440U, 401U,
    }};
    for (std::uint64_t ordinal = 0U; ordinal < expected.size(); ++ordinal) {
        require(
            sample(UINT64_C(0x123456789abcdef0), 0U, ordinal, parameters)
                == expected[static_cast<std::size_t>(ordinal)],
            "clamped-normal insert-length vector changed");
    }
}

void test_fixed_and_zero_dispersion_do_not_change_length()
{
    for (const std::uint64_t ordinal : {
             UINT64_C(0), UINT64_C(1),
             std::numeric_limits<std::uint64_t>::max()}) {
        require(sample(7U, 0U, ordinal, {5U, 5U, 5U, 1000.0}) == 5U,
                "fixed insert length consumed dispersion");
        require(sample(7U, 0U, ordinal, {5U, 8U, 13U, 0.0}) == 8U,
                "zero dispersion changed the mean insert length");
    }
}

void test_extreme_finite_dispersion_clamps_without_integer_overflow()
{
    const Parameters parameters{
        1U,
        UINT32_C(2000000000),
        std::numeric_limits<std::uint32_t>::max(),
        std::numeric_limits<double>::max(),
    };
    bool saw_minimum = false;
    bool saw_maximum = false;
    for (std::uint64_t ordinal = 0U; ordinal < 64U; ++ordinal) {
        const std::uint32_t length = sample(11U, 3U, ordinal, parameters);
        require(length == parameters.minimum || length == parameters.maximum
                    || length == parameters.mean,
                "extreme dispersion escaped its clamp");
        saw_minimum = saw_minimum || length == parameters.minimum;
        saw_maximum = saw_maximum || length == parameters.maximum;
    }
    require(saw_minimum && saw_maximum,
            "extreme-dispersion fixture did not exercise both clamps");
}

void test_addressing_is_order_and_thread_independent()
{
    constexpr Parameters parameters{100U, 400U, 1000U, 25.0};
    std::vector<std::uint32_t> forward(256U);
    for (std::uint64_t ordinal = 0U; ordinal < forward.size(); ++ordinal) {
        forward[static_cast<std::size_t>(ordinal)] =
            sample(19U, 7U, ordinal, parameters);
    }

    std::vector<std::uint32_t> parallel(forward.size());
    std::array<std::future<void>, 4> futures;
    for (std::size_t worker = 0U; worker < futures.size(); ++worker) {
        futures[worker] = std::async(
            std::launch::async,
            [&, worker, parameters] {
                for (std::size_t index = worker;
                     index < parallel.size();
                     index += futures.size()) {
                    parallel[index] = sample(
                        19U,
                        7U,
                        static_cast<std::uint64_t>(index),
                        parameters);
                }
            });
    }
    for (auto &future : futures) {
        future.get();
    }
    require(parallel == forward,
            "thread scheduling changed addressed insert lengths");
    require(sample(19U, 8U, 3U, parameters)
                != sample(19U, 7U, 3U, parameters),
            "contig domain did not affect insert length");
}

void test_invalid_parameters_fail_closed()
{
    for (const Parameters parameters : {
             Parameters{0U, 1U, 1U, 0.0},
             Parameters{5U, 4U, 6U, 1.0},
             Parameters{5U, 7U, 6U, 1.0},
             Parameters{5U, 5U, 6U, -1.0},
             Parameters{5U, 5U, 6U,
                        std::numeric_limits<double>::infinity()},
             Parameters{5U, 5U, 6U,
                        std::numeric_limits<double>::quiet_NaN()},
         }) {
        require_error(
            [&] {sample(0U, 0U, 0U, parameters);},
            "invalid insert-length parameters were accepted");
    }
}

} // namespace

int main()
{
    try {
        test_frozen_default_vector();
        test_fixed_and_zero_dispersion_do_not_change_length();
        test_extreme_finite_dispersion_clamps_without_integer_overflow();
        test_addressing_is_order_and_thread_independent();
        test_invalid_parameters_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "insert_length_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
