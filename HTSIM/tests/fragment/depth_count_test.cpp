#include "fragment.h"

#include <cfenv>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

using htsim::depth_count::DepthCountError;
using htsim::depth_count::read_pairs;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const DepthCountError &) {
        return;
    }
    throw std::runtime_error(message);
}

void test_human_scale_and_pairing()
{
    require(
        read_pairs(20.0, UINT64_C(3100000000), 100U, true)
            == UINT32_C(310000000),
        "human-scale paired-end depth conversion changed");
    require(
        read_pairs(20.0, UINT64_C(3100000000), 100U, false)
            == UINT32_C(620000000),
        "human-scale single-end depth conversion changed");
    require(read_pairs(0.5, 1001U, 100U, false) == 5U,
            "fractional depth did not use frozen floor semantics");
    require(
        read_pairs(
            2.0,
            std::numeric_limits<std::uint32_t>::max(),
            2U,
            false)
            == std::numeric_limits<std::uint32_t>::max(),
        "maximum uint32 result changed");
}

void test_invalid_and_overflow_inputs()
{
    require_error(
        [] {read_pairs(0.0, 100U, 10U, false);},
        "zero depth was accepted");
    require_error(
        [] {
            read_pairs(
                std::numeric_limits<double>::quiet_NaN(),
                100U,
                10U,
                false);
        },
        "NaN depth was accepted");
    require_error(
        [] {read_pairs(1.0, 0U, 10U, false);},
        "zero effective reference length was accepted");
    require_error(
        [] {read_pairs(1.0, 100U, 0U, false);},
        "zero read length was accepted");
    require_error(
        [] {read_pairs(0.001, 10U, 100U, true);},
        "depth-derived zero count was accepted");
    require_error(
        [] {
            read_pairs(
                1.0,
                std::numeric_limits<std::uint64_t>::max(),
                1U,
                false);
        },
        "count above uint32 was accepted");
}

void test_rounding_mode_is_part_of_the_contract()
{
    const int original = std::fegetround();
    require(original != -1, "cannot query floating-point rounding mode");
    require(std::fesetround(FE_DOWNWARD) == 0,
            "cannot set downward rounding mode");
    require_error(
        [] {read_pairs(1.0, 100U, 10U, false);},
        "non-nearest rounding mode was accepted");
    require(std::fesetround(original) == 0,
            "cannot restore floating-point rounding mode");
}

} // namespace

int main()
{
    try {
        test_human_scale_and_pairing();
        test_invalid_and_overflow_inputs();
        test_rounding_mode_is_part_of_the_contract();
    } catch (const std::exception &error) {
        std::cerr << "depth_count_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
