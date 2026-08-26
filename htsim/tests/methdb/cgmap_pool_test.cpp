#include "methdb.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "types.h"
#include "protocol.h"
#include "utilities.h"

namespace {

using htsim::methdb::CgmapPool;
using htsim::methdb::CgmapPoolError;
using htsim::methdb::CgmapRecord;
using htsim::model::MethylationContext;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const CgmapPoolError &) {
        return;
    }
    throw std::runtime_error(message);
}

std::vector<CgmapRecord> records()
{
    return {
        {1U, 0.125F, MethylationContext::cg_c, true, 2U},
        {2U, 0.875F, MethylationContext::cg_g, true, 2U},
        {4U, 0.25F, MethylationContext::chg_c, true, 0U},
        {5U, 0.0F, MethylationContext::chg_g, false, 0U},
        {6U, 0.5F, MethylationContext::chh_c, true, 0U},
        {7U, 0.75F, MethylationContext::chh_g, true, 0U},
    };
}

void test_context_classes_and_exact_addresses()
{
    const CgmapPool pool(records());
    require(pool.size(MethylationContext::cg_c) == 2U
                && pool.size(MethylationContext::cg_g) == 2U
                && pool.size(MethylationContext::chg_c) == 1U
                && pool.size(MethylationContext::chg_g) == 1U
                && pool.size(MethylationContext::chh_c) == 2U
                && pool.size(MethylationContext::chh_g) == 2U,
            "C/G-oriented contexts did not share their typed class pool");

    constexpr std::uint64_t seed = UINT64_C(0x123456789abcdef0);
    constexpr std::uint32_t contig_index = 4U;
    const auto entity0 = htsim::methdb::reference_site_entity(0U);
    const auto entity1 = htsim::methdb::reference_site_entity(1U);
    const auto inserted = htsim::methdb::insertion_site_entity(
        17U,
        2U,
        htsim::model::HaplotypeMask::haplotype_2,
        1U);

    const auto cg0 = pool.sample(
        MethylationContext::cg_c, seed, contig_index, entity0);
    const auto cg1 = pool.sample(
        MethylationContext::cg_g, seed, contig_index, entity1);
    const auto chg = pool.sample(
        MethylationContext::chg_g, seed, contig_index, inserted);
    const auto chh = pool.sample(
        MethylationContext::chh_c, seed, contig_index, inserted);
    require(cg0.has_value() && cg1.has_value() && chg == 0.25F
                && chh.has_value(),
            "CGmap pool failed to sample a populated typed class");

    const std::uint64_t key = htsim::rng::derive_key(
        seed, htsim::rng::Stage::methylation_level, contig_index);
    const auto expected = [&](auto entity, const std::vector<float> &values) {
        const std::uint64_t index = htsim::rng::bounded_integer(
            key,
            entity.value(),
            htsim::methdb::cgmap_pool_local_index,
            values.size());
        return values.at(static_cast<std::size_t>(index));
    };
    require(*cg0 == expected(entity0, {0.125F, 0.875F})
                && *cg1 == expected(entity1, {0.125F, 0.875F})
                && *chh == expected(inserted, {0.5F, 0.75F}),
            "CGmap pool did not use the frozen uint64 address contract");

    // Temporary implementation vector output is intentionally kept tiny; the
    // asserted literal vector below freezes it independently of rng helpers.
    require(*cg0 == 0.125F && *cg1 == 0.125F && *chh == 0.5F,
            "CGmap pool exact selection vector changed");
}

void test_empty_class_falls_back_without_a_draw()
{
    const CgmapPool pool({
        {3U, 0.4F, MethylationContext::cg_c, true, 2U},
    });
    const auto absent = pool.sample(
        MethylationContext::chg_c,
        0U,
        0U,
        htsim::methdb::reference_site_entity(0U));
    require(!absent.has_value(),
            "empty context class did not request Beta fallback");
}

void test_order_and_entity_stability()
{
    const CgmapPool pool(records());
    std::vector<std::uint32_t> positions(128U);
    for (std::uint32_t index = 0U; index < positions.size(); ++index) {
        positions[index] = index;
    }
    const auto sample = [&](std::uint32_t position) {
        return *pool.sample(
            MethylationContext::cg_c,
            UINT64_MAX,
            5U,
            htsim::methdb::reference_site_entity(position));
    };
    std::vector<float> forward;
    forward.reserve(positions.size());
    for (const std::uint32_t position : positions) {
        forward.push_back(sample(position));
    }
    std::reverse(positions.begin(), positions.end());
    for (const std::uint32_t position : positions) {
        require(sample(position) == forward[position],
                "CGmap pool depended on traversal order");
    }

    const auto reference = pool.sample(
        MethylationContext::chh_c,
        9U,
        5U,
        htsim::methdb::reference_site_entity(17U));
    const auto variant = pool.sample(
        MethylationContext::chh_c,
        9U,
        5U,
        htsim::methdb::variant_reference_site_entity(
            17U,
            htsim::model::HaplotypeMask::haplotype_1,
            0U));
    require(reference.has_value() && variant.has_value(),
            "typed site entities did not produce pool samples");
}

void test_invalid_normalized_records_fail_closed()
{
    auto invalid = records();
    std::swap(invalid[0], invalid[1]);
    require_error([&] {(void)CgmapPool(invalid);},
                  "unsorted pool records were accepted");

    invalid = records();
    invalid[1].reference_position = invalid[0].reference_position;
    require_error([&] {(void)CgmapPool(invalid);},
                  "duplicate pool positions were accepted");

    invalid = records();
    invalid[0].context = static_cast<MethylationContext>(0U);
    require_error([&] {(void)CgmapPool(invalid);},
                  "unknown context was accepted");

    invalid = records();
    invalid[0].methylation_probability = 2.0F;
    require_error([&] {(void)CgmapPool(invalid);},
                  "probability above one was accepted");

    invalid = records();
    invalid[3].methylation_probability = 0.25F;
    require_error([&] {(void)CgmapPool(invalid);},
                  "undefined nonzero probability was accepted");

    invalid = records();
    invalid[0].dinucleotide_second = 4U;
    require_error([&] {(void)CgmapPool(invalid);},
                  "invalid dinucleotide encoding was accepted");

    const CgmapPool empty({});
    require(empty.size(MethylationContext::cg_c) == 0U,
            "empty per-contig pool was rejected instead of falling back");
    require_error(
        [&] {empty.size(static_cast<MethylationContext>(0U));},
        "unknown lookup context was accepted");
}

} // namespace

int main()
{
    try {
        test_context_classes_and_exact_addresses();
        test_empty_class_falls_back_without_a_draw();
        test_order_and_entity_stability();
        test_invalid_normalized_records_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
