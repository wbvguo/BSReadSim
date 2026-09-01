#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "methdb.h"
#include "protocol.h"

namespace {

using htsim::methdb::CatalogError;
using htsim::methdb::CgmapRecord;
using htsim::methdb::ContextShapes;
using htsim::methdb::MethylationCatalog;
using htsim::methdb::ShapePair;
using htsim::model::MethylationContext;
using htsim::model::MethylationSource;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

std::uint16_t q(float probability)
{
    return htsim::methdb::probability_to_u16(probability);
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const CatalogError &) {
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

ContextShapes shapes()
{
    return {{2.0, 5.0}, {3.0, 4.0}, {5.0, 2.0}};
}

const ShapePair &shape_for_test(
    MethylationContext context,
    const ContextShapes &value)
{
    switch (context) {
    case MethylationContext::cg_c:
    case MethylationContext::cg_g: return value.cg;
    case MethylationContext::chg_c:
    case MethylationContext::chg_g: return value.chg;
    case MethylationContext::chh_c:
    case MethylationContext::chh_g: return value.chh;
    }
    throw std::runtime_error("unknown context in fixture");
}

void test_catalog_order_shapes_and_reuse()
{
    const auto bases = encode("ACGACAGCAATCGTTGAA");
    const ContextShapes configured = shapes();
    constexpr std::uint32_t contig_index = 3U;
    const MethylationCatalog catalog(
        bases, contig_index, UINT64_C(0x123456789abcdef0), true, configured);
    require(!catalog.sites().empty(), "catalog contains no methylation sites");

    std::uint32_t previous = 0;
    bool first = true;
    bool saw_cg = false;
    bool saw_chg = false;
    bool saw_chh = false;
    for (const auto &site : catalog.sites()) {
        require(first || site.reference_position > previous,
                "catalog sites are not unique and sorted");
        first = false;
        previous = site.reference_position;
        const ShapePair &shape = shape_for_test(site.context, configured);
        const float expected = htsim::beta_sampler::sample_beta(
            UINT64_C(0x123456789abcdef0), contig_index,
            site.reference_position, shape.alpha, shape.beta);
        require(site.probability_u16 == q(expected),
                "catalog used the wrong shape or RNG address");
        saw_cg = saw_cg || site.context == MethylationContext::cg_c
            || site.context == MethylationContext::cg_g;
        saw_chg = saw_chg || site.context == MethylationContext::chg_c
            || site.context == MethylationContext::chg_g;
        saw_chh = saw_chh || site.context == MethylationContext::chh_c
            || site.context == MethylationContext::chh_g;
    }
    require(saw_cg && saw_chg && saw_chh,
            "catalog fixture did not exercise all shape classes");

    const auto range = catalog.sites_in_range(2, 12);
    for (auto site = range.first; site != range.second; ++site) {
        require(site->reference_position >= 2 && site->reference_position < 12,
                "range query escaped its half-open interval");
    }
    require(range.first != range.second, "range query unexpectedly returned no sites");

    const MethylationCatalog repeated(
        bases, contig_index, UINT64_C(0x123456789abcdef0), true, configured);
    require(repeated.sites().size() == catalog.sites().size(),
            "repeated catalog size changed");
    for (std::size_t index = 0; index < catalog.sites().size(); ++index) {
        require(repeated.sites()[index].reference_position
                    == catalog.sites()[index].reference_position
                    && repeated.sites()[index].context
                        == catalog.sites()[index].context
                    && repeated.sites()[index].probability_u16
                        == catalog.sites()[index].probability_u16,
                "repeated catalog changed a genomic site");
    }
}

void test_cpg_filter_and_invalid_inputs()
{
    const auto bases = encode("ACGACAGCAATCGTTGAA");
    const MethylationCatalog cpg_only(bases, 0U, 1, false, shapes());
    require(!cpg_only.sites().empty(), "CpG-only catalog is empty");
    for (const auto &site : cpg_only.sites()) {
        require(site.context == MethylationContext::cg_c
                    || site.context == MethylationContext::cg_g,
                "non-CpG site survived the filter");
    }

    ContextShapes invalid = shapes();
    invalid.chg.alpha = std::numeric_limits<double>::infinity();
    require_error(
        [&] {(void)MethylationCatalog(bases, 0U, 1, true, invalid);},
        "non-finite shape was accepted");
    require_error(
        [&] {(void)MethylationCatalog({0, 5, 1}, 0U, 1, true, shapes());},
        "invalid base encoding was accepted");
    const MethylationCatalog empty({}, 0U, 1, true, shapes());
    require(empty.sites().empty(), "empty contig produced a methylation site");
    require_error(
        [&] {(void)cpg_only.sites_in_range(5, 4);},
        "reversed site range was accepted");
}

void test_cgmap_overlay_and_na_fallback()
{
    const auto bases = encode("ACGACAGCAATCGTTGAA");
    const std::vector<CgmapRecord> records = {
        {1U, q(0.125F), MethylationContext::cg_c, true, 2U},
        {2U, 0U, MethylationContext::cg_g, false, 2U},
        {4U, q(0.75F), MethylationContext::chg_c, true, 0U},
    };
    const MethylationCatalog catalog(
        bases, 2U, 17U, true, shapes(), &records);
    const auto find_position = [&](std::uint32_t position) {
        return std::find_if(
            catalog.sites().begin(), catalog.sites().end(),
            [position](const auto &site) {
                return site.reference_position == position;
            });
    };
    const auto cg_c = find_position(1U);
    const auto cg_g = find_position(2U);
    const auto chg_c = find_position(4U);
    require(
        cg_c != catalog.sites().end()
            && cg_c->methylation_source == MethylationSource::cgmap
            && cg_c->probability_u16 == q(0.125F),
        "defined CGmap value did not replace the Beta level");
    require(
        cg_g != catalog.sites().end()
            && cg_g->methylation_source == MethylationSource::beta,
        "CGmap na did not retain the Beta source");
    require(
        chg_c != catalog.sites().end()
            && chg_c->methylation_source == MethylationSource::cgmap
            && chg_c->probability_u16 == q(0.75F),
        "non-CpG CGmap value was not overlaid");

    const MethylationCatalog cpg_only(
        bases, 2U, 17U, false, shapes(), &records);
    require(
        std::all_of(
            cpg_only.sites().begin(), cpg_only.sites().end(),
            [](const auto &site) {
                return site.context == MethylationContext::cg_c
                    || site.context == MethylationContext::cg_g;
            }),
        "CGmap overlay bypassed the CpG-only filter");

    const std::vector<CgmapRecord> mismatched = {
        {1U, q(0.5F), MethylationContext::chh_c, true, 0U},
    };
    require_error(
        [&] {
            (void)MethylationCatalog(
                bases, 2U, 17U, true, shapes(), &mismatched);
        },
        "reference-inconsistent normalized CGmap record was accepted");
}

void test_cgmap_context_pool_and_beta_fallback()
{
    const auto bases = encode("ACGACAGCAATCGTTGAA");
    const std::vector<CgmapRecord> records = {
        {1U, q(0.125F), MethylationContext::cg_c, true, 2U},
        {2U, q(0.875F), MethylationContext::cg_g, true, 2U},
        {4U, q(0.75F), MethylationContext::chg_c, true, 0U},
    };
    constexpr std::uint64_t seed = 17U;
    constexpr std::uint32_t contig_index = 4U;
    const ContextShapes configured = shapes();
    const htsim::methdb::CgmapPool pool(records);
    const MethylationCatalog catalog(
        bases, contig_index, seed, true, configured, &records, true);
    bool saw_pooled = false;
    bool saw_beta_fallback = false;
    for (const auto &site : catalog.sites()) {
        const auto expected = pool.sample(
            site.context,
            seed,
            contig_index,
            htsim::methdb::reference_site_entity(
                site.reference_position));
        if (expected) {
            require(
                site.methylation_source == MethylationSource::pooled_cgmap
                    && site.probability_u16 == *expected,
                "reference catalog lost a typed CGmap pool draw");
            saw_pooled = true;
        } else {
            const ShapePair &shape = shape_for_test(site.context, configured);
            require(
                site.methylation_source == MethylationSource::beta
                    && site.probability_u16
                        == q(htsim::beta_sampler::sample_beta(
                            seed,
                            contig_index,
                            site.reference_position,
                            shape.alpha,
                            shape.beta)),
                "empty CGmap context class did not retain Beta fallback");
            saw_beta_fallback = true;
        }
    }
    require(saw_pooled && saw_beta_fallback,
            "pool fixture did not exercise pooled and fallback sources");
    require_error(
        [&] {
            (void)MethylationCatalog(
                bases, contig_index, seed, true, shapes(), nullptr, true);
        },
        "catalog accepted CGmap pooling without records");
}

} // namespace

int main()
{
    try {
        test_catalog_order_shapes_and_reuse();
        test_cpg_filter_and_invalid_inputs();
        test_cgmap_overlay_and_na_fallback();
        test_cgmap_context_pool_and_beta_fallback();
    } catch (const std::exception &error) {
        std::cerr << "methdb_catalog_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
