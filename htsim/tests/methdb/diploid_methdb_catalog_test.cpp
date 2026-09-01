#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "methdb.h"
#include "variant.h"
#include "types.h"
#include "protocol.h"
#include "reference.h"

namespace {

using htsim::model::HaplotypeMask;
using htsim::methdb::AsmRecord;
using htsim::methdb::ContextShapes;
using htsim::methdb::CgmapRecord;
using htsim::methdb::DiploidCatalogError;
using htsim::methdb::DiploidMethylationCatalog;
using htsim::methdb::DiploidRuntimeArrays;
using htsim::methdb::MethylationCatalog;
using htsim::methdb::ProbabilityU16;
using htsim::methdb::RuntimeSite;
using htsim::model::Bases;
using htsim::model::MethylationAllele;
using htsim::model::MethylationContext;
using htsim::model::MethylationSource;
using htsim::model::VariantKind;
using htsim::reference::Contig;
using htsim::variant::ContigVariants;
using htsim::variant::Variant;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

ProbabilityU16 quantized(float probability)
{
    return htsim::methdb::probability_to_u16(probability);
}

float quantized_float(float probability)
{
    return htsim::methdb::probability_from_u16(quantized(probability));
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const DiploidCatalogError &) {
        return;
    }
    throw std::runtime_error(message);
}

Bases encode(const std::string &text)
{
    Bases bases;
    for (const char base : text) {
        switch (base) {
        case 'A': bases.push_back(0U); break;
        case 'C': bases.push_back(1U); break;
        case 'G': bases.push_back(2U); break;
        case 'T': bases.push_back(3U); break;
        case 'N': bases.push_back(4U); break;
        default: throw std::runtime_error("invalid test base");
        }
    }
    return bases;
}

Contig make_contig(const std::string &sequence)
{
    Contig contig;
    contig.index = 0;
    contig.name = "chrDiploid";
    contig.bases = encode(sequence);
    contig.length = contig.bases.size();
    return contig;
}

ContextShapes shapes()
{
    return {{2.0, 5.0}, {3.0, 4.0}, {5.0, 2.0}};
}

std::uint32_t insertion_key(
    std::uint32_t event_ordinal,
    std::uint8_t insertion_offset)
{
    return (event_ordinal << 2U) | insertion_offset;
}

const RuntimeSite *find(
    const std::vector<RuntimeSite> &sites,
    std::uint32_t key)
{
    const RuntimeSite query = static_cast<RuntimeSite>(key) << 32U;
    const auto found = std::lower_bound(sites.begin(), sites.end(), query);
    return found == sites.end() || htsim::methdb::runtime_site_key(*found) != key
        ? nullptr : &*found;
}

bool same_arrays(
    const DiploidRuntimeArrays &left,
    const DiploidRuntimeArrays &right)
{
    return left.reference_shared == right.reference_shared
        && left.reference_haplotypes == right.reference_haplotypes
        && left.insertion_shared == right.insertion_shared
        && left.insertion_haplotypes == right.insertion_haplotypes;
}

std::vector<Variant> variant_fixture()
{
    return {
        {0, 3, 4, VariantKind::snv, encode("G"), encode("A"),
         HaplotypeMask::haplotype_1},
        {0, 6, 6, VariantKind::insertion, {}, encode("CG"),
         HaplotypeMask::both},
        {0, 7, 9, VariantKind::deletion, encode("CG"), {},
         HaplotypeMask::haplotype_2},
    };
}

void test_reference_only_parity()
{
    const Contig contig = make_contig("ACGACAGCAATCGTTGAA");
    const ContigVariants variants(contig.bases, {}, contig.index);
    const ContextShapes configured = shapes();
    const MethylationCatalog baseline(
        contig.bases, contig.index, 41U, true, configured);
    const DiploidMethylationCatalog diploid(
        contig, variants, 41U, true, configured);
    const DiploidRuntimeArrays &runtime = diploid.runtime_arrays();
    require(runtime.reference_shared.size() == baseline.sites().size()
                && runtime.insertion_shared.empty()
                && runtime.reference_haplotypes[0].empty()
                && runtime.reference_haplotypes[1].empty()
                && runtime.insertion_haplotypes[0].empty()
                && runtime.insertion_haplotypes[1].empty(),
            "reference-only diploid catalog did not collapse to shared sites");
    for (std::size_t index = 0; index < baseline.sites().size(); ++index) {
        const auto &old_site = baseline.sites()[index];
        const RuntimeSite new_site = runtime.reference_shared[index];
        require(htsim::methdb::runtime_site_key(new_site)
                    == old_site.reference_position
                    && htsim::methdb::runtime_site_context(new_site)
                        == old_site.context
                    && htsim::methdb::runtime_site_allele(new_site)
                        == MethylationAllele::shared
                    && htsim::methdb::runtime_site_probability(new_site)
                        == old_site.probability_u16,
                "reference-only site or Beta address changed");
    }

    const auto projection = htsim::haplotype::project_interval(
        contig,
        variants,
        0U,
        0U,
        static_cast<std::uint32_t>(contig.length));
    const auto projected = diploid.sites_for_projection(projection);
    require(projected.size() == baseline.sites().size(),
            "reference-only site projection lost a site");
    for (std::size_t index = 0; index < projected.size(); ++index) {
        require(projected[index].site_index == index
                    && projected[index].reference_pos
                        == baseline.sites()[index].reference_position
                    && projected[index].context == baseline.sites()[index].context
                    && projected[index].methylation_probability
                        == htsim::methdb::probability_from_u16(
                            baseline.sites()[index].probability_u16),
                "reference-only protocol site projection changed");
    }

    const std::vector<CgmapRecord> records = {
        {1U, quantized(0.125F), MethylationContext::cg_c, true, 2U},
        {2U, 0U, MethylationContext::cg_g, false, 2U},
        {4U, quantized(0.75F), MethylationContext::chg_c, true, 0U},
    };
    const MethylationCatalog cgmap_baseline(
        contig.bases, contig.index, 41U, true, configured, &records);
    const DiploidMethylationCatalog cgmap_diploid(
        contig, variants, 41U, true, configured, &records);
    const DiploidRuntimeArrays &cgmap_runtime =
        cgmap_diploid.runtime_arrays();
    require(
        cgmap_runtime.reference_shared.size()
            == cgmap_baseline.sites().size(),
        "event-free CGmap diploid catalog changed the site count");
    for (std::size_t index = 0U;
         index < cgmap_baseline.sites().size();
        ++index) {
        const auto &reference_site = cgmap_baseline.sites()[index];
        const RuntimeSite diploid_site = cgmap_runtime.reference_shared[index];
        require(
            htsim::methdb::runtime_site_key(diploid_site)
                    == reference_site.reference_position
                && htsim::methdb::runtime_site_context(diploid_site)
                    == reference_site.context
                && htsim::methdb::runtime_site_source(diploid_site)
                    == reference_site.methylation_source
                && htsim::methdb::runtime_site_probability(diploid_site)
                    == reference_site.probability_u16,
            "event-free CGmap reference/diploid catalogs diverged");
    }
}

void test_variant_context_entities_and_alleles()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const auto events = variant_fixture();
    const ContigVariants variants(contig.bases, events, contig.index);
    const ContextShapes configured = shapes();
    constexpr std::uint64_t seed = 73U;
    const DiploidMethylationCatalog catalog(
        contig, variants, seed, true, configured);
    const DiploidRuntimeArrays &runtime = catalog.runtime_arrays();

    const RuntimeSite *alt_context = find(
        runtime.reference_haplotypes[0], 2U);
    const RuntimeSite *ref_context = find(
        runtime.reference_haplotypes[1], 2U);
    require(alt_context && ref_context
                && htsim::methdb::runtime_site_context(*alt_context)
                    == MethylationContext::chh_c
                && htsim::methdb::runtime_site_allele(*alt_context)
                    == MethylationAllele::alternate_haplotype
                && htsim::methdb::runtime_site_context(*ref_context)
                    == MethylationContext::cg_c
                && htsim::methdb::runtime_site_allele(*ref_context)
                    == MethylationAllele::reference_haplotype,
            "heterozygous SNV did not split reference and ALT contexts");
    require(htsim::methdb::runtime_site_probability(*alt_context)
                == quantized(htsim::beta_sampler::sample_beta_for_site(
                    seed,
                    contig.index,
                    htsim::methdb::variant_reference_site_entity(
                        2U, HaplotypeMask::haplotype_1, 0U),
                    configured.chh.alpha,
                    configured.chh.beta))
                && htsim::methdb::runtime_site_probability(*ref_context)
                    == quantized(htsim::beta_sampler::sample_beta_for_site(
                        seed,
                        contig.index,
                        htsim::methdb::reference_site_entity(2U),
                        configured.cg.alpha,
                        configured.cg.beta)),
            "heterozygous context used the wrong 64-bit Beta identity");

    const RuntimeSite *inserted_c = find(
        runtime.insertion_shared, insertion_key(1U, 0U));
    const RuntimeSite *inserted_g = find(
        runtime.insertion_shared, insertion_key(1U, 1U));
    require(inserted_c && inserted_g
                && htsim::methdb::runtime_site_context(*inserted_c)
                    == MethylationContext::cg_c
                && htsim::methdb::runtime_site_context(*inserted_g)
                    == MethylationContext::cg_g
                && htsim::methdb::runtime_site_allele(*inserted_c)
                    == MethylationAllele::shared
                && htsim::methdb::runtime_site_allele(*inserted_g)
                    == MethylationAllele::shared,
            "homozygous inserted CpG was not shared by both haplotypes");
    require(htsim::methdb::runtime_site_probability(*inserted_c)
                == quantized(htsim::beta_sampler::sample_beta_for_site(
                    seed,
                    contig.index,
                    htsim::methdb::insertion_site_entity(
                        1U, 0U, HaplotypeMask::both, 0U),
                    configured.cg.alpha,
                    configured.cg.beta)),
            "inserted CpG used the wrong event/offset Beta identity");

    for (const std::uint32_t position : {7U, 8U}) {
        const RuntimeSite *site = find(
            runtime.reference_haplotypes[0], position);
        require(site
                    && htsim::methdb::runtime_site_allele(*site)
                        == MethylationAllele::reference_haplotype,
                "heterozygous deletion did not retain reference-haplotype ownership");
        require(!find(
                    runtime.reference_haplotypes[1], position),
                "deleted haplotype retained a methylation site");
    }
}

void test_fragment_boundary_uses_complete_haplotype_context()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const ContigVariants variants(
        contig.bases, variant_fixture(), contig.index);
    const DiploidMethylationCatalog catalog(
        contig, variants, 73U, true, shapes());
    const auto haplotype_0 = htsim::haplotype::project_interval(
        contig, variants, 0U, 2U, 3U);
    const auto haplotype_1 = htsim::haplotype::project_interval(
        contig, variants, 1U, 2U, 3U);
    require(haplotype_0.variants.empty()
                && haplotype_1.variants.empty(),
            "boundary fixture unexpectedly included the neighboring SNV");
    const auto alt_sites = catalog.sites_for_projection(haplotype_0);
    const auto ref_sites = catalog.sites_for_projection(haplotype_1);
    require(alt_sites.size() == 1U && ref_sites.size() == 1U
                && alt_sites[0].template_offset == 0U
                && alt_sites[0].context == MethylationContext::chh_c
                && alt_sites[0].allele
                    == MethylationAllele::alternate_haplotype
                && ref_sites[0].context == MethylationContext::cg_c
                && ref_sites[0].allele
                    == MethylationAllele::reference_haplotype,
            "fragment edge was classified from its local slice");
    const RuntimeSite *catalog_site = find(
        catalog.runtime_arrays().reference_haplotypes[0], 2U);
    require(catalog_site
                && alt_sites[0].methylation_probability
                    == htsim::methdb::probability_from_u16(
                        htsim::methdb::runtime_site_probability(*catalog_site)),
            "overlapping fragments did not reuse the same site probability");
}

void test_deletion_joins_previously_separated_context()
{
    const Contig contig = make_contig("CAAG");
    const std::vector<Variant> events = {
        {0, 1, 3, VariantKind::deletion, encode("AA"), {},
         HaplotypeMask::haplotype_1},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const DiploidMethylationCatalog catalog(
        contig, variants, 91U, true, shapes());
    const DiploidRuntimeArrays &runtime = catalog.runtime_arrays();
    const RuntimeSite *alt_c = find(runtime.reference_haplotypes[0], 0U);
    const RuntimeSite *alt_g = find(runtime.reference_haplotypes[0], 3U);
    const RuntimeSite *ref_c = find(runtime.reference_haplotypes[1], 0U);
    const RuntimeSite *ref_g = find(runtime.reference_haplotypes[1], 3U);
    require(alt_c && alt_g && ref_c && ref_g
                && htsim::methdb::runtime_site_context(*alt_c)
                    == MethylationContext::cg_c
                && htsim::methdb::runtime_site_context(*alt_g)
                    == MethylationContext::cg_g
                && htsim::methdb::runtime_site_context(*ref_c)
                    == MethylationContext::chh_c
                && htsim::methdb::runtime_site_context(*ref_g)
                    == MethylationContext::chh_g,
            "deletion did not recompute context across its joined boundary");

    const DiploidMethylationCatalog cpg_only(
        contig, variants, 91U, false, shapes());
    const DiploidRuntimeArrays &cpg_runtime = cpg_only.runtime_arrays();
    require(cpg_runtime.reference_haplotypes[0].size() == 2U
                && cpg_runtime.reference_haplotypes[1].empty()
                && cpg_runtime.reference_shared.empty()
                && cpg_runtime.insertion_shared.empty(),
            "CpG-only filter did not use the final haplotype context");
}

void test_cgmap_overlay_respects_haplotype_equivalence()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const ContigVariants variants(
        contig.bases, variant_fixture(), contig.index);
    const std::vector<CgmapRecord> records = {
        {2U, quantized(0.875F), MethylationContext::cg_c, true, 2U},
        {3U, 0U, MethylationContext::cg_g, false, 2U},
        {7U, quantized(0.625F), MethylationContext::cg_c, true, 2U},
    };
    const DiploidMethylationCatalog catalog(
        contig, variants, 73U, true, shapes(), &records);
    const DiploidRuntimeArrays &runtime = catalog.runtime_arrays();

    const RuntimeSite *alternate = find(
        runtime.reference_haplotypes[0], 2U);
    const RuntimeSite *reference = find(
        runtime.reference_haplotypes[1], 2U);
    require(
        alternate != nullptr && reference != nullptr
            && htsim::methdb::runtime_site_source(*alternate)
                == MethylationSource::beta
            && htsim::methdb::runtime_site_context(*alternate)
                == MethylationContext::chh_c
            && htsim::methdb::runtime_site_source(*reference)
                == MethylationSource::cgmap
            && htsim::methdb::runtime_site_context(*reference)
                == MethylationContext::cg_c
            && htsim::methdb::runtime_site_probability(*reference)
                == quantized(0.875F),
        "CGmap was not restricted to the reference-equivalent haplotype");

    const RuntimeSite *retained_after_deletion = find(
        runtime.reference_haplotypes[0], 7U);
    require(
        retained_after_deletion != nullptr
            && htsim::methdb::runtime_site_source(*retained_after_deletion)
                == MethylationSource::cgmap
            && htsim::methdb::runtime_site_probability(
                   *retained_after_deletion)
                == quantized(0.625F),
        "CGmap did not overlay the retained reference haplotype");

    const RuntimeSite *inserted = find(
        runtime.insertion_shared, insertion_key(1U, 0U));
    require(
        inserted != nullptr
            && htsim::methdb::runtime_site_source(*inserted)
                == MethylationSource::beta,
        "CGmap leaked onto a variant-created methylation site");

    const auto projection = htsim::haplotype::project_interval(
        contig,
        variants,
        1U,
        0U,
        static_cast<std::uint32_t>(contig.length));
    const auto projected = catalog.sites_for_projection(projection);
    const auto projected_reference = std::find_if(
        projected.begin(), projected.end(),
        [](const auto &site) {return site.reference_pos == 2;});
    require(
        projected_reference != projected.end()
            && projected_reference->methylation_source == MethylationSource::cgmap
            && projected_reference->methylation_probability
                == quantized_float(0.875F),
        "protocol projection lost CGmap provenance");
}

void test_cgmap_pool_uses_typed_variant_entities()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const ContigVariants variants(
        contig.bases, variant_fixture(), contig.index);
    const std::vector<CgmapRecord> records = {
        {2U, quantized(0.125F), MethylationContext::cg_c, true, 2U},
        {3U, quantized(0.875F), MethylationContext::cg_g, true, 2U},
    };
    constexpr std::uint64_t seed = 73U;
    const htsim::methdb::CgmapPool pool(records);
    const DiploidMethylationCatalog catalog(
        contig,
        variants,
        seed,
        true,
        shapes(),
        &records,
        nullptr,
        true);
    const DiploidRuntimeArrays &runtime = catalog.runtime_arrays();

    const RuntimeSite *alternate = find(
        runtime.reference_haplotypes[0], 2U);
    const RuntimeSite *reference = find(
        runtime.reference_haplotypes[1], 2U);
    require(
        alternate != nullptr && reference != nullptr
            && htsim::methdb::runtime_site_context(*alternate)
                == MethylationContext::chh_c
            && htsim::methdb::runtime_site_source(*alternate)
                == MethylationSource::beta
            && htsim::methdb::runtime_site_context(*reference)
                == MethylationContext::cg_c
            && htsim::methdb::runtime_site_source(*reference)
                == MethylationSource::pooled_cgmap,
        "CGmap pool ignored the typed context on a split haplotype site");
    const auto expected_reference = pool.sample(
        htsim::methdb::runtime_site_context(*reference),
        seed,
        contig.index,
        htsim::methdb::reference_site_entity(2U));
    require(
        expected_reference.has_value()
            && htsim::methdb::runtime_site_probability(*reference)
                == *expected_reference,
        "reference-equivalent haplotype used the wrong pool address");

    for (const std::uint8_t offset : {std::uint8_t{0U}, std::uint8_t{1U}}) {
        const RuntimeSite *inserted = find(
            runtime.insertion_shared, insertion_key(1U, offset));
        const auto entity = htsim::methdb::insertion_site_entity(
            1U, offset, HaplotypeMask::both, 0U);
        require(inserted != nullptr
                    && htsim::methdb::runtime_site_source(*inserted)
                        == MethylationSource::pooled_cgmap,
                "variant-created CpG did not use the typed context pool");
        const auto expected = pool.sample(
            htsim::methdb::runtime_site_context(*inserted),
            seed,
            contig.index,
            entity);
        require(expected.has_value()
                    && htsim::methdb::runtime_site_probability(*inserted)
                        == *expected,
                "inserted site used the wrong uint64 pool entity");
    }
}

void test_asm_overlay_uses_typed_haplotype_mask()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const std::vector<Variant> events = {
        {0U, 9U, 10U, VariantKind::snv, encode("T"), encode("A"),
         HaplotypeMask::haplotype_1},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const std::vector<CgmapRecord> cgmap = {
        {2U, quantized(0.5F), MethylationContext::cg_c, true, 2U},
    };
    const std::vector<AsmRecord> asm_records = {
        {2U, 9U, quantized(0.2F), quantized(0.8F),
         MethylationContext::cg_c, 2U, 3U, 0U},
    };
    const DiploidMethylationCatalog catalog(
        contig, variants, 73U, true, shapes(), &cgmap, &asm_records);
    DiploidMethylationCatalog reused(
        contig, variants, 73U, true, shapes(), &cgmap, nullptr);
    reused.apply_asm_layer(events, asm_records);
    require(
        same_arrays(reused.runtime_arrays(), catalog.runtime_arrays()),
        "save-and-simulate ASM reuse changed the runtime catalog");
    const std::uint32_t origin = 2U;
    const DiploidRuntimeArrays &runtime = catalog.runtime_arrays();
    require(
        find(runtime.reference_shared, origin) == nullptr,
        "ASM target remained in the shared catalog");
    const RuntimeSite *alternate = find(
        runtime.reference_haplotypes[0], origin);
    const RuntimeSite *reference = find(
        runtime.reference_haplotypes[1], origin);
    require(
        alternate != nullptr && reference != nullptr
            && htsim::methdb::runtime_site_source(*alternate)
                == MethylationSource::asm_source
            && htsim::methdb::runtime_site_allele(*alternate)
                == MethylationAllele::alternate_haplotype
            && htsim::methdb::runtime_site_probability(*alternate)
                == quantized(0.8F)
            && htsim::methdb::runtime_site_source(*reference)
                == MethylationSource::asm_source
            && htsim::methdb::runtime_site_allele(*reference)
                == MethylationAllele::reference_haplotype
            && htsim::methdb::runtime_site_probability(*reference)
                == quantized(0.2F),
        "ASM probabilities did not follow HaplotypeMask value 1");

    for (const std::uint8_t haplotype : {std::uint8_t{0U}, std::uint8_t{1U}}) {
        const auto projection = htsim::haplotype::project_interval(
            contig,
            variants,
            haplotype,
            0U,
            static_cast<std::uint32_t>(contig.length));
        const auto sites = catalog.sites_for_projection(projection);
        const auto target = std::find_if(
            sites.begin(), sites.end(),
            [](const auto &site) {return site.reference_pos == 2;});
        require(
            target != sites.end()
                && target->methylation_source == MethylationSource::asm_source
                && target->methylation_probability
                    == quantized_float(haplotype == 0U ? 0.8F : 0.2F),
            "protocol projection lost ASM allele ownership");
    }

    const std::vector<Variant> reversed_events = {
        {0U, 9U, 10U, VariantKind::snv, encode("T"), encode("A"),
         HaplotypeMask::haplotype_2},
    };
    const ContigVariants reversed_variants(
        contig.bases, reversed_events, contig.index);
    const DiploidMethylationCatalog reversed(
        contig, reversed_variants, 73U, true, shapes(), nullptr, &asm_records);
    const RuntimeSite *reversed_reference = find(
        reversed.runtime_arrays().reference_haplotypes[0], origin);
    const RuntimeSite *reversed_alternate = find(
        reversed.runtime_arrays().reference_haplotypes[1], origin);
    require(
        reversed_reference != nullptr && reversed_alternate != nullptr
            && htsim::methdb::runtime_site_allele(*reversed_reference)
                == MethylationAllele::reference_haplotype
            && htsim::methdb::runtime_site_probability(*reversed_reference)
                == quantized(0.2F)
            && htsim::methdb::runtime_site_allele(*reversed_alternate)
                == MethylationAllele::alternate_haplotype
            && htsim::methdb::runtime_site_probability(*reversed_alternate)
                == quantized(0.8F),
        "ASM probabilities did not follow HaplotypeMask value 2");

    const DiploidMethylationCatalog pooled(
        contig,
        variants,
        73U,
        true,
        shapes(),
        &cgmap,
        &asm_records,
        true);
    const DiploidRuntimeArrays &pooled_runtime = pooled.runtime_arrays();
    require(
        find(pooled_runtime.reference_shared, origin) == nullptr,
        "ASM target remained shared after CGmap pooling");
    const RuntimeSite *pooled_alternate = find(
        pooled_runtime.reference_haplotypes[0], origin);
    const RuntimeSite *pooled_reference = find(
        pooled_runtime.reference_haplotypes[1], origin);
    require(
        pooled_alternate != nullptr && pooled_reference != nullptr
            && htsim::methdb::runtime_site_source(*pooled_alternate)
                == MethylationSource::asm_source
            && htsim::methdb::runtime_site_probability(*pooled_alternate)
                == quantized(0.8F)
            && htsim::methdb::runtime_site_source(*pooled_reference)
                == MethylationSource::asm_source
            && htsim::methdb::runtime_site_probability(*pooled_reference)
                == quantized(0.2F),
        "ASM did not retain precedence over the typed CGmap pool");
}

void test_asm_links_fail_closed()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const std::vector<AsmRecord> records = {
        {2U, 9U, quantized(0.2F), quantized(0.8F),
         MethylationContext::cg_c, 2U, 3U, 0U},
    };
    const ContigVariants no_variants(contig.bases, {}, contig.index);
    require_error(
        [&] {
            (void)DiploidMethylationCatalog(
                contig, no_variants, 73U, true, shapes(), nullptr, &records);
        },
        "ASM row without a linked VCF event was accepted");

    const std::vector<Variant> homozygous_events = {
        {0U, 9U, 10U, VariantKind::snv, encode("T"), encode("A"),
         HaplotypeMask::both},
    };
    const ContigVariants homozygous(
        contig.bases, homozygous_events, contig.index);
    require_error(
        [&] {
            (void)DiploidMethylationCatalog(
                contig, homozygous, 73U, true, shapes(), nullptr, &records);
        },
        "ASM row linked to a homozygous ALT event was accepted");

    const std::vector<Variant> heterozygous_events = {
        {0U, 9U, 10U, VariantKind::snv, encode("T"), encode("A"),
         HaplotypeMask::haplotype_1},
    };
    const ContigVariants heterozygous(
        contig.bases, heterozygous_events, contig.index);
    auto wrong_allele = records;
    wrong_allele.front().linked_alternate_base = 1U;
    require_error(
        [&] {
            (void)DiploidMethylationCatalog(
                contig,
                heterozygous,
                73U,
                true,
                shapes(),
                nullptr,
                &wrong_allele);
        },
        "ASM row with mismatched linked ALT was accepted");

    const Contig context_contig = make_contig("AACGTAACGTT");
    const std::vector<Variant> context_events = {
        {0U, 3U, 4U, VariantKind::snv, encode("G"), encode("A"),
         HaplotypeMask::haplotype_1},
    };
    const ContigVariants context_variants(
        context_contig.bases, context_events, context_contig.index);
    const std::vector<AsmRecord> divergent = {
        {2U, 3U, quantized(0.2F), quantized(0.8F),
         MethylationContext::cg_c, 2U, 2U, 0U},
    };
    require_error(
        [&] {
            (void)DiploidMethylationCatalog(
                context_contig,
                context_variants,
                73U,
                true,
                shapes(),
                nullptr,
                &divergent);
        },
        "ASM context-divergent target was accepted");
}

void test_invalid_inputs_fail_closed()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const ContigVariants variants(
        contig.bases, variant_fixture(), contig.index);
    const DiploidMethylationCatalog catalog(
        contig, variants, 73U, true, shapes());

    DiploidRuntimeArrays invalid_runtime;
    invalid_runtime.reference_shared.push_back(
        htsim::methdb::pack_runtime_site(
            static_cast<std::uint32_t>(contig.length),
            0U,
            MethylationContext::cg_c,
            MethylationSource::beta,
            MethylationAllele::shared,
            true));
    require_error(
        [&] {
            (void)DiploidMethylationCatalog(
                contig.index,
                static_cast<std::uint32_t>(contig.length),
                invalid_runtime);
        },
        "out-of-contig packed runtime site was accepted");

    DiploidRuntimeArrays wrong_ownership;
    wrong_ownership.reference_haplotypes[0].push_back(
        htsim::methdb::pack_runtime_site(
            2U,
            0U,
            MethylationContext::cg_c,
            MethylationSource::beta,
            MethylationAllele::shared,
            true));
    require_error(
        [&] {
            (void)DiploidMethylationCatalog(
                contig.index,
                static_cast<std::uint32_t>(contig.length),
                wrong_ownership);
        },
        "shared allele in a haplotype runtime array was accepted");

    auto projection = htsim::haplotype::project_interval(
        contig,
        variants,
        0U,
        0U,
        static_cast<std::uint32_t>(contig.length));
    auto wrong_haplotype = projection;
    wrong_haplotype.haplotype = 1U;
    require_error(
        [&] {(void)catalog.sites_for_projection(wrong_haplotype);},
        "projection carrying another haplotype's event was accepted");
    const auto inserted = std::find(
        projection.reference_positions.begin(),
        projection.reference_positions.end(),
        -1);
    require(inserted != projection.reference_positions.end(),
            "invalid-input fixture has no inserted base");
    const auto inserted_offset = static_cast<std::size_t>(
        inserted - projection.reference_positions.begin());
    projection.template_bases[inserted_offset] = 3U;
    require_error(
        [&] {(void)catalog.sites_for_projection(projection);},
        "inserted bases inconsistent with their event were accepted");

    auto wrong_contig = htsim::haplotype::project_interval(
        contig,
        variants,
        0U,
        0U,
        static_cast<std::uint32_t>(contig.length));
    wrong_contig.contig_index = 1U;
    require_error(
        [&] {(void)catalog.sites_for_projection(wrong_contig);},
        "projection from another contig was accepted");

    ContextShapes invalid = shapes();
    invalid.cg.alpha = std::numeric_limits<double>::infinity();
    require_error(
        [&] {
            (void)DiploidMethylationCatalog(
                contig, variants, 73U, true, invalid);
        },
        "invalid Beta shapes were accepted");
    Contig wrong = contig;
    wrong.index = 1U;
    require_error(
        [&] {
            (void)DiploidMethylationCatalog(
                wrong, variants, 73U, true, shapes());
        },
        "variant catalog from another contig was accepted");
}

} // namespace

int main()
{
    try {
        test_reference_only_parity();
        test_variant_context_entities_and_alleles();
        test_fragment_boundary_uses_complete_haplotype_context();
        test_deletion_joins_previously_separated_context();
        test_cgmap_overlay_respects_haplotype_equivalence();
        test_cgmap_pool_uses_typed_variant_entities();
        test_asm_overlay_uses_typed_haplotype_mask();
        test_asm_links_fail_closed();
        test_invalid_inputs_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "diploid_methdb_catalog_test: "
                  << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
