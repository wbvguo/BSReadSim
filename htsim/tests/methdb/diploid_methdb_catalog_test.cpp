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
using htsim::methdb::DiploidSite;
using htsim::methdb::MethylationCatalog;
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

const DiploidSite *find(
    const std::vector<DiploidSite> &sites,
    std::uint64_t origin_id)
{
    const auto found = std::find_if(
        sites.begin(), sites.end(),
        [origin_id](const DiploidSite &site) {
            return site.origin_id == origin_id;
        });
    return found == sites.end() ? nullptr : &*found;
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
    require(diploid.shared_sites().size() == baseline.sites().size()
                && diploid.haplotype_sites(0).empty()
                && diploid.haplotype_sites(1).empty(),
            "reference-only diploid catalog did not collapse to shared sites");
    for (std::size_t index = 0; index < baseline.sites().size(); ++index) {
        const auto &old_site = baseline.sites()[index];
        const auto &new_site = diploid.shared_sites()[index];
        require(new_site.origin_id
                    == htsim::methdb::reference_origin_id(
                        old_site.reference_position)
                    && new_site.context == old_site.context
                    && new_site.allele == MethylationAllele::shared
                    && new_site.methylation_probability
                        == old_site.methylation_probability,
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
                        == baseline.sites()[index].methylation_probability,
                "reference-only protocol site projection changed");
    }

    const std::vector<CgmapRecord> records = {
        {1U, 0.125F, MethylationContext::cg_c, true, 2U},
        {2U, 0.0F, MethylationContext::cg_g, false, 2U},
        {4U, 0.75F, MethylationContext::chg_c, true, 0U},
    };
    const MethylationCatalog cgmap_baseline(
        contig.bases, contig.index, 41U, true, configured, &records);
    const DiploidMethylationCatalog cgmap_diploid(
        contig, variants, 41U, true, configured, &records);
    require(
        cgmap_diploid.shared_sites().size() == cgmap_baseline.sites().size(),
        "event-free CGmap diploid catalog changed the site count");
    for (std::size_t index = 0U;
         index < cgmap_baseline.sites().size();
         ++index) {
        const auto &reference_site = cgmap_baseline.sites()[index];
        const auto &diploid_site = cgmap_diploid.shared_sites()[index];
        require(
            diploid_site.origin_id == reference_site.reference_position
                && diploid_site.context == reference_site.context
                && diploid_site.methylation_source == reference_site.methylation_source
                && diploid_site.methylation_probability
                    == reference_site.methylation_probability,
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

    const auto ref2 = htsim::methdb::reference_origin_id(2U);
    const DiploidSite *alt_context = find(catalog.haplotype_sites(0), ref2);
    const DiploidSite *ref_context = find(catalog.haplotype_sites(1), ref2);
    require(alt_context && ref_context
                && alt_context->context == MethylationContext::chh_c
                && alt_context->allele == MethylationAllele::alternate_haplotype
                && ref_context->context == MethylationContext::cg_c
                && ref_context->allele == MethylationAllele::reference_haplotype,
            "heterozygous SNV did not split reference and ALT contexts");
    require(alt_context->methylation_probability
                == htsim::beta_sampler::sample_beta_for_site(
                    seed,
                    contig.index,
                    htsim::methdb::variant_reference_site_entity(
                        2U, HaplotypeMask::haplotype_1, 0U),
                    configured.chh.alpha,
                    configured.chh.beta)
                && ref_context->methylation_probability
                    == htsim::beta_sampler::sample_beta_for_site(
                        seed,
                        contig.index,
                        htsim::methdb::reference_site_entity(2U),
                        configured.cg.alpha,
                        configured.cg.beta),
            "heterozygous context used the wrong 64-bit Beta identity");

    const DiploidSite *inserted_c = find(
        catalog.shared_sites(),
        htsim::methdb::insertion_origin_id(1U, 0U));
    const DiploidSite *inserted_g = find(
        catalog.shared_sites(),
        htsim::methdb::insertion_origin_id(1U, 1U));
    require(inserted_c && inserted_g
                && inserted_c->context == MethylationContext::cg_c
                && inserted_g->context == MethylationContext::cg_g
                && inserted_c->allele == MethylationAllele::shared
                && inserted_g->allele == MethylationAllele::shared,
            "homozygous inserted CpG was not shared by both haplotypes");
    require(inserted_c->methylation_probability
                == htsim::beta_sampler::sample_beta_for_site(
                    seed,
                    contig.index,
                    htsim::methdb::insertion_site_entity(
                        1U, 0U, HaplotypeMask::both, 0U),
                    configured.cg.alpha,
                    configured.cg.beta),
            "inserted CpG used the wrong event/offset Beta identity");

    for (const std::uint32_t position : {7U, 8U}) {
        const DiploidSite *site = find(
            catalog.haplotype_sites(0),
            htsim::methdb::reference_origin_id(position));
        require(site && site->allele == MethylationAllele::reference_haplotype,
                "heterozygous deletion did not retain reference-haplotype ownership");
        require(!find(
                    catalog.haplotype_sites(1),
                    htsim::methdb::reference_origin_id(position)),
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
    const DiploidSite *catalog_site = find(
        catalog.haplotype_sites(0),
        htsim::methdb::reference_origin_id(2U));
    require(catalog_site
                && alt_sites[0].methylation_probability
                    == catalog_site->methylation_probability,
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
    const DiploidSite *alt_c = find(
        catalog.haplotype_sites(0),
        htsim::methdb::reference_origin_id(0U));
    const DiploidSite *alt_g = find(
        catalog.haplotype_sites(0),
        htsim::methdb::reference_origin_id(3U));
    const DiploidSite *ref_c = find(
        catalog.haplotype_sites(1),
        htsim::methdb::reference_origin_id(0U));
    const DiploidSite *ref_g = find(
        catalog.haplotype_sites(1),
        htsim::methdb::reference_origin_id(3U));
    require(alt_c && alt_g && ref_c && ref_g
                && alt_c->context == MethylationContext::cg_c
                && alt_g->context == MethylationContext::cg_g
                && ref_c->context == MethylationContext::chh_c
                && ref_g->context == MethylationContext::chh_g,
            "deletion did not recompute context across its joined boundary");

    const DiploidMethylationCatalog cpg_only(
        contig, variants, 91U, false, shapes());
    require(cpg_only.haplotype_sites(0).size() == 2U
                && cpg_only.haplotype_sites(1).empty()
                && cpg_only.shared_sites().empty(),
            "CpG-only filter did not use the final haplotype context");
}

void test_cgmap_overlay_respects_haplotype_equivalence()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const ContigVariants variants(
        contig.bases, variant_fixture(), contig.index);
    const std::vector<CgmapRecord> records = {
        {2U, 0.875F, MethylationContext::cg_c, true, 2U},
        {3U, 0.0F, MethylationContext::cg_g, false, 2U},
        {7U, 0.625F, MethylationContext::cg_c, true, 2U},
    };
    const DiploidMethylationCatalog catalog(
        contig, variants, 73U, true, shapes(), &records);

    const auto ref2 = htsim::methdb::reference_origin_id(2U);
    const DiploidSite *alternate = find(catalog.haplotype_sites(0), ref2);
    const DiploidSite *reference = find(catalog.haplotype_sites(1), ref2);
    require(
        alternate != nullptr && reference != nullptr
            && alternate->methylation_source == MethylationSource::beta
            && alternate->context == MethylationContext::chh_c
            && reference->methylation_source == MethylationSource::cgmap
            && reference->context == MethylationContext::cg_c
            && reference->methylation_probability == 0.875F,
        "CGmap was not restricted to the reference-equivalent haplotype");

    const DiploidSite *retained_after_deletion = find(
        catalog.haplotype_sites(0),
        htsim::methdb::reference_origin_id(7U));
    require(
        retained_after_deletion != nullptr
            && retained_after_deletion->methylation_source == MethylationSource::cgmap
            && retained_after_deletion->methylation_probability == 0.625F,
        "CGmap did not overlay the retained reference haplotype");

    const DiploidSite *inserted = find(
        catalog.shared_sites(),
        htsim::methdb::insertion_origin_id(1U, 0U));
    require(
        inserted != nullptr && inserted->methylation_source == MethylationSource::beta,
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
            && projected_reference->methylation_probability == 0.875F,
        "protocol projection lost CGmap provenance");
}

void test_cgmap_pool_uses_typed_variant_entities()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const ContigVariants variants(
        contig.bases, variant_fixture(), contig.index);
    const std::vector<CgmapRecord> records = {
        {2U, 0.125F, MethylationContext::cg_c, true, 2U},
        {3U, 0.875F, MethylationContext::cg_g, true, 2U},
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

    const std::uint64_t ref2 =
        htsim::methdb::reference_origin_id(2U);
    const DiploidSite *alternate = find(catalog.haplotype_sites(0U), ref2);
    const DiploidSite *reference = find(catalog.haplotype_sites(1U), ref2);
    require(
        alternate != nullptr && reference != nullptr
            && alternate->context == MethylationContext::chh_c
            && alternate->methylation_source == MethylationSource::beta
            && reference->context == MethylationContext::cg_c
            && reference->methylation_source == MethylationSource::pooled_cgmap,
        "CGmap pool ignored the typed context on a split haplotype site");
    const auto expected_reference = pool.sample(
        reference->context,
        seed,
        contig.index,
        htsim::methdb::reference_site_entity(2U));
    require(
        expected_reference.has_value()
            && reference->methylation_probability == *expected_reference,
        "reference-equivalent haplotype used the wrong pool address");

    for (const std::uint8_t offset : {std::uint8_t{0U}, std::uint8_t{1U}}) {
        const DiploidSite *inserted = find(
            catalog.shared_sites(),
            htsim::methdb::insertion_origin_id(1U, offset));
        const auto entity = htsim::methdb::insertion_site_entity(
            1U, offset, HaplotypeMask::both, 0U);
        require(inserted != nullptr
                    && inserted->methylation_source == MethylationSource::pooled_cgmap,
                "variant-created CpG did not use the typed context pool");
        const auto expected = pool.sample(
            inserted->context, seed, contig.index, entity);
        require(expected.has_value()
                    && inserted->methylation_probability == *expected,
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
        {2U, 0.5F, MethylationContext::cg_c, true, 2U},
    };
    const std::vector<AsmRecord> asm_records = {
        {2U, 9U, 0.2F, 0.8F, MethylationContext::cg_c, 2U, 3U, 0U},
    };
    const DiploidMethylationCatalog catalog(
        contig, variants, 73U, true, shapes(), &cgmap, &asm_records);
    const std::uint64_t origin =
        htsim::methdb::reference_origin_id(2U);
    require(
        find(catalog.shared_sites(), origin) == nullptr,
        "ASM target remained in the shared catalog");
    const DiploidSite *alternate = find(catalog.haplotype_sites(0U), origin);
    const DiploidSite *reference = find(catalog.haplotype_sites(1U), origin);
    require(
        alternate != nullptr && reference != nullptr
            && alternate->methylation_source == MethylationSource::asm_source
            && alternate->allele
                == MethylationAllele::alternate_haplotype
            && alternate->methylation_probability == 0.8F
            && reference->methylation_source == MethylationSource::asm_source
            && reference->allele
                == MethylationAllele::reference_haplotype
            && reference->methylation_probability == 0.2F,
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
                    == (haplotype == 0U ? 0.8F : 0.2F),
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
    const DiploidSite *reversed_reference =
        find(reversed.haplotype_sites(0U), origin);
    const DiploidSite *reversed_alternate =
        find(reversed.haplotype_sites(1U), origin);
    require(
        reversed_reference != nullptr && reversed_alternate != nullptr
            && reversed_reference->allele
                == MethylationAllele::reference_haplotype
            && reversed_reference->methylation_probability == 0.2F
            && reversed_alternate->allele
                == MethylationAllele::alternate_haplotype
            && reversed_alternate->methylation_probability == 0.8F,
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
    require(
        find(pooled.shared_sites(), origin) == nullptr,
        "ASM target remained shared after CGmap pooling");
    const DiploidSite *pooled_alternate =
        find(pooled.haplotype_sites(0U), origin);
    const DiploidSite *pooled_reference =
        find(pooled.haplotype_sites(1U), origin);
    require(
        pooled_alternate != nullptr && pooled_reference != nullptr
            && pooled_alternate->methylation_source == MethylationSource::asm_source
            && pooled_alternate->methylation_probability == 0.8F
            && pooled_reference->methylation_source == MethylationSource::asm_source
            && pooled_reference->methylation_probability == 0.2F,
        "ASM did not retain precedence over the typed CGmap pool");
}

void test_asm_links_fail_closed()
{
    const Contig contig = make_contig("AACGTAACGTT");
    const std::vector<AsmRecord> records = {
        {2U, 9U, 0.2F, 0.8F, MethylationContext::cg_c, 2U, 3U, 0U},
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
        {2U, 3U, 0.2F, 0.8F, MethylationContext::cg_c, 2U, 2U, 0U},
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
    require_error(
        [&] {(void)catalog.haplotype_sites(2U);},
        "invalid haplotype catalog was accepted");
    require_error(
        [] {
            (void)htsim::methdb::insertion_origin_id(
                htsim::model::no_variant_index, 0U);
        },
        "no-event insertion origin was accepted");
    require_error(
        [] {
            (void)htsim::methdb::insertion_origin_id(1U, 4U);
        },
        "fifth insertion origin base was accepted");

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
