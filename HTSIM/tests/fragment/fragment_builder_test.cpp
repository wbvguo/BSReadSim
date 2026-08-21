#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "fragment.h"
#include "variant.h"
#include "types.h"
#include "methdb.h"
#include "protocol.h"
#include "reference.h"

namespace {

using htsim::fragment_builder::FragmentBuilderError;
using htsim::fragment_builder::FragmentDetail;
using htsim::fragment_builder::ReadLayout;
using htsim::methdb::ContextShapes;
using htsim::methdb::MethylationCatalog;
using htsim::methdb::DiploidMethylationCatalog;
using htsim::model::Bases;
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
    } catch (const FragmentBuilderError &) {
        return;
    }
    throw std::runtime_error(message);
}

Bases encode(const std::string &text)
{
    Bases result;
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

Contig contig_for(const std::string &text)
{
    Contig contig;
    contig.index = 0;
    contig.name = "chrBuild";
    contig.bases = encode(text);
    contig.length = contig.bases.size();
    return contig;
}

void require_common_fragment_equal(
    const htsim::model::Fragment &left,
    const htsim::model::Fragment &right)
{
    require(
        left.fragment_ordinal == right.fragment_ordinal
            && left.contig_index == right.contig_index
            && left.haplotype == right.haplotype
            && left.capture_strand == right.capture_strand
            && left.reference_start == right.reference_start
            && left.reference_end == right.reference_end
            && left.template_bases == right.template_bases
            && left.methylation_sites.size()
                == right.methylation_sites.size()
            && left.mates.size() == right.mates.size(),
        "compact construction changed a fragment common column");
    for (std::size_t index = 0U;
         index < left.methylation_sites.size(); ++index) {
        const auto &left_site = left.methylation_sites[index];
        const auto &right_site = right.methylation_sites[index];
        require(
            left_site.site_index == right_site.site_index
                && left_site.template_offset == right_site.template_offset
                && left_site.reference_pos == right_site.reference_pos
                && left_site.context == right_site.context
                && left_site.methylation_source == right_site.methylation_source
                && left_site.allele == right_site.allele
                && left_site.methylation_probability
                    == right_site.methylation_probability,
            "compact construction changed a methylation row");
    }
    for (std::size_t index = 0U; index < left.mates.size(); ++index) {
        const auto &left_mate = left.mates[index];
        const auto &right_mate = right.mates[index];
        require(
            left_mate.mate_index == right_mate.mate_index
                && left_mate.reverse_complement
                    == right_mate.reverse_complement
                && left_mate.template_start == right_mate.template_start
                && left_mate.template_end == right_mate.template_end,
            "compact construction changed a mate common column");
    }
}

void require_compact_annotation_state_empty(
    const htsim::model::Fragment &fragment)
{
    require(
        fragment.reference_positions.empty()
            && fragment.base_variant_indices.empty()
            && fragment.variants.empty(),
        "compact construction retained projection details state");
    for (const auto &mate : fragment.mates) {
        require(
            mate.reference_start == 0U && mate.reference_end == 0U
                && mate.site_refs.empty(),
            "compact construction retained mate details state");
    }
}

ContextShapes shapes()
{
    return {{2.0, 5.0}, {3.0, 4.0}, {5.0, 2.0}};
}

void test_paired_fragment_and_overlap_projection()
{
    const Contig contig = contig_for("ACGCGTACGT");
    const MethylationCatalog catalog(
        contig.bases, contig.index, 17, true, shapes());
    const ReadLayout layout{6, 4, true};
    const auto fragment = htsim::fragment_builder::build_fragment(
        contig,
        catalog,
        0,
        1,
        1,
        htsim::model::CaptureStrand::reverse,
        layout);

    require(fragment.reference_start == 1 && fragment.reference_end == 7
                && fragment.capture_strand
                    == htsim::model::CaptureStrand::reverse
                && fragment.template_bases == encode("CGCGTA"),
            "fragment reference slice is wrong");
    require(fragment.reference_positions
                == std::vector<std::int64_t>({1, 2, 3, 4, 5, 6}),
            "reference positions are wrong");
    require(fragment.base_variant_indices.size() == 6
                && fragment.variants.empty(),
            "variant-free baseline emitted an event");
    for (const auto variant_index : fragment.base_variant_indices) {
        require(variant_index == htsim::model::no_variant_index,
                "variant-free base has an event id");
    }
    require(fragment.mates.size() == 2,
            "paired fragment did not contain two mates");
    require(fragment.mates[0].template_start == 0
                && fragment.mates[0].template_end == 4
                && !fragment.mates[0].reverse_complement
                && fragment.mates[0].reference_start == 1
                && fragment.mates[0].reference_end == 5,
            "R1 projection is wrong");
    require(fragment.mates[1].template_start == 2
                && fragment.mates[1].template_end == 6
                && fragment.mates[1].reverse_complement
                && fragment.mates[1].reference_start == 3
                && fragment.mates[1].reference_end == 7,
            "R2 projection is wrong");

    bool saw_shared_site = false;
    for (const auto &site : fragment.methylation_sites) {
        bool in_r1 = false;
        bool in_r2 = false;
        for (const auto &reference : fragment.mates[0].site_refs) {
            in_r1 = in_r1 || reference.site_index == site.site_index;
        }
        for (const auto &reference : fragment.mates[1].site_refs) {
            in_r2 = in_r2 || reference.site_index == site.site_index;
        }
        saw_shared_site = saw_shared_site || (in_r1 && in_r2);
    }
    require(saw_shared_site,
            "overlapping mates did not share a fragment-level site");
    for (const auto &mate : fragment.mates) {
        for (std::size_t index = 1; index < mate.site_refs.size(); ++index) {
            require(mate.site_refs[index - 1].read_offset
                        < mate.site_refs[index].read_offset,
                    "mate site references are not in oriented read order");
        }
    }

}

void test_single_end_and_payload_boundary()
{
    const Contig contig = contig_for("AACGTTAA");
    const MethylationCatalog catalog(
        contig.bases, contig.index, 19, true, shapes());
    const ReadLayout layout{5, 3, false};
    const auto fragment = htsim::fragment_builder::build_fragment(
        contig,
        catalog,
        7,
        2,
        0,
        htsim::model::CaptureStrand::unknown,
        layout);
    require(fragment.fragment_ordinal == 7 && fragment.mates.size() == 1
                && fragment.mates[0].template_start == 0
                && fragment.mates[0].template_end == 3,
            "single-end fragment is wrong");
    require(htsim::fragment_builder::maximum_payload_bytes(layout)
                == UINT64_C(297),
            "payload upper-bound formula changed");
}

void test_common_columns_fragment_detail()
{
    const Contig contig = contig_for("AACGCGTAACGTT");
    const MethylationCatalog reference_catalog(
        contig.bases, contig.index, 43U, true, shapes());
    const ReadLayout reference_layout{9U, 6U, true};
    const auto full_reference = htsim::fragment_builder::build_fragment(
        contig,
        reference_catalog,
        4U,
        1U,
        1U,
        htsim::model::CaptureStrand::reverse,
        reference_layout,
        FragmentDetail::full);
    const auto compact_reference =
        htsim::fragment_builder::build_fragment(
            contig,
            reference_catalog,
            4U,
            1U,
            1U,
            htsim::model::CaptureStrand::reverse,
            reference_layout,
            FragmentDetail::common_columns);
    require_common_fragment_equal(full_reference, compact_reference);
    require_compact_annotation_state_empty(compact_reference);

    const std::vector<Variant> events = {
        {0U, 3U, 4U, VariantKind::snv, encode("G"), encode("A"),
         htsim::model::HaplotypeMask::haplotype_1},
        {0U, 6U, 6U, VariantKind::insertion, {}, encode("CG"),
         htsim::model::HaplotypeMask::both},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const DiploidMethylationCatalog diploid_catalog(
        contig, variants, 47U, true, shapes());
    const ReadLayout projected_layout{10U, 6U, true};
    const auto full_projected = htsim::fragment_builder::build_fragment(
        htsim::haplotype::project_interval(
            contig, variants, 0U, 1U, 11U),
        diploid_catalog,
        5U,
        htsim::model::CaptureStrand::unknown,
        projected_layout,
        FragmentDetail::full);
    const auto compact_projected =
        htsim::fragment_builder::build_fragment(
            htsim::haplotype::project_interval(
                contig, variants, 0U, 1U, 11U),
            diploid_catalog,
            5U,
            htsim::model::CaptureStrand::unknown,
            projected_layout,
            FragmentDetail::common_columns);
    require_common_fragment_equal(full_projected, compact_projected);
    require_compact_annotation_state_empty(compact_projected);
    require(!full_projected.variants.empty(),
            "variant fixture did not exercise discarded event details");
}

void test_variant_projected_fragment_boundary()
{
    const Contig contig = contig_for("AACGTAACGTT");
    const std::vector<Variant> events = {
        {0U, 3U, 4U, VariantKind::snv, encode("G"), encode("A"),
         htsim::model::HaplotypeMask::haplotype_1},
        {0U, 6U, 6U, VariantKind::insertion, {}, encode("CG"),
         htsim::model::HaplotypeMask::both},
        {0U, 7U, 9U, VariantKind::deletion, encode("CG"), {},
         htsim::model::HaplotypeMask::haplotype_2},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const DiploidMethylationCatalog catalog(
        contig, variants, 29U, true, shapes());
    const ReadLayout layout{9U, 6U, true};
    auto projection = htsim::haplotype::project_interval(
        contig, variants, 0U, 1U, 10U);
    require(
        htsim::fragment_builder::maximum_payload_bytes(projection, layout)
            == UINT64_C(691),
        "variant payload bound omitted projected bases or event bytes");
    const auto fragment = htsim::fragment_builder::build_fragment(
        std::move(projection),
        catalog,
        0U,
        htsim::model::CaptureStrand::unknown,
        layout);

    require(fragment.contig_index == contig.index
                && fragment.haplotype == 0U
                && fragment.reference_start == 1U
                && fragment.reference_end == 10U
                && fragment.template_bases == encode("ACATACGACGT")
                && fragment.variants.size() == 2U,
            "variant projection was not transferred into the fragment");
    require(fragment.mates.size() == 2U
                && fragment.mates[0].template_start == 0U
                && fragment.mates[0].template_end == 6U
                && fragment.mates[0].reference_start == 1U
                && fragment.mates[0].reference_end == 6U
                && fragment.mates[1].template_start == 5U
                && fragment.mates[1].template_end == 11U
                && fragment.mates[1].reference_start == 6U
                && fragment.mates[1].reference_end == 10U,
            "variant-aware mate reference envelopes are wrong");

    const auto inserted_site = std::find_if(
        fragment.methylation_sites.begin(),
        fragment.methylation_sites.end(),
        [](const htsim::model::MethylationSite &site) {
            return site.reference_pos == -1
                && site.context
                    == htsim::model::MethylationContext::cg_c;
        });
    require(inserted_site != fragment.methylation_sites.end(),
            "inserted CpG was not attached to the projected fragment");
    const auto mate_references_site = [&](std::size_t mate_index) {
        return std::any_of(
            fragment.mates[mate_index].site_refs.begin(),
            fragment.mates[mate_index].site_refs.end(),
            [&](const htsim::model::SiteReference &reference) {
                return reference.site_index == inserted_site->site_index;
            });
    };
    require(mate_references_site(0U) && mate_references_site(1U),
            "overlapping mates did not share the inserted methylation site");

}

void test_haplotype_coordinate_insert_length()
{
    const Contig contig = contig_for("AACGTAACGTT");
    const std::vector<Variant> events = {
        {0U, 6U, 6U, VariantKind::insertion, {}, encode("CG"),
         htsim::model::HaplotypeMask::haplotype_1},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const DiploidMethylationCatalog catalog(
        contig, variants, 37U, true, shapes());
    const ReadLayout physical_layout{
        11U,
        6U,
        true,
        ReadLayout::InsertCoordinate::haplotype,
    };
    const auto fragment = htsim::fragment_builder::build_fragment(
        htsim::haplotype::project_interval(
            contig, variants, 0U, 1U, 10U),
        catalog,
        0U,
        htsim::model::CaptureStrand::unknown,
        physical_layout);
    require(fragment.reference_start == 1U
                && fragment.reference_end == 10U
                && fragment.template_bases.size() == 11U,
            "haplotype insert length was confused with reference width");

    require_error(
        [&]() {
            (void)htsim::fragment_builder::build_fragment(
                htsim::haplotype::project_interval(
                    contig, variants, 0U, 1U, 10U),
                catalog,
                0U,
                htsim::model::CaptureStrand::unknown,
                {9U,
                 6U,
                 true,
                 ReadLayout::InsertCoordinate::haplotype});
        },
        "haplotype-coordinate layout accepted the reference-envelope width");
}

void test_insertion_only_mate_envelope()
{
    const Contig contig = contig_for("A");
    const std::vector<Variant> events = {
        {0U, 0U, 1U, VariantKind::deletion, encode("A"), {},
         htsim::model::HaplotypeMask::both},
        {0U, 1U, 1U, VariantKind::insertion, {}, encode("C"),
         htsim::model::HaplotypeMask::both},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const DiploidMethylationCatalog catalog(
        contig, variants, 31U, true, shapes());
    const auto fragment = htsim::fragment_builder::build_fragment(
        htsim::haplotype::project_interval(
            contig, variants, 0U, 0U, 1U),
        catalog,
        0U,
        htsim::model::CaptureStrand::unknown,
        {1U, 1U, false});
    require(fragment.template_bases == encode("C")
                && fragment.reference_positions
                    == std::vector<std::int64_t>({-1})
                && fragment.mates[0].reference_start == 1U
                && fragment.mates[0].reference_end == 1U,
            "insertion-only mate did not use its zero-width anchor envelope");
}

void test_invalid_projected_inputs_fail_closed()
{
    const Contig contig = contig_for("AACGTAACGTT");
    const std::vector<Variant> events = {
        {0U, 3U, 4U, VariantKind::snv, encode("G"), encode("A"),
         htsim::model::HaplotypeMask::haplotype_1},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const DiploidMethylationCatalog catalog(
        contig, variants, 37U, true, shapes());
    const auto good = htsim::haplotype::project_interval(
        contig, variants, 0U, 1U, 8U);
    require_error(
        [&] {
            (void)htsim::fragment_builder::build_fragment(
                good,
                catalog,
                0U,
                htsim::model::CaptureStrand::unknown,
                {6U, 4U, false});
        },
        "projection/reference-span mismatch was accepted");

    auto wrong_haplotype = good;
    wrong_haplotype.haplotype = 1U;
    require_error(
        [&] {
            (void)htsim::fragment_builder::build_fragment(
                std::move(wrong_haplotype),
                catalog,
                0U,
                htsim::model::CaptureStrand::unknown,
                {7U, 4U, false});
        },
        "projection carrying another haplotype's event was accepted");

    auto missing_event_base = good;
    const auto changed = std::find_if(
        missing_event_base.base_variant_indices.begin(),
        missing_event_base.base_variant_indices.end(),
        [](std::uint32_t variant_index) {
            return variant_index != htsim::model::no_variant_index;
        });
    require(changed != missing_event_base.base_variant_indices.end(),
            "invalid projection fixture has no changed base");
    *changed = htsim::model::no_variant_index;
    require_error(
        [&] {
            (void)htsim::fragment_builder::maximum_payload_bytes(
                missing_event_base, {7U, 4U, false});
        },
        "incomplete projected event bases were accepted");

    auto duplicate_position = good;
    duplicate_position.reference_positions[1U] =
        duplicate_position.reference_positions[0U];
    require_error(
        [&] {
            (void)htsim::fragment_builder::maximum_payload_bytes(
                duplicate_position, {7U, 4U, false});
        },
        "non-increasing projected reference coordinates were accepted");

    const Contig short_contig = contig_for("ACGT");
    const ContigVariants shortening(
        short_contig.bases,
        {{0U, 0U, 2U, VariantKind::deletion, encode("AC"), {},
          htsim::model::HaplotypeMask::both}},
        short_contig.index);
    const DiploidMethylationCatalog short_catalog(
        short_contig, shortening, 41U, true, shapes());
    require_error(
        [&] {
            (void)htsim::fragment_builder::build_fragment(
                htsim::haplotype::project_interval(
                    short_contig, shortening, 0U, 0U, 4U),
                short_catalog,
                0U,
                htsim::model::CaptureStrand::unknown,
                {4U, 3U, false});
        },
        "variant-shortened template below read length was accepted");
}

void test_invalid_inputs_fail_closed()
{
    const Contig contig = contig_for("ACGTACGT");
    const MethylationCatalog catalog(
        contig.bases, contig.index, 23, true, shapes());
    require_error(
        [&] {
            (void)htsim::fragment_builder::build_fragment(
                contig,
                catalog,
                0,
                5,
                0,
                htsim::model::CaptureStrand::unknown,
                {4, 2, false});
        },
        "out-of-contig fragment was accepted");
    require_error(
        [&] {
            (void)htsim::fragment_builder::build_fragment(
                contig,
                catalog,
                0,
                0,
                2,
                htsim::model::CaptureStrand::unknown,
                {4, 2, false});
        },
        "invalid haplotype was accepted");
    require_error(
        [&] {
            (void)htsim::fragment_builder::build_fragment(
                contig,
                catalog,
                0,
                0,
                0,
                static_cast<htsim::model::CaptureStrand>(9),
                {4, 2, false});
        },
        "invalid capture strand was accepted");
    require_error(
        [] {
            htsim::fragment_builder::require_payload_fits_protocol(
                {UINT32_C(4294967295), 1, false});
        },
        "oversized protocol fragment was accepted");
    require_error(
        [&] {
            (void)htsim::fragment_builder::build_fragment(
                contig,
                catalog,
                0U,
                0U,
                0U,
                htsim::model::CaptureStrand::unknown,
                {4U, 2U, false},
                static_cast<FragmentDetail>(9U));
        },
        "invalid fragment detail was accepted");
}

} // namespace

int main()
{
    try {
        test_paired_fragment_and_overlap_projection();
        test_single_end_and_payload_boundary();
        test_common_columns_fragment_detail();
        test_variant_projected_fragment_boundary();
        test_haplotype_coordinate_insert_length();
        test_insertion_only_mate_envelope();
        test_invalid_projected_inputs_fail_closed();
        test_invalid_inputs_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "fragment_builder_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
