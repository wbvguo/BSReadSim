#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "variant.h"
#include "types.h"
#include "protocol.h"
#include "reference.h"

namespace {

using htsim::haplotype::HaplotypeLayout;
using htsim::haplotype::HaplotypeLayoutError;
using htsim::haplotype::ProjectionError;
using htsim::model::HaplotypeMask;
using htsim::model::Bases;
using htsim::model::VariantKind;
using htsim::reference::Contig;
using htsim::variant::ContigVariants;
using htsim::variant::Event;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Error, typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const Error &) {
        return;
    }
    throw std::runtime_error(message);
}

Bases encode(const std::string &text)
{
    Bases result;
    for (const char base : text) {
        switch (base) {
        case 'A': result.push_back(0U); break;
        case 'C': result.push_back(1U); break;
        case 'G': result.push_back(2U); break;
        case 'T': result.push_back(3U); break;
        case 'N': result.push_back(4U); break;
        default: throw std::runtime_error("invalid test base");
        }
    }
    return result;
}

Contig make_contig()
{
    Contig contig;
    contig.index = 0U;
    contig.name = "chrLayout";
    contig.bases = encode("ACNGTACG");
    contig.length = contig.bases.size();
    return contig;
}

std::vector<Event> events()
{
    return {
        {0U, 1U, 1U, VariantKind::insertion, {}, encode("TT"),
         HaplotypeMask::haplotype_1},
        {0U, 3U, 5U, VariantKind::deletion, encode("GT"), {},
         HaplotypeMask::haplotype_1},
        {0U, 6U, 7U, VariantKind::snv, encode("C"), encode("G"),
         HaplotypeMask::haplotype_2},
    };
}

void test_materialization_and_two_bit_selection()
{
    const Contig contig = make_contig();
    const ContigVariants variants(contig.bases, events(), contig.index);
    const HaplotypeLayout haplotype_0(contig, variants, 0U, true);
    const HaplotypeLayout haplotype_1(contig, variants, 1U, true);

    require(haplotype_0.length() == 8U
                && haplotype_0.bases() == encode("ATTCNACG"),
            "haplotype-0 materialization ignored its mask bits");
    require(haplotype_1.length() == 8U
                && haplotype_1.bases() == encode("ACNGTAGG"),
            "haplotype-1 materialization ignored its mask bits");
    require(haplotype_0.ambiguous_count(0U, 5U) == 1U
                && haplotype_0.ambiguous_count(5U, 8U) == 0U,
            "haplotype ambiguity rank changed");

    const HaplotypeLayout compact(contig, variants, 0U, false);
    require(!compact.has_materialized_bases() && compact.length() == 8U,
            "compact layout retained the wrong storage contract");
    require_error<HaplotypeLayoutError>(
        [&] {(void)compact.bases();},
        "compact layout exposed absent bases");
}

void test_reference_and_haplotype_boundaries()
{
    const Contig contig = make_contig();
    const ContigVariants variants(contig.bases, events(), contig.index);
    const HaplotypeLayout layout(contig, variants, 0U, false);

    const auto before_insertion = layout.boundary(1U);
    const auto after_insertion = layout.boundary(3U);
    const auto across_deletion = layout.boundary(5U);
    require(before_insertion
                && before_insertion->left_reference_end == 1U
                && before_insertion->right_reference_start == 1U
                && !before_insertion->include_insertion_in_left_fragment
                && before_insertion->include_insertion_in_right_fragment,
            "boundary before insertion changed ownership");
    require(after_insertion
                && after_insertion->left_reference_end == 1U
                && after_insertion->right_reference_start == 1U
                && after_insertion->include_insertion_in_left_fragment
                && !after_insertion->include_insertion_in_right_fragment,
            "boundary after insertion changed ownership");
    require(!layout.boundary(2U),
            "boundary inside an insertion became representable");
    require(across_deletion
                && across_deletion->left_reference_end == 3U
                && across_deletion->right_reference_start == 5U,
            "collapsed deletion boundary lost its two reference sides");

    require(layout.boundary_before_reference(1U) == 1U
                && layout.boundary_before_reference(2U) == 4U
                && !layout.boundary_before_reference(4U)
                && layout.boundary_before_reference(5U) == 5U,
            "reference-to-haplotype boundary projection changed");
}

void test_physical_slices_project_to_typed_truth()
{
    const Contig contig = make_contig();
    const ContigVariants variants(contig.bases, events(), contig.index);
    const HaplotypeLayout layout(contig, variants, 0U, false);

    const auto insertion = layout.project(contig, variants, 1U, 3U);
    require(insertion.reference_start == 1U
                && insertion.reference_end == 1U
                && insertion.template_bases == encode("TT")
                && insertion.reference_positions
                    == std::vector<std::int64_t>({-1, -1})
                && insertion.variant_events.size() == 1U,
            "insertion-only physical slice lost typed truth");

    const auto after_insertion = layout.project(contig, variants, 3U, 5U);
    require(after_insertion.reference_start == 1U
                && after_insertion.reference_end == 3U
                && after_insertion.template_bases == encode("CN")
                && after_insertion.variant_events.empty(),
            "slice after insertion incorrectly retained the insertion");

    const auto after_deletion = layout.project(contig, variants, 5U, 7U);
    require(after_deletion.reference_start == 5U
                && after_deletion.reference_end == 7U
                && after_deletion.template_bases == encode("AC")
                && after_deletion.variant_events.empty(),
            "slice after deletion used the wrong reference envelope");

    const auto spanning = layout.project(contig, variants, 0U, 7U);
    require(spanning.reference_start == 0U
                && spanning.reference_end == 7U
                && spanning.template_bases == encode("ATTCNAC")
                && spanning.variant_events.size() == 2U
                && spanning.variant_events[0].kind == VariantKind::insertion
                && spanning.variant_events[1].kind == VariantKind::deletion,
            "slice spanning insertion and deletion lost its events");

    require_error<ProjectionError>(
        [&] {(void)layout.project(contig, variants, 2U, 5U);},
        "slice beginning inside insertion was accepted");
}

void test_physical_payload_bound_is_window_local()
{
    const Contig contig = make_contig();
    const ContigVariants variants(contig.bases, events(), contig.index);
    const HaplotypeLayout haplotype_0(contig, variants, 0U, false);
    const HaplotypeLayout haplotype_1(contig, variants, 1U, false);

    require(haplotype_0.maximum_variant_event_payload_bytes(variants, 4U)
                == 34U,
            "payload bound summed variants that cannot share one window");
    require(haplotype_0.maximum_variant_event_payload_bytes(variants, 5U)
                == 68U,
            "payload bound missed insertion and deletion in one window");
    require(haplotype_1.maximum_variant_event_payload_bytes(variants, 4U)
                == 34U,
            "payload bound ignored the selected haplotype bit");
    require(haplotype_0.maximum_variant_event_payload_bytes(variants, 9U)
                == 0U,
            "payload bound invented a window longer than its haplotype");
    require_error<HaplotypeLayoutError>(
        [&] {
            (void)haplotype_0.maximum_variant_event_payload_bytes(
                variants, 0U);
        },
        "zero-length payload window was accepted");
}

} // namespace

int main()
{
    try {
        test_materialization_and_two_bit_selection();
        test_reference_and_haplotype_boundaries();
        test_physical_slices_project_to_typed_truth();
        test_physical_payload_bound_is_window_local();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
