#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "variant.h"
#include "types.h"
#include "protocol.h"
#include "reference.h"

namespace {

using htsim::haplotype::ProjectedInterval;
using htsim::haplotype::ProjectionError;
using htsim::haplotype::ProjectionFailure;
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

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const ProjectionError &) {
        return;
    }
    throw std::runtime_error(message);
}

template <typename Operation>
void require_failure(
    Operation operation,
    ProjectionFailure expected,
    const std::string &message)
{
    try {
        operation();
    } catch (const ProjectionError &error) {
        require(error.failure() == expected, message + ": wrong failure kind");
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
    contig.name = "chrProjection";
    contig.bases = encode(sequence);
    contig.length = contig.bases.size();
    return contig;
}

std::vector<Event> mixed_events()
{
    return {
        {0, 1, 2, VariantKind::snv, encode("C"), encode("T"),
         HaplotypeMask::haplotype_1},
        {0, 4, 4, VariantKind::insertion, {}, encode("GG"),
         HaplotypeMask::both},
        {0, 6, 8, VariantKind::deletion, encode("GT"), {},
         HaplotypeMask::haplotype_1},
        {0, 9, 10, VariantKind::snv, encode("C"), encode("A"),
         HaplotypeMask::haplotype_2},
    };
}

void test_haplotype_projection_and_provenance()
{
    const Contig contig = make_contig("ACGTACGTACGT");
    const ContigVariants variants(contig.bases, mixed_events(), contig.index);

    const ProjectedInterval haplotype_0 =
        htsim::haplotype::project_interval(contig, variants, 0, 1, 10);
    require(haplotype_0.contig_index == contig.index
                && haplotype_0.haplotype == 0U
                && haplotype_0.template_bases == encode("TGTGGACAC"),
            "haplotype-0 sequence projection changed");
    require(haplotype_0.reference_positions
                == std::vector<std::int64_t>({1, 2, 3, -1, -1, 4, 5, 8, 9}),
            "haplotype-0 reference positions changed");
    require(haplotype_0.base_event_ids
                == std::vector<std::uint32_t>({
                    0,
                    htsim::model::no_variant_event,
                    htsim::model::no_variant_event,
                    1,
                    1,
                    htsim::model::no_variant_event,
                    htsim::model::no_variant_event,
                    htsim::model::no_variant_event,
                    htsim::model::no_variant_event,
                }),
            "haplotype-0 base provenance changed");
    require(haplotype_0.variant_events.size() == 3U
                && haplotype_0.variant_events[0].event_id == 0U
                && haplotype_0.variant_events[0].phased_haplotype == 0U
                && haplotype_0.variant_events[1].event_id == 1U
                && haplotype_0.variant_events[1].phased_haplotype == 255U
                && haplotype_0.variant_events[2].event_id == 2U
                && haplotype_0.variant_events[2].kind == VariantKind::deletion,
            "haplotype-0 event provenance changed");

    const ProjectedInterval haplotype_1 =
        htsim::haplotype::project_interval(contig, variants, 1, 1, 10);
    require(haplotype_1.contig_index == contig.index
                && haplotype_1.haplotype == 1U
                && haplotype_1.template_bases == encode("CGTGGACGTAA"),
            "haplotype-1 sequence projection changed");
    require(haplotype_1.variant_events.size() == 2U
                && haplotype_1.variant_events[0].event_id == 1U
                && haplotype_1.variant_events[1].event_id == 3U
                && haplotype_1.variant_events[1].phased_haplotype == 1U,
            "sparse per-contig event ids or phase changed");
}

void test_insertion_boundary_ownership()
{
    const Contig contig = make_contig("ACGT");
    const std::vector<Event> events = {
        {0, 0, 0, VariantKind::insertion, {}, encode("T"),
         HaplotypeMask::both},
        {0, 2, 2, VariantKind::insertion, {}, encode("G"),
         HaplotypeMask::both},
        {0, 4, 4, VariantKind::insertion, {}, encode("A"),
         HaplotypeMask::both},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const auto left =
        htsim::haplotype::project_interval(contig, variants, 0, 0, 2);
    const auto right =
        htsim::haplotype::project_interval(contig, variants, 0, 2, 4);
    require(left.template_bases == encode("TAC")
                && left.variant_events.size() == 1U
                && left.variant_events[0].event_id == 0U,
            "start or internal-end insertion ownership changed");
    require(right.template_bases == encode("GGTA")
                && right.reference_positions
                    == std::vector<std::int64_t>({-1, 2, 3, -1})
                && right.variant_events.size() == 2U
                && right.variant_events[0].event_id == 1U
                && right.variant_events[1].event_id == 2U,
            "internal-start or terminal insertion ownership changed");

}

void test_inactive_and_partial_events()
{
    const Contig contig = make_contig("ACGTAC");
    const std::vector<Event> events = {
        {0, 1, 4, VariantKind::deletion, encode("CGT"), {},
         HaplotypeMask::haplotype_1},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const auto inactive =
        htsim::haplotype::project_interval(contig, variants, 1, 2, 5);
    require(inactive.template_bases == encode("GTA")
                && inactive.variant_events.empty(),
            "inactive haplotype event changed the reference sequence");
    require_failure(
        [&] {
            (void)htsim::haplotype::project_interval(
                contig, variants, 0, 2, 5);
        },
        ProjectionFailure::boundary_cut,
        "interval cutting the left side of a deletion was accepted");
    require_failure(
        [&] {
            (void)htsim::haplotype::project_interval(
                contig, variants, 0, 0, 3);
        },
        ProjectionFailure::boundary_cut,
        "interval cutting the right side of a deletion was accepted");
    require_failure(
        [&] {
            (void)htsim::haplotype::project_interval(
                contig, variants, 0, 1, 4);
        },
        ProjectionFailure::empty_projection,
        "fully deleted interval was accepted as a protocol template");
}

void test_invalid_inputs_fail_closed()
{
    const Contig contig = make_contig("ACGT");
    const ContigVariants variants(contig.bases, {}, contig.index);
    require_failure(
        [&] {
            (void)htsim::haplotype::project_interval(
                contig, variants, 2, 0, 2);
        },
        ProjectionFailure::invalid_input,
        "invalid haplotype was accepted");
    require_error(
        [&] {
            (void)htsim::haplotype::project_interval(
                contig, variants, 0, 2, 2);
        },
        "empty interval was accepted");
    require_error(
        [&] {
            (void)htsim::haplotype::project_interval(
                contig, variants, 0, 0, 5);
        },
        "out-of-contig interval was accepted");
    Contig inconsistent = contig;
    inconsistent.length = 3;
    require_error(
        [&] {
            (void)htsim::haplotype::project_interval(
                inconsistent, variants, 0, 0, 2);
        },
        "inconsistent materialized contig was accepted");
    Contig different_index = contig;
    different_index.index = 1;
    const ContigVariants nonempty(
        contig.bases,
        {{0, 1, 2, VariantKind::snv, encode("C"), encode("T"),
          HaplotypeMask::haplotype_2}},
        contig.index);
    require_error(
        [&] {
            (void)htsim::haplotype::project_interval(
                different_index, nonempty, 0, 0, 2);
        },
        "variant catalog from another materialized contig was accepted");
}

} // namespace

int main()
{
    try {
        test_haplotype_projection_and_provenance();
        test_insertion_boundary_ownership();
        test_inactive_and_partial_events();
        test_invalid_inputs_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "haplotype_projector_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
