#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "variant.h"
#include "types.h"
#include "protocol.h"
#include "reference.h"
#include "wgbs.h"

namespace {

using htsim::model::HaplotypeMask;
using htsim::model::Bases;
using htsim::model::VariantKind;
using htsim::reference::Contig;
using htsim::wgbs::FixedFragmentShape;
using htsim::wgbs::HaplotypeStartIndex;
using htsim::wgbs::ValidStartIndex;
using htsim::wgbs::VariantStartIndexError;
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
    } catch (const VariantStartIndexError &) {
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
        default: throw std::runtime_error("invalid fixture base");
        }
    }
    return bases;
}

Contig make_contig(const std::string &sequence)
{
    Contig contig;
    contig.index = 0U;
    contig.name = "chrStarts";
    contig.bases = encode(sequence);
    contig.length = contig.bases.size();
    return contig;
}

void test_reference_only_exact_parity()
{
    const Contig contig = make_contig("ACGTACGTNNACGTACGT");
    const ContigVariants variants(contig.bases, {}, contig.index);
    const FixedFragmentShape shape{5U, 3U, true, 0.0};
    const ValidStartIndex reference(contig.bases, shape);
    const HaplotypeStartIndex diploid(contig, variants, shape);
    require(diploid.possible_start_count() == reference.possible_start_count()
                && diploid.valid_start_count() == reference.valid_start_count(),
            "event-free diploid start domain changed reference eligibility");
    for (std::uint32_t start = 0U;
         start < diploid.possible_start_count();
         ++start) {
        require(diploid.is_valid_start(start) == reference.is_valid_start(start),
                "event-free start validity changed");
        if (diploid.is_valid_start(start)) {
            require(diploid.haplotype_mask(start) == HaplotypeMask::both,
                    "event-free start did not retain both haplotypes");
        }
    }
    for (std::uint32_t rank = 0U; rank < diploid.valid_start_count(); ++rank) {
        require(diploid.start_for_rank(rank) == reference.start_for_rank(rank),
                "event-free rank-select mapping changed");
    }
    require(diploid.sample(contig.index, 17U, 11U, 25U)
                == reference.sample(contig.index, 17U, 11U, 25U),
            "event-free fragment RNG addresses changed");
}

void test_insertion_can_rescue_reference_n_windows()
{
    const Contig contig = make_contig("AANAA");
    const std::vector<Variant> events = {
        {0U, 2U, 2U, VariantKind::insertion, {}, encode("CCC"),
         HaplotypeMask::haplotype_1},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const HaplotypeStartIndex index(
        contig, variants, {4U, 3U, false, 0.0});
    require(index.possible_start_count() == 2U
                && index.valid_start_count() == 2U,
            "insertion-rescued starts were omitted");
    for (std::uint32_t start = 0U; start < 2U; ++start) {
        require(index.haplotype_mask(start) == HaplotypeMask::haplotype_1,
                "insertion rescue was assigned to the wrong haplotype bit");
    }
}

void test_deletion_boundaries_length_and_n_shift()
{
    const Contig contig = make_contig("AAAANAAA");
    const std::vector<Variant> events = {
        {0U, 1U, 3U, VariantKind::deletion, encode("AA"), {},
         HaplotypeMask::haplotype_1},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const HaplotypeStartIndex index(
        contig, variants, {5U, 3U, false, 0.0});
    require(index.possible_start_count() == 4U
                && index.valid_start_count() == 2U,
            "deletion eligibility count changed");
    require(index.haplotype_mask(0U) == HaplotypeMask::haplotype_2
                && index.haplotype_mask(1U) == HaplotypeMask::haplotype_2,
            "deletion-shifted N entered the wrong haplotype domain");
    require(!index.is_valid_start(2U) && !index.is_valid_start(3U),
            "partial deletion boundary or reference N window was accepted");
    require_error(
        [&] {(void)index.haplotype_mask(2U);},
        "zero eligibility bits were converted to a HaplotypeMask");
}

void test_terminal_insertion_after_complete_deletion()
{
    const Contig contig = make_contig("A");
    const std::vector<Variant> events = {
        {0U, 0U, 1U, VariantKind::deletion, encode("A"), {},
         HaplotypeMask::both},
        {0U, 1U, 1U, VariantKind::insertion, {}, encode("C"),
         HaplotypeMask::both},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    const HaplotypeStartIndex index(
        contig, variants, {1U, 1U, true, 0.0});
    require(index.valid_start_count() == 1U
                && index.haplotype_mask(0U) == HaplotypeMask::both,
            "terminal insertion did not rescue a fully deleted interval");
}

std::uint8_t projected_expected_mask(
    const Contig &contig,
    const ContigVariants &variants,
    std::uint32_t start,
    const FixedFragmentShape &shape)
{
    const auto maximum_n = static_cast<std::size_t>(std::floor(
        shape.max_ambiguous_fraction
        * static_cast<double>(shape.read_length)));
    std::uint8_t mask = 0U;
    for (std::uint8_t haplotype = 0U; haplotype < 2U; ++haplotype) {
        try {
            const auto projection = htsim::haplotype::project_interval(
                contig,
                variants,
                haplotype,
                start,
                start + shape.insert_length);
            if (projection.template_bases.size() < shape.read_length) {continue;}
            const auto first_n = static_cast<std::size_t>(std::count(
                projection.template_bases.begin(),
                projection.template_bases.begin()
                    + static_cast<std::ptrdiff_t>(shape.read_length),
                static_cast<std::uint8_t>(4U)));
            const auto last_n = static_cast<std::size_t>(std::count(
                projection.template_bases.end()
                    - static_cast<std::ptrdiff_t>(shape.read_length),
                projection.template_bases.end(),
                static_cast<std::uint8_t>(4U)));
            if (first_n <= maximum_n
                && (!shape.paired_end || last_n <= maximum_n)) {
                mask |= static_cast<std::uint8_t>(1U << haplotype);
            }
        } catch (const htsim::haplotype::ProjectionError &) {
        }
    }
    return mask;
}

void test_exhaustive_small_projection_parity()
{
    const Contig contig = make_contig("ANCGTNAACGTN");
    const std::vector<Variant> events = {
        {0U, 0U, 0U, VariantKind::insertion, {}, encode("T"),
         HaplotypeMask::haplotype_1},
        {0U, 2U, 3U, VariantKind::snv, encode("C"), encode("A"),
         HaplotypeMask::haplotype_2},
        {0U, 3U, 5U, VariantKind::deletion, encode("GT"), {},
         HaplotypeMask::haplotype_1},
        {0U, 6U, 6U, VariantKind::insertion, {}, encode("CG"),
         HaplotypeMask::both},
        {0U, 8U, 10U, VariantKind::deletion, encode("CG"), {},
         HaplotypeMask::haplotype_2},
        {0U, 12U, 12U, VariantKind::insertion, {}, encode("A"),
         HaplotypeMask::haplotype_1},
    };
    const ContigVariants variants(contig.bases, events, contig.index);
    for (const bool paired_end : {false, true}) {
        for (std::uint32_t insert = 1U; insert <= 6U; ++insert) {
            for (std::uint32_t read = 1U; read <= insert; ++read) {
                for (const double ambiguous : {0.0, 0.5, 1.0}) {
                    const FixedFragmentShape shape{
                        insert, read, paired_end, ambiguous};
                    const HaplotypeStartIndex index(contig, variants, shape);
                    std::vector<std::uint32_t> expected_starts;
                    for (std::uint32_t start = 0U;
                         start < index.possible_start_count();
                         ++start) {
                        const std::uint8_t expected = projected_expected_mask(
                            contig, variants, start, shape);
                        require(index.is_valid_start(start) == (expected != 0U),
                                "index validity disagreed with projector");
                        if (expected != 0U) {
                            require(
                                static_cast<std::uint8_t>(
                                    index.haplotype_mask(start)) == expected,
                                "two-bit availability disagreed with projector");
                            expected_starts.push_back(start);
                        }
                    }
                    require(index.valid_start_count() == expected_starts.size(),
                            "rank domain size disagreed with projector");
                    for (std::uint32_t rank = 0U;
                         rank < index.valid_start_count();
                         ++rank) {
                        require(index.start_for_rank(rank)
                                    == expected_starts[rank],
                                "rank order disagreed with reference starts");
                    }
                }
            }
        }
    }
}

void test_invalid_inputs_and_rank_bounds()
{
    const Contig contig = make_contig("ACGTACGT");
    const ContigVariants variants(contig.bases, {}, contig.index);
    const HaplotypeStartIndex index(
        contig, variants, {4U, 2U, false, 0.0});
    require_error(
        [&] {(void)index.start_for_rank(index.valid_start_count());},
        "out-of-range valid-start rank was accepted");
    require_error(
        [&] {
            (void)index.sample(
                contig.index,
                1U,
                std::numeric_limits<std::uint64_t>::max(),
                2U);
        },
        "overflowing candidate ordinal range was accepted");
    require_error(
        [&] {
            (void)HaplotypeStartIndex(
                contig, variants, {4U, 5U, false, 0.0});
        },
        "read longer than insert was accepted");
    Contig wrong = contig;
    wrong.index = 1U;
    require_error(
        [&] {
            (void)HaplotypeStartIndex(
                wrong, variants, {4U, 2U, false, 0.0});
        },
        "variant catalog from another contig was accepted");
}

} // namespace

int main()
{
    try {
        test_reference_only_exact_parity();
        test_insertion_can_rescue_reference_n_windows();
        test_deletion_boundaries_length_and_n_shift();
        test_terminal_insertion_after_complete_deletion();
        test_exhaustive_small_projection_parity();
        test_invalid_inputs_and_rank_bounds();
    } catch (const std::exception &error) {
        std::cerr << "variant_start_index_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
