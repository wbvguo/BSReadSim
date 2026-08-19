#include "wgbs.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

#include "variant.h"
#include "types.h"
#include "protocol.h"
#include "reference.h"

namespace {

using htsim::model::HaplotypeMask;
using htsim::model::Bases;
using htsim::model::VariantKind;
using htsim::reference::Contig;
using htsim::wgbs::VariableHaplotypeCandidate;
using htsim::wgbs::VariableHaplotypeOptions;
using htsim::wgbs::VariableHaplotypeSampler;
using htsim::wgbs::VariableHaplotypeSamplingError;
using htsim::wgbs::VariableWgbsSampler;
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
    } catch (const VariableHaplotypeSamplingError &) {
        return;
    }
    throw std::runtime_error(message);
}

Bases encode(const std::string &text)
{
    Bases bases;
    bases.reserve(text.size());
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

Contig make_contig(const std::string &sequence, std::string name = "chrVariable")
{
    Contig contig;
    contig.index = 0U;
    contig.name = std::move(name);
    contig.bases = encode(sequence);
    contig.length = contig.bases.size();
    return contig;
}

std::optional<HaplotypeMask> projected_mask(
    const Contig &contig,
    const ContigVariants &variants,
    std::uint32_t start,
    std::uint32_t span,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction)
{
    const std::size_t maximum_n = static_cast<std::size_t>(std::floor(
        max_ambiguous_fraction * static_cast<double>(read_length)));
    std::uint8_t mask = 0U;
    for (std::uint8_t haplotype = 0U; haplotype < 2U; ++haplotype) {
        try {
            const auto projection = htsim::haplotype::project_interval(
                contig, variants, haplotype, start, start + span);
            if (projection.template_bases.size() < read_length) {continue;}
            const std::size_t first_n = static_cast<std::size_t>(std::count(
                projection.template_bases.begin(),
                projection.template_bases.begin()
                    + static_cast<std::ptrdiff_t>(read_length),
                static_cast<std::uint8_t>(4U)));
            const std::size_t last_n = static_cast<std::size_t>(std::count(
                projection.template_bases.end()
                    - static_cast<std::ptrdiff_t>(read_length),
                projection.template_bases.end(),
                static_cast<std::uint8_t>(4U)));
            if (first_n <= maximum_n
                && (!paired_end || last_n <= maximum_n)) {
                mask |= static_cast<std::uint8_t>(1U << haplotype);
            }
        } catch (const htsim::haplotype::ProjectionError &) {
        }
    }
    if (mask == 0U) {return std::nullopt;}
    return static_cast<HaplotypeMask>(mask);
}

void test_reference_only_matches_variable_wgbs_exactly()
{
    std::string sequence;
    sequence.reserve(80U);
    constexpr char bases[] = "ACGT";
    for (std::size_t index = 0U; index < 80U; ++index) {
        sequence.push_back(index % 7U == 0U ? 'N' : bases[index % 4U]);
    }
    const Contig contig = make_contig(sequence);
    const ContigVariants variants(contig.bases, {}, contig.index);
    const htsim::insert_length::Parameters parameters{5U, 9U, 16U, 4.0};
    const VariableWgbsSampler reference(
        contig.bases,
        contig.index,
        91U,
        parameters,
        4U,
        true,
        0.0);
    const VariableHaplotypeSampler diploid(
        contig, variants, 91U, parameters, 4U, true, 0.0);
    require(diploid.allocation_weight() == reference.allocation_weight(),
            "event-free allocation weight changed");

    const auto reference_batch = reference.sample(0U, 40U);
    const auto diploid_batch = diploid.sample(0U, 40U);
    require(
        diploid_batch.skipped_count == reference_batch.skipped_count
            && diploid_batch.next_candidate_ordinal
                == reference_batch.next_candidate_ordinal
            && diploid_batch.candidates.size()
                == reference_batch.candidates.size(),
        "event-free candidate accounting changed");
    for (std::size_t index = 0U;
         index < reference_batch.candidates.size();
         ++index) {
        require(
            diploid_batch.candidates[index].reference_start
                    == reference_batch.candidates[index].reference_start
                && diploid_batch.candidates[index].reference_span
                    == reference_batch.candidates[index].insert_length
                && diploid_batch.candidates[index].eligible_haplotypes
                    == HaplotypeMask::both,
            "event-free variable candidate or mask changed");
    }
}

void test_indels_change_only_the_two_bit_eligibility_mask()
{
    const Contig insertion_contig = make_contig("AANAA");
    const ContigVariants insertion_variants(
        insertion_contig.bases,
        {{0U, 2U, 2U, VariantKind::insertion, {}, encode("CCC"),
          HaplotypeMask::haplotype_1}},
        insertion_contig.index);
    const VariableHaplotypeSampler insertion_sampler(
        insertion_contig,
        insertion_variants,
        1U,
        {4U, 4U, 4U, 0.0},
        3U,
        false,
        0.0);
    require(
        insertion_sampler.eligible_haplotypes(0U, 4U)
            == HaplotypeMask::haplotype_1,
        "insertion rescue set the wrong design bit");

    const Contig deletion_contig = make_contig("AAAANAAA");
    const ContigVariants deletion_variants(
        deletion_contig.bases,
        {{0U, 1U, 3U, VariantKind::deletion, encode("AA"), {},
          HaplotypeMask::haplotype_1}},
        deletion_contig.index);
    const VariableHaplotypeSampler deletion_sampler(
        deletion_contig,
        deletion_variants,
        1U,
        {5U, 5U, 5U, 0.0},
        3U,
        false,
        0.0);
    require(
        deletion_sampler.eligible_haplotypes(0U, 5U)
            == HaplotypeMask::haplotype_2,
        "deletion eligibility set the wrong design bit");
    require(
        !deletion_sampler.eligible_haplotypes(2U, 5U),
        "zero eligibility was exposed as a HaplotypeMask");
}

void test_exhaustive_interval_masks_match_projection()
{
    const Contig contig = make_contig("ANCGTNAACGTN");
    const std::vector<Event> events = {
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
        for (std::uint32_t read = 1U; read <= 3U; ++read) {
            for (const double ambiguous : {0.0, 0.5, 1.0}) {
                const VariableHaplotypeSampler sampler(
                    contig,
                    variants,
                    7U,
                    {read, 4U, 6U, 2.0},
                    read,
                    paired_end,
                    ambiguous);
                for (std::uint32_t span = read; span <= 6U; ++span) {
                    for (std::uint32_t start = 0U;
                         start + span <= contig.length;
                         ++start) {
                        require(
                            sampler.eligible_haplotypes(start, span)
                                == projected_mask(
                                    contig,
                                    variants,
                                    start,
                                    span,
                                    read,
                                    paired_end,
                                    ambiguous),
                            "variable interval mask disagreed with projection");
                    }
                }
            }
        }
    }
}

void test_sampling_is_chunk_independent_and_masks_are_nonzero()
{
    const Contig contig = make_contig("ANCGTNAACGTNACGTACGTACGT");
    const ContigVariants variants(
        contig.bases,
        {{0U, 3U, 5U, VariantKind::deletion, encode("GT"), {},
          HaplotypeMask::haplotype_1},
         {0U, 8U, 8U, VariantKind::insertion, {}, encode("CC"),
          HaplotypeMask::haplotype_2}},
        contig.index);
    const VariableHaplotypeSampler sampler(
        contig, variants, 43U, {3U, 6U, 9U, 3.0}, 3U, true, 0.0);
    const auto whole = sampler.sample(0U, 40U);
    std::vector<VariableHaplotypeCandidate> addressed;
    std::uint64_t addressed_ordinal = 0U;
    std::uint64_t addressed_skips = 0U;
    while (addressed.size() < 40U) {
        const auto candidate = sampler.candidate_at(addressed_ordinal);
        ++addressed_ordinal;
        if (candidate) {
            addressed.push_back(*candidate);
        } else {
            ++addressed_skips;
        }
    }
    require(
        addressed_skips == whole.skipped_count
            && addressed_ordinal == whole.next_candidate_ordinal
            && addressed.size() == whole.candidates.size(),
        "single-proposal traversal changed haplotype accounting");
    for (std::size_t index = 0U; index < addressed.size(); ++index) {
        require(
            addressed[index].reference_start
                    == whole.candidates[index].reference_start
                && addressed[index].reference_span
                    == whole.candidates[index].reference_span
                && addressed[index].eligible_haplotypes
                    == whole.candidates[index].eligible_haplotypes,
            "single-proposal interface changed a haplotype candidate");
    }
    const auto first = sampler.sample(0U, 13U);
    const auto second = sampler.sample(first.next_candidate_ordinal, 27U);
    std::vector<VariableHaplotypeCandidate> combined = first.candidates;
    combined.insert(
        combined.end(), second.candidates.begin(), second.candidates.end());
    require(
        whole.skipped_count > 0U
            && whole.skipped_count
                == first.skipped_count + second.skipped_count
            && whole.next_candidate_ordinal
                == second.next_candidate_ordinal
            && whole.candidates.size() == combined.size(),
        "chunking changed haplotype candidate accounting");
    for (std::size_t index = 0U; index < combined.size(); ++index) {
        require(
            whole.candidates[index].reference_start
                    == combined[index].reference_start
                && whole.candidates[index].reference_span
                    == combined[index].reference_span
                && whole.candidates[index].eligible_haplotypes
                    == combined[index].eligible_haplotypes
                && htsim::model::is_haplotype_mask(
                    static_cast<std::uint8_t>(
                        whole.candidates[index].eligible_haplotypes)),
            "chunking changed a variable haplotype candidate");
    }
}

void test_attempt_cap_and_invalid_domains_fail_closed()
{
    std::string sequence(20U, 'N');
    for (const std::size_t index : {0U, 1U, 2U, 7U, 8U, 9U}) {
        sequence[index] = 'A';
    }
    const Contig contig = make_contig(sequence, "chrCap");
    const ContigVariants variants(contig.bases, {}, contig.index);
    const VariableHaplotypeSampler impossible(
        contig,
        variants,
        3U,
        {10U, 10U, 10U, 0.0},
        3U,
        true,
        0.0,
        VariableHaplotypeOptions{1U});
    require(impossible.allocation_weight() > 0U,
            "attempt-cap fixture lost its maximum-span gate");
    require_error(
        [&] {(void)impossible.sample(0U, 1U);},
        "attempt-cap exhaustion was accepted");
    require_error(
        [&] {(void)impossible.eligible_haplotypes(0U, 2U);},
        "span below insert_min was accepted");
    require_error(
        [&] {(void)impossible.eligible_haplotypes(11U, 10U);},
        "interval beyond the contig was accepted");
    require_error(
        [&] {
            (void)VariableHaplotypeSampler(
                contig,
                variants,
                0U,
                {3U, 3U, 10U, 0.0},
                4U,
                true,
                0.0);
        },
        "read longer than insert_min was accepted");
    require_error(
        [&] {
            (void)VariableHaplotypeSampler(
                contig,
                variants,
                0U,
                {3U, 3U, 10U, 0.0},
                3U,
                true,
                0.0,
                VariableHaplotypeOptions{0U});
        },
        "zero attempt cap was accepted");

    const auto empty = impossible.sample(
        std::numeric_limits<std::uint64_t>::max(), 0U);
    require(empty.next_candidate_ordinal
                == std::numeric_limits<std::uint64_t>::max()
                && empty.candidates.empty(),
            "zero-count sample changed the ordinal");
}

} // namespace

int main()
{
    try {
        test_reference_only_matches_variable_wgbs_exactly();
        test_indels_change_only_the_two_bit_eligibility_mask();
        test_exhaustive_interval_masks_match_projection();
        test_sampling_is_chunk_independent_and_masks_are_nonzero();
        test_attempt_cap_and_invalid_domains_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "variable_haplotype_sampler_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
