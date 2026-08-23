#include "variant.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "types.h"
#include "protocol.h"
#include "reference.h"

namespace {

using htsim::model::HaplotypeMask;
using htsim::model::Bases;
using htsim::model::VariantKind;
using htsim::reference::Contig;
using htsim::variant::ContigVariants;
using htsim::variant::Variant;
using htsim::variant::MutationCatalogError;
using htsim::variant::MutationParameters;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const MutationCatalogError &) {
        return;
    }
    throw std::runtime_error(message);
}

Bases encode(const std::string &sequence)
{
    Bases result;
    result.reserve(sequence.size());
    for (const char base : sequence) {
        switch (base) {
        case 'A': result.push_back(0U); break;
        case 'C': result.push_back(1U); break;
        case 'G': result.push_back(2U); break;
        case 'T': result.push_back(3U); break;
        case 'N': result.push_back(4U); break;
        default: throw std::runtime_error("invalid fixture base");
        }
    }
    return result;
}

Contig make_contig(
    const std::string &sequence,
    const std::string &name = "chrMutation")
{
    Contig contig;
    contig.index = 7U;
    contig.name = name;
    contig.bases = encode(sequence);
    contig.length = contig.bases.size();
    return contig;
}

char base_character(std::uint8_t base)
{
    return "ACGTN"[base];
}

std::string bases_text(const Bases &bases)
{
    std::string result;
    for (const std::uint8_t base : bases) {
        result.push_back(base_character(base));
    }
    return result.empty() ? "-" : result;
}

std::string event_summary(const std::vector<Variant> &events)
{
    std::ostringstream summary;
    for (const Variant &event : events) {
        const char kind = event.kind == VariantKind::snv
            ? 'S'
            : event.kind == VariantKind::insertion ? 'I' : 'D';
        summary << kind << ':' << event.reference_start << '-'
                << event.reference_end << ':' << bases_text(event.ref_bases)
                << '>' << bases_text(event.alt_bases) << ':'
                << static_cast<unsigned>(event.alt_haplotypes) << ';';
    }
    return summary.str();
}

void test_snv_only_and_unresolved_reference()
{
    const Contig contig = make_contig("ACGTNACGT");
    const MutationParameters parameters{1.0, 0.0, 1.0, true};
    const auto events = htsim::variant::generate_de_novo_events(
        contig, UINT64_C(123), parameters);
    require(events.size() == 8U, "SNV-only run mutated N or missed a base");
    for (const Variant &event : events) {
        require(
            event.kind == VariantKind::snv
                && event.reference_end == event.reference_start + 1U
                && event.ref_bases.size() == 1U
                && event.alt_bases.size() == 1U
                && event.ref_bases != event.alt_bases
                && event.alt_haplotypes == HaplotypeMask::both,
            "SNV-only event shape or homozygous mask changed");
        require(
            event.reference_start != 4U,
            "unresolved reference position became a mutation event");
    }
}

void test_indel_shape_extension_and_canonical_order()
{
    const Contig contig = make_contig("ACGTACGTACGTACGTACGTACGT");
    const MutationParameters one_base{1.0, 1.0, 0.0, false};
    const auto short_events = htsim::variant::generate_de_novo_events(
        contig, UINT64_C(9), one_base);
    bool saw_insertion = false;
    bool saw_deletion = false;
    for (const Variant &event : short_events) {
        require(event.kind != VariantKind::snv, "indel-only run emitted an SNV");
        saw_insertion = saw_insertion || event.kind == VariantKind::insertion;
        saw_deletion = saw_deletion || event.kind == VariantKind::deletion;
        const std::size_t length = event.kind == VariantKind::insertion
            ? event.alt_bases.size()
            : event.ref_bases.size();
        require(length == 1U, "zero extension did not produce one-base indels");
        require(
            htsim::model::is_haplotype_mask(
                static_cast<std::uint8_t>(event.alt_haplotypes)),
            "de novo indel has an invalid haplotype mask");
    }
    require(saw_insertion && saw_deletion,
            "indel fixture did not cover insertion and deletion");
    (void)ContigVariants(contig.bases, short_events, contig.index);

    const MutationParameters extended{1.0, 1.0, 1.0, true};
    const auto long_events = htsim::variant::generate_de_novo_events(
        contig, UINT64_C(9), extended);
    bool saw_full_length = false;
    for (const Variant &event : long_events) {
        const std::size_t length = event.kind == VariantKind::insertion
            ? event.alt_bases.size()
            : event.ref_bases.size();
        require(length >= 1U && length <= 4U,
                "extended indel escaped the typed four-base boundary");
        saw_full_length = saw_full_length || length == 4U;
    }
    require(saw_full_length, "extension fixture did not reach four bases");
    (void)ContigVariants(contig.bases, long_events, contig.index);
}

void test_addressed_golden_and_domain_separation()
{
    const Contig contig = make_contig(
        "ACGTACGTACGTACGTACGTACGTACGTACGT");
    const MutationParameters parameters{0.8, 0.65, 0.6, false};
    const auto first = htsim::variant::generate_de_novo_events(
        contig, UINT64_C(0x123456789abcdef0), parameters);
    const auto repeated = htsim::variant::generate_de_novo_events(
        contig, UINT64_C(0x123456789abcdef0), parameters);
    (void)ContigVariants(contig.bases, first, contig.index);
    const std::string summary = event_summary(first);
    require(summary == event_summary(repeated),
            "same mutation address did not reproduce its event stream");
    require(
        summary ==
            "S:0-1:A>C:1;I:2-2:->CT:3;S:3-4:T>A:3;"
            "D:6-8:GT>-:3;I:9-9:->G:3;D:10-11:G>-:3;"
            "S:12-13:A>C:1;S:13-14:C>T:2;I:15-15:->TG:1;"
            "D:16-19:ACG>-:2;D:19-23:TACG>-:1;"
            "D:23-27:TACG>-:3;D:27-30:TAC>-:3;"
            "D:30-31:G>-:3;S:31-32:T>C:2;",
        "de novo mutation golden changed: " + summary);

    const auto other_seed = htsim::variant::generate_de_novo_events(
        contig, UINT64_C(0x123456789abcdef1), parameters);
    Contig other_contig = make_contig(
        "ACGTACGTACGTACGTACGTACGTACGTACGT",
        "chrOther");
    other_contig.index = 8U;
    const auto other_index = htsim::variant::generate_de_novo_events(
        other_contig, UINT64_C(0x123456789abcdef0), parameters);
    require(summary != event_summary(other_seed)
                && summary != event_summary(other_index),
            "mutation seed or contig domain did not isolate its stream");
}

void test_snv_distribution_sanity()
{
    std::string sequence;
    sequence.reserve(100000U);
    for (std::size_t index = 0U; index < 25000U; ++index) {
        sequence += "ACGT";
    }
    const Contig contig = make_contig(sequence, "chrDistribution");
    const MutationParameters parameters{0.2, 0.0, 0.0, false};
    const auto events = htsim::variant::generate_de_novo_events(
        contig, UINT64_C(20260812), parameters);
    require(events.size() > 19000U && events.size() < 21000U,
            "SNV mutation-rate sanity check escaped its tolerance");

    std::uint32_t mask_counts[3] = {0U, 0U, 0U};
    std::uint32_t alternate_offsets[3] = {0U, 0U, 0U};
    for (const Variant &event : events) {
        require(event.kind == VariantKind::snv,
                "SNV distribution fixture emitted an indel");
        const std::uint8_t mask =
            static_cast<std::uint8_t>(event.alt_haplotypes);
        require(mask >= 1U && mask <= 3U,
                "SNV distribution fixture emitted an invalid mask");
        ++mask_counts[mask - 1U];
        const std::uint8_t offset = static_cast<std::uint8_t>(
            (event.alt_bases.front() + 4U - event.ref_bases.front()) & 3U);
        require(offset >= 1U && offset <= 3U,
                "SNV alternate base did not differ from REF");
        ++alternate_offsets[offset - 1U];
    }
    const double expected = static_cast<double>(events.size()) / 3.0;
    for (std::size_t index = 0U; index < 3U; ++index) {
        require(
            std::abs(static_cast<double>(mask_counts[index]) - expected)
                    < expected * 0.05
                && std::abs(
                       static_cast<double>(alternate_offsets[index]) - expected)
                    < expected * 0.05,
            "SNV allele or haplotype distribution escaped its tolerance");
    }
}

void test_invalid_inputs_fail_closed()
{
    Contig contig = make_contig("ACGT");
    const MutationParameters valid{0.1, 0.2, 0.3, false};
    contig.length += 1U;
    require_error(
        [&] {
            (void)htsim::variant::generate_de_novo_events(
                contig, 1U, valid);
        },
        "contig length mismatch was accepted");

    contig = make_contig("ACGT");
    contig.bases[2] = 5U;
    require_error(
        [&] {
            (void)htsim::variant::generate_de_novo_events(
                contig, 1U, valid);
        },
        "invalid base code was accepted");

    contig = make_contig("ACGT");
    for (const MutationParameters &invalid : {
             MutationParameters{-0.1, 0.2, 0.3, false},
             MutationParameters{0.1, 1.1, 0.3, false},
             MutationParameters{
                 0.1,
                 0.2,
                 std::numeric_limits<double>::quiet_NaN(),
                 false},
         }) {
        require_error(
            [&] {
                (void)htsim::variant::generate_de_novo_events(
                    contig, 1U, invalid);
            },
            "invalid mutation probability was accepted");
    }
}

} // namespace

int main()
{
    try {
        test_snv_only_and_unresolved_reference();
        test_indel_shape_extension_and_canonical_order();
        test_addressed_golden_and_domain_separation();
        test_snv_distribution_sanity();
        test_invalid_inputs_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "mutation_catalog_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
