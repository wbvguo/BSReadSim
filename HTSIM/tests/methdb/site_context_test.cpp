#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include "methdb.h"
#include "protocol.h"

namespace {

using htsim::model::Bases;
using htsim::model::MethylationContext;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
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

void require_context(
    const std::string &sequence,
    std::uint64_t position,
    MethylationContext expected)
{
    require(
        htsim::methdb::classify_context(
            encode(sequence), position, true)
            == std::optional<MethylationContext>(expected),
        "wrong context for " + sequence);
}

void test_six_contexts()
{
    require_context("CG", 0, MethylationContext::cg_c);
    require_context("CAG", 0, MethylationContext::chg_c);
    require_context("CTG", 0, MethylationContext::chg_c);
    require_context("CAA", 0, MethylationContext::chh_c);
    require_context("CCT", 0, MethylationContext::chh_c);
    require_context("CG", 1, MethylationContext::cg_g);
    require_context("CAG", 2, MethylationContext::chg_g);
    require_context("AAG", 2, MethylationContext::chh_g);
    require_context("TTG", 2, MethylationContext::chh_g);
}

void test_boundaries_and_ambiguous_flanks()
{
    using htsim::methdb::classify_context;
    const auto none = std::optional<MethylationContext>{};
    require(classify_context(encode("C"), 0, true) == none,
            "terminal C was classified without a flank");
    require(classify_context(encode("CA"), 0, true) == none,
            "C was classified without its second flank");
    require(classify_context(encode("G"), 0, true) == none,
            "initial G was classified without a flank");
    require(classify_context(encode("AG"), 1, true) == none,
            "G was classified without its second flank");
    for (const auto &fixture : {
             std::pair<std::string, std::uint64_t>{"CN", 0},
             {"CAN", 0},
             {"NG", 1},
             {"NAG", 2},
         }) {
        require(classify_context(encode(fixture.first), fixture.second, true) == none,
                "N in a required flank was classified");
    }
    require(classify_context(encode("ACGTA"), 1, true)
                == MethylationContext::cg_c,
            "context was not evaluated from the complete contig");
    require(classify_context(encode("CGN"), 0, true)
                == MethylationContext::cg_c,
            "unused third CpG flank incorrectly blocked classification");
    require(classify_context(encode("NCG"), 2, true)
                == MethylationContext::cg_g,
            "unused third reverse-CpG flank incorrectly blocked classification");
}

void test_cpg_only_and_invalid_values()
{
    using htsim::methdb::classify_context;
    require(classify_context(encode("CG"), 0, false)
                == MethylationContext::cg_c,
            "CpG C was removed by the CpG-only filter");
    require(classify_context(encode("CG"), 1, false)
                == MethylationContext::cg_g,
            "CpG G was removed by the CpG-only filter");
    require(!classify_context(encode("CAG"), 0, false).has_value(),
            "CHG C survived the CpG-only filter");
    require(!classify_context(encode("CAG"), 2, false).has_value(),
            "CHG G survived the CpG-only filter");
    require(!classify_context(encode("A"), 0, true).has_value(),
            "non-C/G center was classified");

    try {
        (void)classify_context(encode("A"), 1, true);
        throw std::runtime_error("out-of-range position was accepted");
    } catch (const htsim::methdb::ContextError &) {
    }
    try {
        (void)classify_context({}, UINT64_MAX, true);
        throw std::runtime_error("maximum position on an empty contig was accepted");
    } catch (const htsim::methdb::ContextError &) {
    }
    try {
        (void)classify_context(Bases({1, 5}), 0, true);
        throw std::runtime_error("invalid protocol base was accepted");
    } catch (const htsim::methdb::ContextError &) {
    }
    try {
        (void)classify_context(Bases({5, 2}), 1, true);
        throw std::runtime_error("invalid reverse flank base was accepted");
    } catch (const htsim::methdb::ContextError &) {
    }
}

void test_explicit_haplotype_neighborhood()
{
    using htsim::methdb::ContextNeighborhood;
    using htsim::methdb::classify_context;
    require(classify_context(
                ContextNeighborhood{0U, 3U, 2U, std::nullopt, std::nullopt},
                true)
                == MethylationContext::chh_g,
            "explicit upstream haplotype neighborhood was misclassified");
    require(classify_context(
                ContextNeighborhood{
                    std::nullopt, std::nullopt, 1U, 0U, 2U},
                true)
                == MethylationContext::chg_c,
            "explicit downstream haplotype neighborhood was misclassified");
    require(!classify_context(
                 ContextNeighborhood{
                     std::nullopt, std::nullopt, 1U, 0U, std::nullopt},
                 true).has_value(),
            "missing true haplotype flank was accepted");
    try {
        (void)classify_context(
            ContextNeighborhood{
                std::nullopt, std::nullopt, 1U, 5U, 0U},
            true);
        throw std::runtime_error("invalid explicit flank was accepted");
    } catch (const htsim::methdb::ContextError &) {
    }
}

} // namespace

int main()
{
    try {
        test_six_contexts();
        test_boundaries_and_ambiguous_flanks();
        test_cpg_only_and_invalid_values();
        test_explicit_haplotype_neighborhood();
    } catch (const std::exception &error) {
        std::cerr << "site_context_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
