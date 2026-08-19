#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

#include "protocol.h"
#include "rrbs.h"

namespace {

using htsim::model::Bases;
using htsim::rrbs::CandidateCatalog;
using htsim::rrbs::CandidateBed;
using htsim::rrbs::DiploidCandidateCatalog;
using htsim::rrbs::ProfileSampler;
using htsim::rrbs::RrbsCatalogError;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const RrbsCatalogError &) {
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

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-rrbs-candidates-XXXXXX";
        const int descriptor = mkstemp(pattern);
        if (descriptor < 0) {throw std::runtime_error("mkstemp failed");}
        path_ = pattern;
        if (close(descriptor) != 0) {
            throw std::runtime_error("temporary descriptor close failed");
        }
    }
    ~TempFile() {if (!path_.empty()) {(void)unlink(path_.c_str());}}
    const std::string &path() const noexcept {return path_;}

private:
    std::string path_;
};

void write_text(const std::string &path, const std::string &text)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output << text;
    if (!output) {throw std::runtime_error("temporary BED write failed");}
}

void test_parse_and_wildcard_matching()
{
    const auto sites = htsim::rrbs::parse_cut_sites(
        {"C|CGG", "CCTN|AGG"});
    require(sites.size() == 2
                && sites[0].motif == encode("CCGG")
                && sites[0].cut_offset == 1
                && sites[1].motif == encode("CCTNAGG")
                && sites[1].cut_offset == 4,
            "cut declarations were not parsed exactly");

    const auto cuts = htsim::rrbs::find_cut_positions(
        encode("CCGGACCTAAGGNCCTNAGG"), sites);
    require(cuts.size() == 2
                && cuts[0].position == 1
                && cuts[0].recognition_count == 1
                && cuts[1].position == 9
                && cuts[1].recognition_count == 1,
            "wildcard recognition or reference-N rejection changed");

    const auto shared_cut = htsim::rrbs::find_cut_positions(
        encode("CC"), htsim::rrbs::parse_cut_sites({"|C", "|CC"}));
    require(shared_cut.size() == 2
                && shared_cut[0].position == 0
                && shared_cut[0].recognition_count == 2
                && shared_cut[1].position == 1
                && shared_cut[1].recognition_count == 1,
            "motifs sharing one cut coordinate were not coalesced");

    require_error(
        [] {(void)htsim::rrbs::parse_cut_sites({"CCGG"});},
        "cut declaration without | was accepted");
    require_error(
        [] {(void)htsim::rrbs::parse_cut_sites({"CCGG|"});},
        "cut declaration with an empty right side was accepted");
    require_error(
        [] {(void)htsim::rrbs::parse_cut_sites({"C|CGG", "C|CGG"});},
        "duplicate cut declaration was accepted");
    require_error(
        [] {
            (void)htsim::rrbs::parse_cut_sites(
                {std::string(1025, 'A') + "|A"});
        },
        "oversized cut declaration was accepted");
}

void test_all_candidate_pairs_and_site_counts()
{
    const Bases bases = encode("CAACAAC");
    const auto sites = htsim::rrbs::parse_cut_sites({"|C"});
    const CandidateCatalog catalog(bases, sites, 3, 6, 2, true, 0.0);
    require(catalog.candidate_count() == 3,
            "all valid restriction-site pairs were not generated");
    require(catalog.allocation_weight() == 6U,
            "reference RRBS candidates did not count both physical copies");
    const auto &first = catalog.candidate(0);
    const auto &second = catalog.candidate(1);
    const auto &third = catalog.candidate(2);
    require(first.reference_start == 0 && first.reference_end == 3
                && first.gc_count == 1U
                && first.restriction_site_count == 2,
            "first RRBS candidate is wrong");
    require(second.reference_start == 0 && second.reference_end == 6
                && second.restriction_site_count == 3,
            "candidate did not include its internal restriction site");
    require(third.reference_start == 3 && third.reference_end == 6
                && third.restriction_site_count == 2,
            "last RRBS candidate is wrong");
    for (std::uint32_t index = 0; index < catalog.candidate_count(); ++index) {
        require(catalog.candidate(index).haplotype_mask
                    == htsim::model::HaplotypeMask::both,
                "reference-only candidate did not preserve the both-haplotype mask");
    }
}

void test_short_ids_bed_roundtrip_and_profile_sampling()
{
    const std::vector<htsim::rrbs::Candidate> candidates = {
        {1U, 5U, 4U, 2U, 2U,
         htsim::model::HaplotypeMask::haplotype_1, true, false},
        {1U, 5U, 4U, 2U, 2U,
         htsim::model::HaplotypeMask::haplotype_2, true, false},
        {6U, 10U, 4U, 1U, 3U,
         htsim::model::HaplotypeMask::both, true, false},
    };
    const auto ids = htsim::rrbs::candidate_ids("chrR", candidates);
    require(ids == std::vector<std::string>{
                "chrR:1-5~0", "chrR:1-5~1", "chrR:6-10"},
            "RRBS candidate IDs are not minimally disambiguated");

    std::ostringstream exported;
    htsim::rrbs::write_candidate_bed_header(exported);
    htsim::rrbs::write_candidate_bed_contig(
        exported, "chrR", candidates);
    const std::string text = exported.str();
    require(text.find("sha") == std::string::npos
                && text.find("fingerprint") == std::string::npos
                && text.find("@h") == std::string::npos,
            "RRBS exchange BED unexpectedly contains a hash or haplotype ID");

    const std::vector<htsim::reference::ContigMetadata> reference = {
        {"chrR", 20U, {}},
    };
    TempFile exported_file;
    write_text(exported_file.path(), text);
    const CandidateBed roundtrip(exported_file.path(), reference);
    require(roundtrip.row_count() == 3U
                && roundtrip.match_scores(
                    0U, "chrR", candidates, true)
                    == std::vector<double>{1.0, 1.0, 1.0},
            "exported RRBS candidates did not round-trip exactly");

    // Row order is intentionally different.  Only candidate_id controls the
    // join; a missing score is accepted for uniform matching and rejected for
    // profile matching.
    const std::string reordered_missing_score =
        "#chrom\tstart\tend\tcandidate_id\tscore\tstrand\thaplotype_mask"
        "\ttemplate_length\tgc_count\trestriction_site_count\n"
        "chrR\t6\t10\tchrR:6-10\t9\t.\t3\t4\t1\t3\n"
        "chrR\t1\t5\tchrR:1-5~1\t.\t.\t2\t4\t2\t2\n"
        "chrR\t1\t5\tchrR:1-5~0\t0\t.\t1\t4\t2\t2\n";
    TempFile missing_file;
    write_text(missing_file.path(), reordered_missing_score);
    const CandidateBed missing(missing_file.path(), reference);
    (void)missing.match_scores(0U, "chrR", candidates, false);
    require_error(
        [&] {(void)missing.match_scores(0U, "chrR", candidates, true);},
        "RRBS profile accepted a missing candidate score");

    const std::string scored =
        "chrR\t6\t10\tchrR:6-10\t0\t.\t3\t4\t1\t3\n"
        "chrR\t1\t5\tchrR:1-5~1\t5\t.\t2\t4\t2\t2\n"
        "chrR\t1\t5\tchrR:1-5~0\t0\t.\t1\t4\t2\t2\n";
    TempFile scored_file;
    write_text(scored_file.path(), scored);
    const CandidateBed scored_bed(scored_file.path(), reference);
    const auto scores = scored_bed.match_scores(
        0U, "chrR", candidates, true);
    const ProfileSampler sampler(candidates, scores);
    require(sampler.allocation_weight() == 5.0,
            "RRBS profile allocation did not retain relative score mass");
    const auto all = sampler.sample_indices(0U, 73U, 0U, 23U);
    const auto first = sampler.sample_indices(0U, 73U, 0U, 7U);
    const auto second = sampler.sample_indices(0U, 73U, 7U, 16U);
    std::vector<std::uint32_t> joined = first;
    joined.insert(joined.end(), second.begin(), second.end());
    require(all == joined
                && std::all_of(all.begin(), all.end(),
                    [](std::uint32_t index) {return index == 1U;}),
            "RRBS profile sampling changed across chunks or selected zero weight");
    const ProfileSampler both_copy_sampler(
        candidates, std::vector<double>{0.0, 0.0, 3.0});
    require(both_copy_sampler.allocation_weight() == 6.0,
            "mask-3 RRBS profile score did not carry two-copy mass");

    const std::string changed_fixed_field =
        "chrR\t1\t5\tchrR:1-5~0\t1\t.\t1\t4\t3\t2\n"
        "chrR\t1\t5\tchrR:1-5~1\t1\t.\t2\t4\t2\t2\n"
        "chrR\t6\t10\tchrR:6-10\t1\t.\t3\t4\t1\t3\n";
    TempFile changed_file;
    write_text(changed_file.path(), changed_fixed_field);
    const CandidateBed changed(changed_file.path(), reference);
    require_error(
        [&] {(void)changed.match_scores(0U, "chrR", candidates, true);},
        "RRBS candidate BED accepted a changed fixed feature");

    const std::string missing_row =
        "chrR\t1\t5\tchrR:1-5~0\t1\t.\t1\t4\t2\t2\n"
        "chrR\t1\t5\tchrR:1-5~1\t1\t.\t2\t4\t2\t2\n";
    TempFile missing_row_file;
    write_text(missing_row_file.path(), missing_row);
    const CandidateBed incomplete(missing_row_file.path(), reference);
    require_error(
        [&] {(void)incomplete.match_scores(0U, "chrR", candidates, true);},
        "RRBS candidate BED accepted a missing candidate row");

    const std::string duplicate_id =
        "chrR\t1\t5\tchrR:1-5~0\t1\t.\t1\t4\t2\t2\n"
        "chrR\t1\t5\tchrR:1-5~0\t1\t.\t2\t4\t2\t2\n";
    TempFile duplicate_file;
    write_text(duplicate_file.path(), duplicate_id);
    require_error(
        [&] {(void)CandidateBed(duplicate_file.path(), reference);},
        "RRBS candidate BED accepted a duplicate candidate ID");
}

void test_ambiguity_filter_and_coordinate_edges()
{
    const auto sites = htsim::rrbs::parse_cut_sites({"|C"});
    const CandidateCatalog rejected(
        encode("CNANC"), sites, 4, 4, 2, true, 0.0);
    require(rejected.candidate_count() == 0,
            "candidate with ambiguous mate bases was accepted");

    const CandidateCatalog accepted(
        encode("CNANC"), sites, 4, 4, 2, true, 0.5);
    require(accepted.candidate_count() == 1
                && accepted.candidate(0).reference_start == 0
                && accepted.candidate(0).reference_end == 4,
            "fraction boundary or cut-at-zero coordinate is wrong");
}

void test_sampling_is_chunk_independent()
{
    const auto sites = htsim::rrbs::parse_cut_sites({"|C"});
    const CandidateCatalog catalog(
        encode("CAACAACAAC"), sites, 3, 9, 2, false, 0.0);
    const auto all = catalog.sample_indices(0U, 17, 0, 9);
    const auto first = catalog.sample_indices(0U, 17, 0, 4);
    const auto second = catalog.sample_indices(0U, 17, 4, 5);
    std::vector<std::uint32_t> joined = first;
    joined.insert(joined.end(), second.begin(), second.end());
    require(all == joined, "RRBS sampling changed across chunk boundaries");
    require(all == catalog.sample_indices(0U, 17, 0, 9),
            "RRBS sampling is not deterministic");
    for (const std::uint32_t index : all) {
        require(index < catalog.candidate_count(),
                "RRBS sample index is outside its catalog");
    }

    require_error(
        [&] {(void)catalog.candidate(catalog.candidate_count());},
        "out-of-range candidate lookup was accepted");
    const CandidateCatalog empty(
        encode("AAAA"), sites, 2, 4, 1, false, 0.0);
    require_error(
        [&] {(void)empty.sample_indices(0U, 17, 0, 1);},
        "sampling an empty RRBS catalog was accepted");
    require_error(
        [&] {
            (void)catalog.sample_indices(
                0U, 17, std::numeric_limits<std::uint64_t>::max(), 2);
        },
        "overflowing candidate ordinal range was accepted");
}

htsim::reference::Contig make_contig(const std::string &sequence)
{
    htsim::reference::Contig contig;
    contig.index = 0U;
    contig.name = "chrDiploidRrbs";
    contig.bases = encode(sequence);
    contig.length = contig.bases.size();
    return contig;
}

void test_diploid_motifs_follow_haplotype_bits()
{
    const auto contig = make_contig("AACCGGAAACCGGAA");
    const std::vector<htsim::variant::Event> events = {
        {0U, 10U, 11U, htsim::model::VariantKind::snv,
         encode("C"), encode("T"),
         htsim::model::HaplotypeMask::haplotype_1},
    };
    const htsim::variant::ContigVariants variants(
        contig.bases, events, contig.index);
    const DiploidCandidateCatalog catalog(
        contig,
        variants,
        htsim::rrbs::parse_cut_sites({"C|CGG"}),
        7U,
        7U,
        2U,
        true,
        0.0);
    require(catalog.candidate_count() == 1U,
            "haplotype-specific motif loss did not change the RRBS domain");
    const auto &candidate = catalog.candidate(0U);
    require(candidate.reference_start == 3U
                && candidate.reference_end == 10U
                && candidate.template_length == 7U
                && candidate.haplotype_mask
                    == htsim::model::HaplotypeMask::haplotype_2,
            "RRBS motif availability used the wrong two-bit haplotype mask");
}

void test_diploid_indel_lengths_are_physical()
{
    const auto contig = make_contig("AACCGGAAACCGGAA");
    const std::vector<htsim::variant::Event> events = {
        {0U, 6U, 6U, htsim::model::VariantKind::insertion,
         {}, encode("TT"),
         htsim::model::HaplotypeMask::haplotype_1},
    };
    const htsim::variant::ContigVariants variants(
        contig.bases, events, contig.index);
    const DiploidCandidateCatalog catalog(
        contig,
        variants,
        htsim::rrbs::parse_cut_sites({"C|CGG"}),
        7U,
        9U,
        2U,
        true,
        0.0);
    require(catalog.candidate_count() == 2U,
            "diploid RRBS did not retain one physical candidate per haplotype");
    require(catalog.allocation_weight() == 2U,
            "diploid RRBS allocation did not count physical copies fairly");
    const auto &haplotype_0 = catalog.candidate(0U);
    const auto &haplotype_1 = catalog.candidate(1U);
    require(haplotype_0.reference_start == 3U
                && haplotype_0.reference_end == 10U
                && haplotype_0.template_length == 9U
                && haplotype_0.haplotype_mask
                    == htsim::model::HaplotypeMask::haplotype_1,
            "RRBS insertion did not extend the selected haplotype fragment");
    require(haplotype_1.reference_start == 3U
                && haplotype_1.reference_end == 10U
                && haplotype_1.template_length == 7U
                && haplotype_1.haplotype_mask
                    == htsim::model::HaplotypeMask::haplotype_2,
            "RRBS reference haplotype length changed beside an insertion");
}

} // namespace

int main()
{
    try {
        test_parse_and_wildcard_matching();
        test_all_candidate_pairs_and_site_counts();
        test_short_ids_bed_roundtrip_and_profile_sampling();
        test_ambiguity_filter_and_coordinate_edges();
        test_sampling_is_chunk_independent();
        test_diploid_motifs_follow_haplotype_bits();
        test_diploid_indel_lengths_are_physical();
    } catch (const std::exception &error) {
        std::cerr << "rrbs_catalog_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
