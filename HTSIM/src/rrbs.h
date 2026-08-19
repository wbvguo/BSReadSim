#ifndef HTSIM_RRBS_H
#define HTSIM_RRBS_H

#include <cstdint>
#include <iosfwd>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "types.h"
#include "reference.h"
#include "variant.h"

// ---- catalog --------------------------------------------------------

namespace htsim::rrbs {

inline constexpr std::uint32_t maximum_motif_length = 1024;

class RrbsCatalogError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct CutSite {
    model::Bases motif;
    std::uint32_t cut_offset = 0;
};

struct CutPosition {
    // A cut is a boundary between bases and may equal zero or contig length.
    std::uint32_t position = 0;
    // Multiple configured motifs may recognize the same cut boundary.
    std::uint32_t recognition_count = 0;
};

struct Candidate {
    std::uint32_t reference_start = 0;
    std::uint32_t reference_end = 0;
    // Physical length on the selected haplotype.  It equals the reference span
    // only when no active indel changes the fragment.
    std::uint32_t template_length = 0;
    // Number of physical C/G bases in the enzyme-bounded template.  For a
    // haplotype catalog this is measured on that constructed haplotype, not
    // on the reference envelope.
    std::uint32_t gc_count = 0;
    std::uint32_t restriction_site_count = 0;
    // Typed view of the frozen availability bits. It is not a sampled
    // zero-based protocol haplotype.
    model::HaplotypeMask haplotype_mask =
        model::HaplotypeMask::both;
    bool include_start_anchor_insertion = true;
    bool include_end_anchor_insertion = false;
};

// Plain, BED-compatible RRBS exchange rows.  The score is deliberately the
// sole mutable field: an external model may replace it, while C++ regenerates
// and verifies every candidate field before sampling.  No row or file digest
// is part of this contract.
struct CandidateBedRow {
    std::uint32_t reference_start = 0;
    std::uint32_t reference_end = 0;
    std::string candidate_id;
    std::optional<double> score;
    model::HaplotypeMask haplotype_mask = model::HaplotypeMask::both;
    std::uint32_t template_length = 0;
    std::uint32_t gc_count = 0;
    std::uint32_t restriction_site_count = 0;
};

// A short, deterministic, globally unique ID for each row in one contig.
// The common case is "contig:start-end".  When multiple physical candidates
// share that reference envelope, all members receive a local base-36 suffix
// such as "~0" and "~1".  The suffix is only a disambiguator; haplotype_mask
// remains the authoritative biological field.
std::vector<std::string> candidate_ids(
    std::string_view contig_name,
    const std::vector<Candidate> &candidates);

void write_candidate_bed_header(std::ostream &output);
void write_candidate_bed_contig(
    std::ostream &output,
    std::string_view contig_name,
    const std::vector<Candidate> &candidates);

// Parsed external-score table grouped by reference contig.  Construction
// validates syntax and unique IDs but performs no hash check.  match_scores()
// then requires an exact one-to-one match with the regenerated C++ catalog;
// row order is irrelevant and only score may differ.
class CandidateBed {
public:
    CandidateBed(
        const std::string &path,
        const std::vector<reference::ContigMetadata> &reference_catalog);

    std::vector<double> match_scores(
        std::uint32_t contig_index,
        std::string_view contig_name,
        const std::vector<Candidate> &candidates,
        bool require_scores) const;

    std::uint64_t row_count() const noexcept {return row_count_;}

private:
    std::vector<std::vector<CandidateBedRow>> rows_by_contig_;
    std::uint64_t row_count_ = 0;
};

// Reusable categorical sampler for model scores.  A score is a non-negative
// relative weight per physical haplotype copy, so a row with mask=3 carries
// twice the mass of an otherwise identical row with mask=1 or mask=2.
class ProfileSampler {
public:
    ProfileSampler(
        const std::vector<Candidate> &candidates,
        const std::vector<double> &scores);

    double allocation_weight() const noexcept {return allocation_weight_;}
    std::vector<std::uint32_t> sample_indices(
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count) const;

private:
    std::vector<long double> cumulative_weights_;
    long double total_weight_ = 0.0L;
    double allocation_weight_ = 0.0;
};

std::vector<CutSite> parse_cut_sites(
    const std::vector<std::string> &declarations);

// Motif N matches any concrete A/C/G/T reference base, but never reference N.
// Overlapping motif matches are retained. Matches sharing a cut coordinate are
// coalesced and counted in CutPosition::recognition_count.
std::vector<CutPosition> find_cut_positions(
    const model::Bases &contig_bases,
    const std::vector<CutSite> &cut_sites);

class CandidateCatalog {
public:
    CandidateCatalog(
        const model::Bases &contig_bases,
        const std::vector<CutSite> &cut_sites,
        std::uint32_t minimum_insert_length,
        std::uint32_t maximum_insert_length,
        std::uint32_t read_length,
        bool paired_end,
        double max_ambiguous_fraction);

    std::uint32_t candidate_count() const noexcept;
    std::uint32_t allocation_weight() const;
    const Candidate &candidate(std::uint32_t index) const;
    const std::vector<Candidate> &candidates() const noexcept {
        return candidates_;
    }

    // Uniform sampling with replacement. The RNG address is the shared
    // fragment-stage contract: contig key, per-contig candidate ordinal, and
    // local index 1. Returned indices remain in candidate-ordinal order.
    std::vector<std::uint32_t> sample_indices(
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count) const;

private:
    std::vector<Candidate> candidates_;
};

// Variant-aware RRBS catalog.  Restriction motifs and mate ambiguity are
// evaluated independently on the two already constructed haplotypes.  Each
// candidate carries exactly one of the frozen two mask bits; sampling converts
// that mask to protocol haplotype 0/1 only after validation.
class DiploidCandidateCatalog {
public:
    DiploidCandidateCatalog(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        const std::vector<CutSite> &cut_sites,
        std::uint32_t minimum_insert_length,
        std::uint32_t maximum_insert_length,
        std::uint32_t read_length,
        bool paired_end,
        double max_ambiguous_fraction);

    std::uint32_t candidate_count() const noexcept;
    std::uint32_t allocation_weight() const;
    const Candidate &candidate(std::uint32_t index) const;
    const std::vector<Candidate> &candidates() const noexcept {
        return candidates_;
    }
    std::vector<std::uint32_t> sample_indices(
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count) const;

private:
    std::vector<Candidate> candidates_;
};

} // namespace htsim::rrbs

#endif // HTSIM_RRBS_H
