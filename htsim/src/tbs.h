#ifndef HTSIM_TBS_H
#define HTSIM_TBS_H

#include <array>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "types.h"
#include "variant.h"
#include "reference.h"
#include "utilities.h"

// ---- catalog --------------------------------------------------------

namespace htsim::tbs {

inline constexpr std::uint32_t maximum_attempts_per_fragment = 100000;

enum class SamplingMode : std::uint8_t {
    uniform = 0,
    output_weight = 1,
};

class TbsCatalogError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct Target {
    std::uint32_t contig_index = 0;
    std::uint32_t interval_start = 0;
    std::uint32_t interval_end = 0;
    std::string name;
    double score = 0.0;
    model::CaptureStrand capture_strand =
        model::CaptureStrand::unknown;
};

// A fully verified and canonically ordered projection of BED6 target rows.
// Targets are grouped by reference FASTA index and sorted by coordinate,
// strand, name, and score. Exact coordinate+strand duplicates are rejected.
class TargetFile {
public:
    TargetFile(
        const std::string &path,
        const std::vector<reference::ContigMetadata> &reference_catalog);

    const std::vector<Target> &targets(std::uint32_t contig_index) const;
    std::uint64_t target_count() const noexcept {return target_count_;}
    const crypto::Sha256Digest &file_sha256() const noexcept {return file_sha256_;}

private:
    crypto::Sha256Digest file_sha256_ = {};
    std::vector<std::vector<Target>> targets_by_contig_;
    std::uint64_t target_count_ = 0;
};

struct Candidate {
    std::uint32_t reference_start = 0;
    std::uint32_t reference_end = 0;
    std::uint32_t template_length = 0;
    std::uint32_t target_start = 0;
    std::uint32_t target_end = 0;
    std::uint32_t target_ordinal = 0;
    double output_weight = 0.0;
    model::CaptureStrand capture_strand =
        model::CaptureStrand::unknown;
    model::HaplotypeMask haplotype_mask =
        model::HaplotypeMask::both;
    bool include_start_anchor_insertion = true;
    bool include_end_anchor_insertion = false;
};

struct SampleBatch {
    std::vector<Candidate> candidates;
    std::uint64_t skipped_count = 0;
};

// Target-first diploid TBS sampler. A BED row is one logical capture target and
// contributes exactly one uniform choice or one target-level output weight.
// After selecting that target, the sampler chooses uniformly among its
// eligible haplotypes and then materializes a physical haplotype slice.
class DiploidCandidateCatalog {
public:
    DiploidCandidateCatalog(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        const std::vector<Target> &targets,
        double center_stddev,
        std::uint32_t insert_length,
        std::uint32_t read_length,
        bool paired_end,
        double max_ambiguous_fraction,
        SamplingMode sampling_mode = SamplingMode::uniform);

    std::uint32_t choice_count() const noexcept;
    std::uint32_t allocation_weight() const noexcept;
    SampleBatch sample(
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        std::uint64_t first_candidate_ordinal,
        std::uint64_t first_fragment_ordinal,
        std::uint32_t output_count) const;

private:
    struct TargetChoice {
        std::uint32_t target_start = 0;
        std::uint32_t target_end = 0;
        std::uint32_t target_ordinal = 0;
        std::array<std::uint32_t, 2> haplotype_centers = {};
        double output_weight = 0.0;
        model::CaptureStrand capture_strand =
            model::CaptureStrand::unknown;
        model::HaplotypeMask eligible_haplotypes =
            model::HaplotypeMask::both;
    };

    bool project(
        const TargetChoice &target,
        std::uint8_t haplotype,
        std::int64_t center_displacement,
        Candidate &candidate) const;
    bool sequenceable(
        const haplotype::HaplotypeLayout &layout,
        std::uint32_t fragment_start,
        std::uint32_t fragment_end) const;

    std::array<std::unique_ptr<haplotype::HaplotypeLayout>, 2> layouts_;
    std::vector<TargetChoice> target_choices_;
    std::vector<std::uint32_t> cumulative_weights_;
    double center_stddev_ = 0.0;
    std::uint32_t insert_length_ = 0;
    std::uint32_t read_length_ = 0;
    std::uint32_t maximum_ambiguous_count_ = 0;
    std::uint32_t total_weight_ = 0;
    bool paired_end_ = false;
    SamplingMode sampling_mode_ = SamplingMode::uniform;
};

// Per-contig fixed-insert target sampler. At center_stddev=0, one eligible BED
// target produces one fixed choice. At center_stddev>0, every target on a
// contig with at least one sequenceable start is a choice; each output first
// selects one target, then retries normal center displacements for that same
// target until it obtains a valid fragment or exhausts the fixed attempt cap.
class CandidateCatalog {
public:
    CandidateCatalog(
        const model::Bases &contig_bases,
        const std::vector<Target> &targets,
        std::uint32_t contig_index,
        double center_stddev,
        std::uint32_t insert_length,
        std::uint32_t read_length,
        bool paired_end,
        double max_ambiguous_fraction,
        SamplingMode sampling_mode = SamplingMode::uniform);

    std::uint32_t choice_count() const noexcept;
    std::uint32_t allocation_weight() const noexcept;

    SampleBatch sample(
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count) const;

private:
    struct TargetChoice {
        std::uint32_t target_start = 0;
        std::uint32_t target_end = 0;
        std::uint32_t target_ordinal = 0;
        double output_weight = 0.0;
        model::CaptureStrand capture_strand =
            model::CaptureStrand::unknown;
    };

    bool project(
        const TargetChoice &target,
        std::int64_t center_displacement,
        Candidate &candidate) const;
    std::uint32_t ambiguous_count(
        std::uint32_t begin,
        std::uint32_t end) const;
    bool sequenceable(
        std::uint32_t fragment_start,
        std::uint32_t fragment_end) const;

    std::vector<std::uint64_t> ambiguous_words_;
    std::vector<std::uint32_t> ambiguous_prefix_;
    std::vector<TargetChoice> target_choices_;
    std::vector<Candidate> fixed_candidates_;
    std::vector<std::uint32_t> cumulative_weights_;
    double center_stddev_ = 0.0;
    std::uint32_t contig_length_ = 0;
    std::uint32_t insert_length_ = 0;
    std::uint32_t read_length_ = 0;
    std::uint32_t maximum_ambiguous_count_ = 0;
    std::uint32_t total_weight_ = 0;
    bool paired_end_ = false;
    SamplingMode sampling_mode_ = SamplingMode::uniform;
};

} // namespace htsim::tbs

#endif // HTSIM_TBS_H
