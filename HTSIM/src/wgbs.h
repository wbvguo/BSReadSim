#ifndef HTSIM_WGBS_H
#define HTSIM_WGBS_H

#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <vector>
#include <string>
#include <optional>

#include "types.h"
#include "utilities.h"
#include "fragment.h"
#include "reference.h"
#include "variant.h"

// ---- uniform_sampler --------------------------------------------------------

namespace htsim::wgbs {

class UniformSamplingError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct FixedFragmentShape {
    std::uint32_t insert_length = 0;
    std::uint32_t read_length = 0;
    bool paired_end = false;
    double max_ambiguous_fraction = 0.0;
};

// Immutable rank-select index for one contig and one fixed fragment shape.
// Eligibility uses one bit per possible start plus a sparse prefix table; it
// does not retain the contig bases or a uint32_t value for every valid start.
// The index can serve arbitrarily many chunks without rescanning the contig.
class ValidStartIndex {
public:
    ValidStartIndex(
        const model::Bases &contig_bases,
        const FixedFragmentShape &shape);

    std::uint32_t possible_start_count() const noexcept;
    std::uint32_t valid_start_count() const noexcept;
    bool is_valid_start(std::uint32_t zero_based_start) const noexcept;
    std::uint32_t start_for_rank(std::uint32_t zero_based_rank) const;

    std::vector<std::uint32_t> sample(
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count) const;

private:
    std::vector<std::uint64_t> valid_words_;
    std::vector<std::uint32_t> superblock_prefix_;
    std::uint32_t possible_start_count_ = 0;
    std::uint32_t valid_start_count_ = 0;
};

// A start is eligible when a complete fixed insert fits and each emitted mate
// contains at most floor(max_ambiguous_fraction * read_length) N bases.  Bases
// in the unsequenced insert interior do not affect eligibility.
std::uint32_t count_valid_starts(
    const model::Bases &contig_bases,
    const FixedFragmentShape &shape);

// Draw output_count starts uniformly with replacement from the valid starts.
// RNG stage is Stage::fragment, keyed by contig_index. Each output uses its per-contig
// candidate ordinal as entity_ordinal and local_index 1 for the valid-start
// rank. Direct valid-rank sampling means every candidate is accepted, so the
// caller advances first_candidate_ordinal by output_count across chunks and
// reports zero skipped fragments for this baseline. Requested ranks are sorted
// internally and resolved against a rank-select index; returned starts remain
// in candidate order and are chunk independent. Global fragment ordinals are
// assigned later. This convenience function builds a temporary index; the core
// constructs ValidStartIndex once per contig and reuses it across chunks.
std::vector<std::uint32_t> sample_valid_starts(
    const model::Bases &contig_bases,
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count,
    const FixedFragmentShape &shape,
    std::uint32_t expected_valid_start_count);

} // namespace htsim::wgbs

// ---- coverage_profile --------------------------------------------------------

namespace htsim::wgbs {

inline constexpr std::string_view wgbs_gc_format = "tsv";
inline constexpr std::string_view wgbs_gc_version = "wgbs-gc-target-v1";
inline constexpr std::uint32_t maximum_attempts_per_fragment = 100000;

class CoverageProfileError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// A verified line-oriented WGBS GC target distribution. Every physical line
// is exactly one probability; its zero-based line index is its bin. The
// probabilities must sum to one. Empty lines, comments, and extra fields are
// invalid. Plain and gzip input are distinguished by their bytes and therefore
// by their required digest.
class WgbsGcProfile {
public:
    WgbsGcProfile(
        const std::string &path,
        const crypto::Sha256Digest &expected_file_sha256);

    std::uint32_t bin_count() const noexcept;
    std::uint32_t bin_for_counts(
        std::uint32_t gc_count,
        std::uint32_t fragment_length) const;
    double target_probability_for_counts(
        std::uint32_t gc_count,
        std::uint32_t fragment_length) const;
    const std::vector<double> &target_probabilities() const noexcept;
    const crypto::Sha256Digest &file_sha256() const noexcept;

private:
    crypto::Sha256Digest file_sha256_ = {};
    std::vector<double> probabilities_;
};

struct SampleBatch {
    std::vector<std::uint32_t> starts;
    std::uint64_t skipped_count = 0;
};

struct WgbsGcTargetCalibration {
    std::vector<double> acceptance_probabilities;
    std::vector<double> contig_allocation_weights;
};

// Calibrate a global target distribution against fixed-insert opportunity
// counts. If N_i is the global eligible-start count and p_i the requested
// output probability, acceptance_i is proportional to p_i/N_i. Contig
// allocation weights are sum_i N_ci*acceptance_i, which preserves p_i after
// mixing contigs. A positive target on an unreachable bin fails closed.
WgbsGcTargetCalibration calibrate_gc_target(
    const WgbsGcProfile &profile,
    const std::vector<std::vector<std::uint32_t>> &contig_bin_counts);

// Target-profile WGBS sampling over the same valid-start domain as uniform
// WGBS. Each output entity repeatedly draws a valid-start rank and an
// independent calibrated acceptance variate. Rejections never advance another
// entity's RNG address.
class WgbsGcSampler {
public:
    WgbsGcSampler(
        const model::Bases &contig_bases,
        const wgbs::FixedFragmentShape &shape,
        const WgbsGcProfile &profile);

    std::uint32_t valid_start_count() const noexcept;
    const std::vector<std::uint32_t> &bin_opportunity_counts() const noexcept;

    SampleBatch sample(
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count,
        const std::vector<double> &acceptance_probabilities) const;

private:
    std::uint32_t gc_count(std::uint32_t begin, std::uint32_t end) const;
    std::uint32_t bin(std::uint32_t start) const;

    wgbs::FixedFragmentShape shape_;
    wgbs::ValidStartIndex valid_starts_;
    WgbsGcProfile profile_;
    std::vector<std::uint64_t> gc_words_;
    std::vector<std::uint32_t> gc_prefix_;
    std::vector<std::uint32_t> bin_opportunity_counts_;
    std::uint32_t contig_length_ = 0;
};

} // namespace htsim::wgbs

// ---- variable_wgbs_sampler --------------------------------------------------------

namespace htsim::wgbs {

class VariableWgbsSamplingError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct VariableWgbsOptions {
    std::uint32_t maximum_attempts_per_fragment = 100000U;
};

struct VariableWgbsCandidate {
    std::uint32_t reference_start = 0;
    std::uint32_t insert_length = 0;
};

struct VariableWgbsBatch {
    std::vector<VariableWgbsCandidate> candidates;
    std::uint64_t skipped_count = 0;
    std::uint64_t next_candidate_ordinal = 0;
};

class VariableWgbsSampler {
public:
    VariableWgbsSampler(
        const model::Bases &contig_bases,
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        const insert_length::Parameters &insert_parameters,
        std::uint32_t read_length,
        bool paired_end,
        double max_ambiguous_fraction,
        const VariableWgbsOptions &options = {});

    std::uint32_t allocation_weight() const noexcept;

    // Evaluate exactly one addressed proposal. The insert draw uses local
    // index 0 and the start draw uses local index 1. An invalid span or mate
    // ambiguity returns nullopt; no hidden retry or ordinal increment occurs.
    std::optional<VariableWgbsCandidate> candidate_at(
        std::uint64_t candidate_ordinal) const;

    VariableWgbsBatch sample(
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count) const;

private:
    bool valid_candidate(
        std::uint32_t reference_start,
        std::uint32_t insert_length) const noexcept;

    const model::Bases *bases_ = nullptr;
    insert_length::Sampler insert_sampler_;
    std::uint64_t fragment_key_ = 0;
    std::uint64_t maximum_n_count_ = 0;
    std::uint32_t read_length_ = 0;
    std::uint32_t allocation_weight_ = 0;
    std::uint32_t maximum_attempts_ = 0;
    bool paired_end_ = false;
};

} // namespace htsim::wgbs

// ---- variable_haplotype_sampler --------------------------------------------------------

namespace htsim::wgbs {

class VariableHaplotypeSamplingError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct VariableHaplotypeOptions {
    std::uint32_t maximum_attempts_per_fragment = 100000U;
};

struct VariableHaplotypeCandidate {
    std::uint32_t reference_start = 0;
    std::uint32_t reference_span = 0;
    model::HaplotypeMask eligible_haplotypes =
        model::HaplotypeMask::both;
};

struct VariableHaplotypeBatch {
    std::vector<VariableHaplotypeCandidate> candidates;
    std::uint64_t skipped_count = 0;
    std::uint64_t next_candidate_ordinal = 0;
};

class VariableHaplotypeSampler {
public:
    VariableHaplotypeSampler(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        std::uint64_t master_seed,
        const insert_length::Parameters &insert_parameters,
        std::uint32_t read_length,
        bool paired_end,
        double max_ambiguous_fraction,
        const VariableHaplotypeOptions &options = {});

    std::uint32_t allocation_weight() const noexcept;

    // Return design-deck HaplotypeMask values 1, 2, or 3. nullopt is the
    // ineligible zero-bit state; it is never exposed as a HaplotypeMask.
    std::optional<model::HaplotypeMask> eligible_haplotypes(
        std::uint32_t reference_start,
        std::uint32_t reference_span) const;

    // Resolve exactly one addressed proposal. The candidate ordinal owns the
    // insert draw at fragment/local 0 and start draw at fragment/local 1.
    // nullopt means that proposal is ineligible on both haplotypes.
    std::optional<VariableHaplotypeCandidate> candidate_at(
        std::uint64_t candidate_ordinal) const;

    VariableHaplotypeBatch sample(
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count) const;

private:
    const reference::Contig *contig_ = nullptr;
    const variant::ContigVariants *variants_ = nullptr;
    insert_length::Sampler insert_sampler_;
    insert_length::Parameters insert_parameters_;
    std::uint64_t fragment_key_ = 0;
    std::uint32_t maximum_n_count_ = 0;
    std::uint32_t read_length_ = 0;
    std::uint32_t allocation_weight_ = 0;
    std::uint32_t maximum_attempts_ = 0;
    bool paired_end_ = false;
};

} // namespace htsim::wgbs

// ---- haplotype_start_index --------------------------------------------------------

namespace htsim::wgbs {

class VariantStartIndexError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Compact rank-select index over original-reference starts. Every start owns
// two HaplotypeMask bits: bit 0 is zero-based haplotype 0 and bit 1 is
// haplotype 1.
// A zero pair is ineligible and is never returned by rank sampling.
class HaplotypeStartIndex {
public:
    HaplotypeStartIndex(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        const FixedFragmentShape &shape);

    std::uint32_t possible_start_count() const noexcept;
    std::uint32_t valid_start_count() const noexcept;
    bool is_valid_start(std::uint32_t zero_based_start) const noexcept;
    model::HaplotypeMask haplotype_mask(
        std::uint32_t zero_based_start) const;
    std::uint32_t start_for_rank(std::uint32_t zero_based_rank) const;

    // Identical fragment-stage RNG address to ValidStartIndex: per-contig
    // candidate ordinal with local_index=1. Haplotype choice remains a later
    // global-fragment-ordinal draw constrained by haplotype_mask(start).
    std::vector<std::uint32_t> sample(
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        std::uint64_t first_candidate_ordinal,
        std::uint32_t output_count) const;

private:
    std::uint8_t mask_bits(std::uint32_t zero_based_start) const noexcept;
    void set_eligible(
        std::uint32_t zero_based_start,
        std::uint8_t zero_based_haplotype);
    void build_rank_index();

    // Two adjacent bits per reference start, 32 starts per uint64 word.
    std::vector<std::uint64_t> haplotype_words_;
    std::vector<std::uint32_t> superblock_prefix_;
    std::uint32_t possible_start_count_ = 0;
    std::uint32_t valid_start_count_ = 0;
};

} // namespace htsim::wgbs

#endif // HTSIM_WGBS_H
