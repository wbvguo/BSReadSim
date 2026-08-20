#include "wgbs.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <charconv>
#include <string>
#include <cstdint>
#include <optional>
#include <utility>
#include <array>

#include "utilities.h"
#include "variant.h"

// ---- uniform_sampler --------------------------------------------------------

namespace htsim::wgbs {
namespace {

struct ValidatedShape {
    std::size_t insert_length;
    std::size_t read_length;
    std::uint64_t maximum_n_count;
    bool paired_end;
};

ValidatedShape validate_shape(
    const model::Bases &bases,
    const FixedFragmentShape &shape)
{
    if (bases.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw UniformSamplingError("contig length exceeds uint32");
    }
    if (shape.insert_length == 0 || shape.read_length == 0) {
        throw UniformSamplingError("insert and read lengths must be positive");
    }
    if (shape.read_length > shape.insert_length) {
        throw UniformSamplingError("read length must not exceed insert length");
    }
    if (!std::isfinite(shape.max_ambiguous_fraction)
        || shape.max_ambiguous_fraction < 0.0
        || shape.max_ambiguous_fraction > 1.0) {
        throw UniformSamplingError(
            "maximum ambiguous fraction must be finite and in [0, 1]");
    }
    for (const std::uint8_t base : bases) {
        if (base > 4) {
            throw UniformSamplingError(
                "contig contains a base outside protocol encoding");
        }
    }
    const double maximum_n = std::floor(
        shape.max_ambiguous_fraction
        * static_cast<double>(shape.read_length));
    return {
        static_cast<std::size_t>(shape.insert_length),
        static_cast<std::size_t>(shape.read_length),
        static_cast<std::uint64_t>(maximum_n),
        shape.paired_end,
    };
}

std::uint64_t n_count(
    const model::Bases &bases,
    std::size_t start,
    std::size_t length) noexcept
{
    return static_cast<std::uint64_t>(std::count(
        bases.begin() + static_cast<std::ptrdiff_t>(start),
        bases.begin() + static_cast<std::ptrdiff_t>(start + length),
        static_cast<std::uint8_t>(4)));
}

class ValidStartCursor {
public:
    ValidStartCursor(
        const model::Bases &bases,
        const ValidatedShape &shape)
        : bases_(&bases), shape_(shape)
    {
        if (shape_.insert_length > bases.size()) {return;}
        has_candidate_ = true;
        first_n_ = n_count(bases, 0, shape_.read_length);
        if (shape_.paired_end) {
            second_n_ = n_count(
                bases,
                shape_.insert_length - shape_.read_length,
                shape_.read_length);
        }
    }

    bool next(std::uint64_t &start)
    {
        while (has_candidate_) {
            const std::size_t candidate = position_;
            const bool valid = first_n_ <= shape_.maximum_n_count
                && (!shape_.paired_end
                    || second_n_ <= shape_.maximum_n_count);
            advance();
            if (valid) {
                start = static_cast<std::uint64_t>(candidate);
                return true;
            }
        }
        return false;
    }

private:
    void advance()
    {
        const std::size_t last = bases_->size() - shape_.insert_length;
        if (position_ == last) {
            has_candidate_ = false;
            return;
        }
        if ((*bases_)[position_] == 4) {--first_n_;}
        if ((*bases_)[position_ + shape_.read_length] == 4) {++first_n_;}
        if (shape_.paired_end) {
            const std::size_t second_start =
                position_ + shape_.insert_length - shape_.read_length;
            if ((*bases_)[second_start] == 4) {--second_n_;}
            if ((*bases_)[second_start + shape_.read_length] == 4) {++second_n_;}
        }
        ++position_;
    }

    const model::Bases *bases_;
    ValidatedShape shape_;
    std::size_t position_ = 0;
    std::uint64_t first_n_ = 0;
    std::uint64_t second_n_ = 0;
    bool has_candidate_ = false;
};

constexpr std::size_t words_per_superblock = 16;

std::uint32_t population_count(std::uint64_t value) noexcept
{
    std::uint32_t count = 0;
    while (value != 0) {
        value &= value - 1;
        ++count;
    }
    return count;
}

std::uint32_t least_set_bit_index(std::uint64_t value) noexcept
{
    std::uint32_t index = 0;
    while ((value & UINT64_C(1)) == 0) {
        value >>= 1U;
        ++index;
    }
    return index;
}

} // namespace

ValidStartIndex::ValidStartIndex(
    const model::Bases &contig_bases,
    const FixedFragmentShape &shape)
{
    const ValidatedShape validated = validate_shape(contig_bases, shape);
    if (validated.insert_length > contig_bases.size()) {return;}
    const std::size_t possible =
        contig_bases.size() - validated.insert_length + 1U;
    if (possible > std::numeric_limits<std::uint32_t>::max()) {
        throw UniformSamplingError("possible-start count exceeds uint32");
    }
    possible_start_count_ = static_cast<std::uint32_t>(possible);
    const std::size_t word_count = (possible + 63U) / 64U;
    valid_words_.assign(word_count, 0);

    ValidStartCursor cursor(contig_bases, validated);
    std::uint64_t start = 0;
    while (cursor.next(start)) {
        if (start >= possible) {
            throw UniformSamplingError("internal valid-start position overflow");
        }
        const std::size_t word = static_cast<std::size_t>(start / 64U);
        const unsigned bit = static_cast<unsigned>(start % 64U);
        valid_words_[word] |= UINT64_C(1) << bit;
        if (valid_start_count_ == std::numeric_limits<std::uint32_t>::max()) {
            throw UniformSamplingError("valid-start count exceeds uint32");
        }
        ++valid_start_count_;
    }

    const std::size_t superblock_count =
        (valid_words_.size() + words_per_superblock - 1U)
        / words_per_superblock;
    superblock_prefix_.reserve(superblock_count + 1U);
    std::uint32_t prefix = 0;
    for (std::size_t block = 0; block < superblock_count; ++block) {
        superblock_prefix_.push_back(prefix);
        const std::size_t begin = block * words_per_superblock;
        const std::size_t end = std::min(
            begin + words_per_superblock, valid_words_.size());
        for (std::size_t word = begin; word < end; ++word) {
            const std::uint32_t count = population_count(valid_words_[word]);
            if (count > std::numeric_limits<std::uint32_t>::max() - prefix) {
                throw UniformSamplingError("valid-start prefix exceeds uint32");
            }
            prefix += count;
        }
    }
    superblock_prefix_.push_back(prefix);
    if (prefix != valid_start_count_) {
        throw UniformSamplingError("valid-start rank index is inconsistent");
    }
}

std::uint32_t ValidStartIndex::possible_start_count() const noexcept
{
    return possible_start_count_;
}

std::uint32_t ValidStartIndex::valid_start_count() const noexcept
{
    return valid_start_count_;
}

bool ValidStartIndex::is_valid_start(
    std::uint32_t zero_based_start) const noexcept
{
    if (zero_based_start >= possible_start_count_) {return false;}
    const std::size_t word = static_cast<std::size_t>(zero_based_start / 64U);
    const unsigned bit = static_cast<unsigned>(zero_based_start % 64U);
    return (valid_words_[word] & (UINT64_C(1) << bit)) != 0U;
}

std::uint32_t ValidStartIndex::start_for_rank(
    std::uint32_t zero_based_rank) const
{
    if (zero_based_rank >= valid_start_count_) {
        throw UniformSamplingError("valid-start rank is out of range");
    }
    const auto upper = std::upper_bound(
        superblock_prefix_.begin(), superblock_prefix_.end(), zero_based_rank);
    if (upper == superblock_prefix_.begin()) {
        throw UniformSamplingError("internal valid-start prefix underflow");
    }
    const std::size_t block = static_cast<std::size_t>(
        (upper - superblock_prefix_.begin()) - 1);
    std::uint32_t observed = superblock_prefix_[block];
    const std::size_t begin = block * words_per_superblock;
    const std::size_t end = std::min(
        begin + words_per_superblock, valid_words_.size());
    for (std::size_t word_index = begin; word_index < end; ++word_index) {
        std::uint64_t word = valid_words_[word_index];
        const std::uint32_t count = population_count(word);
        if (zero_based_rank < observed + count) {
            std::uint32_t within_word = zero_based_rank - observed;
            while (within_word != 0) {
                word &= word - 1;
                --within_word;
            }
            const std::uint64_t start =
                static_cast<std::uint64_t>(word_index) * 64U
                + least_set_bit_index(word);
            if (start >= possible_start_count_) {
                throw UniformSamplingError("rank index resolved outside the contig");
            }
            return static_cast<std::uint32_t>(start);
        }
        observed += count;
    }
    throw UniformSamplingError("valid-start rank was not resolved");
}

std::vector<std::uint32_t> ValidStartIndex::sample(
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count) const
{
    if (output_count == 0) {return {};}
    if (valid_start_count_ == 0) {
        throw UniformSamplingError("cannot sample from zero valid starts");
    }
    if (static_cast<std::uint64_t>(output_count) - 1U
        > std::numeric_limits<std::uint64_t>::max() - first_candidate_ordinal) {
        throw UniformSamplingError("candidate ordinal range exceeds uint64");
    }
    const std::uint64_t key =
        rng::derive_key(master_seed, rng::Stage::fragment, contig_index);
    std::vector<std::uint32_t> starts;
    starts.reserve(output_count);
    for (std::uint32_t index = 0; index < output_count; ++index) {
        const std::uint64_t candidate = first_candidate_ordinal + index;
        const std::uint64_t rank = rng::bounded_integer(
            key, candidate, UINT64_C(1), valid_start_count_);
        starts.push_back(start_for_rank(static_cast<std::uint32_t>(rank)));
    }
    return starts;
}

std::uint32_t count_valid_starts(
    const model::Bases &contig_bases,
    const FixedFragmentShape &shape)
{
    const ValidatedShape validated = validate_shape(contig_bases, shape);
    ValidStartCursor cursor(contig_bases, validated);
    std::uint32_t count = 0;
    std::uint64_t ignored_start = 0;
    while (cursor.next(ignored_start)) {
        if (count == std::numeric_limits<std::uint32_t>::max()) {
            throw UniformSamplingError("valid-start count exceeds uint32");
        }
        ++count;
    }
    return count;
}

std::vector<std::uint32_t> sample_valid_starts(
    const model::Bases &contig_bases,
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count,
    const FixedFragmentShape &shape,
    std::uint32_t expected_valid_start_count)
{
    const ValidStartIndex index(contig_bases, shape);
    if (index.valid_start_count() != expected_valid_start_count) {
        throw UniformSamplingError(
            "contig changed or valid-start count disagrees with its plan");
    }
    return index.sample(
        contig_index, master_seed, first_candidate_ordinal, output_count);
}

} // namespace htsim::wgbs

// ---- coverage_profile --------------------------------------------------------

namespace htsim::wgbs {
namespace {

double parse_probability(std::string_view text)
{
    if (text.empty()) {
        throw CoverageProfileError("coverage probability is empty");
    }
    double value = 0.0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()
        || !std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw CoverageProfileError(
            "coverage probability must be finite and in [0, 1]");
    }
    return value;
}

} // namespace

WgbsGcProfile::WgbsGcProfile(
    const std::string &path,
    const crypto::Sha256Digest &expected_file_sha256)
    : file_sha256_(expected_file_sha256)
{
    try {
        text::TextSnapshot snapshot(path, expected_file_sha256);
        snapshot.visit_lines([&](std::string_view line, std::uint64_t line_number) {
            try {
                if (probabilities_.size()
                    == std::numeric_limits<std::uint32_t>::max()) {
                    throw CoverageProfileError("coverage bin count exceeds uint32");
                }
                probabilities_.push_back(parse_probability(line));
            } catch (const CoverageProfileError &error) {
                throw CoverageProfileError(
                    "coverage profile line " + std::to_string(line_number)
                    + ": " + error.what());
            }
        });
        if (probabilities_.size() < 2U) {
            throw CoverageProfileError(
                "coverage profile requires at least two GC bins");
        }
        double sum = 0.0;
        double compensation = 0.0;
        for (const double probability : probabilities_) {
            const double corrected = probability - compensation;
            const double next = sum + corrected;
            compensation = (next - sum) - corrected;
            sum = next;
        }
        if (!std::isfinite(sum) || std::abs(sum - 1.0) > 1e-9) {
            throw CoverageProfileError(
                "target profile probabilities must sum to one");
        }
        if (snapshot.file_sha256() != file_sha256_) {
            throw CoverageProfileError("coverage profile digest changed");
        }
    } catch (const CoverageProfileError &) {
        throw;
    } catch (const std::exception &error) {
        throw CoverageProfileError(error.what());
    }
}

std::uint32_t WgbsGcProfile::bin_count() const noexcept
{
    return static_cast<std::uint32_t>(probabilities_.size());
}

std::uint32_t WgbsGcProfile::bin_for_counts(
    std::uint32_t gc_count,
    std::uint32_t fragment_length) const
{
    if (fragment_length == 0U || gc_count > fragment_length) {
        throw CoverageProfileError("GC count is outside the fragment length");
    }
    const std::uint64_t last_bin = bin_count() - 1U;
    // Both factors are uint32_t, so their product fits uint64_t exactly. This
    // is round-half-up of gc_count/fragment_length onto [0,last_bin] without a
    // floating-point boundary.
    const std::uint64_t product =
        static_cast<std::uint64_t>(gc_count) * last_bin;
    std::uint64_t result = product / fragment_length;
    const std::uint64_t remainder = product % fragment_length;
    if (UINT64_C(2) * remainder >= fragment_length) {++result;}
    if (result > last_bin) {
        throw CoverageProfileError("GC bin mapping exceeded the profile");
    }
    return static_cast<std::uint32_t>(result);
}

double WgbsGcProfile::target_probability_for_counts(
    std::uint32_t gc_count,
    std::uint32_t fragment_length) const
{
    return probabilities_.at(bin_for_counts(gc_count, fragment_length));
}

const std::vector<double> &WgbsGcProfile::target_probabilities() const noexcept
{
    return probabilities_;
}

const crypto::Sha256Digest &WgbsGcProfile::file_sha256() const noexcept
{
    return file_sha256_;
}

GcRankIndex::GcRankIndex(const model::Bases &contig_bases)
{
    if (contig_bases.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw CoverageProfileError("contig length exceeds uint32");
    }
    length_ = static_cast<std::uint32_t>(contig_bases.size());
    gc_words_.assign((contig_bases.size() + 63U) / 64U, 0U);
    for (std::size_t index = 0U; index < contig_bases.size(); ++index) {
        if (contig_bases[index] == 1U || contig_bases[index] == 2U) {
            gc_words_[index / 64U] |= UINT64_C(1) << (index % 64U);
        }
    }
    gc_prefix_.reserve(gc_words_.size() + 1U);
    gc_prefix_.push_back(0U);
    for (const std::uint64_t word : gc_words_) {
        const std::uint32_t word_count = population_count(word);
        if (word_count > std::numeric_limits<std::uint32_t>::max()
                - gc_prefix_.back()) {
            throw CoverageProfileError("contig GC count exceeds uint32");
        }
        gc_prefix_.push_back(gc_prefix_.back() + word_count);
    }
}

std::uint32_t GcRankIndex::length() const noexcept
{
    return length_;
}

std::uint32_t GcRankIndex::count(
    std::uint32_t begin,
    std::uint32_t end) const
{
    if (begin > end || end > length_) {
        throw CoverageProfileError("GC interval is outside the contig");
    }
    const auto rank = [&](std::uint32_t offset) {
        const std::size_t word = static_cast<std::size_t>(offset / 64U);
        const unsigned bit = static_cast<unsigned>(offset % 64U);
        std::uint32_t result = gc_prefix_.at(word);
        if (bit != 0U) {
            const std::uint64_t mask = (UINT64_C(1) << bit) - 1U;
            result += population_count(gc_words_.at(word) & mask);
        }
        return result;
    };
    return rank(end) - rank(begin);
}

WgbsGcTargetCalibration calibrate_gc_target(
    const WgbsGcProfile &profile,
    const std::vector<std::vector<std::uint32_t>> &contig_bin_counts,
    UnreachableTargetPolicy unreachable_policy)
{
    if (contig_bin_counts.empty()) {
        throw CoverageProfileError(
            "target calibration requires at least one contig");
    }
    const std::size_t bins = profile.bin_count();
    std::vector<std::uint64_t> global_counts(bins, 0U);
    for (const auto &counts : contig_bin_counts) {
        if (counts.size() != bins) {
            throw CoverageProfileError(
                "target calibration bin count disagrees with its profile");
        }
        for (std::size_t bin = 0U; bin < bins; ++bin) {
            if (counts[bin]
                > std::numeric_limits<std::uint64_t>::max()
                    - global_counts[bin]) {
                throw CoverageProfileError(
                    "global GC opportunity count exceeds uint64");
            }
            global_counts[bin] += counts[bin];
        }
    }

    const auto &targets = profile.target_probabilities();
    if (unreachable_policy != UnreachableTargetPolicy::reject
        && unreachable_policy
            != UnreachableTargetPolicy::drop_and_renormalize) {
        throw CoverageProfileError(
            "target calibration unreachable-bin policy is invalid");
    }
    WgbsGcTargetCalibration result;
    std::vector<double> ratios(bins, 0.0);
    double maximum_ratio = 0.0;
    for (std::size_t bin = 0U; bin < bins; ++bin) {
        if (targets[bin] == 0.0) {continue;}
        if (global_counts[bin] == 0U) {
            if (unreachable_policy == UnreachableTargetPolicy::reject) {
                throw CoverageProfileError(
                    "target profile bin " + std::to_string(bin)
                    + " has positive mass but no eligible fragment start");
            }
            result.dropped_target_probability += targets[bin];
            continue;
        }
        ratios[bin] = targets[bin]
            / static_cast<double>(global_counts[bin]);
        if (!std::isfinite(ratios[bin]) || ratios[bin] <= 0.0) {
            throw CoverageProfileError(
                "target calibration ratio is outside binary64 range");
        }
        maximum_ratio = std::max(maximum_ratio, ratios[bin]);
    }
    if (maximum_ratio <= 0.0) {
        throw CoverageProfileError(
            "target profile has no reachable positive probability");
    }

    result.acceptance_probabilities.resize(bins, 0.0);
    for (std::size_t bin = 0U; bin < bins; ++bin) {
        result.acceptance_probabilities[bin] = ratios[bin] / maximum_ratio;
    }
    result.contig_allocation_weights.reserve(contig_bin_counts.size());
    bool has_positive_contig = false;
    for (const auto &counts : contig_bin_counts) {
        long double weight = 0.0L;
        for (std::size_t bin = 0U; bin < bins; ++bin) {
            weight += static_cast<long double>(counts[bin])
                * result.acceptance_probabilities[bin];
        }
        const double stored = static_cast<double>(weight);
        if (!std::isfinite(stored) || stored < 0.0) {
            throw CoverageProfileError(
                "target contig allocation weight is outside binary64 range");
        }
        has_positive_contig = has_positive_contig || stored > 0.0;
        result.contig_allocation_weights.push_back(stored);
    }
    if (!has_positive_contig) {
        throw CoverageProfileError(
            "target profile has no eligible contig");
    }
    return result;
}

WgbsGcSampler::WgbsGcSampler(
    const model::Bases &contig_bases,
    const wgbs::FixedFragmentShape &shape,
    const WgbsGcProfile &profile)
    : shape_(shape),
      valid_starts_(contig_bases, shape),
      profile_(profile),
      gc_index_(contig_bases)
{
    bin_opportunity_counts_.assign(profile_.bin_count(), 0U);
    for (std::uint32_t start = 0U;
         start < valid_starts_.possible_start_count();
         ++start) {
        if (!valid_starts_.is_valid_start(start)) {continue;}
        std::uint32_t &count = bin_opportunity_counts_.at(bin(start));
        if (count == std::numeric_limits<std::uint32_t>::max()) {
            throw CoverageProfileError(
                "contig GC-bin opportunity count exceeds uint32");
        }
        ++count;
    }
}

std::uint32_t WgbsGcSampler::valid_start_count() const noexcept
{
    return valid_starts_.valid_start_count();
}

const std::vector<std::uint32_t> &
WgbsGcSampler::bin_opportunity_counts() const noexcept
{
    return bin_opportunity_counts_;
}

std::uint32_t WgbsGcSampler::bin(std::uint32_t start) const
{
    const std::uint64_t end =
        static_cast<std::uint64_t>(start) + shape_.insert_length;
    if (end > gc_index_.length()) {
        throw CoverageProfileError("profiled fragment exceeds the contig");
    }
    return profile_.bin_for_counts(
        gc_index_.count(start, static_cast<std::uint32_t>(end)),
        shape_.insert_length);
}

SampleBatch WgbsGcSampler::sample(
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count,
    const std::vector<double> &acceptance_probabilities) const
{
    if (output_count == 0U) {return {};}
    if (acceptance_probabilities.size() != profile_.bin_count()) {
        throw CoverageProfileError(
            "calibrated acceptance bin count disagrees with its profile");
    }
    bool has_positive_opportunity = false;
    for (std::size_t bin_index = 0U;
         bin_index < acceptance_probabilities.size();
         ++bin_index) {
        const double probability = acceptance_probabilities[bin_index];
        if (!std::isfinite(probability)
            || probability < 0.0 || probability > 1.0) {
            throw CoverageProfileError(
                "calibrated acceptance probability is outside [0,1]");
        }
        has_positive_opportunity = has_positive_opportunity
            || (probability > 0.0
                && bin_opportunity_counts_[bin_index] > 0U);
    }
    if (!has_positive_opportunity) {
        throw CoverageProfileError(
            "cannot sample a contig without positive calibrated opportunity");
    }
    if (static_cast<std::uint64_t>(output_count) - 1U
        > std::numeric_limits<std::uint64_t>::max()
            - first_candidate_ordinal) {
        throw CoverageProfileError("coverage candidate ordinal exceeds uint64");
    }
    const std::uint64_t key =
        rng::derive_key(master_seed, rng::Stage::fragment, contig_index);
    SampleBatch result;
    result.starts.reserve(output_count);
    for (std::uint32_t output_index = 0;
         output_index < output_count;
         ++output_index) {
        const std::uint64_t ordinal = first_candidate_ordinal + output_index;
        bool accepted = false;
        for (std::uint32_t attempt = 0;
             attempt < maximum_attempts_per_fragment;
             ++attempt) {
            const std::uint64_t start_local = UINT64_C(1)
                + UINT64_C(2) * attempt;
            const std::uint64_t rank = rng::bounded_integer(
                key,
                ordinal,
                start_local,
                valid_starts_.valid_start_count());
            const std::uint32_t start = valid_starts_.start_for_rank(
                static_cast<std::uint32_t>(rank));
            const double draw = rng::uniform01(
                key, ordinal, start_local + 1U);
            if (draw < acceptance_probabilities.at(bin(start))) {
                result.starts.push_back(start);
                accepted = true;
                break;
            }
            if (result.skipped_count
                == std::numeric_limits<std::uint64_t>::max()) {
                throw CoverageProfileError("coverage skipped count exceeds uint64");
            }
            ++result.skipped_count;
        }
        if (!accepted) {
            throw CoverageProfileError(
                "WGBS GC sampler exhausted its attempt cap");
        }
    }
    return result;
}

} // namespace htsim::wgbs

// ---- variable_wgbs_sampler --------------------------------------------------------

namespace htsim::wgbs {
namespace {

std::uint64_t count_n(
    const model::Bases &bases,
    std::size_t start,
    std::size_t length) noexcept
{
    return static_cast<std::uint64_t>(std::count(
        bases.begin() + static_cast<std::ptrdiff_t>(start),
        bases.begin() + static_cast<std::ptrdiff_t>(start + length),
        static_cast<std::uint8_t>(4U)));
}

} // namespace

VariableWgbsSampler::VariableWgbsSampler(
    const model::Bases &contig_bases,
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    const insert_length::Parameters &insert_parameters,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction,
    const VariableWgbsOptions &options)
    : bases_(&contig_bases),
      insert_sampler_(master_seed, contig_index, insert_parameters),
      fragment_key_(rng::derive_key(
          master_seed, rng::Stage::fragment, contig_index)),
      read_length_(read_length),
      maximum_attempts_(options.maximum_attempts_per_fragment),
      paired_end_(paired_end)
{
    if (contig_bases.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw VariableWgbsSamplingError("WGBS contig length exceeds uint32");
    }
    if (read_length == 0U || read_length > insert_parameters.minimum) {
        throw VariableWgbsSamplingError(
            "WGBS read length must fit the minimum insert length");
    }
    if (!std::isfinite(max_ambiguous_fraction)
        || max_ambiguous_fraction < 0.0
        || max_ambiguous_fraction > 1.0) {
        throw VariableWgbsSamplingError(
            "maximum ambiguous fraction must be finite and in [0, 1]");
    }
    if (maximum_attempts_ == 0U) {
        throw VariableWgbsSamplingError(
            "maximum attempts per fragment must be positive");
    }
    for (const std::uint8_t base : contig_bases) {
        if (base > 4U) {
            throw VariableWgbsSamplingError(
                "WGBS contig contains a base outside protocol encoding");
        }
    }

    maximum_n_count_ = static_cast<std::uint64_t>(std::floor(
        max_ambiguous_fraction * static_cast<double>(read_length)));
    try {
        allocation_weight_ = count_valid_starts(
            contig_bases,
            {insert_parameters.maximum,
             read_length,
             paired_end,
             max_ambiguous_fraction});
    } catch (const UniformSamplingError &error) {
        throw VariableWgbsSamplingError(error.what());
    }
}

std::uint32_t VariableWgbsSampler::allocation_weight() const noexcept
{
    return allocation_weight_;
}

std::optional<VariableWgbsCandidate> VariableWgbsSampler::candidate_at(
    std::uint64_t candidate_ordinal) const
{
    const std::uint32_t insert = insert_sampler_.sample(candidate_ordinal);
    if (insert > bases_->size()) {return std::nullopt;}
    const std::uint64_t possible =
        static_cast<std::uint64_t>(bases_->size()) - insert + 1U;
    const std::uint64_t selected = rng::bounded_integer(
        fragment_key_, candidate_ordinal, UINT64_C(1), possible);
    const std::uint32_t start = static_cast<std::uint32_t>(selected);
    if (!valid_candidate(start, insert)) {return std::nullopt;}
    return VariableWgbsCandidate{start, insert};
}

VariableWgbsBatch VariableWgbsSampler::sample(
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count) const
{
    VariableWgbsBatch batch;
    batch.next_candidate_ordinal = first_candidate_ordinal;
    if (output_count == 0U) {return batch;}
    if (allocation_weight_ == 0U) {
        throw VariableWgbsSamplingError(
            "cannot sample variable inserts from zero valid maximum-span starts");
    }
    batch.candidates.reserve(output_count);

    std::uint64_t candidate_ordinal = first_candidate_ordinal;
    for (std::uint32_t output_index = 0U;
         output_index < output_count;
         ++output_index) {
        bool accepted = false;
        for (std::uint32_t attempt = 0U;
             attempt < maximum_attempts_;
             ++attempt) {
            const auto candidate = candidate_at(candidate_ordinal);
            if (candidate) {
                batch.candidates.push_back(*candidate);
                accepted = true;
            }
            if (candidate_ordinal
                == std::numeric_limits<std::uint64_t>::max()) {
                throw VariableWgbsSamplingError(
                    "variable-insert candidate ordinal exceeds uint64");
            }
            ++candidate_ordinal;
            if (accepted) {break;}
            if (batch.skipped_count
                == std::numeric_limits<std::uint64_t>::max()) {
                throw VariableWgbsSamplingError(
                    "variable-insert skipped count exceeds uint64");
            }
            ++batch.skipped_count;
        }
        if (!accepted) {
            throw VariableWgbsSamplingError(
                "variable-insert sampler exhausted its attempt cap");
        }
    }
    batch.next_candidate_ordinal = candidate_ordinal;
    return batch;
}

bool VariableWgbsSampler::valid_candidate(
    std::uint32_t reference_start,
    std::uint32_t insert_length) const noexcept
{
    const std::size_t start = reference_start;
    const std::size_t insert = insert_length;
    const std::size_t read = read_length_;
    if (read > insert || start > bases_->size()
        || insert > bases_->size() - start) {
        return false;
    }
    if (count_n(*bases_, start, read) > maximum_n_count_) {return false;}
    return !paired_end_
        || count_n(*bases_, start + insert - read, read) <= maximum_n_count_;
}

VariableWgbsGcSampler::VariableWgbsGcSampler(
    const model::Bases &contig_bases,
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    const insert_length::Parameters &insert_parameters,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction,
    const WgbsGcProfile &profile,
    const VariableWgbsOptions &options)
    : proposals_(
          contig_bases,
          contig_index,
          master_seed,
          insert_parameters,
          read_length,
          paired_end,
          max_ambiguous_fraction,
          options),
      profile_(profile),
      gc_index_(contig_bases),
      fragment_key_(rng::derive_key(
          master_seed, rng::Stage::fragment, contig_index)),
      maximum_attempts_(options.maximum_attempts_per_fragment)
{
}

std::uint32_t VariableWgbsGcSampler::allocation_weight() const noexcept
{
    return proposals_.allocation_weight();
}

std::uint32_t VariableWgbsGcSampler::bin(
    const VariableWgbsCandidate &candidate) const
{
    const std::uint64_t end =
        static_cast<std::uint64_t>(candidate.reference_start)
        + candidate.insert_length;
    if (end > gc_index_.length()) {
        throw CoverageProfileError("profiled fragment exceeds the contig");
    }
    return profile_.bin_for_counts(
        gc_index_.count(
            candidate.reference_start, static_cast<std::uint32_t>(end)),
        candidate.insert_length);
}

VariableWgbsBatch VariableWgbsGcSampler::sample(
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count,
    const std::vector<double> &acceptance_probabilities) const
{
    if (acceptance_probabilities.size() != profile_.bin_count()) {
        throw CoverageProfileError(
            "calibrated acceptance bin count disagrees with its profile");
    }
    bool has_positive_acceptance = false;
    for (const double probability : acceptance_probabilities) {
        if (!std::isfinite(probability)
            || probability < 0.0 || probability > 1.0) {
            throw CoverageProfileError(
                "calibrated acceptance probability is outside [0,1]");
        }
        has_positive_acceptance = has_positive_acceptance || probability > 0.0;
    }
    if (!has_positive_acceptance && output_count != 0U) {
        throw CoverageProfileError(
            "cannot sample without a positive GC acceptance probability");
    }

    VariableWgbsBatch result;
    result.next_candidate_ordinal = first_candidate_ordinal;
    result.candidates.reserve(output_count);
    std::uint64_t candidate_ordinal = first_candidate_ordinal;
    for (std::uint32_t output_index = 0U;
         output_index < output_count;
         ++output_index) {
        bool accepted = false;
        for (std::uint32_t attempt = 0U;
             attempt < maximum_attempts_;
             ++attempt) {
            const auto candidate = proposals_.candidate_at(candidate_ordinal);
            if (candidate
                && rng::uniform01(fragment_key_, candidate_ordinal, UINT64_C(2))
                    < acceptance_probabilities.at(bin(*candidate))) {
                result.candidates.push_back(*candidate);
                accepted = true;
            }
            if (candidate_ordinal == std::numeric_limits<std::uint64_t>::max()) {
                throw VariableWgbsSamplingError(
                    "profiled variable-insert candidate ordinal exceeds uint64");
            }
            ++candidate_ordinal;
            if (accepted) {break;}
            if (result.skipped_count
                == std::numeric_limits<std::uint64_t>::max()) {
                throw VariableWgbsSamplingError(
                    "profiled variable-insert skipped count exceeds uint64");
            }
            ++result.skipped_count;
        }
        if (!accepted) {
            throw VariableWgbsSamplingError(
                "profiled variable-insert sampler exhausted its attempt cap");
        }
    }
    result.next_candidate_ordinal = candidate_ordinal;
    return result;
}

} // namespace htsim::wgbs

// ---- variable_haplotype_sampler --------------------------------------------------------

namespace htsim::wgbs {
namespace {

std::uint32_t count_projected_n(
    const model::Bases &bases,
    std::size_t begin,
    std::size_t end)
{
    if (begin > end || end > bases.size()) {
        throw VariableHaplotypeSamplingError(
            "projected mate interval is outside its template");
    }
    return static_cast<std::uint32_t>(std::count(
        bases.begin() + static_cast<std::ptrdiff_t>(begin),
        bases.begin() + static_cast<std::ptrdiff_t>(end),
        static_cast<std::uint8_t>(4U)));
}

} // namespace

VariableHaplotypeSampler::VariableHaplotypeSampler(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint64_t master_seed,
    const insert_length::Parameters &insert_parameters,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction,
    const VariableHaplotypeOptions &options)
    : contig_(&contig),
      variants_(&variants),
      insert_sampler_(master_seed, contig.index, insert_parameters),
      insert_parameters_(insert_parameters),
      fragment_key_(rng::derive_key(
          master_seed, rng::Stage::fragment, contig.index)),
      read_length_(read_length),
      maximum_attempts_(options.maximum_attempts_per_fragment),
      paired_end_(paired_end)
{
    if (read_length == 0U || read_length > insert_parameters.minimum) {
        throw VariableHaplotypeSamplingError(
            "haplotype read length must fit the minimum insert length");
    }
    if (!std::isfinite(max_ambiguous_fraction)
        || max_ambiguous_fraction < 0.0
        || max_ambiguous_fraction > 1.0) {
        throw VariableHaplotypeSamplingError(
            "maximum ambiguous fraction must be finite and in [0, 1]");
    }
    if (maximum_attempts_ == 0U) {
        throw VariableHaplotypeSamplingError(
            "maximum attempts per fragment must be positive");
    }
    maximum_n_count_ = static_cast<std::uint32_t>(std::floor(
        max_ambiguous_fraction * static_cast<double>(read_length)));
    try {
        allocation_weight_ = HaplotypeStartIndex(
            contig,
            variants,
            {insert_parameters.maximum,
             read_length,
             paired_end,
             max_ambiguous_fraction}).valid_start_count();
    } catch (const VariantStartIndexError &error) {
        throw VariableHaplotypeSamplingError(error.what());
    }
}

std::uint32_t VariableHaplotypeSampler::allocation_weight() const noexcept
{
    return allocation_weight_;
}

std::optional<model::HaplotypeMask>
VariableHaplotypeSampler::eligible_haplotypes(
    std::uint32_t reference_start,
    std::uint32_t reference_span) const
{
    if (reference_span < insert_parameters_.minimum
        || reference_span > insert_parameters_.maximum
        || reference_start > contig_->length
        || reference_span > contig_->length - reference_start) {
        throw VariableHaplotypeSamplingError(
            "candidate reference interval is outside its configured domain");
    }
    const std::uint32_t reference_end = reference_start + reference_span;
    std::uint8_t mask = 0U;
    for (std::uint8_t haplotype = 0U; haplotype < 2U; ++haplotype) {
        try {
            const auto projection = haplotype::project_interval(
                *contig_,
                *variants_,
                haplotype,
                reference_start,
                reference_end);
            const std::size_t template_length = projection.template_bases.size();
            if (template_length < read_length_) {continue;}
            if (count_projected_n(
                    projection.template_bases, 0U, read_length_)
                > maximum_n_count_) {
                continue;
            }
            if (paired_end_
                && count_projected_n(
                       projection.template_bases,
                       template_length - read_length_,
                       template_length)
                    > maximum_n_count_) {
                continue;
            }
            mask |= static_cast<std::uint8_t>(1U << haplotype);
        } catch (const haplotype::ProjectionError &error) {
            if (error.failure() != haplotype::ProjectionFailure::boundary_cut
                && error.failure()
                    != haplotype::ProjectionFailure::empty_projection) {
                throw VariableHaplotypeSamplingError(error.what());
            }
            // A candidate that cuts an active deletion or is fully deleted is
            // ineligible for that haplotype, not a malformed typed catalog.
        }
    }
    if (mask == 0U) {return std::nullopt;}
    if (!model::is_haplotype_mask(mask)) {
        throw VariableHaplotypeSamplingError(
            "candidate produced an invalid HaplotypeMask");
    }
    return static_cast<model::HaplotypeMask>(mask);
}

std::optional<VariableHaplotypeCandidate>
VariableHaplotypeSampler::candidate_at(
    std::uint64_t candidate_ordinal) const
{
    const std::uint32_t span = insert_sampler_.sample(candidate_ordinal);
    if (span > contig_->length) {return std::nullopt;}
    const std::uint64_t possible = contig_->length - span + 1U;
    const std::uint64_t selected = rng::bounded_integer(
        fragment_key_, candidate_ordinal, UINT64_C(1), possible);
    const std::uint32_t start = static_cast<std::uint32_t>(selected);
    const auto mask = eligible_haplotypes(start, span);
    if (!mask) {return std::nullopt;}
    return VariableHaplotypeCandidate{start, span, *mask};
}

VariableHaplotypeBatch VariableHaplotypeSampler::sample(
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count) const
{
    VariableHaplotypeBatch batch;
    batch.next_candidate_ordinal = first_candidate_ordinal;
    if (output_count == 0U) {return batch;}
    if (allocation_weight_ == 0U) {
        throw VariableHaplotypeSamplingError(
            "cannot sample from zero valid maximum-span haplotype starts");
    }
    batch.candidates.reserve(output_count);

    std::uint64_t candidate_ordinal = first_candidate_ordinal;
    for (std::uint32_t output_index = 0U;
         output_index < output_count;
         ++output_index) {
        bool accepted = false;
        for (std::uint32_t attempt = 0U;
             attempt < maximum_attempts_;
             ++attempt) {
            const auto candidate = candidate_at(candidate_ordinal);
            if (candidate) {
                batch.candidates.push_back(*candidate);
                accepted = true;
            }
            if (candidate_ordinal
                == std::numeric_limits<std::uint64_t>::max()) {
                throw VariableHaplotypeSamplingError(
                    "haplotype candidate ordinal exceeds uint64");
            }
            ++candidate_ordinal;
            if (accepted) {break;}
            if (batch.skipped_count
                == std::numeric_limits<std::uint64_t>::max()) {
                throw VariableHaplotypeSamplingError(
                    "haplotype skipped count exceeds uint64");
            }
            ++batch.skipped_count;
        }
        if (!accepted) {
            throw VariableHaplotypeSamplingError(
                "variable haplotype sampler exhausted its attempt cap");
        }
    }
    batch.next_candidate_ordinal = candidate_ordinal;
    return batch;
}

} // namespace htsim::wgbs

// ---- haplotype_start_index --------------------------------------------------------

namespace htsim::wgbs {
namespace {

inline constexpr std::size_t starts_per_word = 32U;
inline constexpr std::uint64_t low_pair_bits = UINT64_C(0x5555555555555555);

std::uint64_t valid_start_bits(std::uint64_t haplotype_pairs) noexcept
{
    return (haplotype_pairs | (haplotype_pairs >> 1U)) & low_pair_bits;
}

struct ValidatedHaplotypeShape {
    std::uint32_t insert_length;
    std::uint32_t read_length;
    std::uint32_t maximum_n_count;
    bool paired_end;
};

ValidatedHaplotypeShape validate_inputs(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    const FixedFragmentShape &shape)
{
    if (contig.length != contig.bases.size()
        || contig.length > std::numeric_limits<std::uint32_t>::max()
        || variants.contig_index() != contig.index
        || variants.reference_length() != contig.length) {
        throw VariantStartIndexError(
            "variant start index inputs do not identify one uint32 contig");
    }
    if (shape.insert_length == 0U || shape.read_length == 0U
        || shape.read_length > shape.insert_length) {
        throw VariantStartIndexError(
            "variant start index requires positive fitting read/insert lengths");
    }
    if (!std::isfinite(shape.max_ambiguous_fraction)
        || shape.max_ambiguous_fraction < 0.0
        || shape.max_ambiguous_fraction > 1.0) {
        throw VariantStartIndexError(
            "maximum ambiguous fraction must be finite and in [0, 1]");
    }
    if (!std::all_of(
            contig.bases.begin(), contig.bases.end(),
            [](std::uint8_t base) {return base <= 4U;})) {
        throw VariantStartIndexError(
            "contig contains a base outside protocol encoding");
    }
    const double maximum_n = std::floor(
        shape.max_ambiguous_fraction
        * static_cast<double>(shape.read_length));
    return {
        shape.insert_length,
        shape.read_length,
        static_cast<std::uint32_t>(maximum_n),
        shape.paired_end,
    };
}

class NRankIndex {
public:
    void append(bool is_n)
    {
        if (length_ == std::numeric_limits<std::uint64_t>::max()) {
            throw VariantStartIndexError("haplotype length exceeds uint64");
        }
        const std::uint64_t word_index = length_ / 64U;
        if (word_index > std::numeric_limits<std::size_t>::max()) {
            throw VariantStartIndexError("haplotype bit index exceeds size_t");
        }
        if (static_cast<std::size_t>(word_index) == words_.size()) {
            words_.push_back(0U);
        }
        if (is_n) {
            words_.back() |= UINT64_C(1) << static_cast<unsigned>(length_ % 64U);
        }
        ++length_;
    }

    void finish()
    {
        prefix_.reserve(words_.size() + 1U);
        prefix_.push_back(0U);
        std::uint32_t count = 0U;
        for (const std::uint64_t word : words_) {
            const std::uint32_t increment = population_count(word);
            if (increment > std::numeric_limits<std::uint32_t>::max() - count) {
                throw VariantStartIndexError("haplotype N count exceeds uint32");
            }
            count += increment;
            prefix_.push_back(count);
        }
    }

    std::uint64_t length() const noexcept {return length_;}

    std::uint32_t count(std::uint64_t begin, std::uint64_t end) const
    {
        if (begin > end || end > length_ || prefix_.size() != words_.size() + 1U) {
            throw VariantStartIndexError("haplotype N-rank interval is invalid");
        }
        return rank(end) - rank(begin);
    }

private:
    std::uint32_t rank(std::uint64_t end) const
    {
        const std::uint64_t word_index = end / 64U;
        const unsigned offset = static_cast<unsigned>(end % 64U);
        if (word_index > words_.size()) {
            throw VariantStartIndexError("haplotype N-rank position is invalid");
        }
        const std::size_t index = static_cast<std::size_t>(word_index);
        std::uint32_t result = prefix_.at(index);
        if (offset != 0U) {
            if (index >= words_.size()) {
                throw VariantStartIndexError(
                    "haplotype N-rank partial word is absent");
            }
            const std::uint64_t mask = (UINT64_C(1) << offset) - 1U;
            result += population_count(words_[index] & mask);
        }
        return result;
    }

    std::vector<std::uint64_t> words_;
    std::vector<std::uint32_t> prefix_;
    std::uint64_t length_ = 0U;
};

void append_reference_n(
    NRankIndex &index,
    const model::Bases &reference,
    std::uint32_t begin,
    std::uint32_t end)
{
    if (begin > end || end > reference.size()) {
        throw VariantStartIndexError("haplotype reference cursor is invalid");
    }
    for (std::uint32_t position = begin; position < end; ++position) {
        index.append(reference[position] == 4U);
    }
}

NRankIndex build_n_rank(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint8_t haplotype)
{
    NRankIndex index;
    std::uint32_t cursor = 0U;
    for (const variant::Event &event : variants.events()) {
        if (!model::mask_contains(event.alt_haplotypes, haplotype)) {
            continue;
        }
        append_reference_n(index, contig.bases, cursor, event.reference_start);
        switch (event.kind) {
        case model::VariantKind::insertion:
            for (std::size_t offset = 0U; offset < event.alt_bases.size(); ++offset) {
                index.append(false);
            }
            cursor = event.reference_start;
            break;
        case model::VariantKind::snv:
            index.append(false);
            cursor = event.reference_end;
            break;
        case model::VariantKind::deletion:
            cursor = event.reference_end;
            break;
        default:
            throw VariantStartIndexError("variant kind is outside the typed catalog");
        }
    }
    append_reference_n(
        index,
        contig.bases,
        cursor,
        static_cast<std::uint32_t>(contig.length));
    index.finish();
    return index;
}

class BoundaryCursor {
public:
    BoundaryCursor(
        const variant::ContigVariants &variants,
        std::uint8_t haplotype)
        : events_(&variants.events()), haplotype_(haplotype)
    {}

    std::optional<std::uint64_t> before(std::uint32_t target)
    {
        if (target < previous_target_) {
            throw VariantStartIndexError(
                "haplotype boundary cursor moved backwards");
        }
        previous_target_ = target;
        while (reference_cursor_ < target) {
            if (event_index_ < events_->size()) {
                const variant::Event &event = events_->at(event_index_);
                if (event.reference_start < reference_cursor_) {
                    throw VariantStartIndexError(
                        "variant boundary cursor passed an event");
                }
                if (event.reference_start == reference_cursor_) {
                    ++event_index_;
                    if (!model::mask_contains(
                            event.alt_haplotypes, haplotype_)) {
                        continue;
                    }
                    if (event.kind == model::VariantKind::insertion) {
                        add_haplotype_length(event.alt_bases.size());
                        continue;
                    }
                    if (event.kind == model::VariantKind::snv) {
                        add_haplotype_length(event.alt_bases.size());
                        reference_cursor_ = event.reference_end;
                        continue;
                    }
                    if (event.kind == model::VariantKind::deletion) {
                        reference_cursor_ = event.reference_end;
                        continue;
                    }
                    throw VariantStartIndexError(
                        "variant kind is outside the typed catalog");
                }
            }
            ++reference_cursor_;
            add_haplotype_length(1U);
        }
        if (reference_cursor_ > target) {return std::nullopt;}
        return haplotype_cursor_;
    }

private:
    void add_haplotype_length(std::size_t increment)
    {
        if (increment
            > std::numeric_limits<std::uint64_t>::max() - haplotype_cursor_) {
            throw VariantStartIndexError("haplotype boundary exceeds uint64");
        }
        haplotype_cursor_ += increment;
    }

    const std::vector<variant::Event> *events_;
    std::uint8_t haplotype_ = 0U;
    std::uint32_t reference_cursor_ = 0U;
    std::uint32_t previous_target_ = 0U;
    std::uint64_t haplotype_cursor_ = 0U;
    std::size_t event_index_ = 0U;
};

std::uint32_t terminal_insertion_length(
    const variant::ContigVariants &variants,
    std::uint8_t haplotype)
{
    if (variants.events().empty()) {return 0U;}
    const variant::Event &event = variants.events().back();
    if (event.kind != model::VariantKind::insertion
        || event.reference_start != variants.reference_length()
        || !model::mask_contains(event.alt_haplotypes, haplotype)) {
        return 0U;
    }
    return static_cast<std::uint32_t>(event.alt_bases.size());
}

} // namespace

HaplotypeStartIndex::HaplotypeStartIndex(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    const FixedFragmentShape &shape)
{
    const ValidatedHaplotypeShape validated =
        validate_inputs(contig, variants, shape);
    if (validated.insert_length > contig.length) {return;}
    const std::uint64_t possible =
        contig.length - validated.insert_length + 1U;
    if (possible > std::numeric_limits<std::uint32_t>::max()) {
        throw VariantStartIndexError("possible start count exceeds uint32");
    }
    possible_start_count_ = static_cast<std::uint32_t>(possible);
    const std::uint64_t pair_bits = possible * 2U;
    const std::uint64_t word_count = (pair_bits + 63U) / 64U;
    if (word_count > std::numeric_limits<std::size_t>::max()) {
        throw VariantStartIndexError("haplotype start bitset exceeds size_t");
    }
    haplotype_words_.assign(static_cast<std::size_t>(word_count), 0U);

    for (std::uint8_t haplotype = 0U; haplotype < 2U; ++haplotype) {
        const NRankIndex n_rank = build_n_rank(contig, variants, haplotype);
        BoundaryCursor starts(variants, haplotype);
        BoundaryCursor ends(variants, haplotype);
        const std::uint64_t terminal =
            terminal_insertion_length(variants, haplotype);
        for (std::uint32_t start = 0U;
             start < possible_start_count_;
             ++start) {
            const std::uint32_t end = start + validated.insert_length;
            const auto projected_start = starts.before(start);
            auto projected_end = ends.before(end);
            if (!projected_start || !projected_end) {continue;}
            if (end == variants.reference_length()) {
                *projected_end += terminal;
            }
            if (*projected_end < *projected_start
                || *projected_end - *projected_start
                    < validated.read_length) {
                continue;
            }
            const std::uint64_t first_end =
                *projected_start + validated.read_length;
            if (n_rank.count(*projected_start, first_end)
                > validated.maximum_n_count) {
                continue;
            }
            if (validated.paired_end
                && n_rank.count(
                       *projected_end - validated.read_length,
                       *projected_end)
                    > validated.maximum_n_count) {
                continue;
            }
            set_eligible(start, haplotype);
        }
        const auto final_boundary = ends.before(variants.reference_length());
        if (!final_boundary
            || terminal
                > std::numeric_limits<std::uint64_t>::max() - *final_boundary
            || n_rank.length() != *final_boundary + terminal) {
            throw VariantStartIndexError(
                "haplotype boundary and N-rank lengths disagree");
        }
    }
    build_rank_index();
}

std::uint32_t HaplotypeStartIndex::possible_start_count() const noexcept
{
    return possible_start_count_;
}

std::uint32_t HaplotypeStartIndex::valid_start_count() const noexcept
{
    return valid_start_count_;
}

std::uint8_t HaplotypeStartIndex::mask_bits(
    std::uint32_t zero_based_start) const noexcept
{
    if (zero_based_start >= possible_start_count_) {return 0U;}
    const std::uint64_t bit_index =
        static_cast<std::uint64_t>(zero_based_start) * 2U;
    const std::size_t word = static_cast<std::size_t>(bit_index / 64U);
    const unsigned shift = static_cast<unsigned>(bit_index % 64U);
    return static_cast<std::uint8_t>(
        (haplotype_words_[word] >> shift) & UINT64_C(3));
}

bool HaplotypeStartIndex::is_valid_start(
    std::uint32_t zero_based_start) const noexcept
{
    return mask_bits(zero_based_start) != 0U;
}

model::HaplotypeMask HaplotypeStartIndex::haplotype_mask(
    std::uint32_t zero_based_start) const
{
    const std::uint8_t mask = mask_bits(zero_based_start);
    if (!model::is_haplotype_mask(mask)) {
        throw VariantStartIndexError(
            "start is out of range or has no eligible haplotype");
    }
    return static_cast<model::HaplotypeMask>(mask);
}

void HaplotypeStartIndex::set_eligible(
    std::uint32_t zero_based_start,
    std::uint8_t zero_based_haplotype)
{
    if (zero_based_start >= possible_start_count_
        || zero_based_haplotype > 1U) {
        throw VariantStartIndexError("eligible start bit address is invalid");
    }
    const std::uint64_t bit_index =
        static_cast<std::uint64_t>(zero_based_start) * 2U
        + zero_based_haplotype;
    const std::size_t word = static_cast<std::size_t>(bit_index / 64U);
    const unsigned shift = static_cast<unsigned>(bit_index % 64U);
    haplotype_words_[word] |= UINT64_C(1) << shift;
}

void HaplotypeStartIndex::build_rank_index()
{
    const std::size_t superblock_count =
        (haplotype_words_.size() + words_per_superblock - 1U)
        / words_per_superblock;
    superblock_prefix_.reserve(superblock_count + 1U);
    std::uint32_t prefix = 0U;
    for (std::size_t block = 0U; block < superblock_count; ++block) {
        superblock_prefix_.push_back(prefix);
        const std::size_t begin = block * words_per_superblock;
        const std::size_t end = std::min(
            begin + words_per_superblock, haplotype_words_.size());
        for (std::size_t word = begin; word < end; ++word) {
            const std::uint32_t increment =
                population_count(valid_start_bits(haplotype_words_[word]));
            if (increment > std::numeric_limits<std::uint32_t>::max() - prefix) {
                throw VariantStartIndexError("valid start count exceeds uint32");
            }
            prefix += increment;
        }
    }
    superblock_prefix_.push_back(prefix);
    valid_start_count_ = prefix;
}

std::uint32_t HaplotypeStartIndex::start_for_rank(
    std::uint32_t zero_based_rank) const
{
    if (zero_based_rank >= valid_start_count_) {
        throw VariantStartIndexError("valid start rank is out of range");
    }
    const auto upper = std::upper_bound(
        superblock_prefix_.begin(), superblock_prefix_.end(), zero_based_rank);
    if (upper == superblock_prefix_.begin()) {
        throw VariantStartIndexError("valid start prefix underflow");
    }
    const std::size_t block = static_cast<std::size_t>(
        (upper - superblock_prefix_.begin()) - 1);
    std::uint32_t observed = superblock_prefix_[block];
    const std::size_t begin = block * words_per_superblock;
    const std::size_t end = std::min(
        begin + words_per_superblock, haplotype_words_.size());
    for (std::size_t word_index = begin; word_index < end; ++word_index) {
        std::uint64_t word = valid_start_bits(haplotype_words_[word_index]);
        const std::uint32_t count = population_count(word);
        if (zero_based_rank < observed + count) {
            std::uint32_t within_word = zero_based_rank - observed;
            while (within_word != 0U) {
                word &= word - 1U;
                --within_word;
            }
            const std::uint64_t start =
                static_cast<std::uint64_t>(word_index) * starts_per_word
                + least_set_bit_index(word) / 2U;
            if (start >= possible_start_count_) {
                throw VariantStartIndexError(
                    "valid start rank resolved outside the contig");
            }
            return static_cast<std::uint32_t>(start);
        }
        observed += count;
    }
    throw VariantStartIndexError("valid start rank was not resolved");
}

std::vector<std::uint32_t> HaplotypeStartIndex::sample(
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count) const
{
    if (output_count == 0U) {return {};}
    if (valid_start_count_ == 0U) {
        throw VariantStartIndexError("cannot sample from zero valid starts");
    }
    if (static_cast<std::uint64_t>(output_count) - 1U
        > std::numeric_limits<std::uint64_t>::max() - first_candidate_ordinal) {
        throw VariantStartIndexError("candidate ordinal range exceeds uint64");
    }
    const std::uint64_t key =
        rng::derive_key(master_seed, rng::Stage::fragment, contig_index);
    std::vector<std::uint32_t> starts;
    starts.reserve(output_count);
    for (std::uint32_t index = 0U; index < output_count; ++index) {
        const std::uint64_t rank = rng::bounded_integer(
            key,
            first_candidate_ordinal + index,
            UINT64_C(1),
            valid_start_count_);
        starts.push_back(start_for_rank(static_cast<std::uint32_t>(rank)));
    }
    return starts;
}

} // namespace htsim::wgbs
