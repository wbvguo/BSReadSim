#include "tbs.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <tuple>
#include <unordered_map>

#include "utilities.h"

// ---- catalog --------------------------------------------------------

namespace htsim::tbs {
namespace {

bool metadata_line(std::string_view line) noexcept
{
    return line.empty() || line.front() == '#'
        || line == "track" || line.substr(0, 6) == "track "
        || line == "browser" || line.substr(0, 8) == "browser ";
}

bool target_less(const Target &left, const Target &right) noexcept
{
    return std::tie(
               left.interval_start,
               left.interval_end,
               left.capture_strand,
               left.name,
               left.score)
        < std::tie(
               right.interval_start,
               right.interval_end,
               right.capture_strand,
               right.name,
               right.score);
}

std::vector<std::string_view> split_bed6(std::string_view line)
{
    std::vector<std::string_view> fields;
    fields.reserve(6);
    std::size_t begin = 0;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(
            begin,
            end == std::string_view::npos ? line.size() - begin : end - begin));
        if (end == std::string_view::npos) {break;}
        begin = end + 1U;
    }
    if (fields.size() != 6U) {
        throw TbsCatalogError("TBS target row must contain exactly six BED fields");
    }
    return fields;
}

std::uint32_t parse_coordinate(std::string_view text, const char *field)
{
    if (text.empty()) {
        throw TbsCatalogError(std::string("TBS ") + field + " is empty");
    }
    std::uint32_t value = 0;
    const auto parsed = std::from_chars(text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()) {
        throw TbsCatalogError(
            std::string("TBS ") + field + " is not a uint32 decimal");
    }
    return value;
}

double parse_score(std::string_view text)
{
    if (text.empty()) {throw TbsCatalogError("TBS score is empty");}
    double value = 0.0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()
        || !std::isfinite(value) || value < 0.0) {
        throw TbsCatalogError("TBS score must be a finite non-negative number");
    }
    return value;
}

model::CaptureStrand parse_strand(std::string_view text)
{
    if (text == "+") {return model::CaptureStrand::forward;}
    if (text == "-") {return model::CaptureStrand::reverse;}
    if (text == ".") {return model::CaptureStrand::unknown;}
    throw TbsCatalogError("TBS strand must be +, -, or .");
}

void validate_name(std::string_view name)
{
    if (name.empty()) {throw TbsCatalogError("TBS target name is empty");}
    for (const char value : name) {
        const auto byte = static_cast<unsigned char>(value);
        if (byte == 0U || byte < 0x20U || byte == 0x7fU) {
            throw TbsCatalogError("TBS target name contains a control byte");
        }
    }
}

std::uint32_t maximum_ambiguous_count(
    std::uint32_t read_length,
    double max_ambiguous_fraction)
{
    if (!std::isfinite(max_ambiguous_fraction)
        || max_ambiguous_fraction < 0.0
        || max_ambiguous_fraction > 1.0) {
        throw TbsCatalogError(
            "maximum ambiguous fraction must be finite and in [0, 1]");
    }
    return static_cast<std::uint32_t>(std::floor(
        max_ambiguous_fraction * static_cast<double>(read_length)));
}

std::uint32_t population_count(std::uint64_t value) noexcept
{
    std::uint32_t count = 0;
    while (value != 0U) {
        value &= value - 1U;
        ++count;
    }
    return count;
}

std::uint32_t exact_output_weight(double score)
{
    if (!std::isfinite(score) || score < 0.0
        || score > static_cast<double>(
            std::numeric_limits<std::uint32_t>::max())
        || std::trunc(score) != score) {
        throw TbsCatalogError(
            "TBS output weights must be exact uint32 integers");
    }
    return static_cast<std::uint32_t>(score);
}

} // namespace

TargetFile::TargetFile(
    const std::string &path,
    const std::vector<reference::ContigMetadata> &reference_catalog)
    : targets_by_contig_(reference_catalog.size())
{
    try {
        std::unordered_map<std::string, std::uint32_t> contig_indices;
        contig_indices.reserve(reference_catalog.size());
        for (std::size_t index = 0; index < reference_catalog.size(); ++index) {
            if (index > std::numeric_limits<std::uint32_t>::max()) {
                throw TbsCatalogError("reference contig index exceeds uint32");
            }
            if (!contig_indices.emplace(
                    reference_catalog[index].name,
                    static_cast<std::uint32_t>(index)).second) {
                throw TbsCatalogError("reference catalog contains duplicate names");
            }
        }

        text::TextSnapshot snapshot(path);
        file_sha256_ = snapshot.file_sha256();
        snapshot.visit_lines([&](std::string_view line, std::uint64_t line_number) {
            if (metadata_line(line)) {return;}
            try {
                const auto fields = split_bed6(line);
                const auto found = contig_indices.find(std::string(fields[0]));
                if (found == contig_indices.end()) {
                    throw TbsCatalogError("TBS target names an unknown contig");
                }
                const std::uint32_t start = parse_coordinate(fields[1], "start");
                const std::uint32_t end = parse_coordinate(fields[2], "end");
                if (start >= end) {
                    throw TbsCatalogError(
                        "TBS target must satisfy 0 <= start < end");
                }
                const auto &metadata = reference_catalog[found->second];
                if (end > metadata.length) {
                    throw TbsCatalogError("TBS target exceeds its reference contig");
                }
                validate_name(fields[3]);
                const double score = parse_score(fields[4]);
                const model::CaptureStrand strand = parse_strand(fields[5]);
                auto &contig_targets = targets_by_contig_[found->second];
                if (contig_targets.size()
                    == std::numeric_limits<std::uint32_t>::max()) {
                    throw TbsCatalogError(
                        "TBS target count for one contig exceeds uint32");
                }
                if (target_count_ == std::numeric_limits<std::uint64_t>::max()) {
                    throw TbsCatalogError("TBS target count exceeds uint64");
                }
                contig_targets.push_back(Target{
                    found->second,
                    start,
                    end,
                    std::string(fields[3]),
                    score,
                    strand,
                });
                ++target_count_;
            } catch (const TbsCatalogError &error) {
                throw TbsCatalogError(
                    "TBS BED line " + std::to_string(line_number)
                    + ": " + error.what());
            }
        });

        if (target_count_ == 0U) {
            throw TbsCatalogError("TBS BED contains no target rows");
        }
        for (auto &targets : targets_by_contig_) {
            std::sort(targets.begin(), targets.end(), target_less);
            for (std::size_t index = 1; index < targets.size(); ++index) {
                const Target &previous = targets[index - 1U];
                const Target &current = targets[index];
                if (previous.interval_start == current.interval_start
                    && previous.interval_end == current.interval_end
                    && previous.capture_strand == current.capture_strand) {
                    throw TbsCatalogError(
                        "TBS BED contains duplicate coordinate+strand targets");
                }
            }
        }
    } catch (const TbsCatalogError &) {
        throw;
    } catch (const std::exception &error) {
        throw TbsCatalogError(error.what());
    }
}

const std::vector<Target> &TargetFile::targets(std::uint32_t contig_index) const
{
    if (contig_index >= targets_by_contig_.size()) {
        throw TbsCatalogError("TBS contig index is out of range");
    }
    return targets_by_contig_[contig_index];
}

CandidateCatalog::CandidateCatalog(
    const model::Bases &contig_bases,
    const std::vector<Target> &targets,
    std::uint32_t contig_index,
    double center_sd,
    std::uint32_t insert_length,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction,
    SamplingMode sampling_mode)
    : CandidateCatalog(
          contig_bases,
          targets,
          contig_index,
          center_sd,
          {insert_length, insert_length, insert_length, 0.0},
          read_length,
          paired_end,
          max_ambiguous_fraction,
          sampling_mode)
{}

CandidateCatalog::CandidateCatalog(
    const model::Bases &contig_bases,
    const std::vector<Target> &targets,
    std::uint32_t contig_index,
    double center_sd,
    const insert_length::Parameters &insert_parameters,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction,
    SamplingMode sampling_mode)
    : insert_parameters_(insert_parameters)
{
    if (contig_bases.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw TbsCatalogError("contig length exceeds uint32");
    }
    if (!std::isfinite(center_sd) || center_sd < 0.0) {
        throw TbsCatalogError("TBS center SD is invalid");
    }
    if (insert_parameters_.minimum == 0U
        || insert_parameters_.mean < insert_parameters_.minimum
        || insert_parameters_.maximum < insert_parameters_.mean
        || !std::isfinite(insert_parameters_.standard_deviation)
        || insert_parameters_.standard_deviation < 0.0) {
        throw TbsCatalogError("TBS insert-length parameters are invalid");
    }
    variable_insert_ = insert_parameters_.standard_deviation != 0.0
        && insert_parameters_.minimum != insert_parameters_.maximum;
    const std::uint32_t read_boundary = variable_insert_
        ? insert_parameters_.minimum
        : insert_parameters_.mean;
    if (read_length == 0U || read_length > read_boundary) {
        throw TbsCatalogError("TBS insert/read lengths are invalid");
    }
    for (const std::uint8_t base : contig_bases) {
        if (base > 4U) {
            throw TbsCatalogError(
                "contig contains a base outside protocol encoding");
        }
    }
    if (sampling_mode != SamplingMode::uniform
        && sampling_mode != SamplingMode::output_weight) {
        throw TbsCatalogError("TBS sampling mode is invalid");
    }
    center_sd_ = center_sd;
    sampling_mode_ = sampling_mode;
    read_length_ = read_length;
    paired_end_ = paired_end;
    maximum_ambiguous_count_ = maximum_ambiguous_count(
        read_length, max_ambiguous_fraction);
    contig_length_ = static_cast<std::uint32_t>(contig_bases.size());

    ambiguous_words_.assign((contig_bases.size() + 63U) / 64U, 0U);
    for (std::size_t index = 0; index < contig_bases.size(); ++index) {
        if (contig_bases[index] == 4U) {
            ambiguous_words_[index / 64U] |=
                UINT64_C(1) << (index % 64U);
        }
    }
    ambiguous_prefix_.reserve(ambiguous_words_.size() + 1U);
    ambiguous_prefix_.push_back(0U);
    for (const std::uint64_t word : ambiguous_words_) {
        const std::uint32_t count = population_count(word);
        if (count > std::numeric_limits<std::uint32_t>::max()
                - ambiguous_prefix_.back()) {
            throw TbsCatalogError("ambiguous-base prefix exceeds uint32");
        }
        ambiguous_prefix_.push_back(ambiguous_prefix_.back() + count);
    }

    bool has_sequenceable_start = false;
    const std::uint32_t calibration_insert = insert_parameters_.mean;
    if (center_sd_ > 0.0 && calibration_insert <= contig_length_) {
        const std::uint32_t last = contig_length_ - calibration_insert;
        for (std::uint32_t start = 0;; ++start) {
            if (sequenceable(
                    start,
                    start + calibration_insert,
                    calibration_insert)) {
                has_sequenceable_start = true;
                break;
            }
            if (start == last) {break;}
        }
    }

    std::vector<std::uint32_t> choice_weights;
    choice_weights.reserve(targets.size());
    for (std::size_t index = 0; index < targets.size(); ++index) {
        const Target &target = targets[index];
        if (target.contig_index != contig_index) {
            throw TbsCatalogError("TBS target belongs to a different contig");
        }
        if (index != 0U && target_less(target, targets[index - 1U])) {
            throw TbsCatalogError("TBS targets are not in canonical order");
        }
        if (index != 0U) {
            const Target &previous = targets[index - 1U];
            if (previous.interval_start == target.interval_start
                && previous.interval_end == target.interval_end
                && previous.capture_strand == target.capture_strand) {
                throw TbsCatalogError("TBS target coordinates are duplicated");
            }
        }
        if (target.interval_start >= target.interval_end
            || target.interval_end > contig_length_) {
            throw TbsCatalogError("TBS target interval is outside the contig");
        }
        if (target.capture_strand != model::CaptureStrand::unknown
            && target.capture_strand != model::CaptureStrand::forward
            && target.capture_strand != model::CaptureStrand::reverse) {
            throw TbsCatalogError("TBS target has an invalid capture strand");
        }
        if (!std::isfinite(target.score) || target.score < 0.0) {
            throw TbsCatalogError("TBS target has an invalid score");
        }
        const std::uint32_t output_weight =
            sampling_mode_ == SamplingMode::output_weight
            ? exact_output_weight(target.score)
            : 0U;
        if (index >= std::numeric_limits<std::uint32_t>::max()) {
            throw TbsCatalogError("TBS target count exceeds uint32");
        }
        TargetChoice choice{
            target.interval_start,
            target.interval_end,
            static_cast<std::uint32_t>(index),
            target.score,
            target.capture_strand,
        };
        if (center_sd_ == 0.0) {
            Candidate candidate;
            if (!project(choice, 0, calibration_insert, candidate)) {continue;}
            if (variable_insert_) {
                target_choices_.push_back(std::move(choice));
            } else {
                fixed_candidates_.push_back(candidate);
            }
            if (sampling_mode_ == SamplingMode::output_weight) {
                choice_weights.push_back(output_weight);
            }
        } else if (has_sequenceable_start) {
            // A displaced center can land at any complete fragment start. If
            // no sequenceable start exists on this contig, no target can ever
            // be accepted and the contig must carry zero allocation weight.
            target_choices_.push_back(std::move(choice));
            if (sampling_mode_ == SamplingMode::output_weight) {
                choice_weights.push_back(output_weight);
            }
        }
    }
    if (sampling_mode_ == SamplingMode::output_weight) {
        cumulative_weights_.reserve(choice_weights.size());
        for (const std::uint32_t weight : choice_weights) {
            if (weight > std::numeric_limits<std::uint32_t>::max()
                    - total_weight_) {
                throw TbsCatalogError(
                    "TBS eligible target weight sum exceeds uint32");
            }
            total_weight_ += weight;
            cumulative_weights_.push_back(total_weight_);
        }
    }
}

std::uint32_t CandidateCatalog::choice_count() const noexcept
{
    return static_cast<std::uint32_t>(
        center_sd_ == 0.0 && !variable_insert_
            ? fixed_candidates_.size()
            : target_choices_.size());
}

std::uint32_t CandidateCatalog::allocation_weight() const noexcept
{
    return sampling_mode_ == SamplingMode::output_weight
        ? total_weight_
        : choice_count();
}

std::uint32_t CandidateCatalog::ambiguous_count(
    std::uint32_t begin,
    std::uint32_t end) const
{
    if (begin > end || end > contig_length_) {
        throw TbsCatalogError("TBS mate slice is outside the contig");
    }
    const auto rank = [&](std::uint32_t offset) {
        const std::size_t word = static_cast<std::size_t>(offset / 64U);
        const unsigned bit = static_cast<unsigned>(offset % 64U);
        std::uint32_t result = ambiguous_prefix_.at(word);
        if (bit != 0U) {
            const std::uint64_t mask = (UINT64_C(1) << bit) - 1U;
            result += population_count(ambiguous_words_.at(word) & mask);
        }
        return result;
    };
    return rank(end) - rank(begin);
}

bool CandidateCatalog::sequenceable(
    std::uint32_t fragment_start,
    std::uint32_t fragment_end,
    std::uint32_t insert_length) const
{
    if (fragment_start > fragment_end
        || fragment_end > contig_length_
        || fragment_end - fragment_start != insert_length
        || insert_length < read_length_) {
        return false;
    }
    if (ambiguous_count(
            fragment_start, fragment_start + read_length_)
        > maximum_ambiguous_count_) {
        return false;
    }
    return !paired_end_ || ambiguous_count(
        fragment_end - read_length_, fragment_end)
        <= maximum_ambiguous_count_;
}

bool CandidateCatalog::project(
    const TargetChoice &target,
    std::int64_t center_displacement,
    std::uint32_t insert_length,
    Candidate &candidate) const
{
    const std::uint32_t base_center = target.target_start
        + (target.target_end - target.target_start) / 2U;
    const std::int64_t contig_limit = static_cast<std::int64_t>(contig_length_);
    if (center_displacement < -static_cast<std::int64_t>(base_center)
        || center_displacement > contig_limit) {
        return false;
    }
    const std::int64_t center = static_cast<std::int64_t>(base_center)
        + center_displacement;
    const std::int64_t fragment_start = center
        - static_cast<std::int64_t>(insert_length / 2U);
    if (fragment_start < 0) {return false;}
    const std::uint64_t start = static_cast<std::uint64_t>(fragment_start);
    const std::uint64_t end = start + insert_length;
    if (end > contig_length_) {return false;}
    const auto start_u32 = static_cast<std::uint32_t>(start);
    const auto end_u32 = static_cast<std::uint32_t>(end);
    if (!sequenceable(start_u32, end_u32, insert_length)) {return false;}
    candidate = {
        start_u32,
        end_u32,
        insert_length,
        target.target_start,
        target.target_end,
        target.target_ordinal,
        target.output_weight,
        target.capture_strand,
        model::HaplotypeMask::both,
        true,
        false,
    };
    return true;
}

SampleBatch CandidateCatalog::sample(
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count) const
{
    if (output_count == 0U) {return {};}
    if (choice_count() == 0U) {
        throw TbsCatalogError("cannot sample from zero TBS candidates");
    }
    const std::uint64_t key = rng::derive_key(
        master_seed, rng::Stage::fragment, contig_index);
    const insert_length::Sampler insert_sampler(
        master_seed, contig_index, insert_parameters_);
    SampleBatch result;
    result.candidates.reserve(output_count);
    if (static_cast<std::uint64_t>(output_count) - 1U
        > std::numeric_limits<std::uint64_t>::max()
            - first_candidate_ordinal) {
        throw TbsCatalogError("TBS candidate ordinal exceeds uint64");
    }
    for (std::uint32_t output_index = 0;
         output_index < output_count;
         ++output_index) {
        const std::uint64_t ordinal = first_candidate_ordinal + output_index;
        std::uint32_t selected = 0U;
        if (sampling_mode_ == SamplingMode::uniform) {
            selected = static_cast<std::uint32_t>(rng::bounded_integer(
                key, ordinal, UINT64_C(1), choice_count()));
        } else {
            if (total_weight_ == 0U
                || cumulative_weights_.size() != choice_count()) {
                throw TbsCatalogError(
                    "cannot sample from zero TBS output weight");
            }
            const auto draw = static_cast<std::uint32_t>(
                rng::bounded_integer(
                    key, ordinal, UINT64_C(1), total_weight_));
            const auto found = std::upper_bound(
                cumulative_weights_.begin(), cumulative_weights_.end(), draw);
            if (found == cumulative_weights_.end()) {
                throw TbsCatalogError(
                    "TBS output-weight selection exceeded its cumulative weight");
            }
            selected = static_cast<std::uint32_t>(
                found - cumulative_weights_.begin());
        }
        Candidate candidate;
        if (center_sd_ == 0.0 && !variable_insert_) {
            candidate = fixed_candidates_.at(
                static_cast<std::size_t>(selected));
        } else {
            const TargetChoice &target = target_choices_.at(
                static_cast<std::size_t>(selected));
            bool accepted = false;
            for (std::uint32_t attempt = 0;
                 attempt < maximum_attempts_per_fragment;
                 ++attempt) {
                std::uint32_t insert_length = insert_parameters_.mean;
                if (variable_insert_) {
                    const std::uint64_t length_local = center_sd_ == 0.0
                        ? UINT64_C(2) + attempt
                        : UINT64_C(2) + UINT64_C(2) * attempt;
                    insert_length = insert_sampler.sample(
                        ordinal, length_local);
                }
                std::int64_t center_displacement = 0;
                if (center_sd_ > 0.0) {
                    const std::uint64_t center_local = variable_insert_
                        ? UINT64_C(3) + UINT64_C(2) * attempt
                        : UINT64_C(2) + attempt;
                    const double normal = normal_sampler::standard_normal(
                        key, ordinal, center_local);
                    const double displacement = center_sd_ * normal;
                    const long double extended_displacement = displacement;
                    if (!std::isfinite(displacement)
                        || extended_displacement
                            < static_cast<long double>(INT64_MIN)
                        || extended_displacement
                            > static_cast<long double>(INT64_MAX)) {
                        throw TbsCatalogError(
                            "TBS center displacement exceeds int64");
                    }
                    // Truncate the continuous normal displacement toward zero
                    // before applying it to the BED center.
                    center_displacement = static_cast<std::int64_t>(
                        displacement);
                }
                if (project(
                        target,
                        center_displacement,
                        insert_length,
                        candidate)) {
                    accepted = true;
                    break;
                }
                if (result.skipped_count
                    == std::numeric_limits<std::uint64_t>::max()) {
                    throw TbsCatalogError("TBS skipped count exceeds uint64");
                }
                ++result.skipped_count;
            }
            if (!accepted) {
                throw TbsCatalogError(
                    "TBS fragment sampler exhausted its attempt cap");
            }
        }
        result.candidates.push_back(candidate);
    }
    return result;
}

DiploidCandidateCatalog::DiploidCandidateCatalog(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    const std::vector<Target> &targets,
    double center_sd,
    std::uint32_t insert_length,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction,
    SamplingMode sampling_mode)
    : DiploidCandidateCatalog(
          contig,
          variants,
          targets,
          center_sd,
          {insert_length, insert_length, insert_length, 0.0},
          read_length,
          paired_end,
          max_ambiguous_fraction,
          sampling_mode)
{}

DiploidCandidateCatalog::DiploidCandidateCatalog(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    const std::vector<Target> &targets,
    double center_sd,
    const insert_length::Parameters &insert_parameters,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction,
    SamplingMode sampling_mode)
    : insert_parameters_(insert_parameters),
      center_sd_(center_sd),
      read_length_(read_length),
      maximum_ambiguous_count_(maximum_ambiguous_count(
          read_length, max_ambiguous_fraction)),
      paired_end_(paired_end),
      variable_insert_(insert_parameters.standard_deviation != 0.0
          && insert_parameters.minimum != insert_parameters.maximum),
      sampling_mode_(sampling_mode)
{
    if (!std::isfinite(center_sd_) || center_sd_ < 0.0) {
        throw TbsCatalogError("TBS center SD is invalid");
    }
    if (insert_parameters_.minimum == 0U
        || insert_parameters_.mean < insert_parameters_.minimum
        || insert_parameters_.maximum < insert_parameters_.mean
        || !std::isfinite(insert_parameters_.standard_deviation)
        || insert_parameters_.standard_deviation < 0.0) {
        throw TbsCatalogError("TBS insert-length parameters are invalid");
    }
    const std::uint32_t read_boundary = variable_insert_
        ? insert_parameters_.minimum
        : insert_parameters_.mean;
    if (read_length_ == 0U || read_length_ > read_boundary) {
        throw TbsCatalogError("TBS insert/read lengths are invalid");
    }
    if (sampling_mode_ != SamplingMode::uniform
        && sampling_mode_ != SamplingMode::output_weight) {
        throw TbsCatalogError("TBS sampling mode is invalid");
    }
    for (std::uint8_t haplotype_index = 0U;
         haplotype_index < 2U;
         ++haplotype_index) {
        layouts_[haplotype_index] =
            std::make_unique<haplotype::HaplotypeLayout>(
                contig, variants, haplotype_index, false);
    }

    std::array<bool, 2> has_sequenceable_start = {false, false};
    const std::uint32_t calibration_insert = insert_parameters_.mean;
    if (center_sd_ > 0.0) {
        for (std::uint8_t haplotype_index = 0U;
             haplotype_index < 2U;
            ++haplotype_index) {
            const auto &layout = *layouts_[haplotype_index];
            if (calibration_insert > layout.length()) {continue;}
            const std::uint32_t last = layout.length() - calibration_insert;
            for (std::uint32_t start = 0U;; ++start) {
                if (sequenceable(
                        layout,
                        start,
                        start + calibration_insert,
                        calibration_insert)
                    && layout.boundary(start)
                    && layout.boundary(start + calibration_insert)) {
                    has_sequenceable_start[haplotype_index] = true;
                    break;
                }
                if (start == last) {break;}
            }
        }
    }

    std::vector<std::uint32_t> choice_weights;
    choice_weights.reserve(targets.size());
    target_choices_.reserve(targets.size());
    if (targets.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw TbsCatalogError("TBS target count exceeds uint32");
    }
    for (std::size_t index = 0U; index < targets.size(); ++index) {
        const Target &target = targets[index];
        if (target.contig_index != contig.index) {
            throw TbsCatalogError("TBS target belongs to a different contig");
        }
        if (index != 0U && target_less(target, targets[index - 1U])) {
            throw TbsCatalogError("TBS targets are not in canonical order");
        }
        if (index != 0U) {
            const Target &previous = targets[index - 1U];
            if (previous.interval_start == target.interval_start
                && previous.interval_end == target.interval_end
                && previous.capture_strand == target.capture_strand) {
                throw TbsCatalogError("TBS target coordinates are duplicated");
            }
        }
        if (target.interval_start >= target.interval_end
            || target.interval_end > contig.length) {
            throw TbsCatalogError("TBS target interval is outside the contig");
        }
        if (target.capture_strand != model::CaptureStrand::unknown
            && target.capture_strand != model::CaptureStrand::forward
            && target.capture_strand != model::CaptureStrand::reverse) {
            throw TbsCatalogError("TBS target has an invalid capture strand");
        }
        if (!std::isfinite(target.score) || target.score < 0.0) {
            throw TbsCatalogError("TBS target has an invalid score");
        }
        const std::uint32_t output_weight =
            sampling_mode_ == SamplingMode::output_weight
            ? exact_output_weight(target.score)
            : 0U;
        const std::uint32_t reference_center = target.interval_start
            + (target.interval_end - target.interval_start) / 2U;
        TargetChoice choice{
            target.interval_start,
            target.interval_end,
            static_cast<std::uint32_t>(index),
            {},
            target.score,
            target.capture_strand,
            model::HaplotypeMask::both,
        };
        std::uint8_t eligible_bits = 0U;
        for (std::uint8_t haplotype_index = 0U;
             haplotype_index < 2U;
             ++haplotype_index) {
            const auto center = layouts_[haplotype_index]
                ->boundary_before_reference(reference_center);
            if (!center) {continue;}
            choice.haplotype_centers[haplotype_index] = *center;
            bool accepted = false;
            if (center_sd_ == 0.0) {
                Candidate candidate;
                if (project(
                        choice,
                        haplotype_index,
                        0,
                        calibration_insert,
                        candidate)) {
                    accepted = true;
                }
            } else if (has_sequenceable_start[haplotype_index]) {
                accepted = true;
            }
            if (accepted) {
                eligible_bits = static_cast<std::uint8_t>(
                    eligible_bits | (UINT8_C(1) << haplotype_index));
            }
        }
        if (eligible_bits == 0U) {continue;}
        if (!model::is_haplotype_mask(eligible_bits)) {
            throw TbsCatalogError(
                "diploid TBS target produced an invalid haplotype mask");
        }
        choice.eligible_haplotypes =
            static_cast<model::HaplotypeMask>(eligible_bits);
        if (target_choices_.size()
            == std::numeric_limits<std::uint32_t>::max()) {
            throw TbsCatalogError(
                "diploid TBS target choice count exceeds uint32");
        }
        target_choices_.push_back(std::move(choice));
        if (sampling_mode_ == SamplingMode::output_weight) {
            // A BED score is one target-level output weight. It is never
            // multiplied by the number of eligible haplotypes.
            choice_weights.push_back(output_weight);
        }
    }

    if (sampling_mode_ == SamplingMode::output_weight) {
        cumulative_weights_.reserve(choice_weights.size());
        for (const std::uint32_t weight : choice_weights) {
            if (weight > std::numeric_limits<std::uint32_t>::max()
                    - total_weight_) {
                throw TbsCatalogError(
                    "TBS eligible target weight sum exceeds uint32");
            }
            total_weight_ += weight;
            cumulative_weights_.push_back(total_weight_);
        }
    }
}

std::uint32_t DiploidCandidateCatalog::choice_count() const noexcept
{
    return static_cast<std::uint32_t>(target_choices_.size());
}

std::uint32_t DiploidCandidateCatalog::allocation_weight() const noexcept
{
    return sampling_mode_ == SamplingMode::output_weight
        ? total_weight_
        : choice_count();
}

bool DiploidCandidateCatalog::sequenceable(
    const haplotype::HaplotypeLayout &layout,
    std::uint32_t fragment_start,
    std::uint32_t fragment_end,
    std::uint32_t insert_length) const
{
    if (fragment_start > fragment_end
        || fragment_end > layout.length()
        || fragment_end - fragment_start != insert_length
        || insert_length < read_length_) {
        return false;
    }
    if (layout.ambiguous_count(
            fragment_start, fragment_start + read_length_)
        > maximum_ambiguous_count_) {
        return false;
    }
    return !paired_end_ || layout.ambiguous_count(
        fragment_end - read_length_, fragment_end)
        <= maximum_ambiguous_count_;
}

bool DiploidCandidateCatalog::project(
    const TargetChoice &target,
    std::uint8_t haplotype,
    std::int64_t center_displacement,
    std::uint32_t insert_length,
    Candidate &candidate) const
{
    if (haplotype > 1U
        || !model::is_haplotype_mask(
            static_cast<std::uint8_t>(target.eligible_haplotypes))
        || !model::mask_contains(target.eligible_haplotypes, haplotype)) {
        throw TbsCatalogError("TBS target haplotype is invalid");
    }
    const auto &layout = *layouts_[haplotype];
    const std::uint32_t haplotype_center =
        target.haplotype_centers[haplotype];
    if (center_displacement
            < -static_cast<std::int64_t>(haplotype_center)
        || center_displacement
            > static_cast<std::int64_t>(layout.length())) {
        return false;
    }
    const std::int64_t center =
        static_cast<std::int64_t>(haplotype_center)
        + center_displacement;
    const std::int64_t fragment_start = center
        - static_cast<std::int64_t>(insert_length / 2U);
    if (fragment_start < 0) {return false;}
    const std::uint64_t start = static_cast<std::uint64_t>(fragment_start);
    const std::uint64_t end = start + insert_length;
    if (end > layout.length()) {return false;}
    const auto start_u32 = static_cast<std::uint32_t>(start);
    const auto end_u32 = static_cast<std::uint32_t>(end);
    if (!sequenceable(layout, start_u32, end_u32, insert_length)) {
        return false;
    }
    const auto left = layout.boundary(start_u32);
    const auto right = layout.boundary(end_u32);
    if (!left || !right
        || left->right_reference_start > right->left_reference_end) {
        return false;
    }
    candidate = {
        left->right_reference_start,
        right->left_reference_end,
        insert_length,
        target.target_start,
        target.target_end,
        target.target_ordinal,
        target.output_weight,
        target.capture_strand,
        haplotype == 0U
            ? model::HaplotypeMask::haplotype_1
            : model::HaplotypeMask::haplotype_2,
        left->include_insertion_in_right_fragment,
        right->include_insertion_in_left_fragment,
    };
    return true;
}

SampleBatch DiploidCandidateCatalog::sample(
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint64_t first_fragment_ordinal,
    std::uint32_t output_count) const
{
    if (output_count == 0U) {return {};}
    if (choice_count() == 0U) {
        throw TbsCatalogError("cannot sample from zero TBS candidates");
    }
    if (static_cast<std::uint64_t>(output_count) - 1U
        > std::numeric_limits<std::uint64_t>::max()
            - first_candidate_ordinal) {
        throw TbsCatalogError("TBS candidate ordinal exceeds uint64");
    }
    if (static_cast<std::uint64_t>(output_count) - 1U
        > std::numeric_limits<std::uint64_t>::max()
            - first_fragment_ordinal) {
        throw TbsCatalogError("TBS fragment ordinal exceeds uint64");
    }
    const std::uint64_t fragment_key = rng::derive_key(
        master_seed, rng::Stage::fragment, contig_index);
    const insert_length::Sampler insert_sampler(
        master_seed, contig_index, insert_parameters_);
    const std::uint64_t haplotype_key = rng::derive_key(
        master_seed, rng::Stage::haplotype, contig_index);
    SampleBatch result;
    result.candidates.reserve(output_count);
    for (std::uint32_t output_index = 0U;
         output_index < output_count;
         ++output_index) {
        const std::uint64_t ordinal = first_candidate_ordinal + output_index;
        const std::uint64_t fragment_ordinal =
            first_fragment_ordinal + output_index;
        std::uint32_t selected = 0U;
        if (sampling_mode_ == SamplingMode::uniform) {
            selected = static_cast<std::uint32_t>(rng::bounded_integer(
                fragment_key, ordinal, UINT64_C(1), choice_count()));
        } else {
            if (total_weight_ == 0U
                || cumulative_weights_.size() != choice_count()) {
                throw TbsCatalogError(
                    "cannot sample from zero TBS output weight");
            }
            const auto draw = static_cast<std::uint32_t>(
                rng::bounded_integer(
                    fragment_key, ordinal, UINT64_C(1), total_weight_));
            const auto found = std::upper_bound(
                cumulative_weights_.begin(), cumulative_weights_.end(), draw);
            if (found == cumulative_weights_.end()) {
                throw TbsCatalogError(
                    "TBS output-weight selection exceeded its cumulative weight");
            }
            selected = static_cast<std::uint32_t>(
                found - cumulative_weights_.begin());
        }

        const TargetChoice &target = target_choices_.at(
            static_cast<std::size_t>(selected));
        std::uint8_t selected_haplotype = 0U;
        switch (target.eligible_haplotypes) {
        case model::HaplotypeMask::haplotype_1:
            selected_haplotype = 0U;
            break;
        case model::HaplotypeMask::haplotype_2:
            selected_haplotype = 1U;
            break;
        case model::HaplotypeMask::both:
            selected_haplotype = static_cast<std::uint8_t>(
                rng::bernoulli(
                    haplotype_key,
                    fragment_ordinal,
                    UINT64_C(0),
                    0.5)
                ? 1U : 0U);
            break;
        default:
            throw TbsCatalogError(
                "TBS target has an invalid eligible-haplotype mask");
        }

        Candidate candidate;
        if (center_sd_ == 0.0 && !variable_insert_) {
            if (!project(
                    target,
                    selected_haplotype,
                    0,
                    insert_parameters_.mean,
                    candidate)) {
                throw TbsCatalogError(
                    "fixed diploid TBS target changed eligibility");
            }
        } else {
            bool accepted = false;
            for (std::uint32_t attempt = 0U;
                 attempt < maximum_attempts_per_fragment;
                 ++attempt) {
                std::uint32_t insert_length = insert_parameters_.mean;
                if (variable_insert_) {
                    const std::uint64_t length_local = center_sd_ == 0.0
                        ? UINT64_C(2) + attempt
                        : UINT64_C(2) + UINT64_C(2) * attempt;
                    insert_length = insert_sampler.sample(
                        ordinal, length_local);
                }
                std::int64_t center_displacement = 0;
                if (center_sd_ > 0.0) {
                    const std::uint64_t center_local = variable_insert_
                        ? UINT64_C(3) + UINT64_C(2) * attempt
                        : UINT64_C(2) + attempt;
                    const double normal = normal_sampler::standard_normal(
                        fragment_key, ordinal, center_local);
                    const double displacement = center_sd_ * normal;
                    const long double extended_displacement = displacement;
                    if (!std::isfinite(displacement)
                        || extended_displacement
                            < static_cast<long double>(INT64_MIN)
                        || extended_displacement
                            > static_cast<long double>(INT64_MAX)) {
                        throw TbsCatalogError(
                            "TBS center displacement exceeds int64");
                    }
                    center_displacement = static_cast<std::int64_t>(
                        displacement);
                }
                if (project(
                        target,
                        selected_haplotype,
                        center_displacement,
                        insert_length,
                        candidate)) {
                    accepted = true;
                    break;
                }
                if (result.skipped_count
                    == std::numeric_limits<std::uint64_t>::max()) {
                    throw TbsCatalogError("TBS skipped count exceeds uint64");
                }
                ++result.skipped_count;
            }
            if (!accepted) {
                throw TbsCatalogError(
                    "TBS fragment sampler exhausted its attempt cap");
            }
        }
        result.candidates.push_back(candidate);
    }
    return result;
}

} // namespace htsim::tbs
