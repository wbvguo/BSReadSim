#include "fragment.h"

#include <algorithm>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <cstdint>
#include <optional>
#include <unordered_map>
#include <utility>

#include "utilities.h"
#include "protocol.h"

// ---- allocation --------------------------------------------------------

namespace htsim::allocation {
namespace {

struct RemainderSeat {
    std::size_t index;
    std::uint64_t remainder;
};

struct RealRemainderSeat {
    std::size_t index;
    long double remainder;
};

} // namespace

std::vector<std::uint32_t> largest_remainder(
    const std::vector<std::uint32_t> &weights,
    std::uint32_t total_fragments)
{
    if (weights.empty()) {
        throw AllocationError("allocation requires at least one contig");
    }
    if (total_fragments == 0) {
        throw AllocationError("allocation total must be positive");
    }

    std::uint64_t total_weight = 0;
    for (const std::uint32_t weight : weights) {
        if (weight > std::numeric_limits<std::uint64_t>::max() - total_weight) {
            throw AllocationError("allocation weight sum exceeds uint64");
        }
        total_weight += weight;
    }
    if (total_weight == 0) {
        throw AllocationError("no contig has an eligible fragment start");
    }

    std::vector<std::uint32_t> result(weights.size(), 0);
    std::vector<RemainderSeat> seats;
    seats.reserve(weights.size());
    std::uint64_t assigned = 0;
    for (std::size_t index = 0; index < weights.size(); ++index) {
        if (weights[index] == 0) {continue;}
        // Both factors are bounded by UINT32_MAX.  Their largest product is
        // 18446744065119617025, which is UINT64_MAX - 8589934590 and therefore
        // fits exactly in uint64_t.  The sum of weights uses uint64_t because
        // multiple human-sized contigs can exceed UINT32_MAX in aggregate.
        const std::uint64_t numerator =
            static_cast<std::uint64_t>(total_fragments) * weights[index];
        const std::uint64_t quotient = numerator / total_weight;
        const std::uint64_t remainder = numerator % total_weight;
        if (quotient > std::numeric_limits<std::uint32_t>::max()) {
            throw AllocationError("internal allocation quotient overflow");
        }
        result[index] = static_cast<std::uint32_t>(quotient);
        if (quotient > total_fragments - assigned) {
            throw AllocationError("internal allocation total overflow");
        }
        assigned += quotient;
        seats.push_back({index, remainder});
    }

    const std::uint64_t unassigned = total_fragments - assigned;
    if (unassigned > seats.size()) {
        throw AllocationError("internal largest-remainder invariant failed");
    }
    std::sort(
        seats.begin(), seats.end(),
        [](const RemainderSeat &left, const RemainderSeat &right) {
            return left.remainder != right.remainder
                ? left.remainder > right.remainder
                : left.index < right.index;
        });
    for (std::size_t offset = 0;
         offset < static_cast<std::size_t>(unassigned);
         ++offset) {
        ++result[seats[offset].index];
    }

    const std::uint64_t observed = std::accumulate(
        result.begin(), result.end(), UINT64_C(0));
    if (observed != total_fragments) {
        throw AllocationError("allocation did not preserve the requested total");
    }
    return result;
}

std::vector<std::uint32_t> largest_remainder_real(
    const std::vector<double> &weights,
    std::uint32_t total_fragments)
{
    if (weights.empty()) {
        throw AllocationError("allocation requires at least one contig");
    }
    if (total_fragments == 0U) {
        throw AllocationError("allocation total must be positive");
    }

    long double total_weight = 0.0L;
    for (const double weight : weights) {
        if (!std::isfinite(weight) || weight < 0.0) {
            throw AllocationError(
                "real allocation weights must be finite and non-negative");
        }
        total_weight += static_cast<long double>(weight);
    }
    if (!std::isfinite(total_weight) || total_weight <= 0.0L) {
        throw AllocationError("no contig has positive real allocation weight");
    }

    std::vector<std::uint32_t> result(weights.size(), 0U);
    std::vector<RealRemainderSeat> seats;
    seats.reserve(weights.size());
    std::uint64_t assigned = 0U;
    for (std::size_t index = 0U; index < weights.size(); ++index) {
        if (weights[index] == 0.0) {continue;}
        const long double quota =
            static_cast<long double>(total_fragments)
            * static_cast<long double>(weights[index]) / total_weight;
        const long double floored = std::floor(quota);
        if (floored < 0.0L
            || floored > std::numeric_limits<std::uint32_t>::max()) {
            throw AllocationError("internal real allocation quotient overflow");
        }
        const std::uint32_t count = static_cast<std::uint32_t>(floored);
        if (count > total_fragments - assigned) {
            throw AllocationError("internal real allocation total overflow");
        }
        result[index] = count;
        assigned += count;
        seats.push_back({index, quota - floored});
    }

    const std::uint64_t unassigned = total_fragments - assigned;
    if (unassigned > seats.size()) {
        throw AllocationError("internal real largest-remainder invariant failed");
    }
    std::sort(
        seats.begin(), seats.end(),
        [](const RealRemainderSeat &left, const RealRemainderSeat &right) {
            return left.remainder != right.remainder
                ? left.remainder > right.remainder
                : left.index < right.index;
        });
    for (std::size_t offset = 0U;
         offset < static_cast<std::size_t>(unassigned);
         ++offset) {
        ++result[seats[offset].index];
    }
    if (std::accumulate(result.begin(), result.end(), UINT64_C(0))
        != total_fragments) {
        throw AllocationError(
            "real allocation did not preserve the requested total");
    }
    return result;
}

} // namespace htsim::allocation

namespace htsim::depth_count {

std::uint32_t fragments(
    double depth,
    std::uint64_t effective_reference_bases,
    std::uint32_t read_length,
    bool paired_end)
{
    if (std::fegetround() != FE_TONEAREST) {
        throw DepthCountError(
            "depth conversion requires round-to-nearest floating point");
    }
    if (!std::isfinite(depth) || depth <= 0.0) {
        throw DepthCountError("depth must be finite and positive");
    }
    if (effective_reference_bases == 0U) {
        throw DepthCountError(
            "depth conversion requires positive effective reference length");
    }
    if (read_length == 0U) {
        throw DepthCountError("depth conversion requires positive read length");
    }

    const std::uint64_t sequenced_bases_per_fragment =
        static_cast<std::uint64_t>(read_length)
        * (paired_end ? UINT64_C(2) : UINT64_C(1));
    const double numerator =
        static_cast<double>(effective_reference_bases) * depth;
    const double raw_count = numerator
        / static_cast<double>(sequenced_bases_per_fragment);
    if (!std::isfinite(raw_count)) {
        throw DepthCountError("depth-derived read-pair count is not finite");
    }
    const double floored = std::ceil(raw_count);
    if (floored < 1.0) {
        throw DepthCountError(
            "depth derives zero read pairs for the effective reference");
    }
    if (floored
        > static_cast<double>(std::numeric_limits<std::uint32_t>::max())) {
        throw DepthCountError("depth-derived read-pair count exceeds uint32");
    }
    return static_cast<std::uint32_t>(floored);
}

} // namespace htsim::depth_count

// ---- insert_length --------------------------------------------------------

namespace htsim::insert_length {

Sampler::Sampler(
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    const Parameters &parameters)
    : parameters_(parameters)
{
    if (parameters.minimum == 0U
        || parameters.mean < parameters.minimum
        || parameters.maximum < parameters.mean) {
        throw InsertLengthError(
            "insert lengths must satisfy 0 < minimum <= mean <= maximum");
    }
    if (!std::isfinite(parameters.standard_deviation)
        || parameters.standard_deviation < 0.0) {
        throw InsertLengthError(
            "insert-length standard deviation must be finite and non-negative");
    }

    key_ = rng::derive_key(master_seed, rng::Stage::fragment, contig_index);
}

std::uint32_t Sampler::sample(std::uint64_t candidate_ordinal) const
{
    if (parameters_.standard_deviation == 0.0
        || parameters_.minimum == parameters_.maximum) {
        return parameters_.mean;
    }

    const double normal = normal_sampler::standard_normal(
        key_, candidate_ordinal, UINT64_C(0));
    if (normal == 0.0) {return parameters_.mean;}

    const std::int64_t lower_delta =
        static_cast<std::int64_t>(parameters_.minimum)
        - static_cast<std::int64_t>(parameters_.mean);
    const std::int64_t upper_delta =
        static_cast<std::int64_t>(parameters_.maximum)
        - static_cast<std::int64_t>(parameters_.mean);
    const double lower_normal = static_cast<double>(lower_delta)
        / parameters_.standard_deviation;
    const double upper_normal = static_cast<double>(upper_delta)
        / parameters_.standard_deviation;
    if (normal <= lower_normal) {return parameters_.minimum;}
    if (normal >= upper_normal) {return parameters_.maximum;}

    const double scaled = parameters_.standard_deviation * normal;
    if (!std::isfinite(scaled)) {
        throw InsertLengthError(
            "insert-length deviation is unexpectedly non-finite");
    }
    const std::int64_t deviation = static_cast<std::int64_t>(scaled);
    const std::int64_t length =
        static_cast<std::int64_t>(parameters_.mean) + deviation;
    if (length < static_cast<std::int64_t>(parameters_.minimum)
        || length > static_cast<std::int64_t>(parameters_.maximum)) {
        throw InsertLengthError(
            "insert-length clamp invariant was violated");
    }
    return static_cast<std::uint32_t>(length);
}

std::uint32_t sample(
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    std::uint64_t candidate_ordinal,
    const Parameters &parameters)
{
    return Sampler(master_seed, contig_index, parameters)
        .sample(candidate_ordinal);
}

} // namespace htsim::insert_length

// ---- fragment_builder --------------------------------------------------------

namespace htsim::fragment_builder {
namespace {

void validate_layout(const ReadLayout &layout)
{
    if (layout.insert_length == 0 || layout.read_length == 0) {
        throw FragmentBuilderError("insert and read lengths must be positive");
    }
    if (layout.read_length > layout.insert_length) {
        throw FragmentBuilderError("read length must not exceed insert length");
    }
    if (layout.insert_coordinate != ReadLayout::InsertCoordinate::reference
        && layout.insert_coordinate
            != ReadLayout::InsertCoordinate::haplotype) {
        throw FragmentBuilderError("insert coordinate system is invalid");
    }
}

void validate_capture_strand(model::CaptureStrand capture_strand)
{
    if (capture_strand != model::CaptureStrand::unknown
        && capture_strand != model::CaptureStrand::forward
        && capture_strand != model::CaptureStrand::reverse) {
        throw FragmentBuilderError("capture strand is invalid");
    }
}

void validate_detail(FragmentDetail detail)
{
    if (detail != FragmentDetail::full
        && detail != FragmentDetail::common_columns) {
        throw FragmentBuilderError("fragment detail is invalid");
    }
}

void add_payload_bytes(
    std::uint64_t &total,
    std::uint64_t increment)
{
    if (increment > std::numeric_limits<std::uint64_t>::max() - total) {
        throw FragmentBuilderError("fragment payload byte bound exceeds uint64");
    }
    total += increment;
}

std::uint64_t payload_base_bound(
    std::uint64_t template_length,
    std::uint64_t read_length,
    std::uint64_t mate_count)
{
    // 56 bytes covers fixed fields and all vector length prefixes. Per
    // template base: 1 base + 8 reference position + 4 event id + 24 site.
    // Per mate: 32 fixed bytes plus at most one 8-byte site reference per read
    // base. PE overlap legitimately references the shared fragment site twice.
    return UINT64_C(56) + UINT64_C(37) * template_length
        + UINT64_C(32) * mate_count
        + UINT64_C(8) * mate_count * read_length;
}

bool valid_event_base(const model::Bases &bases) noexcept
{
    return std::all_of(
        bases.begin(), bases.end(),
        [](std::uint8_t base) {return base <= 3U;});
}

void validate_projected_interval(
    const haplotype::ProjectedInterval &projection,
    const ReadLayout &layout)
{
    validate_layout(layout);
    if (projection.haplotype > 1U) {
        throw FragmentBuilderError("projected haplotype must be zero or one");
    }
    if (projection.reference_start > projection.reference_end) {
        throw FragmentBuilderError("projection reference interval is reversed");
    }
    if (layout.insert_coordinate == ReadLayout::InsertCoordinate::reference
        && projection.reference_end - projection.reference_start
            != layout.insert_length) {
        throw FragmentBuilderError(
            "projection reference span disagrees with insert length");
    }
    if (projection.template_bases.empty()
        || projection.template_bases.size()
            > std::numeric_limits<std::uint32_t>::max()
        || projection.template_bases.size()
            != projection.reference_positions.size()
        || projection.template_bases.size()
            != projection.base_variant_indices.size()
        || projection.template_bases.size() < layout.read_length
        || projection.variants.size()
            > std::numeric_limits<std::uint32_t>::max()) {
        throw FragmentBuilderError("projected fragment shape is invalid");
    }
    if (layout.insert_coordinate == ReadLayout::InsertCoordinate::haplotype
        && projection.template_bases.size() != layout.insert_length) {
        throw FragmentBuilderError(
            "projection haplotype span disagrees with insert length");
    }

    std::unordered_map<std::uint32_t, std::size_t> event_indices;
    std::vector<std::uint32_t> observed_event_bases(
        projection.variants.size(), 0U);
    std::vector<std::size_t> last_variant_offsets(
        projection.variants.size(),
        std::numeric_limits<std::size_t>::max());
    for (std::size_t index = 0; index < projection.variants.size(); ++index) {
        const model::Variant &event = projection.variants[index];
        if (event.index == model::no_variant_index
            || !event_indices.emplace(event.index, index).second
            || (index != 0U
                && (projection.variants[index - 1U].index
                        >= event.index
                    || projection.variants[index - 1U].reference_start
                        > event.reference_start))
            || (event.phased_haplotype != 255U
                && event.phased_haplotype != projection.haplotype)
            || event.reference_start < projection.reference_start
            || event.reference_end > projection.reference_end
            || event.reference_start > event.reference_end
            || !valid_event_base(event.ref_bases)
            || !valid_event_base(event.alt_bases)) {
            throw FragmentBuilderError("projected variant event is invalid");
        }
        switch (event.kind) {
        case model::VariantKind::snv:
            if (event.reference_end - event.reference_start != 1U
                || event.ref_bases.size() != 1U
                || event.alt_bases.size() != 1U) {
                throw FragmentBuilderError("projected SNV shape is invalid");
            }
            break;
        case model::VariantKind::insertion:
            if (event.reference_start != event.reference_end
                || !event.ref_bases.empty()
                || event.alt_bases.empty()
                || event.alt_bases.size() > 4U) {
                throw FragmentBuilderError("projected insertion shape is invalid");
            }
            break;
        case model::VariantKind::deletion:
            if (event.reference_start == event.reference_end
                || event.ref_bases.empty()
                || event.ref_bases.size() > 4U
                || !event.alt_bases.empty()
                || event.ref_bases.size()
                    != event.reference_end - event.reference_start) {
                throw FragmentBuilderError("projected deletion shape is invalid");
            }
            break;
        default:
            throw FragmentBuilderError("projected variant kind is invalid");
        }
    }

    std::optional<std::int64_t> previous_mapped_position;
    for (std::size_t offset = 0; offset < projection.template_bases.size(); ++offset) {
        const std::uint8_t base = projection.template_bases[offset];
        const std::int64_t position = projection.reference_positions[offset];
        const std::uint32_t variant_index = projection.base_variant_indices[offset];
        if (base > 4U || position < -1
            || (position >= 0
                && (static_cast<std::uint64_t>(position)
                        < projection.reference_start
                    || static_cast<std::uint64_t>(position)
                        >= projection.reference_end))) {
            throw FragmentBuilderError("projected base or coordinate is invalid");
        }
        if (position >= 0) {
            if (previous_mapped_position
                && position <= *previous_mapped_position) {
                throw FragmentBuilderError(
                    "projected mapped positions are not strictly increasing");
            }
            previous_mapped_position = position;
        }
        if (variant_index == model::no_variant_index) {
            if (position == -1) {
                throw FragmentBuilderError(
                    "inserted projected base lacks an event id");
            }
            continue;
        }
        const auto found = event_indices.find(variant_index);
        if (found == event_indices.end()) {
            throw FragmentBuilderError(
                "projected base refers to an unknown event");
        }
        const model::Variant &event =
            projection.variants[found->second];
        std::uint32_t &observed = observed_event_bases[found->second];
        std::size_t &last_offset = last_variant_offsets[found->second];
        if (event.kind == model::VariantKind::deletion
            || observed >= event.alt_bases.size()
            || base != event.alt_bases[observed]
            || (observed != 0U && last_offset + 1U != offset)) {
            throw FragmentBuilderError(
                "projected bases disagree with their variant event");
        }
        if (event.kind == model::VariantKind::insertion) {
            if (position != -1) {
                throw FragmentBuilderError(
                    "inserted event base has a mapped reference position");
            }
        } else if (position != static_cast<std::int64_t>(
                       event.reference_start + observed)) {
            throw FragmentBuilderError(
                "mapped event base has the wrong reference position");
        }
        last_offset = offset;
        ++observed;
    }
    for (std::size_t index = 0; index < projection.variants.size(); ++index) {
        const model::Variant &event = projection.variants[index];
        const std::uint32_t expected = event.kind == model::VariantKind::deletion
            ? 0U
            : static_cast<std::uint32_t>(event.alt_bases.size());
        if (observed_event_bases[index] != expected) {
            throw FragmentBuilderError(
                "projected event bases are incomplete or duplicated");
        }
    }
}

const model::Variant &insertion_event(
    const model::Fragment &fragment,
    std::uint32_t variant_index)
{
    const auto event = std::find_if(
        fragment.variants.begin(), fragment.variants.end(),
        [variant_index](const model::Variant &candidate) {
            return candidate.index == variant_index;
        });
    if (event == fragment.variants.end()
        || event->kind != model::VariantKind::insertion) {
        throw FragmentBuilderError(
            "insertion-only mate refers to an invalid event");
    }
    return *event;
}

model::Mate build_mate(
    const model::Fragment &fragment,
    std::uint8_t mate_index,
    bool reverse_complement,
    std::uint32_t template_start,
    std::uint32_t template_end)
{
    model::Mate mate;
    mate.mate_index = mate_index;
    mate.reverse_complement = reverse_complement;
    mate.template_start = template_start;
    mate.template_end = template_end;
    std::optional<std::uint64_t> insertion_anchor;
    bool multiple_insertion_anchors = false;
    std::optional<std::uint64_t> minimum_position;
    std::optional<std::uint64_t> maximum_position;
    for (std::uint32_t offset = template_start; offset < template_end; ++offset) {
        const std::int64_t position = fragment.reference_positions[offset];
        if (position >= 0) {
            const auto mapped = static_cast<std::uint64_t>(position);
            minimum_position = minimum_position
                ? std::min(*minimum_position, mapped)
                : mapped;
            maximum_position = maximum_position
                ? std::max(*maximum_position, mapped)
                : mapped;
            continue;
        }
        const auto &event = insertion_event(
            fragment, fragment.base_variant_indices[offset]);
        multiple_insertion_anchors = multiple_insertion_anchors
            || (insertion_anchor
                && *insertion_anchor != event.reference_start);
        insertion_anchor = event.reference_start;
    }
    if (minimum_position) {
        mate.reference_start = *minimum_position;
        mate.reference_end = *maximum_position + 1U;
    } else if (insertion_anchor) {
        if (multiple_insertion_anchors) {
            throw FragmentBuilderError(
                "an insertion-only mate spans multiple reference anchors");
        }
        mate.reference_start = *insertion_anchor;
        mate.reference_end = *insertion_anchor;
    } else {
        throw FragmentBuilderError("mate slice has no projected bases");
    }
    for (const model::MethylationSite &site : fragment.methylation_sites) {
        if (site.template_offset < template_start
            || site.template_offset >= template_end) {
            continue;
        }
        const std::uint32_t read_offset = reverse_complement
            ? template_end - 1U - site.template_offset
            : site.template_offset - template_start;
        mate.site_refs.push_back({read_offset, site.site_index});
    }
    std::sort(
        mate.site_refs.begin(), mate.site_refs.end(),
        [](const model::SiteReference &left,
           const model::SiteReference &right) {
            return left.read_offset < right.read_offset;
        });
    return mate;
}

model::Mate build_common_mate(
    std::uint8_t mate_index,
    bool reverse_complement,
    std::uint32_t template_start,
    std::uint32_t template_end)
{
    model::Mate mate;
    mate.mate_index = mate_index;
    mate.reverse_complement = reverse_complement;
    mate.template_start = template_start;
    mate.template_end = template_end;
    return mate;
}

} // namespace

std::uint64_t maximum_payload_bytes(const ReadLayout &layout)
{
    validate_layout(layout);
    const std::uint64_t insert = layout.insert_length;
    const std::uint64_t read = layout.read_length;
    const std::uint64_t mate_count = layout.paired_end ? 2U : 1U;

    return payload_base_bound(insert, read, mate_count);
}

void require_payload_fits_protocol(const ReadLayout &layout)
{
    if (maximum_payload_bytes(layout) > protocol::maximum_frame_payload) {
        throw FragmentBuilderError(
            "maximum fragment payload exceeds the protocol frame limit");
    }
}

std::uint64_t maximum_payload_bytes(
    const haplotype::ProjectedInterval &projection,
    const ReadLayout &layout)
{
    validate_projected_interval(projection, layout);
    const std::uint64_t mate_count = layout.paired_end ? 2U : 1U;
    std::uint64_t total = payload_base_bound(
        projection.template_bases.size(), layout.read_length, mate_count);
    for (const model::Variant &event : projection.variants) {
        add_payload_bytes(total, 32U);
        add_payload_bytes(total, event.ref_bases.size());
        add_payload_bytes(total, event.alt_bases.size());
    }
    return total;
}

void require_payload_fits_protocol(
    const haplotype::ProjectedInterval &projection,
    const ReadLayout &layout)
{
    if (maximum_payload_bytes(projection, layout)
        > protocol::maximum_frame_payload) {
        throw FragmentBuilderError(
            "maximum projected fragment payload exceeds the protocol frame limit");
    }
}

model::Fragment build_fragment(
    const reference::Contig &contig,
    const methdb::MethylationCatalog &catalog,
    std::uint64_t fragment_ordinal,
    std::uint32_t reference_start,
    std::uint8_t haplotype,
    model::CaptureStrand capture_strand,
    const ReadLayout &layout,
    FragmentDetail detail)
{
    validate_layout(layout);
    require_payload_fits_protocol(layout);
    validate_detail(detail);
    if (haplotype > 1) {
        throw FragmentBuilderError("haplotype must be zero or one");
    }
    validate_capture_strand(capture_strand);
    if (contig.length != contig.bases.size()
        || contig.length > std::numeric_limits<std::uint32_t>::max()) {
        throw FragmentBuilderError("materialized contig length is inconsistent");
    }
    const std::uint64_t end =
        static_cast<std::uint64_t>(reference_start) + layout.insert_length;
    if (end > contig.length) {
        throw FragmentBuilderError("fragment exceeds its contig");
    }
    const std::uint32_t reference_end = static_cast<std::uint32_t>(end);

    model::Fragment fragment;
    fragment.fragment_ordinal = fragment_ordinal;
    fragment.contig_index = contig.index;
    fragment.haplotype = haplotype;
    fragment.capture_strand = capture_strand;
    fragment.reference_start = reference_start;
    fragment.reference_end = reference_end;

    const auto first = contig.bases.begin()
        + static_cast<std::ptrdiff_t>(reference_start);
    const auto last = first + static_cast<std::ptrdiff_t>(layout.insert_length);
    fragment.template_bases.assign(first, last);
    if (detail == FragmentDetail::full) {
        fragment.reference_positions.reserve(layout.insert_length);
        fragment.base_variant_indices.assign(
            layout.insert_length, model::no_variant_index);
        for (std::uint32_t offset = 0; offset < layout.insert_length; ++offset) {
            fragment.reference_positions.push_back(
                static_cast<std::int64_t>(reference_start) + offset);
        }
    }

    const auto sites = catalog.sites_in_range(reference_start, reference_end);
    for (auto site = sites.first; site != sites.second; ++site) {
        const std::uint32_t site_index = static_cast<std::uint32_t>(
            fragment.methylation_sites.size());
        fragment.methylation_sites.push_back(model::MethylationSite{
            site_index,
            site->reference_position - reference_start,
            static_cast<std::int64_t>(site->reference_position),
            site->context,
            site->methylation_source,
            model::MethylationAllele::shared,
            site->methylation_probability,
        });
    }

    const auto append_mate = [&](std::uint8_t mate_index,
                                 bool reverse_complement,
                                 std::uint32_t template_start,
                                 std::uint32_t template_end) {
        fragment.mates.push_back(detail == FragmentDetail::full
            ? build_mate(
                  fragment,
                  mate_index,
                  reverse_complement,
                  template_start,
                  template_end)
            : build_common_mate(
                  mate_index,
                  reverse_complement,
                  template_start,
                  template_end));
    };
    append_mate(0U, false, 0U, layout.read_length);
    if (layout.paired_end) {
        const std::uint32_t second_start =
            layout.insert_length - layout.read_length;
        append_mate(1U, true, second_start, layout.insert_length);
    }
    return fragment;
}


model::Fragment build_fragment(
    haplotype::ProjectedInterval projection,
    const methdb::DiploidMethylationCatalog &catalog,
    std::uint64_t fragment_ordinal,
    model::CaptureStrand capture_strand,
    const ReadLayout &layout,
    FragmentDetail detail)
{
    validate_capture_strand(capture_strand);
    require_payload_fits_protocol(projection, layout);
    validate_detail(detail);

    model::Fragment fragment;
    fragment.fragment_ordinal = fragment_ordinal;
    fragment.contig_index = projection.contig_index;
    fragment.haplotype = projection.haplotype;
    fragment.capture_strand = capture_strand;
    fragment.reference_start = projection.reference_start;
    fragment.reference_end = projection.reference_end;
    fragment.methylation_sites = catalog.sites_for_projection(projection);
    fragment.template_bases = std::move(projection.template_bases);
    if (detail == FragmentDetail::full) {
        fragment.reference_positions = std::move(projection.reference_positions);
        fragment.base_variant_indices = std::move(projection.base_variant_indices);
        fragment.variants = std::move(projection.variants);
    }

    const auto template_length = static_cast<std::uint32_t>(
        fragment.template_bases.size());
    const auto append_mate = [&](std::uint8_t mate_index,
                                 bool reverse_complement,
                                 std::uint32_t template_start,
                                 std::uint32_t template_end) {
        fragment.mates.push_back(detail == FragmentDetail::full
            ? build_mate(
                  fragment,
                  mate_index,
                  reverse_complement,
                  template_start,
                  template_end)
            : build_common_mate(
                  mate_index,
                  reverse_complement,
                  template_start,
                  template_end));
    };
    append_mate(0U, false, 0U, layout.read_length);
    if (layout.paired_end) {
        append_mate(
            1U,
            true,
            template_length - layout.read_length,
            template_length);
    }
    return fragment;
}

} // namespace htsim::fragment_builder
