#include "methdb.h"

#include <cstring>
#include <fstream>
#include <iomanip>

#include <cstddef>
#include <limits>
#include <array>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>
#include <sys/types.h>
#include <cfenv>
#include <algorithm>
#include <memory>
#include <tuple>

#include <zlib.h>

#include "utilities.h"
#include "types.h"

// ---- site --------------------------------------------------------

namespace htsim::methdb {
namespace {

constexpr std::uint8_t base_a = 0;
constexpr std::uint8_t base_c = 1;
constexpr std::uint8_t base_g = 2;
constexpr std::uint8_t base_t = 3;
constexpr std::uint8_t base_n = 4;

std::uint8_t checked_base(const model::Bases &bases, std::size_t position)
{
    const std::uint8_t base = bases[position];
    if (base > base_n) {
        throw ContextError("contig contains a base outside protocol encoding");
    }
    return base;
}

std::uint8_t checked_base(std::uint8_t base)
{
    if (base > base_n) {
        throw ContextError("context contains a base outside protocol encoding");
    }
    return base;
}

bool is_resolved_non_c(std::uint8_t base) noexcept
{
    return base == base_a || base == base_g || base == base_t;
}

bool is_resolved_non_g(std::uint8_t base) noexcept
{
    return base == base_a || base == base_c || base == base_t;
}

} // namespace

std::optional<model::MethylationContext> classify_context(
    const ContextNeighborhood &neighborhood,
    bool collect_non_cpg)
{
    const std::uint8_t center = checked_base(neighborhood.center);
    if (center == base_c) {
        if (!neighborhood.downstream_first) {return std::nullopt;}
        const std::uint8_t next = checked_base(*neighborhood.downstream_first);
        if (next == base_g) {return model::MethylationContext::cg_c;}
        if (!collect_non_cpg) {return std::nullopt;}
        if (!is_resolved_non_g(next) || !neighborhood.downstream_second) {
            return std::nullopt;
        }
        const std::uint8_t second = checked_base(*neighborhood.downstream_second);
        if (second == base_g) {return model::MethylationContext::chg_c;}
        if (is_resolved_non_g(second)) {
            return model::MethylationContext::chh_c;
        }
        return std::nullopt;
    }

    if (center == base_g) {
        if (!neighborhood.upstream_first) {return std::nullopt;}
        const std::uint8_t previous = checked_base(*neighborhood.upstream_first);
        if (previous == base_c) {return model::MethylationContext::cg_g;}
        if (!collect_non_cpg) {return std::nullopt;}
        if (!is_resolved_non_c(previous) || !neighborhood.upstream_second) {
            return std::nullopt;
        }
        const std::uint8_t second = checked_base(*neighborhood.upstream_second);
        if (second == base_c) {return model::MethylationContext::chg_g;}
        if (is_resolved_non_c(second)) {
            return model::MethylationContext::chh_g;
        }
        return std::nullopt;
    }
    return std::nullopt;
}

std::optional<model::MethylationContext> classify_context(
    const model::Bases &contig_bases,
    std::uint64_t reference_position,
    bool collect_non_cpg)
{
    if (reference_position > std::numeric_limits<std::size_t>::max()
        || static_cast<std::size_t>(reference_position) >= contig_bases.size()) {
        throw ContextError("methylation position is outside its contig");
    }
    const std::size_t position = static_cast<std::size_t>(reference_position);
    return classify_context(
        ContextNeighborhood{
            position >= 2U
                ? std::optional<std::uint8_t>(
                      contig_bases[position - 2U])
                : std::nullopt,
            position >= 1U
                ? std::optional<std::uint8_t>(
                      contig_bases[position - 1U])
                : std::nullopt,
            checked_base(contig_bases, position),
            position + 1U < contig_bases.size()
                ? std::optional<std::uint8_t>(
                      contig_bases[position + 1U])
                : std::nullopt,
            position + 2U < contig_bases.size()
                ? std::optional<std::uint8_t>(
                      contig_bases[position + 2U])
                : std::nullopt,
        },
        collect_non_cpg);
}

} // namespace htsim::methdb

namespace htsim::methdb {
namespace {

inline constexpr unsigned kind_shift = 60U;
inline constexpr std::uint64_t payload_mask = UINT64_C(0x0fffffffffffffff);
inline constexpr std::uint64_t reference_payload_mask =
    UINT64_C(0x00000000ffffffff);
inline constexpr std::uint64_t insertion_payload_mask =
    UINT64_C(0x00000003ffffffff);
inline constexpr unsigned insertion_event_shift = 2U;

std::uint64_t encode(SiteEntityKind kind, std::uint64_t payload) noexcept
{
    return (static_cast<std::uint64_t>(kind) << kind_shift) | payload;
}

void require_selected_haplotype(
    model::HaplotypeMask mask,
    std::uint8_t zero_based_haplotype)
{
    if (!model::is_haplotype_mask(static_cast<std::uint8_t>(mask))) {
        throw EntityError("methylation site has an invalid haplotype mask");
    }
    if (zero_based_haplotype > 1U
        || !model::mask_contains(mask, zero_based_haplotype)) {
        throw EntityError(
            "methylation site is not present on the selected haplotype");
    }
}

SiteEntityKind reference_variant_kind(
    model::HaplotypeMask mask,
    std::uint8_t zero_based_haplotype)
{
    require_selected_haplotype(mask, zero_based_haplotype);
    if (mask == model::HaplotypeMask::both) {
        return SiteEntityKind::variant_reference_shared;
    }
    return zero_based_haplotype == 0U
        ? SiteEntityKind::variant_reference_haplotype_0
        : SiteEntityKind::variant_reference_haplotype_1;
}

SiteEntityKind insertion_kind(
    model::HaplotypeMask mask,
    std::uint8_t zero_based_haplotype)
{
    require_selected_haplotype(mask, zero_based_haplotype);
    if (mask == model::HaplotypeMask::both) {
        return SiteEntityKind::insertion_shared;
    }
    return zero_based_haplotype == 0U
        ? SiteEntityKind::insertion_haplotype_0
        : SiteEntityKind::insertion_haplotype_1;
}

bool valid_kind(std::uint8_t kind) noexcept
{
    return kind <= static_cast<std::uint8_t>(
        SiteEntityKind::insertion_haplotype_1);
}

} // namespace

SiteEntity reference_site_entity(std::uint32_t reference_position)
{
    return SiteEntity(reference_position);
}

SiteEntity variant_reference_site_entity(
    std::uint32_t reference_position,
    model::HaplotypeMask alt_haplotypes,
    std::uint8_t zero_based_haplotype)
{
    return SiteEntity(encode(
        reference_variant_kind(alt_haplotypes, zero_based_haplotype),
        reference_position));
}

SiteEntity insertion_site_entity(
    std::uint32_t event_ordinal,
    std::uint8_t insertion_offset,
    model::HaplotypeMask alt_haplotypes,
    std::uint8_t zero_based_haplotype)
{
    require_selected_haplotype(alt_haplotypes, zero_based_haplotype);
    if (event_ordinal == model::no_variant_index) {
        throw EntityError("insertion event ordinal uses the no-event sentinel");
    }
    if (insertion_offset >= model::maximum_insertion_bases) {
        throw EntityError("insertion offset must be in [0, 3]");
    }
    const std::uint64_t payload =
        (static_cast<std::uint64_t>(event_ordinal) << insertion_event_shift)
        | insertion_offset;
    return SiteEntity(encode(
        insertion_kind(alt_haplotypes, zero_based_haplotype), payload));
}

DecodedSiteEntity decode_site_entity(SiteEntity entity) noexcept
{
    const std::uint64_t encoded = entity.value();
    const auto kind = static_cast<SiteEntityKind>(encoded >> kind_shift);
    const std::uint64_t payload = encoded & payload_mask;
    if (kind == SiteEntityKind::reference_baseline
        || kind == SiteEntityKind::variant_reference_shared
        || kind == SiteEntityKind::variant_reference_haplotype_0
        || kind == SiteEntityKind::variant_reference_haplotype_1) {
        return {kind, static_cast<std::uint32_t>(payload), 0U, 0U};
    }
    return {
        kind,
        0U,
        static_cast<std::uint32_t>(payload >> insertion_event_shift),
        static_cast<std::uint8_t>(payload & UINT64_C(0x3)),
    };
}

SiteEntity decode_site_entity(std::uint64_t encoded)
{
    const std::uint8_t kind_value = static_cast<std::uint8_t>(
        encoded >> kind_shift);
    if (!valid_kind(kind_value)) {
        throw EntityError("methylation site entity uses a reserved kind tag");
    }
    const auto kind = static_cast<SiteEntityKind>(kind_value);
    const std::uint64_t payload = encoded & payload_mask;
    const bool insertion =
        kind == SiteEntityKind::insertion_shared
        || kind == SiteEntityKind::insertion_haplotype_0
        || kind == SiteEntityKind::insertion_haplotype_1;
    if (!insertion) {
        if ((payload & ~reference_payload_mask) != 0U) {
            throw EntityError(
                "reference site entity has non-zero reserved bits");
        }
    } else {
        if ((payload & ~insertion_payload_mask) != 0U) {
            throw EntityError(
                "insertion site entity has non-zero reserved bits");
        }
        const std::uint32_t event_ordinal = static_cast<std::uint32_t>(
            payload >> insertion_event_shift);
        if (event_ordinal == model::no_variant_index) {
            throw EntityError(
                "insertion site entity uses the no-event sentinel");
        }
    }
    return SiteEntity(encoded);
}

bool entity_is_insertion(SiteEntity entity) noexcept
{
    const SiteEntityKind kind = decode_site_entity(entity).kind;
    return kind == SiteEntityKind::insertion_shared
        || kind == SiteEntityKind::insertion_haplotype_0
        || kind == SiteEntityKind::insertion_haplotype_1;
}

} // namespace htsim::methdb

// ---- cgmap_profile --------------------------------------------------------

namespace htsim::methdb {

struct CgmapFileCloser {
    void operator()(std::FILE *file) const noexcept
    {
        if (file != nullptr) {(void)std::fclose(file);}
    }
};

using CgmapFilePointer = std::unique_ptr<std::FILE, CgmapFileCloser>;

namespace {

inline constexpr std::size_t cgmap_spool_record_bytes = 12U;
inline constexpr std::uint8_t pending_bed_dinucleotide = 0xffU;

static_assert(sizeof(float) == sizeof(std::uint32_t));
static_assert(std::numeric_limits<float>::is_iec559);
static_assert(sizeof(off_t) >= sizeof(std::int64_t));
static_assert(std::numeric_limits<off_t>::is_signed);

std::vector<std::string_view> cgmap_split_fields(std::string_view line)
{
    std::vector<std::string_view> fields;
    fields.reserve(8U);
    std::size_t begin = 0U;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(
            begin,
            end == std::string_view::npos ? line.size() - begin : end - begin));
        if (end == std::string_view::npos) {break;}
        begin = end + 1U;
    }
    if (fields.size() != 8U) {
        throw CgmapProfileError(
            "CGmap row must contain exactly eight tab-separated fields");
    }
    return fields;
}

std::vector<std::string_view> bed_methyl_split_fields(std::string_view line)
{
    std::vector<std::string_view> fields;
    fields.reserve(18U);
    std::size_t begin = 0U;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(
            begin,
            end == std::string_view::npos ? line.size() - begin : end - begin));
        if (end == std::string_view::npos) {break;}
        begin = end + 1U;
    }
    if (fields.size() != 11U && fields.size() != 18U) {
        throw CgmapProfileError(
            "bedMethyl row must contain exactly eleven or eighteen tab-separated fields");
    }
    return fields;
}

bool bed_methyl_header(std::string_view line) noexcept
{
    return line.rfind("track ", 0U) == 0U
        || line.rfind("browser ", 0U) == 0U;
}

std::uint32_t bed_methyl_parse_u32(
    std::string_view text,
    const char *field)
{
    if (text.empty()) {
        throw CgmapProfileError(
            std::string("bedMethyl ") + field + " is empty");
    }
    std::uint32_t value = 0U;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()) {
        throw CgmapProfileError(
            std::string("bedMethyl ") + field + " is not a uint32 decimal");
    }
    return value;
}

float bed_methyl_parse_percent(std::string_view text)
{
    if (text.empty()) {
        throw CgmapProfileError("bedMethyl percent modified is empty");
    }
    double value = 0.0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()
        || !std::isfinite(value) || value < 0.0 || value > 100.0) {
        throw CgmapProfileError(
            "bedMethyl percent modified must be finite and in [0, 100]");
    }
    if (value == 0.0) {return 0.0F;}
    return static_cast<float>(value / 100.0);
}

void bed_methyl_validate_color(std::string_view text)
{
    if (text == "0") {return;}
    std::size_t begin = 0U;
    for (unsigned component = 0U; component < 3U; ++component) {
        const std::size_t end = text.find(',', begin);
        if ((component < 2U && end == std::string_view::npos)
            || (component == 2U && end != std::string_view::npos)) {
            throw CgmapProfileError(
                "bedMethyl itemRgb must be 0 or an R,G,B triple");
        }
        const std::string_view value = text.substr(
            begin,
            end == std::string_view::npos ? text.size() - begin : end - begin);
        const std::uint32_t parsed = bed_methyl_parse_u32(value, "itemRgb component");
        if (parsed > 255U) {
            throw CgmapProfileError(
                "bedMethyl itemRgb components must be in [0, 255]");
        }
        begin = end == std::string_view::npos ? text.size() : end + 1U;
    }
}

model::MethylationContext bed_methyl_parse_strand(std::string_view strand)
{
    if (strand == "+") {return model::MethylationContext::cg_c;}
    if (strand == "-") {return model::MethylationContext::cg_g;}
    if (strand == ".") {
        return static_cast<model::MethylationContext>(0U);
    }
    throw CgmapProfileError("bedMethyl strand must be +, -, or .");
}

void bed_methyl_validate_extended_counts(
    const std::vector<std::string_view> &fields,
    std::uint32_t coverage)
{
    if (fields.size() != 18U) {return;}
    std::array<std::uint32_t, 7U> counts = {};
    for (std::size_t index = 0U; index < counts.size(); ++index) {
        counts[index] = bed_methyl_parse_u32(
            fields[11U + index], "extended count");
    }
    const std::uint64_t valid = static_cast<std::uint64_t>(counts[0])
        + counts[1] + counts[2];
    if (valid != coverage) {
        throw CgmapProfileError(
            "bedMethyl modified, canonical, and other counts must sum to valid coverage");
    }
}

std::uint32_t cgmap_parse_u32(std::string_view text, const char *field)
{
    if (text.empty()) {
        throw CgmapProfileError(std::string("CGmap ") + field + " is empty");
    }
    std::uint32_t value = 0U;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()) {
        throw CgmapProfileError(
            std::string("CGmap ") + field + " is not a uint32 decimal");
    }
    return value;
}

std::optional<float> cgmap_parse_probability(std::string_view text)
{
    if (text == "na") {return std::nullopt;}
    if (text.empty()) {
        throw CgmapProfileError("CGmap methylation level is empty");
    }
    double value = 0.0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()
        || !std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw CgmapProfileError(
            "CGmap methylation level must be finite and in [0, 1], or na");
    }
    if (value == 0.0) {return 0.0F;}
    return static_cast<float>(value);
}

model::MethylationContext cgmap_parse_context(
    std::string_view nucleotide,
    std::string_view context)
{
    std::uint8_t value = 0U;
    if (context == "CG") {
        value = 1U;
    } else if (context == "CHG") {
        value = 3U;
    } else if (context == "CHH") {
        value = 7U;
    } else {
        throw CgmapProfileError("CGmap context must be CG, CHG, or CHH");
    }
    if (nucleotide == "G") {
        value = static_cast<std::uint8_t>(value | 8U);
    } else if (nucleotide != "C") {
        throw CgmapProfileError("CGmap nucleotide must be C or G");
    }
    return static_cast<model::MethylationContext>(value);
}

bool cgmap_valid_context(model::MethylationContext context) noexcept
{
    switch (context) {
    case model::MethylationContext::cg_c:
    case model::MethylationContext::chg_c:
    case model::MethylationContext::chh_c:
    case model::MethylationContext::cg_g:
    case model::MethylationContext::chg_g:
    case model::MethylationContext::chh_g:
        return true;
    }
    return false;
}

std::uint8_t cgmap_parse_dinucleotide(
    std::string_view value,
    model::MethylationContext context)
{
    std::uint8_t second = 0U;
    if (value == "CA") {
        second = 0U;
    } else if (value == "CC") {
        second = 1U;
    } else if (value == "CG") {
        second = 2U;
    } else if (value == "CT") {
        second = 3U;
    } else {
        throw CgmapProfileError(
            "CGmap dinucleotide must be CA, CC, CG, or CT");
    }
    const bool is_cpg =
        context == model::MethylationContext::cg_c
        || context == model::MethylationContext::cg_g;
    if ((second == 2U) != is_cpg) {
        throw CgmapProfileError(
            "CGmap dinucleotide disagrees with its context class");
    }
    return second;
}

void cgmap_put_u32(
    std::array<std::uint8_t, cgmap_spool_record_bytes> &bytes,
    std::size_t offset,
    std::uint32_t value) noexcept
{
    for (unsigned shift = 0U; shift < 32U; shift += 8U) {
        bytes[offset++] = static_cast<std::uint8_t>(value >> shift);
    }
}

std::uint32_t cgmap_get_u32(
    const std::array<std::uint8_t, cgmap_spool_record_bytes> &bytes,
    std::size_t offset) noexcept
{
    std::uint32_t value = 0U;
    for (unsigned shift = 0U; shift < 32U; shift += 8U) {
        value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
    }
    return value;
}

std::array<std::uint8_t, cgmap_spool_record_bytes> cgmap_encode_record(
    const CgmapRecord &record)
{
    std::array<std::uint8_t, cgmap_spool_record_bytes> bytes = {};
    cgmap_put_u32(bytes, 0U, record.reference_position);
    bytes[4] = static_cast<std::uint8_t>(record.context);
    bytes[5] = record.has_probability ? 1U : 0U;
    bytes[6] = record.dinucleotide_second;
    std::uint32_t probability_bits = 0U;
    std::memcpy(
        &probability_bits,
        &record.methylation_probability,
        sizeof(probability_bits));
    cgmap_put_u32(bytes, 8U, probability_bits);
    return bytes;
}

CgmapRecord cgmap_decode_record(
    const std::array<std::uint8_t, cgmap_spool_record_bytes> &bytes,
    MethylationProfileFormat format)
{
    const bool bed_methyl = format == MethylationProfileFormat::bed_methyl;
    const bool valid_context = bed_methyl
        ? (bytes[4] == 0U
            || bytes[4] == static_cast<std::uint8_t>(
                model::MethylationContext::cg_c)
            || bytes[4] == static_cast<std::uint8_t>(
                model::MethylationContext::cg_g))
        : cgmap_valid_context(
            static_cast<model::MethylationContext>(bytes[4]));
    if (bytes[5] > 1U || bytes[7] != 0U || !valid_context
        || (bed_methyl
            ? (bytes[5] != 1U || bytes[6] != pending_bed_dinucleotide)
            : bytes[6] > 3U)) {
        throw CgmapProfileError("CGmap spool record flags are corrupt");
    }
    CgmapRecord record;
    record.reference_position = cgmap_get_u32(bytes, 0U);
    record.context = static_cast<model::MethylationContext>(bytes[4]);
    record.has_probability = bytes[5] == 1U;
    record.dinucleotide_second = bytes[6];
    const std::uint32_t probability_bits = cgmap_get_u32(bytes, 8U);
    std::memcpy(
        &record.methylation_probability,
        &probability_bits,
        sizeof(probability_bits));
    if (!std::isfinite(record.methylation_probability)
        || record.methylation_probability < 0.0F
        || record.methylation_probability > 1.0F
        || (!record.has_probability && record.methylation_probability != 0.0F)) {
        throw CgmapProfileError("CGmap spool record is corrupt");
    }
    return record;
}

CgmapRecord normalize_bed_methyl_record(
    const model::Bases &contig_bases,
    CgmapRecord record)
{
    if (record.reference_position >= contig_bases.size()
        || record.dinucleotide_second != pending_bed_dinucleotide) {
        throw CgmapProfileError(
            "bedMethyl record is outside the contig or not normalized");
    }
    const auto observed = classify_context(
        contig_bases, record.reference_position, true);
    if (!observed) {
        throw CgmapProfileError(
            "bedMethyl target is not a resolved reference C/G context");
    }
    const std::uint8_t strand_hint = static_cast<std::uint8_t>(record.context);
    const bool observed_is_g = static_cast<std::uint8_t>(*observed) >= 8U;
    if ((strand_hint
            == static_cast<std::uint8_t>(model::MethylationContext::cg_c)
            && observed_is_g)
        || (strand_hint
            == static_cast<std::uint8_t>(model::MethylationContext::cg_g)
            && !observed_is_g)) {
        throw CgmapProfileError(
            "bedMethyl strand disagrees with the reference C/G base");
    }
    record.context = *observed;
    record.dinucleotide_second = observed_is_g
        ? static_cast<std::uint8_t>(
              3U - contig_bases[record.reference_position - 1U])
        : contig_bases[record.reference_position + 1U];
    return record;
}

std::string cgmap_io_error(const char *operation)
{
    return std::string("cannot ") + operation + " CGmap spool: "
        + std::strerror(errno);
}

} // namespace

void validate_cgmap_records(
    const model::Bases &contig_bases,
    const std::vector<CgmapRecord> &records)
{
    bool first = true;
    std::uint32_t previous = 0U;
    for (const CgmapRecord &record : records) {
        if ((!first && record.reference_position <= previous)
            || record.reference_position >= contig_bases.size()) {
            throw CgmapProfileError(
                "CGmap records are unsorted, duplicate, or outside the contig");
        }
        first = false;
        previous = record.reference_position;
        if (!cgmap_valid_context(record.context)
            || !std::isfinite(record.methylation_probability)
            || record.methylation_probability < 0.0F
            || record.methylation_probability > 1.0F
            || record.dinucleotide_second > 3U
            || (!record.has_probability
                && record.methylation_probability != 0.0F)) {
            throw CgmapProfileError("CGmap record value is invalid");
        }
        const auto observed = classify_context(
            contig_bases, record.reference_position, true);
        if (!observed || *observed != record.context) {
            throw CgmapProfileError(
                "CGmap nucleotide/context disagrees with the reference");
        }
        const bool is_g = static_cast<std::uint8_t>(record.context) >= 8U;
        const std::uint8_t expected_dinucleotide_second = is_g
            ? static_cast<std::uint8_t>(
                  3U - contig_bases[record.reference_position - 1U])
            : contig_bases[record.reference_position + 1U];
        if (record.dinucleotide_second != expected_dinucleotide_second) {
            throw CgmapProfileError(
                "CGmap dinucleotide disagrees with the reference");
        }
    }
}

class CgmapProfile::Impl {
public:
    Impl(
        const std::string &path,
        const std::vector<reference::ContigMetadata> &reference_catalog,
        MethylationProfileFormat format)
        : reference_catalog_(reference_catalog),
          row_counts_(reference_catalog.size(), 0U),
          format_(format)
    {
        try {
            if (format_ != MethylationProfileFormat::cgmap
                && format_ != MethylationProfileFormat::bed_methyl) {
                throw CgmapProfileError(
                    "unsupported methylation profile format");
            }
            const char *const label =
                format_ == MethylationProfileFormat::bed_methyl
                ? "bedMethyl"
                : "CGmap";
            if (reference_catalog_.empty()) {
                throw CgmapProfileError(
                    std::string(label) + " requires a non-empty reference catalog");
            }
            if (reference_catalog_.size()
                > std::numeric_limits<std::uint32_t>::max()) {
                throw CgmapProfileError(
                    "CGmap reference contig count exceeds uint32");
            }
            std::unordered_map<std::string, std::uint32_t> contig_indices;
            contig_indices.reserve(reference_catalog_.size());
            for (std::size_t index = 0U;
                 index < reference_catalog_.size();
                 ++index) {
                if (reference_catalog_[index].length
                        > std::numeric_limits<std::uint32_t>::max()) {
                    throw CgmapProfileError(
                        "CGmap reference catalog exceeds uint32 boundaries");
                }
                if (!contig_indices.emplace(
                        reference_catalog_[index].name,
                        static_cast<std::uint32_t>(index)).second) {
                    throw CgmapProfileError(
                        "CGmap reference catalog contains duplicate names");
                }
            }

            spool_.reset(std::tmpfile());
            if (!spool_) {throw CgmapProfileError(cgmap_io_error("create"));}

            text::TextSnapshot snapshot(path);
            file_sha256_ = snapshot.file_sha256();
            bool have_previous = false;
            std::uint32_t previous_contig = 0U;
            std::uint32_t previous_position = 0U;
            std::optional<std::size_t> bed_field_count;
            snapshot.visit_lines(
                [&](std::string_view line, std::uint64_t line_number) {
                    if (line.empty() || line.front() == '#'
                        || (format_ == MethylationProfileFormat::bed_methyl
                            && bed_methyl_header(line))) {
                        return;
                    }
                    try {
                        const auto fields =
                            format_ == MethylationProfileFormat::bed_methyl
                            ? bed_methyl_split_fields(line)
                            : cgmap_split_fields(line);
                        if (format_ == MethylationProfileFormat::bed_methyl) {
                            if (bed_field_count
                                && *bed_field_count != fields.size()) {
                                throw CgmapProfileError(
                                    "bedMethyl rows must use one consistent field count");
                            }
                            bed_field_count = fields.size();
                        }
                        const auto found = contig_indices.find(std::string(fields[0]));
                        if (found == contig_indices.end()) {
                            throw CgmapProfileError(
                                std::string(label) + " row names an unknown contig");
                        }
                        const std::uint32_t contig_index = found->second;
                        std::uint32_t position = 0U;
                        CgmapRecord record;
                        if (format_ == MethylationProfileFormat::bed_methyl) {
                            const std::uint32_t start = bed_methyl_parse_u32(
                                fields[1], "chromStart");
                            const std::uint32_t end = bed_methyl_parse_u32(
                                fields[2], "chromEnd");
                            if (start >= reference_catalog_[contig_index].length
                                || end != static_cast<std::uint64_t>(start) + 1U
                                || end > reference_catalog_[contig_index].length) {
                                throw CgmapProfileError(
                                    "bedMethyl target must be one in-range half-open base");
                            }
                            if (fields[3].empty()
                                || fields[3].find('\0') != std::string_view::npos) {
                                throw CgmapProfileError(
                                    "bedMethyl name must not be empty");
                            }
                            (void)bed_methyl_parse_u32(fields[4], "score");
                            const auto strand = bed_methyl_parse_strand(fields[5]);
                            const std::uint32_t thick_start =
                                bed_methyl_parse_u32(fields[6], "thickStart");
                            const std::uint32_t thick_end =
                                bed_methyl_parse_u32(fields[7], "thickEnd");
                            if (thick_start != start || thick_end != end) {
                                throw CgmapProfileError(
                                    "bedMethyl thickStart/thickEnd must match its target interval");
                            }
                            bed_methyl_validate_color(fields[8]);
                            const std::uint32_t coverage =
                                bed_methyl_parse_u32(fields[9], "valid coverage");
                            const float probability =
                                bed_methyl_parse_percent(fields[10]);
                            bed_methyl_validate_extended_counts(fields, coverage);
                            position = start;
                            record = CgmapRecord{
                                position,
                                probability,
                                strand,
                                true,
                                pending_bed_dinucleotide,
                            };
                        } else {
                            const std::uint32_t one_based_position =
                                cgmap_parse_u32(fields[2], "position");
                            if (one_based_position == 0U
                                || one_based_position
                                    > reference_catalog_[contig_index].length) {
                                throw CgmapProfileError(
                                    "CGmap position is outside its reference contig");
                            }
                            position = one_based_position - 1U;
                            const auto context =
                                cgmap_parse_context(fields[1], fields[3]);
                            const std::uint8_t dinucleotide_second =
                                cgmap_parse_dinucleotide(fields[4], context);
                            const auto probability =
                                cgmap_parse_probability(fields[5]);
                            const std::uint32_t methylated_count =
                                cgmap_parse_u32(fields[6], "methylated count");
                            const std::uint32_t total_count =
                                cgmap_parse_u32(fields[7], "total count");
                            if (methylated_count > total_count) {
                                throw CgmapProfileError(
                                    "CGmap methylated count exceeds total count");
                            }
                            record = CgmapRecord{
                                position,
                                probability.value_or(0.0F),
                                context,
                                probability.has_value(),
                                dinucleotide_second,
                            };
                        }
                        if (have_previous
                            && (contig_index < previous_contig
                                || (contig_index == previous_contig
                                    && position <= previous_position))) {
                            throw CgmapProfileError(
                                std::string(label)
                                + " rows must follow FASTA order with unique positions");
                        }
                        append(contig_index, record);
                        have_previous = true;
                        previous_contig = contig_index;
                        previous_position = position;
                    } catch (const CgmapProfileError &error) {
                        throw CgmapProfileError(
                            std::string(label) + " line "
                            + std::to_string(line_number)
                            + ": " + error.what());
                    }
                });
            if (row_count_ == 0U) {
                throw CgmapProfileError(
                    std::string(label) + " contains no data rows");
            }
            if (std::fflush(spool_.get()) != 0) {
                throw CgmapProfileError(cgmap_io_error("flush"));
            }
            first_record_.reserve(row_counts_.size() + 1U);
            first_record_.push_back(0U);
            for (const std::uint32_t count : row_counts_) {
                first_record_.push_back(first_record_.back() + count);
            }
            if (first_record_.back() != row_count_
                || snapshot.file_sha256() != file_sha256_) {
                throw CgmapProfileError("CGmap spool accounting or digest changed");
            }
        } catch (const CgmapProfileError &) {
            throw;
        } catch (const std::exception &error) {
            throw CgmapProfileError(error.what());
        }
    }

    const crypto::Sha256Digest &file_sha256() const noexcept {return file_sha256_;}
    std::uint64_t row_count() const noexcept {return row_count_;}
    std::uint64_t defined_probability_count() const noexcept
    {
        return defined_probability_count_;
    }

    std::vector<CgmapRecord> records(std::uint32_t contig_index) const
    {
        if (contig_index >= row_counts_.size()) {
            throw CgmapProfileError("CGmap contig index is out of range");
        }
        const std::uint64_t first = first_record_[contig_index];
        if (first > static_cast<std::uint64_t>(
                        std::numeric_limits<off_t>::max())
                / cgmap_spool_record_bytes) {
            throw CgmapProfileError("CGmap spool offset exceeds off_t");
        }
        const std::uint64_t byte_offset = first * cgmap_spool_record_bytes;
        if (::fseeko(
                spool_.get(), static_cast<off_t>(byte_offset), SEEK_SET) != 0) {
            throw CgmapProfileError(cgmap_io_error("seek"));
        }
        const std::uint32_t count = row_counts_[contig_index];
        std::vector<CgmapRecord> result;
        result.reserve(count);
        for (std::uint32_t index = 0U; index < count; ++index) {
            std::array<std::uint8_t, cgmap_spool_record_bytes> bytes = {};
            if (std::fread(bytes.data(), bytes.size(), 1U, spool_.get()) != 1U) {
                throw CgmapProfileError(
                    "CGmap spool is truncated or unreadable");
            }
            result.push_back(cgmap_decode_record(bytes, format_));
        }
        return result;
    }

    std::vector<CgmapRecord> records(const reference::Contig &contig) const
    {
        if (contig.index >= reference_catalog_.size()) {
            throw CgmapProfileError("CGmap contig index is out of range");
        }
        const reference::ContigMetadata &expected =
            reference_catalog_[contig.index];
        if (contig.name != expected.name || contig.length != expected.length
            || contig.reference_sha256 != expected.reference_sha256
            || contig.length != contig.bases.size()) {
            throw CgmapProfileError(
                "CGmap validation contig disagrees with the reference catalog");
        }
        std::vector<CgmapRecord> result = records(contig.index);
        if (format_ == MethylationProfileFormat::bed_methyl) {
            for (CgmapRecord &record : result) {
                record = normalize_bed_methyl_record(contig.bases, record);
            }
        }
        validate_cgmap_records(contig.bases, result);
        return result;
    }

private:
    void append(std::uint32_t contig_index, const CgmapRecord &record)
    {
        if (row_counts_[contig_index]
            == std::numeric_limits<std::uint32_t>::max()) {
            throw CgmapProfileError("CGmap per-contig row count exceeds uint32");
        }
        const std::uint64_t maximum_records =
            static_cast<std::uint64_t>(std::numeric_limits<off_t>::max())
            / cgmap_spool_record_bytes;
        if (row_count_ >= maximum_records) {
            throw CgmapProfileError("CGmap spool size exceeds off_t");
        }
        const auto bytes = cgmap_encode_record(record);
        if (std::fwrite(bytes.data(), bytes.size(), 1U, spool_.get()) != 1U) {
            throw CgmapProfileError(cgmap_io_error("write"));
        }
        ++row_counts_[contig_index];
        ++row_count_;
        if (record.has_probability) {++defined_probability_count_;}
    }

    crypto::Sha256Digest file_sha256_ = {};
    std::vector<reference::ContigMetadata> reference_catalog_;
    std::vector<std::uint32_t> row_counts_;
    std::vector<std::uint64_t> first_record_;
    CgmapFilePointer spool_;
    std::uint64_t row_count_ = 0U;
    std::uint64_t defined_probability_count_ = 0U;
    MethylationProfileFormat format_ = MethylationProfileFormat::cgmap;
};

CgmapProfile::CgmapProfile(
    const std::string &path,
    const std::vector<reference::ContigMetadata> &reference_catalog,
    MethylationProfileFormat format)
    : impl_(std::make_unique<Impl>(path, reference_catalog, format))
{}

CgmapProfile::~CgmapProfile() = default;

const crypto::Sha256Digest &CgmapProfile::file_sha256() const noexcept
{
    return impl_->file_sha256();
}

std::uint64_t CgmapProfile::row_count() const noexcept
{
    return impl_->row_count();
}

std::uint64_t CgmapProfile::defined_probability_count() const noexcept
{
    return impl_->defined_probability_count();
}

std::vector<CgmapRecord> CgmapProfile::records(
    const reference::Contig &contig) const
{
    return impl_->records(contig);
}

void CgmapProfile::validate_contig(const reference::Contig &contig) const
{
    (void)impl_->records(contig);
}

} // namespace htsim::methdb

// ---- asm_profile --------------------------------------------------------

namespace htsim::methdb {

struct AsmFileCloser {
    void operator()(std::FILE *file) const noexcept
    {
        if (file != nullptr) {(void)std::fclose(file);}
    }
};

using AsmFilePointer = std::unique_ptr<std::FILE, AsmFileCloser>;

namespace {

inline constexpr std::size_t asm_spool_record_bytes = 20U;

static_assert(sizeof(float) == sizeof(std::uint32_t));
static_assert(std::numeric_limits<float>::is_iec559);
static_assert(sizeof(off_t) >= sizeof(std::int64_t));
static_assert(std::numeric_limits<off_t>::is_signed);

std::vector<std::string_view> asm_split_fields(std::string_view line)
{
    std::vector<std::string_view> fields;
    fields.reserve(14U);
    std::size_t begin = 0U;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(
            begin,
            end == std::string_view::npos ? line.size() - begin : end - begin));
        if (end == std::string_view::npos) {break;}
        begin = end + 1U;
    }
    if (fields.size() != 14U) {
        throw AsmProfileError(
            "ASM row must contain exactly fourteen tab-separated fields");
    }
    return fields;
}

std::vector<std::string_view> asm_bed_split_fields(std::string_view line)
{
    std::vector<std::string_view> fields;
    fields.reserve(12U);
    std::size_t begin = 0U;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(
            begin,
            end == std::string_view::npos ? line.size() - begin : end - begin));
        if (end == std::string_view::npos) {break;}
        begin = end + 1U;
    }
    if (fields.size() != 12U) {
        throw AsmProfileError(
            "ASM BED row must contain exactly twelve tab-separated fields");
    }
    return fields;
}

bool asm_bed_header(std::string_view line) noexcept
{
    return line.rfind("track ", 0U) == 0U
        || line.rfind("browser ", 0U) == 0U;
}

std::uint32_t asm_parse_u32(std::string_view text, const char *field)
{
    if (text.empty()) {
        throw AsmProfileError(std::string("ASM BED ") + field + " is empty");
    }
    std::uint32_t value = 0U;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()) {
        throw AsmProfileError(
            std::string("ASM BED ") + field + " is not a uint32 decimal");
    }
    return value;
}

model::MethylationContext asm_bed_parse_strand(std::string_view strand)
{
    if (strand == "+") {return model::MethylationContext::cg_c;}
    if (strand == "-") {return model::MethylationContext::cg_g;}
    if (strand == ".") {
        return static_cast<model::MethylationContext>(0U);
    }
    throw AsmProfileError("ASM BED strand must be +, -, or .");
}

std::uint32_t asm_parse_positive_u32(std::string_view text, const char *field)
{
    if (text.empty()) {
        throw AsmProfileError(std::string("ASM ") + field + " is empty");
    }
    std::uint32_t value = 0U;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()
        || value == 0U) {
        throw AsmProfileError(
            std::string("ASM ") + field
            + " must be a positive uint32 decimal");
    }
    return value;
}

std::optional<float> asm_parse_probability(
    std::string_view text,
    const char *field,
    bool allow_na)
{
    if (allow_na && text == "na") {return std::nullopt;}
    if (text.empty()) {
        throw AsmProfileError(std::string("ASM ") + field + " is empty");
    }
    double value = 0.0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()
        || !std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw AsmProfileError(
            std::string("ASM ") + field
            + " must be finite and in [0, 1]"
            + (allow_na ? ", or na" : ""));
    }
    if (value == 0.0) {return 0.0F;}
    return static_cast<float>(value);
}

void asm_validate_optional_finite(std::string_view text, const char *field)
{
    if (text == "na") {return;}
    if (text.empty()) {
        throw AsmProfileError(std::string("ASM ") + field + " is empty");
    }
    double value = 0.0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()
        || !std::isfinite(value)) {
        throw AsmProfileError(
            std::string("ASM ") + field + " must be finite, or na");
    }
}

std::uint8_t asm_parse_base(std::string_view text, const char *field)
{
    if (text == "A") {return 0U;}
    if (text == "C") {return 1U;}
    if (text == "G") {return 2U;}
    if (text == "T") {return 3U;}
    throw AsmProfileError(
        std::string("ASM ") + field + " must be one uppercase A/C/G/T base");
}

model::MethylationContext asm_parse_context(
    std::string_view nucleotide,
    std::string_view context)
{
    std::uint8_t value = 0U;
    if (context == "CG") {
        value = 1U;
    } else if (context == "CHG") {
        value = 3U;
    } else if (context == "CHH") {
        value = 7U;
    } else {
        throw AsmProfileError("ASM context must be CG, CHG, or CHH");
    }
    if (nucleotide == "G") {
        value = static_cast<std::uint8_t>(value | 8U);
    } else if (nucleotide != "C") {
        throw AsmProfileError("ASM nucleotide must be C or G");
    }
    return static_cast<model::MethylationContext>(value);
}

bool asm_valid_context(model::MethylationContext context) noexcept
{
    switch (context) {
    case model::MethylationContext::cg_c:
    case model::MethylationContext::chg_c:
    case model::MethylationContext::chh_c:
    case model::MethylationContext::cg_g:
    case model::MethylationContext::chg_g:
    case model::MethylationContext::chh_g:
        return true;
    }
    return false;
}

std::uint8_t asm_parse_dinucleotide(
    std::string_view value,
    model::MethylationContext context)
{
    std::uint8_t second = 0U;
    if (value == "CA") {
        second = 0U;
    } else if (value == "CC") {
        second = 1U;
    } else if (value == "CG") {
        second = 2U;
    } else if (value == "CT") {
        second = 3U;
    } else {
        throw AsmProfileError(
            "ASM dinucleotide must be CA, CC, CG, or CT");
    }
    const bool is_cpg =
        context == model::MethylationContext::cg_c
        || context == model::MethylationContext::cg_g;
    if ((second == 2U) != is_cpg) {
        throw AsmProfileError(
            "ASM dinucleotide disagrees with its context class");
    }
    return second;
}

void asm_put_u32(
    std::array<std::uint8_t, asm_spool_record_bytes> &bytes,
    std::size_t offset,
    std::uint32_t value) noexcept
{
    for (unsigned shift = 0U; shift < 32U; shift += 8U) {
        bytes[offset++] = static_cast<std::uint8_t>(value >> shift);
    }
}

std::uint32_t asm_get_u32(
    const std::array<std::uint8_t, asm_spool_record_bytes> &bytes,
    std::size_t offset) noexcept
{
    std::uint32_t value = 0U;
    for (unsigned shift = 0U; shift < 32U; shift += 8U) {
        value |= static_cast<std::uint32_t>(bytes[offset++]) << shift;
    }
    return value;
}

std::array<std::uint8_t, asm_spool_record_bytes> asm_encode_record(
    const AsmRecord &record)
{
    std::array<std::uint8_t, asm_spool_record_bytes> bytes = {};
    asm_put_u32(bytes, 0U, record.target_reference_position);
    asm_put_u32(bytes, 4U, record.linked_variant_position);
    std::uint32_t reference_probability_bits = 0U;
    std::uint32_t alternate_probability_bits = 0U;
    std::memcpy(
        &reference_probability_bits,
        &record.reference_methylation_probability,
        sizeof(reference_probability_bits));
    std::memcpy(
        &alternate_probability_bits,
        &record.alternate_methylation_probability,
        sizeof(alternate_probability_bits));
    asm_put_u32(bytes, 8U, reference_probability_bits);
    asm_put_u32(bytes, 12U, alternate_probability_bits);
    bytes[16] = static_cast<std::uint8_t>(record.context);
    bytes[17] = record.dinucleotide_second;
    bytes[18] = record.linked_reference_base;
    bytes[19] = record.linked_alternate_base;
    return bytes;
}

AsmRecord asm_decode_record(
    const std::array<std::uint8_t, asm_spool_record_bytes> &bytes,
    AsmProfileFormat format)
{
    AsmRecord record;
    record.target_reference_position = asm_get_u32(bytes, 0U);
    record.linked_variant_position = asm_get_u32(bytes, 4U);
    const std::uint32_t reference_probability_bits = asm_get_u32(bytes, 8U);
    const std::uint32_t alternate_probability_bits = asm_get_u32(bytes, 12U);
    std::memcpy(
        &record.reference_methylation_probability,
        &reference_probability_bits,
        sizeof(reference_probability_bits));
    std::memcpy(
        &record.alternate_methylation_probability,
        &alternate_probability_bits,
        sizeof(alternate_probability_bits));
    record.context = static_cast<model::MethylationContext>(bytes[16]);
    record.dinucleotide_second = bytes[17];
    record.linked_reference_base = bytes[18];
    record.linked_alternate_base = bytes[19];
    const bool bed = format == AsmProfileFormat::bed;
    const std::uint8_t context_value = static_cast<std::uint8_t>(record.context);
    const bool valid_context = bed
        ? (context_value == 0U
            || context_value == static_cast<std::uint8_t>(
                model::MethylationContext::cg_c)
            || context_value == static_cast<std::uint8_t>(
                model::MethylationContext::cg_g))
        : asm_valid_context(record.context);
    if (!valid_context
        || !std::isfinite(record.reference_methylation_probability)
        || !std::isfinite(record.alternate_methylation_probability)
        || record.reference_methylation_probability < 0.0F
        || record.reference_methylation_probability > 1.0F
        || record.alternate_methylation_probability < 0.0F
        || record.alternate_methylation_probability > 1.0F
        || (bed
            ? record.dinucleotide_second != pending_bed_dinucleotide
            : record.dinucleotide_second > 3U)
        || record.linked_reference_base > 3U
        || record.linked_alternate_base > 3U
        || record.linked_reference_base == record.linked_alternate_base) {
        throw AsmProfileError("ASM spool record is corrupt");
    }
    return record;
}

AsmRecord normalize_asm_bed_record(
    const model::Bases &contig_bases,
    AsmRecord record)
{
    if (record.target_reference_position >= contig_bases.size()
        || record.linked_variant_position >= contig_bases.size()
        || record.dinucleotide_second != pending_bed_dinucleotide) {
        throw AsmProfileError(
            "ASM BED record is outside the contig or not normalized");
    }
    const auto observed = classify_context(
        contig_bases, record.target_reference_position, true);
    if (!observed) {
        throw AsmProfileError(
            "ASM BED target is not a resolved reference C/G context");
    }
    const std::uint8_t strand_hint = static_cast<std::uint8_t>(record.context);
    const bool observed_is_g = static_cast<std::uint8_t>(*observed) >= 8U;
    if ((strand_hint
            == static_cast<std::uint8_t>(model::MethylationContext::cg_c)
            && observed_is_g)
        || (strand_hint
            == static_cast<std::uint8_t>(model::MethylationContext::cg_g)
            && !observed_is_g)) {
        throw AsmProfileError(
            "ASM BED strand disagrees with the reference C/G base");
    }
    record.context = *observed;
    record.dinucleotide_second = observed_is_g
        ? static_cast<std::uint8_t>(
              3U - contig_bases[record.target_reference_position - 1U])
        : contig_bases[record.target_reference_position + 1U];
    return record;
}

std::string asm_io_error(const char *operation)
{
    return std::string("cannot ") + operation + " ASM spool: "
        + std::strerror(errno);
}

} // namespace

void validate_asm_records(
    const model::Bases &contig_bases,
    const std::vector<AsmRecord> &records)
{
    bool first = true;
    std::uint32_t previous = 0U;
    for (const AsmRecord &record : records) {
        if ((!first && record.target_reference_position <= previous)
            || record.target_reference_position >= contig_bases.size()
            || record.linked_variant_position >= contig_bases.size()) {
            throw AsmProfileError(
                "ASM records are unsorted, duplicate, or outside the contig");
        }
        first = false;
        previous = record.target_reference_position;
        if (!asm_valid_context(record.context)
            || !std::isfinite(record.reference_methylation_probability)
            || !std::isfinite(record.alternate_methylation_probability)
            || record.reference_methylation_probability < 0.0F
            || record.reference_methylation_probability > 1.0F
            || record.alternate_methylation_probability < 0.0F
            || record.alternate_methylation_probability > 1.0F
            || record.dinucleotide_second > 3U
            || record.linked_reference_base > 3U
            || record.linked_alternate_base > 3U
            || record.linked_reference_base == record.linked_alternate_base) {
            throw AsmProfileError("ASM record value is invalid");
        }
        const auto observed = classify_context(
            contig_bases, record.target_reference_position, true);
        if (!observed || *observed != record.context) {
            throw AsmProfileError(
                "ASM nucleotide/context disagrees with the reference");
        }
        const bool is_g = static_cast<std::uint8_t>(record.context) >= 8U;
        const std::uint8_t expected_dinucleotide_second = is_g
            ? static_cast<std::uint8_t>(
                  3U - contig_bases[record.target_reference_position - 1U])
            : contig_bases[record.target_reference_position + 1U];
        if (record.dinucleotide_second != expected_dinucleotide_second) {
            throw AsmProfileError(
                "ASM dinucleotide disagrees with the reference");
        }
        if (record.linked_reference_base
            != contig_bases[record.linked_variant_position]) {
            throw AsmProfileError(
                "ASM linked REF base disagrees with the reference");
        }
    }
}

class AsmProfile::Impl {
public:
    Impl(
        const std::string &path,
        const std::vector<reference::ContigMetadata> &reference_catalog,
        AsmProfileFormat format)
        : reference_catalog_(reference_catalog),
          row_counts_(reference_catalog.size(), 0U),
          format_(format)
    {
        try {
            if (format_ != AsmProfileFormat::htsim
                && format_ != AsmProfileFormat::bed) {
                throw AsmProfileError("unsupported ASM profile format");
            }
            const char *const label = format_ == AsmProfileFormat::bed
                ? "ASM BED"
                : "ASM";
            if (reference_catalog_.empty()) {
                throw AsmProfileError(
                    std::string(label) + " requires a non-empty reference catalog");
            }
            if (reference_catalog_.size()
                > std::numeric_limits<std::uint32_t>::max()) {
                throw AsmProfileError(
                    "ASM reference contig count exceeds uint32");
            }
            std::unordered_map<std::string, std::uint32_t> contig_indices;
            contig_indices.reserve(reference_catalog_.size());
            for (std::size_t index = 0U;
                 index < reference_catalog_.size();
                 ++index) {
                if (reference_catalog_[index].length
                        > std::numeric_limits<std::uint32_t>::max()) {
                    throw AsmProfileError(
                        "ASM reference catalog exceeds uint32 boundaries");
                }
                if (!contig_indices.emplace(
                        reference_catalog_[index].name,
                        static_cast<std::uint32_t>(index)).second) {
                    throw AsmProfileError(
                        "ASM reference catalog contains duplicate names");
                }
            }

            spool_.reset(std::tmpfile());
            if (!spool_) {throw AsmProfileError(asm_io_error("create"));}

            text::TextSnapshot snapshot(path);
            file_sha256_ = snapshot.file_sha256();
            bool have_previous = false;
            std::uint32_t previous_contig = 0U;
            std::uint32_t previous_position = 0U;
            snapshot.visit_lines(
                [&](std::string_view line, std::uint64_t line_number) {
                    if (line.empty() || line.front() == '#'
                        || (format_ == AsmProfileFormat::bed
                            && asm_bed_header(line))) {
                        return;
                    }
                    try {
                        const auto fields = format_ == AsmProfileFormat::bed
                            ? asm_bed_split_fields(line)
                            : asm_split_fields(line);
                        const auto found = contig_indices.find(std::string(fields[0]));
                        if (found == contig_indices.end()) {
                            throw AsmProfileError(
                                std::string(label) + " row names an unknown contig");
                        }
                        const std::uint32_t contig_index = found->second;
                        std::uint32_t target_position = 0U;
                        std::uint32_t variant_position = 0U;
                        AsmRecord record;
                        if (format_ == AsmProfileFormat::bed) {
                            const std::uint32_t target_start =
                                asm_parse_u32(fields[1], "chromStart");
                            const std::uint32_t target_end =
                                asm_parse_u32(fields[2], "chromEnd");
                            const std::uint32_t variant_start =
                                asm_parse_u32(fields[6], "linkedStart");
                            const std::uint32_t variant_end =
                                asm_parse_u32(fields[7], "linkedEnd");
                            const std::uint64_t contig_length =
                                reference_catalog_[contig_index].length;
                            if (target_start >= contig_length
                                || target_end
                                    != static_cast<std::uint64_t>(target_start) + 1U
                                || target_end > contig_length
                                || variant_start >= contig_length
                                || variant_end
                                    != static_cast<std::uint64_t>(variant_start) + 1U
                                || variant_end > contig_length) {
                                throw AsmProfileError(
                                    "ASM BED target and linked SNV must be in-range one-base intervals");
                            }
                            if (fields[3].empty()
                                || fields[3].find('\0') != std::string_view::npos) {
                                throw AsmProfileError(
                                    "ASM BED name must not be empty");
                            }
                            const std::uint32_t score =
                                asm_parse_u32(fields[4], "score");
                            if (score > 1000U) {
                                throw AsmProfileError(
                                    "ASM BED score must be in [0, 1000]");
                            }
                            const auto context =
                                asm_bed_parse_strand(fields[5]);
                            const std::uint8_t linked_reference_base =
                                asm_parse_base(fields[8], "linked REF");
                            const std::uint8_t linked_alternate_base =
                                asm_parse_base(fields[9], "linked ALT");
                            if (linked_reference_base == linked_alternate_base) {
                                throw AsmProfileError(
                                    "ASM BED linked REF and ALT bases must differ");
                            }
                            const float reference_probability =
                                *asm_parse_probability(
                                    fields[10],
                                    "reference methylation level",
                                    false);
                            const float alternate_probability =
                                *asm_parse_probability(
                                    fields[11],
                                    "alternate methylation level",
                                    false);
                            target_position = target_start;
                            variant_position = variant_start;
                            record = AsmRecord{
                                target_position,
                                variant_position,
                                reference_probability,
                                alternate_probability,
                                context,
                                pending_bed_dinucleotide,
                                linked_reference_base,
                                linked_alternate_base,
                            };
                        } else {
                            const std::uint32_t one_based_target =
                                asm_parse_positive_u32(
                                    fields[2], "target position");
                            const std::uint32_t one_based_variant =
                                asm_parse_positive_u32(
                                    fields[6], "linked variant position");
                            if (one_based_target
                                    > reference_catalog_[contig_index].length
                                || one_based_variant
                                    > reference_catalog_[contig_index].length) {
                                throw AsmProfileError(
                                    "ASM position is outside its reference contig");
                            }
                            target_position = one_based_target - 1U;
                            variant_position = one_based_variant - 1U;
                            const auto context =
                                asm_parse_context(fields[1], fields[3]);
                            const std::uint8_t dinucleotide_second =
                                asm_parse_dinucleotide(fields[4], context);
                            (void)asm_parse_probability(
                                fields[5], "total methylation level", true);
                            const std::uint8_t linked_reference_base =
                                asm_parse_base(fields[7], "linked REF");
                            const std::uint8_t linked_alternate_base =
                                asm_parse_base(fields[8], "linked ALT");
                            if (linked_reference_base == linked_alternate_base) {
                                throw AsmProfileError(
                                    "ASM linked REF and ALT bases must differ");
                            }
                            const float reference_probability =
                                *asm_parse_probability(
                                    fields[9],
                                    "reference methylation level",
                                    false);
                            const float alternate_probability =
                                *asm_parse_probability(
                                    fields[10],
                                    "alternate methylation level",
                                    false);
                            asm_validate_optional_finite(
                                fields[11], "fold change");
                            (void)asm_parse_probability(
                                fields[12], "p-value", true);
                            if (fields[13].empty()
                                || fields[13].find('\0')
                                    != std::string_view::npos) {
                                throw AsmProfileError(
                                    "ASM comment must not be empty");
                            }
                            record = AsmRecord{
                                target_position,
                                variant_position,
                                reference_probability,
                                alternate_probability,
                                context,
                                dinucleotide_second,
                                linked_reference_base,
                                linked_alternate_base,
                            };
                        }
                        if (have_previous
                            && (contig_index < previous_contig
                                || (contig_index == previous_contig
                                    && target_position <= previous_position))) {
                            throw AsmProfileError(
                                std::string(label)
                                + " rows must follow FASTA order with unique target positions");
                        }
                        append(contig_index, record);
                        have_previous = true;
                        previous_contig = contig_index;
                        previous_position = target_position;
                    } catch (const AsmProfileError &error) {
                        throw AsmProfileError(
                            std::string(label) + " line "
                            + std::to_string(line_number)
                            + ": " + error.what());
                    }
                });
            if (row_count_ == 0U) {
                throw AsmProfileError(
                    std::string(label) + " contains no data rows");
            }
            if (std::fflush(spool_.get()) != 0) {
                throw AsmProfileError(asm_io_error("flush"));
            }
            first_record_.reserve(row_counts_.size() + 1U);
            first_record_.push_back(0U);
            for (const std::uint32_t count : row_counts_) {
                first_record_.push_back(first_record_.back() + count);
            }
            if (first_record_.back() != row_count_
                || snapshot.file_sha256() != file_sha256_) {
                throw AsmProfileError("ASM spool accounting or digest changed");
            }
        } catch (const AsmProfileError &) {
            throw;
        } catch (const std::exception &error) {
            throw AsmProfileError(error.what());
        }
    }

    const crypto::Sha256Digest &file_sha256() const noexcept {return file_sha256_;}
    std::uint64_t row_count() const noexcept {return row_count_;}

    std::vector<AsmRecord> records(std::uint32_t contig_index) const
    {
        if (contig_index >= row_counts_.size()) {
            throw AsmProfileError("ASM contig index is out of range");
        }
        const std::uint64_t first = first_record_[contig_index];
        if (first > static_cast<std::uint64_t>(
                        std::numeric_limits<off_t>::max())
                / asm_spool_record_bytes) {
            throw AsmProfileError("ASM spool offset exceeds off_t");
        }
        const std::uint64_t byte_offset = first * asm_spool_record_bytes;
        if (::fseeko(
                spool_.get(), static_cast<off_t>(byte_offset), SEEK_SET) != 0) {
            throw AsmProfileError(asm_io_error("seek"));
        }
        const std::uint32_t count = row_counts_[contig_index];
        std::vector<AsmRecord> result;
        result.reserve(count);
        for (std::uint32_t index = 0U; index < count; ++index) {
            std::array<std::uint8_t, asm_spool_record_bytes> bytes = {};
            if (std::fread(bytes.data(), bytes.size(), 1U, spool_.get()) != 1U) {
                throw AsmProfileError("ASM spool is truncated or unreadable");
            }
            result.push_back(asm_decode_record(bytes, format_));
        }
        return result;
    }

    std::vector<AsmRecord> records(const reference::Contig &contig) const
    {
        if (contig.index >= reference_catalog_.size()) {
            throw AsmProfileError("ASM contig index is out of range");
        }
        const reference::ContigMetadata &expected =
            reference_catalog_[contig.index];
        if (contig.name != expected.name || contig.length != expected.length
            || contig.reference_sha256 != expected.reference_sha256
            || contig.length != contig.bases.size()) {
            throw AsmProfileError(
                "ASM validation contig disagrees with the reference catalog");
        }
        std::vector<AsmRecord> result = records(contig.index);
        if (format_ == AsmProfileFormat::bed) {
            for (AsmRecord &record : result) {
                record = normalize_asm_bed_record(contig.bases, record);
            }
        }
        validate_asm_records(contig.bases, result);
        return result;
    }

private:
    void append(std::uint32_t contig_index, const AsmRecord &record)
    {
        if (row_counts_[contig_index]
            == std::numeric_limits<std::uint32_t>::max()) {
            throw AsmProfileError("ASM per-contig row count exceeds uint32");
        }
        const std::uint64_t maximum_records =
            static_cast<std::uint64_t>(std::numeric_limits<off_t>::max())
            / asm_spool_record_bytes;
        if (row_count_ >= maximum_records) {
            throw AsmProfileError("ASM spool size exceeds off_t");
        }
        const auto bytes = asm_encode_record(record);
        if (std::fwrite(bytes.data(), bytes.size(), 1U, spool_.get()) != 1U) {
            throw AsmProfileError(asm_io_error("write"));
        }
        ++row_counts_[contig_index];
        ++row_count_;
    }

    crypto::Sha256Digest file_sha256_ = {};
    std::vector<reference::ContigMetadata> reference_catalog_;
    std::vector<std::uint32_t> row_counts_;
    std::vector<std::uint64_t> first_record_;
    AsmFilePointer spool_;
    std::uint64_t row_count_ = 0U;
    AsmProfileFormat format_ = AsmProfileFormat::htsim;
};

AsmProfile::AsmProfile(
    const std::string &path,
    const std::vector<reference::ContigMetadata> &reference_catalog,
    AsmProfileFormat format)
    : impl_(std::make_unique<Impl>(path, reference_catalog, format))
{}

AsmProfile::~AsmProfile() = default;

const crypto::Sha256Digest &AsmProfile::file_sha256() const noexcept
{
    return impl_->file_sha256();
}

std::uint64_t AsmProfile::row_count() const noexcept
{
    return impl_->row_count();
}

std::vector<AsmRecord> AsmProfile::records(
    const reference::Contig &contig) const
{
    return impl_->records(contig);
}

void AsmProfile::validate_contig(const reference::Contig &contig) const
{
    (void)impl_->records(contig);
}

} // namespace htsim::methdb

// ---- beta_sampler --------------------------------------------------------

namespace htsim::beta_sampler {
namespace {

static_assert(std::numeric_limits<double>::is_iec559
                  && std::numeric_limits<double>::radix == 2
                  && std::numeric_limits<double>::digits == 53,
              "beta sampler requires IEEE-754 binary64");
static_assert(std::numeric_limits<float>::is_iec559
                  && std::numeric_limits<float>::radix == 2
                  && std::numeric_limits<float>::digits == 24,
              "beta sampler requires IEEE-754 binary32 output");

constexpr std::uint64_t role_shift = 62;
constexpr std::uint64_t role_offset_mask =
    (UINT64_C(1) << role_shift) - UINT64_C(1);
constexpr double inverse_uniform_scale = 0x1p-53;

enum class Role : std::uint64_t {
    alpha_candidate = 0,
    alpha_boost = 1,
    beta_candidate = 2,
    beta_boost = 3,
};

std::uint64_t local_index(Role role, std::uint64_t offset)
{
    if (offset > role_offset_mask) {
        throw SamplingError("beta sampler counter offset overflow");
    }
    return (static_cast<std::uint64_t>(role) << role_shift) | offset;
}

double uniform_open_closed(
    std::uint64_t key,
    std::uint64_t entity,
    Role role,
    std::uint64_t offset,
    std::size_t pair)
{
    const std::uint64_t significand =
        rng::u64(key, entity, local_index(role, offset), pair) >> 11;
    return static_cast<double>(significand + UINT64_C(1))
        * inverse_uniform_scale;
}

double gamma_at_least_one(
    std::uint64_t key,
    std::uint64_t entity,
    double shape,
    Role candidate_role,
    std::uint32_t maximum_attempts)
{
    const double d = shape - (1.0 / 3.0);
    // sqrt(d) avoids an otherwise overflowing intermediate at shapes > DBL_MAX/9.
    const double c = (1.0 / 3.0) / std::sqrt(d);
    if (!(d > 0.0) || !std::isfinite(d) || !std::isfinite(c)) {
        throw SamplingError("Gamma setup produced a non-finite value");
    }

    for (std::uint32_t attempt = 0; attempt < maximum_attempts; ++attempt) {
        const std::uint64_t offset = static_cast<std::uint64_t>(attempt) * 2U;
        const double normal = normal_sampler::standard_normal(
            key, entity, local_index(candidate_role, offset));
        const double transformed = 1.0 + c * normal;
        if (!(transformed > 0.0)) {continue;}
        const double transformed_squared = transformed * transformed;
        const double v = transformed_squared * transformed;
        if (!(v > 0.0) || !std::isfinite(v)) {continue;}

        const double acceptance_uniform =
            uniform_open_closed(key, entity, candidate_role, offset + 1U, 0);
        const double normal_squared = normal * normal;
        const double normal_fourth = normal_squared * normal_squared;
        bool accepted =
            acceptance_uniform < 1.0 - 0.0331 * normal_fourth;
        if (!accepted) {
            const double log_threshold =
                0.5 * normal_squared
                + d * ((1.0 - v) + std::log(v));
            accepted = std::log(acceptance_uniform) < log_threshold;
        }
        if (!accepted) {continue;}

        const double result = d * v;
        if (!std::isfinite(result) || !(result > 0.0)) {
            throw SamplingError("accepted Gamma draw is not finite and positive");
        }
        return result;
    }
    throw SamplingError("Gamma rejection sampler exhausted its attempt cap");
}

double gamma_draw(
    std::uint64_t key,
    std::uint64_t entity,
    double shape,
    Role candidate_role,
    Role boost_role,
    std::uint32_t maximum_attempts)
{
    if (shape >= 1.0) {
        return gamma_at_least_one(
            key, entity, shape, candidate_role, maximum_attempts);
    }

    const double gamma = gamma_at_least_one(
        key, entity, shape + 1.0, candidate_role, maximum_attempts);
    const double boost_uniform =
        uniform_open_closed(key, entity, boost_role, 0, 0);
    const double boost = std::pow(boost_uniform, 1.0 / shape);
    if (!std::isfinite(boost) || boost < 0.0 || boost > 1.0) {
        throw SamplingError("Gamma shape boost produced an invalid value");
    }
    const double result = gamma * boost;
    if (!std::isfinite(result) || result < 0.0) {
        throw SamplingError("Gamma draw produced an invalid value");
    }
    return result;
}

double stable_beta_ratio(double alpha_gamma, double beta_gamma)
{
    if (!std::isfinite(alpha_gamma) || !std::isfinite(beta_gamma)
        || alpha_gamma < 0.0 || beta_gamma < 0.0) {
        throw SamplingError("Beta Gamma inputs are outside their finite domain");
    }
    if (alpha_gamma == 0.0 && beta_gamma == 0.0) {
        throw SamplingError("Beta ratio is indeterminate after Gamma underflow");
    }
    if (alpha_gamma == 0.0) {return 0.0;}
    if (beta_gamma == 0.0) {return 1.0;}

    if (alpha_gamma > beta_gamma) {
        return 1.0 / (1.0 + beta_gamma / alpha_gamma);
    }
    const double ratio = alpha_gamma / beta_gamma;
    return ratio / (1.0 + ratio);
}

void validate_shape(double shape, const char *name)
{
    if (!std::isfinite(shape) || !(shape > 0.0)) {
        throw SamplingError(std::string(name) + " must be finite and positive");
    }
}

float sample_beta_at_entity(
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    std::uint64_t entity,
    double alpha,
    double beta,
    Options options)
{
    validate_shape(alpha, "alpha");
    validate_shape(beta, "beta");
    if (options.max_gamma_attempts == 0
        || options.max_gamma_attempts > default_max_gamma_attempts) {
        throw SamplingError("Gamma attempt cap must be in [1, 1024]");
    }
    if (std::fegetround() != FE_TONEAREST) {
        throw SamplingError("beta sampler requires round-to-nearest mode");
    }
    const std::uint64_t key = rng::derive_key(
        master_seed, rng::Stage::methylation_level, contig_index);

    const double alpha_gamma = gamma_draw(
        key,
        entity,
        alpha,
        Role::alpha_candidate,
        Role::alpha_boost,
        options.max_gamma_attempts);
    const double beta_gamma = gamma_draw(
        key,
        entity,
        beta,
        Role::beta_candidate,
        Role::beta_boost,
        options.max_gamma_attempts);
    const double result = stable_beta_ratio(alpha_gamma, beta_gamma);
    const float output = static_cast<float>(result);
    if (!std::isfinite(output) || output < 0.0F || output > 1.0F) {
        throw SamplingError("Beta draw is outside [0, 1]");
    }
    return output;
}

} // namespace

float sample_beta(
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    std::uint64_t reference_position,
    double alpha,
    double beta,
    Options options)
{
    return sample_beta_at_entity(
        master_seed,
        contig_index,
        reference_position,
        alpha,
        beta,
        options);
}

float sample_beta_for_site(
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    methdb::SiteEntity entity,
    double alpha,
    double beta,
    Options options)
{
    return sample_beta_at_entity(
        master_seed, contig_index, entity.value(), alpha, beta, options);
}

} // namespace htsim::beta_sampler

// ---- catalog --------------------------------------------------------

namespace htsim::methdb {
namespace {

void validate_shape(const ShapePair &shape, const char *name)
{
    if (!std::isfinite(shape.alpha) || !std::isfinite(shape.beta)
        || shape.alpha <= 0.0 || shape.beta <= 0.0) {
        throw CatalogError(std::string(name) + " Beta shapes must be finite and positive");
    }
}

} // namespace

void validate_context_shapes(const ContextShapes &shapes)
{
    validate_shape(shapes.cg, "CG");
    validate_shape(shapes.chg, "CHG");
    validate_shape(shapes.chh, "CHH");
}

const ShapePair &shape_for_context(
    model::MethylationContext context,
    const ContextShapes &shapes)
{
    switch (context) {
    case model::MethylationContext::cg_c:
    case model::MethylationContext::cg_g:
        return shapes.cg;
    case model::MethylationContext::chg_c:
    case model::MethylationContext::chg_g:
        return shapes.chg;
    case model::MethylationContext::chh_c:
    case model::MethylationContext::chh_g:
        return shapes.chh;
    }
    throw CatalogError("unknown methylation context");
}

MethylationCatalog::MethylationCatalog(
    std::uint32_t reference_length,
    std::vector<CatalogSite> sites)
    : sites_(std::move(sites))
{
    std::uint32_t previous = 0U;
    bool first = true;
    for (const CatalogSite &site : sites_) {
        if (site.reference_position >= reference_length
            || (!first && site.reference_position <= previous)
            || !std::isfinite(site.methylation_probability)
            || site.methylation_probability < 0.0F
            || site.methylation_probability > 1.0F) {
            throw CatalogError("loaded methylation profile row is invalid");
        }
        switch (site.context) {
        case model::MethylationContext::cg_c:
        case model::MethylationContext::cg_g:
        case model::MethylationContext::chg_c:
        case model::MethylationContext::chg_g:
        case model::MethylationContext::chh_c:
        case model::MethylationContext::chh_g:
            break;
        default:
            throw CatalogError("loaded methylation context is invalid");
        }
        switch (site.methylation_source) {
        case model::MethylationSource::beta:
        case model::MethylationSource::cgmap:
        case model::MethylationSource::asm_source:
        case model::MethylationSource::pooled_cgmap:
            break;
        default:
            throw CatalogError("loaded methylation source is invalid");
        }
        previous = site.reference_position;
        first = false;
    }
}

MethylationCatalog::MethylationCatalog(
    const model::Bases &contig_bases,
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    bool collect_non_cpg,
    const ContextShapes &shapes,
    const std::vector<CgmapRecord> *cgmap_records,
    bool pool_cgmap)
{
    if (contig_bases.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw CatalogError("contig length exceeds uint32");
    }
    validate_context_shapes(shapes);

    try {
        if (pool_cgmap && cgmap_records == nullptr) {
            throw CatalogError("CGmap pooling requires normalized CGmap records");
        }
        std::unique_ptr<CgmapPool> context_pool;
        if (cgmap_records != nullptr) {
            validate_cgmap_records(contig_bases, *cgmap_records);
            if (pool_cgmap) {
                context_pool = std::make_unique<CgmapPool>(*cgmap_records);
            }
        }
        for (std::size_t position = 0; position < contig_bases.size(); ++position) {
            const auto context = classify_context(
                contig_bases, static_cast<std::uint64_t>(position), collect_non_cpg);
            if (!context.has_value()) {continue;}
            const ShapePair &shape = shape_for_context(*context, shapes);
            const SiteEntity entity = reference_site_entity(
                static_cast<std::uint32_t>(position));
            const std::optional<float> pooled = context_pool
                ? context_pool->sample(
                      *context, master_seed, contig_index, entity)
                : std::nullopt;
            const float probability = pooled
                ? *pooled
                : beta_sampler::sample_beta_for_site(
                      master_seed,
                      contig_index,
                      entity,
                      shape.alpha,
                      shape.beta);
            sites_.push_back(CatalogSite{
                static_cast<std::uint32_t>(position),
                probability,
                *context,
                pooled
                    ? model::MethylationSource::pooled_cgmap
                    : model::MethylationSource::beta,
            });
        }
        if (cgmap_records != nullptr && !pool_cgmap) {
            auto site = sites_.begin();
            for (const CgmapRecord &record : *cgmap_records) {
                while (site != sites_.end()
                       && site->reference_position
                           < record.reference_position) {
                    ++site;
                }
                if (!record.has_probability || site == sites_.end()
                    || site->reference_position != record.reference_position) {
                    continue;
                }
                if (site->context != record.context) {
                    throw CatalogError(
                        "CGmap context disagrees with the prepared "
                        "methylation profile");
                }
                site->methylation_probability =
                    record.methylation_probability;
                site->methylation_source = model::MethylationSource::cgmap;
            }
        }
    } catch (const CatalogError &) {
        throw;
    } catch (const std::exception &error) {
        throw CatalogError(
            std::string("methylation profile preparation failed: ")
            + error.what());
    }
}

const std::vector<CatalogSite> &MethylationCatalog::sites() const noexcept
{
    return sites_;
}

std::pair<MethylationCatalog::const_iterator, MethylationCatalog::const_iterator>
MethylationCatalog::sites_in_range(
    std::uint32_t begin,
    std::uint32_t end) const
{
    if (begin > end) {throw CatalogError("methylation range is reversed");}
    const auto first = std::lower_bound(
        sites_.begin(), sites_.end(), begin,
        [](const CatalogSite &site, std::uint32_t position) {
            return site.reference_position < position;
        });
    const auto last = std::lower_bound(
        first, sites_.end(), end,
        [](const CatalogSite &site, std::uint32_t position) {
            return site.reference_position < position;
        });
    return {first, last};
}

} // namespace htsim::methdb

// ---- cgmap_pool --------------------------------------------------------

namespace htsim::methdb {
namespace {

std::size_t pool_index(model::MethylationContext context)
{
    switch (context) {
    case model::MethylationContext::cg_c:
    case model::MethylationContext::cg_g:
        return 0U;
    case model::MethylationContext::chg_c:
    case model::MethylationContext::chg_g:
        return 1U;
    case model::MethylationContext::chh_c:
    case model::MethylationContext::chh_g:
        return 2U;
    }
    throw CgmapPoolError("CGmap pool received an unknown methylation context");
}

void validate_record(const CgmapRecord &record)
{
    (void)pool_index(record.context);
    if (!std::isfinite(record.methylation_probability)
        || record.methylation_probability < 0.0F
        || record.methylation_probability > 1.0F
        || record.dinucleotide_second > 3U
        || (!record.has_probability
            && record.methylation_probability != 0.0F)) {
        throw CgmapPoolError("CGmap pool record is not normalized");
    }
}

} // namespace

CgmapPool::CgmapPool(const std::vector<CgmapRecord> &records)
{
    bool first = true;
    std::uint32_t previous_position = 0U;
    for (const CgmapRecord &record : records) {
        validate_record(record);
        if (!first && record.reference_position <= previous_position) {
            throw CgmapPoolError(
                "CGmap pool records must have unique increasing positions");
        }
        first = false;
        previous_position = record.reference_position;
        if (!record.has_probability) {continue;}
        auto &values = values_[pool_index(record.context)];
        if (values.size() >= std::numeric_limits<std::uint32_t>::max()) {
            throw CgmapPoolError("CGmap context pool exceeds uint32");
        }
        values.push_back(record.methylation_probability);
    }
}

std::uint32_t CgmapPool::size(model::MethylationContext context) const
{
    const std::size_t observed = values_[pool_index(context)].size();
    if (observed > std::numeric_limits<std::uint32_t>::max()) {
        throw CgmapPoolError("CGmap context pool exceeds uint32");
    }
    return static_cast<std::uint32_t>(observed);
}

std::optional<float> CgmapPool::sample(
    model::MethylationContext context,
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    SiteEntity entity) const
{
    const auto &values = values_[pool_index(context)];
    if (values.empty()) {return std::nullopt;}
    if (values.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw CgmapPoolError("CGmap context pool exceeds uint32");
    }
    try {
        const std::uint64_t key = rng::derive_key(
            master_seed, rng::Stage::methylation_level, contig_index);
        const std::uint64_t selected = rng::bounded_integer(
            key,
            entity.value(),
            cgmap_pool_local_index,
            static_cast<std::uint64_t>(values.size()));
        return values.at(static_cast<std::size_t>(selected));
    } catch (const CgmapPoolError &) {
        throw;
    } catch (const std::exception &error) {
        throw CgmapPoolError(
            std::string("CGmap pool sampling failed: ") + error.what());
    }
}

} // namespace htsim::methdb

// ---- diploid_catalog --------------------------------------------------------

namespace htsim::methdb {
namespace {

inline constexpr std::uint64_t insertion_origin_flag = UINT64_C(1) << 63U;

struct BaseOrigin {
    bool insertion = false;
    std::uint32_t reference_anchor = 0;
    std::uint32_t event_ordinal = 0;
    std::uint8_t insertion_offset = 0;

    std::uint64_t id() const
    {
        return insertion
            ? insertion_origin_id(event_ordinal, insertion_offset)
            : reference_origin_id(reference_anchor);
    }
};

bool origin_less(const BaseOrigin &left, const BaseOrigin &right) noexcept
{
    // Insertions are emitted immediately before reference[anchor].
    return std::make_tuple(
               left.reference_anchor,
               left.insertion ? 0U : 1U,
               left.event_ordinal,
               left.insertion_offset)
        < std::make_tuple(
               right.reference_anchor,
               right.insertion ? 0U : 1U,
               right.event_ordinal,
               right.insertion_offset);
}

struct BaseRecord {
    std::uint8_t base = 0;
    BaseOrigin origin;
};

class HaplotypeBaseStream {
public:
    HaplotypeBaseStream(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        std::uint8_t zero_based_haplotype)
        : contig_(contig),
          variants_(variants.variants()),
          zero_based_haplotype_(zero_based_haplotype)
    {
        if (zero_based_haplotype > 1U) {
            throw DiploidCatalogError("haplotype must be zero or one");
        }
        if (contig.length != contig.bases.size()
            || contig.length > std::numeric_limits<std::uint32_t>::max()
            || variants.contig_index() != contig.index
            || variants.reference_length() != contig.length) {
            throw DiploidCatalogError(
                "variant set does not match the materialized contig");
        }
        contig_length_ = static_cast<std::uint32_t>(contig.length);
    }

    std::optional<BaseRecord> next()
    {
        if (emitting_insertion_) {
            return emit_insertion_base();
        }
        while (true) {
            if (variant_index_ < variants_.size()) {
                const variant::Variant &event = variants_[variant_index_];
                if (event.reference_start < reference_cursor_) {
                    throw DiploidCatalogError(
                        "variant stream moved behind the reference cursor");
                }
                if (event.reference_start == reference_cursor_) {
                    if (!model::mask_contains(
                            event.alt_haplotypes, zero_based_haplotype_)) {
                        ++variant_index_;
                        continue;
                    }
                    if (event.kind == model::VariantKind::insertion) {
                        emitting_insertion_ = true;
                        insertion_offset_ = 0;
                        return emit_insertion_base();
                    }
                    if (event.kind == model::VariantKind::snv) {
                        const BaseRecord record{
                            event.alt_bases.front(),
                            {false, reference_cursor_, 0U, 0U},
                        };
                        reference_cursor_ = event.reference_end;
                        ++variant_index_;
                        return record;
                    }
                    if (event.kind == model::VariantKind::deletion) {
                        reference_cursor_ = event.reference_end;
                        ++variant_index_;
                        continue;
                    }
                    throw DiploidCatalogError("variant stream has an invalid event kind");
                }
            }
            if (reference_cursor_ < contig_length_) {
                const std::uint32_t position = reference_cursor_++;
                return BaseRecord{
                    contig_.bases[position],
                    {false, position, 0U, 0U},
                };
            }
            if (variant_index_ != variants_.size()) {
                throw DiploidCatalogError(
                    "variant stream has an event beyond the contig");
            }
            return std::nullopt;
        }
    }

private:
    std::optional<BaseRecord> emit_insertion_base()
    {
        const variant::Variant &event = variants_.at(variant_index_);
        if (variant_index_
                >= static_cast<std::size_t>(model::no_variant_index)
            || insertion_offset_ >= event.alt_bases.size()
            || insertion_offset_ >= static_cast<std::size_t>(
                model::maximum_insertion_bases)) {
            throw DiploidCatalogError("insertion stream state is invalid");
        }
        const auto offset = static_cast<std::uint8_t>(insertion_offset_);
        const BaseRecord record{
            event.alt_bases[insertion_offset_],
            {true,
             event.reference_start,
             static_cast<std::uint32_t>(variant_index_),
             offset},
        };
        ++insertion_offset_;
        if (insertion_offset_ == event.alt_bases.size()) {
            emitting_insertion_ = false;
            insertion_offset_ = 0;
            ++variant_index_;
        }
        return record;
    }

    const reference::Contig &contig_;
    const std::vector<variant::Variant> &variants_;
    std::uint8_t zero_based_haplotype_ = 0;
    std::uint32_t contig_length_ = 0;
    std::uint32_t reference_cursor_ = 0;
    std::size_t variant_index_ = 0;
    std::size_t insertion_offset_ = 0;
    bool emitting_insertion_ = false;
};

struct RawSite {
    BaseOrigin origin;
    model::MethylationContext context = model::MethylationContext::cg_c;
    bool equals_reference = false;
};

class ClassifiedSiteStream {
public:
    ClassifiedSiteStream(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        std::uint8_t zero_based_haplotype,
        bool collect_non_cpg)
        : contig_(contig),
          bases_(contig, variants, zero_based_haplotype),
          collect_non_cpg_(collect_non_cpg)
    {
        current_ = bases_.next();
        downstream_first_ = bases_.next();
        downstream_second_ = bases_.next();
    }

    std::optional<RawSite> next()
    {
        while (current_) {
            const auto context = classify_context(
                ContextNeighborhood{
                    base(upstream_second_),
                    base(upstream_first_),
                    current_->base,
                    base(downstream_first_),
                    base(downstream_second_),
                },
                collect_non_cpg_);
            std::optional<RawSite> result;
            if (context) {
                bool equals_reference = false;
                if (!current_->origin.insertion) {
                    const std::uint32_t position =
                        current_->origin.reference_anchor;
                    const auto reference_context = classify_context(
                        contig_.bases, position, collect_non_cpg_);
                    equals_reference =
                        current_->base == contig_.bases[position]
                        && reference_context == context;
                }
                result = RawSite{current_->origin, *context, equals_reference};
            }
            advance();
            if (result) {return result;}
        }
        return std::nullopt;
    }

private:
    static std::optional<std::uint8_t> base(
        const std::optional<BaseRecord> &record) noexcept
    {
        return record
            ? std::optional<std::uint8_t>(record->base)
            : std::nullopt;
    }

    void advance()
    {
        upstream_second_ = upstream_first_;
        upstream_first_ = current_;
        current_ = downstream_first_;
        downstream_first_ = downstream_second_;
        downstream_second_ = bases_.next();
    }

    const reference::Contig &contig_;
    HaplotypeBaseStream bases_;
    bool collect_non_cpg_ = false;
    std::optional<BaseRecord> upstream_second_;
    std::optional<BaseRecord> upstream_first_;
    std::optional<BaseRecord> current_;
    std::optional<BaseRecord> downstream_first_;
    std::optional<BaseRecord> downstream_second_;
};

model::HaplotypeMask output_mask(
    bool shared,
    std::uint8_t zero_based_haplotype)
{
    if (shared) {return model::HaplotypeMask::both;}
    return zero_based_haplotype == 0U
        ? model::HaplotypeMask::haplotype_1
        : model::HaplotypeMask::haplotype_2;
}

SiteEntity site_entity(
    const RawSite &site,
    bool shared,
    std::uint8_t zero_based_haplotype)
{
    if (site.equals_reference) {
        return reference_site_entity(site.origin.reference_anchor);
    }
    const auto mask = output_mask(shared, zero_based_haplotype);
    return site.origin.insertion
        ? insertion_site_entity(
              site.origin.event_ordinal,
              site.origin.insertion_offset,
              mask,
              zero_based_haplotype)
        : variant_reference_site_entity(
              site.origin.reference_anchor, mask, zero_based_haplotype);
}

DiploidSite make_site(
    const RawSite &raw,
    bool shared,
    std::uint8_t zero_based_haplotype,
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    const ContextShapes &shapes,
    const CgmapPool *context_pool)
{
    const SiteEntity entity = site_entity(raw, shared, zero_based_haplotype);
    const ShapePair &shape = shape_for_context(raw.context, shapes);
    const model::MethylationAllele allele = shared
        ? model::MethylationAllele::shared
        : raw.equals_reference
            ? model::MethylationAllele::reference_haplotype
            : model::MethylationAllele::alternate_haplotype;
    const std::optional<float> pooled = context_pool
        ? context_pool->sample(
              raw.context, master_seed, contig_index, entity)
        : std::nullopt;
    return {
        raw.origin.id(),
        raw.context,
        pooled
            ? model::MethylationSource::pooled_cgmap
            : model::MethylationSource::beta,
        allele,
        pooled
            ? *pooled
            : beta_sampler::sample_beta_for_site(
                  master_seed,
                  contig_index,
                  entity,
                  shape.alpha,
                  shape.beta),
    };
}

void sort_and_validate(std::vector<DiploidSite> &sites)
{
    std::sort(
        sites.begin(), sites.end(),
        [](const DiploidSite &left, const DiploidSite &right) {
            return left.origin_id < right.origin_id;
        });
    for (std::size_t index = 1; index < sites.size(); ++index) {
        if (sites[index - 1U].origin_id == sites[index].origin_id) {
            throw DiploidCatalogError("methylation profile has duplicate origins");
        }
    }
}

const DiploidSite *find_site(
    const std::vector<DiploidSite> &sites,
    std::uint64_t origin_id) noexcept
{
    const auto found = std::lower_bound(
        sites.begin(), sites.end(), origin_id,
        [](const DiploidSite &site, std::uint64_t value) {
            return site.origin_id < value;
        });
    return found != sites.end() && found->origin_id == origin_id
        ? &*found
        : nullptr;
}

void overlay_cgmap(
    std::vector<DiploidSite> &sites,
    const std::vector<CgmapRecord> &records)
{
    for (const CgmapRecord &record : records) {
        if (!record.has_probability) {continue;}
        const std::uint64_t origin =
            reference_origin_id(record.reference_position);
        const auto found = std::lower_bound(
            sites.begin(), sites.end(), origin,
            [](const DiploidSite &site, std::uint64_t value) {
                return site.origin_id < value;
            });
        if (found == sites.end() || found->origin_id != origin
            || found->context != record.context
            || found->allele
                == model::MethylationAllele::alternate_haplotype) {
            continue;
        }
        found->methylation_probability = record.methylation_probability;
        found->methylation_source = model::MethylationSource::cgmap;
    }
}

const variant::Variant &resolve_asm_snv(
    const std::vector<variant::Variant> &variants,
    const AsmRecord &record)
{
    const auto found = std::lower_bound(
        variants.begin(), variants.end(), record.linked_variant_position,
        [](const variant::Variant &event, std::uint32_t position) {
            return event.reference_start < position;
        });
    if (found == variants.end()
        || found->reference_start != record.linked_variant_position
        || found->reference_end != record.linked_variant_position + 1U
        || found->kind != model::VariantKind::snv
        || found->ref_bases.size() != 1U
        || found->alt_bases.size() != 1U
        || found->ref_bases.front() != record.linked_reference_base
        || found->alt_bases.front() != record.linked_alternate_base) {
        throw DiploidCatalogError(
            "ASM row does not resolve to its exact typed VCF SNV");
    }
    if (found->alt_haplotypes == model::HaplotypeMask::both) {
        throw DiploidCatalogError(
            "ASM linked VCF SNV must retain one reference haplotype");
    }
    if (found->alt_haplotypes != model::HaplotypeMask::haplotype_1
        && found->alt_haplotypes
            != model::HaplotypeMask::haplotype_2) {
        throw DiploidCatalogError("ASM linked VCF SNV has an invalid mask");
    }
    return *found;
}

std::uint8_t alternate_haplotype(
    model::HaplotypeMask mask)
{
    if (mask == model::HaplotypeMask::haplotype_1) {return 0U;}
    if (mask == model::HaplotypeMask::haplotype_2) {return 1U;}
    throw DiploidCatalogError(
        "ASM linked VCF SNV must be heterozygous");
}

void overlay_asm(
    std::vector<DiploidSite> &shared_sites,
    std::array<std::vector<DiploidSite>, 2> &haplotype_sites,
    const std::vector<variant::Variant> &variants,
    const std::vector<AsmRecord> &records)
{
    std::vector<DiploidSite> retained_shared;
    retained_shared.reserve(shared_sites.size());
    std::array<std::vector<DiploidSite>, 2> additions;
    additions[0].reserve(records.size());
    additions[1].reserve(records.size());

    std::size_t shared_index = 0U;
    for (const AsmRecord &record : records) {
        const std::uint64_t target =
            reference_origin_id(record.target_reference_position);
        while (shared_index < shared_sites.size()
               && shared_sites[shared_index].origin_id < target) {
            retained_shared.push_back(shared_sites[shared_index++]);
        }
        if (shared_index == shared_sites.size()
            || shared_sites[shared_index].origin_id != target
            || shared_sites[shared_index].context != record.context
            || shared_sites[shared_index].allele
                != model::MethylationAllele::shared
            || find_site(haplotype_sites[0], target) != nullptr
            || find_site(haplotype_sites[1], target) != nullptr) {
            throw DiploidCatalogError(
                "ASM target is not one shared reference-equivalent diploid site");
        }

        const variant::Variant &event = resolve_asm_snv(variants, record);
        const std::uint8_t alt_haplotype =
            alternate_haplotype(event.alt_haplotypes);
        const std::uint8_t ref_haplotype =
            static_cast<std::uint8_t>(1U - alt_haplotype);

        DiploidSite reference_site = shared_sites[shared_index];
        reference_site.methylation_source = model::MethylationSource::asm_source;
        reference_site.allele =
            model::MethylationAllele::reference_haplotype;
        reference_site.methylation_probability =
            record.reference_methylation_probability;
        additions[ref_haplotype].push_back(reference_site);

        DiploidSite alternate_site = shared_sites[shared_index];
        alternate_site.methylation_source = model::MethylationSource::asm_source;
        alternate_site.allele =
            model::MethylationAllele::alternate_haplotype;
        alternate_site.methylation_probability =
            record.alternate_methylation_probability;
        additions[alt_haplotype].push_back(alternate_site);
        ++shared_index;
    }
    retained_shared.insert(
        retained_shared.end(),
        shared_sites.begin() + static_cast<std::ptrdiff_t>(shared_index),
        shared_sites.end());
    shared_sites.swap(retained_shared);
    for (std::size_t haplotype = 0U; haplotype < additions.size(); ++haplotype) {
        haplotype_sites[haplotype].insert(
            haplotype_sites[haplotype].end(),
            additions[haplotype].begin(),
            additions[haplotype].end());
        sort_and_validate(haplotype_sites[haplotype]);
    }
}

bool context_matches_base(
    model::MethylationContext context,
    std::uint8_t base) noexcept
{
    const auto value = static_cast<std::uint8_t>(context);
    return (value < 8U && base == 1U) || (value >= 8U && base == 2U);
}

} // namespace

std::uint64_t reference_origin_id(std::uint32_t reference_position) noexcept
{
    return reference_position;
}

std::uint64_t insertion_origin_id(
    std::uint32_t event_ordinal,
    std::uint8_t insertion_offset)
{
    if (event_ordinal == model::no_variant_index) {
        throw DiploidCatalogError("insertion origin uses the no-event sentinel");
    }
    if (insertion_offset >= model::maximum_insertion_bases) {
        throw DiploidCatalogError("insertion origin offset must be in [0, 3]");
    }
    return insertion_origin_flag
        | (static_cast<std::uint64_t>(event_ordinal) << 2U)
        | insertion_offset;
}

DiploidMethylationCatalog::DiploidMethylationCatalog(
    std::uint32_t contig_index,
    std::uint32_t reference_length,
    std::vector<DiploidSite> shared_sites,
    std::array<std::vector<DiploidSite>, 2> haplotype_sites)
    : contig_index_(contig_index),
      reference_length_(reference_length),
      shared_sites_(std::move(shared_sites)),
      haplotype_sites_(std::move(haplotype_sites))
{
    const auto validate = [&](const std::vector<DiploidSite> &sites) {
        std::uint64_t previous = 0U;
        bool first = true;
        for (const DiploidSite &site : sites) {
            if ((!first && site.origin_id <= previous)
                || !std::isfinite(site.methylation_probability)
                || site.methylation_probability < 0.0F
                || site.methylation_probability > 1.0F
                || static_cast<std::uint8_t>(site.allele) > 2U) {
                throw DiploidCatalogError(
                    "loaded diploid methylation profile row is invalid");
            }
            if ((site.origin_id >> 63U) == 0U
                && site.origin_id >= reference_length_) {
                throw DiploidCatalogError(
                    "loaded diploid reference origin is outside the contig");
            }
            previous = site.origin_id;
            first = false;
        }
    };
    validate(shared_sites_);
    validate(haplotype_sites_[0]);
    validate(haplotype_sites_[1]);
}

DiploidMethylationCatalog::DiploidMethylationCatalog(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint64_t master_seed,
    bool collect_non_cpg,
    const ContextShapes &shapes,
    const std::vector<CgmapRecord> *cgmap_records,
    const std::vector<AsmRecord> *asm_records,
    bool pool_cgmap)
{
    try {
        if (contig.length != contig.bases.size()
            || contig.length > std::numeric_limits<std::uint32_t>::max()
            || variants.contig_index() != contig.index
            || variants.reference_length() != contig.length) {
            throw DiploidCatalogError(
                "variant set does not match the materialized contig");
        }
        contig_index_ = contig.index;
        reference_length_ = static_cast<std::uint32_t>(contig.length);
        validate_context_shapes(shapes);
        if (pool_cgmap && cgmap_records == nullptr) {
            throw DiploidCatalogError(
                "CGmap pooling requires normalized CGmap records");
        }
        std::unique_ptr<CgmapPool> context_pool;
        if (cgmap_records != nullptr) {
            validate_cgmap_records(contig.bases, *cgmap_records);
            if (pool_cgmap) {
                context_pool = std::make_unique<CgmapPool>(*cgmap_records);
            }
        }

        ClassifiedSiteStream haplotype_0(
            contig, variants, 0U, collect_non_cpg);
        ClassifiedSiteStream haplotype_1(
            contig, variants, 1U, collect_non_cpg);
        auto left = haplotype_0.next();
        auto right = haplotype_1.next();
        while (left || right) {
            if (left && right && left->origin.id() == right->origin.id()) {
                if (left->context == right->context) {
                    if (left->equals_reference != right->equals_reference) {
                        throw DiploidCatalogError(
                            "identical diploid sites disagree with the reference");
                    }
                    shared_sites_.push_back(make_site(
                        *left,
                        true,
                        0U,
                        master_seed,
                        contig.index,
                        shapes,
                        context_pool.get()));
                } else {
                    haplotype_sites_[0].push_back(make_site(
                        *left,
                        false,
                        0U,
                        master_seed,
                        contig.index,
                        shapes,
                        context_pool.get()));
                    haplotype_sites_[1].push_back(make_site(
                        *right,
                        false,
                        1U,
                        master_seed,
                        contig.index,
                        shapes,
                        context_pool.get()));
                }
                left = haplotype_0.next();
                right = haplotype_1.next();
            } else if (left
                       && (!right || origin_less(left->origin, right->origin))) {
                haplotype_sites_[0].push_back(make_site(
                    *left,
                    false,
                    0U,
                    master_seed,
                    contig.index,
                    shapes,
                    context_pool.get()));
                left = haplotype_0.next();
            } else {
                haplotype_sites_[1].push_back(make_site(
                    *right,
                    false,
                    1U,
                    master_seed,
                    contig.index,
                    shapes,
                    context_pool.get()));
                right = haplotype_1.next();
            }
        }
        sort_and_validate(shared_sites_);
        sort_and_validate(haplotype_sites_[0]);
        sort_and_validate(haplotype_sites_[1]);
        if (cgmap_records != nullptr && !pool_cgmap) {
            overlay_cgmap(shared_sites_, *cgmap_records);
            overlay_cgmap(haplotype_sites_[0], *cgmap_records);
            overlay_cgmap(haplotype_sites_[1], *cgmap_records);
        }
        if (asm_records != nullptr) {
            validate_asm_records(contig.bases, *asm_records);
            overlay_asm(
                shared_sites_, haplotype_sites_, variants.variants(), *asm_records);
        }
    } catch (const DiploidCatalogError &) {
        throw;
    } catch (const std::exception &error) {
        throw DiploidCatalogError(error.what());
    }
}

const std::vector<DiploidSite> &
DiploidMethylationCatalog::shared_sites() const noexcept
{
    return shared_sites_;
}

const std::vector<DiploidSite> &DiploidMethylationCatalog::haplotype_sites(
    std::uint8_t zero_based_haplotype) const
{
    if (zero_based_haplotype > 1U) {
        throw DiploidCatalogError("haplotype must be zero or one");
    }
    return haplotype_sites_[zero_based_haplotype];
}

std::vector<model::MethylationSite>
DiploidMethylationCatalog::sites_for_projection(
    const haplotype::ProjectedInterval &projection) const
{
    const std::uint8_t zero_based_haplotype = projection.haplotype;
    if (zero_based_haplotype > 1U) {
        throw DiploidCatalogError("haplotype must be zero or one");
    }
    if (projection.contig_index != contig_index_
        || projection.reference_start > projection.reference_end
        || projection.reference_end > reference_length_
        || projection.template_bases.size()
            != projection.reference_positions.size()
        || projection.template_bases.size()
            != projection.base_variant_indices.size()
        || projection.template_bases.size()
            > std::numeric_limits<std::uint32_t>::max()) {
        throw DiploidCatalogError(
            "haplotype projection identity or shape is invalid");
    }

    std::unordered_map<std::uint32_t, const model::Variant *> variants;
    for (const model::Variant &event : projection.variants) {
        if (!variants.emplace(event.index, &event).second) {
            throw DiploidCatalogError("haplotype projection has duplicate variants");
        }
        if (event.phased_haplotype != 255U
            && event.phased_haplotype != zero_based_haplotype) {
            throw DiploidCatalogError(
                "haplotype projection contains an event from another haplotype");
        }
    }
    std::unordered_map<std::uint32_t, std::uint8_t> insertion_offsets;
    std::vector<model::MethylationSite> result;
    for (std::size_t offset = 0; offset < projection.template_bases.size(); ++offset) {
        const std::int64_t reference_position =
            projection.reference_positions[offset];
        std::uint64_t origin_id = 0;
        if (reference_position >= 0) {
            if (static_cast<std::uint64_t>(reference_position)
                    > std::numeric_limits<std::uint32_t>::max()
                || static_cast<std::uint64_t>(reference_position)
                    < projection.reference_start
                || static_cast<std::uint64_t>(reference_position)
                    >= projection.reference_end) {
                throw DiploidCatalogError(
                    "mapped projection base is outside its reference interval");
            }
            origin_id = reference_origin_id(
                static_cast<std::uint32_t>(reference_position));
        } else if (reference_position == -1) {
            const std::uint32_t variant_index = projection.base_variant_indices[offset];
            const auto found = variants.find(variant_index);
            if (found == variants.end()
                || found->second->kind != model::VariantKind::insertion) {
                throw DiploidCatalogError(
                    "inserted projection base has no insertion event");
            }
            std::uint8_t &insertion_offset = insertion_offsets[variant_index];
            if (insertion_offset >= found->second->alt_bases.size()
                || projection.template_bases[offset]
                    != found->second->alt_bases[insertion_offset]) {
                throw DiploidCatalogError(
                    "inserted projection bases disagree with their event");
            }
            origin_id = insertion_origin_id(variant_index, insertion_offset);
            ++insertion_offset;
        } else {
            throw DiploidCatalogError(
                "projection reference position must be non-negative or -1");
        }

        const DiploidSite *shared = find_site(shared_sites_, origin_id);
        const DiploidSite *specific = find_site(
            haplotype_sites_[zero_based_haplotype], origin_id);
        if (shared && specific) {
            throw DiploidCatalogError(
                "shared and haplotype profile sites overlap at one origin");
        }
        const DiploidSite *site = shared ? shared : specific;
        if (!site) {continue;}
        if (!context_matches_base(site->context, projection.template_bases[offset])) {
            throw DiploidCatalogError(
                "methylation context disagrees with projected center base");
        }
        result.push_back(model::MethylationSite{
            static_cast<std::uint32_t>(result.size()),
            static_cast<std::uint32_t>(offset),
            reference_position,
            site->context,
            site->methylation_source,
            site->allele,
            site->methylation_probability,
        });
    }
    for (const auto &entry : variants) {
        const model::Variant &event = *entry.second;
        if (event.kind != model::VariantKind::insertion) {continue;}
        const auto count = insertion_offsets.find(entry.first);
        if (count == insertion_offsets.end()
            || count->second != event.alt_bases.size()) {
            throw DiploidCatalogError(
                "projection did not contain every inserted event base");
        }
    }
    return result;
}

} // namespace htsim::methdb

// ---- fixed snapshot --------------------------------------------------------

namespace htsim::methdb {
class SnapshotWriterImpl {
public:
    explicit SnapshotWriterImpl(std::ostream &output) : output_(output)
    {
        if (deflateInit(&stream_, Z_DEFAULT_COMPRESSION) != Z_OK) {
            throw SnapshotError("cannot initialize MethDB compression");
        }
        active_ = true;
    }

    ~SnapshotWriterImpl()
    {
        if (active_) {(void)deflateEnd(&stream_);}
    }

    void write(const void *data, std::size_t size)
    {
        const auto *cursor = static_cast<const std::uint8_t *>(data);
        while (size != 0U) {
            const std::size_t chunk = std::min<std::size_t>(
                size, std::numeric_limits<uInt>::max());
            stream_.next_in = const_cast<Bytef *>(cursor);
            stream_.avail_in = static_cast<uInt>(chunk);
            while (stream_.avail_in != 0U) {pump(Z_NO_FLUSH);}
            cursor += chunk;
            size -= chunk;
        }
    }

    void finish()
    {
        if (!active_) {throw SnapshotError("MethDB compression was finished twice");}
        int status = Z_OK;
        while (status != Z_STREAM_END) {status = pump(Z_FINISH);}
        if (deflateEnd(&stream_) != Z_OK) {
            active_ = false;
            throw SnapshotError("cannot finalize MethDB compression");
        }
        active_ = false;
        if (!output_) {throw SnapshotError("failed while writing MethDB snapshot");}
    }

private:
    int pump(int flush)
    {
        stream_.next_out = output_buffer_.data();
        stream_.avail_out = static_cast<uInt>(output_buffer_.size());
        const int status = deflate(&stream_, flush);
        if (status != Z_OK && status != Z_STREAM_END) {
            throw SnapshotError("MethDB compression failed");
        }
        const std::size_t produced = output_buffer_.size() - stream_.avail_out;
        output_.write(
            reinterpret_cast<const char *>(output_buffer_.data()),
            static_cast<std::streamsize>(produced));
        if (!output_) {throw SnapshotError("failed while writing MethDB snapshot");}
        return status;
    }

    std::ostream &output_;
    z_stream stream_ = {};
    std::array<std::uint8_t, 65536> output_buffer_ = {};
    bool active_ = false;
};

namespace {

void snapshot_write_raw(
    std::ostream &output, const void *data, std::size_t size)
{
    output.write(static_cast<const char *>(data), static_cast<std::streamsize>(size));
    if (!output) {throw SnapshotError("failed while writing MethDB snapshot");}
}

void snapshot_write_bytes(
    SnapshotWriterImpl &output, const void *data, std::size_t size)
{
    output.write(data, size);
}

void snapshot_write_u8(SnapshotWriterImpl &output, std::uint8_t value)
{
    snapshot_write_bytes(output, &value, 1U);
}

void snapshot_write_u32(SnapshotWriterImpl &output, std::uint32_t value)
{
    std::uint8_t bytes[4];
    for (unsigned index = 0U; index < 4U; ++index) {
        bytes[index] = static_cast<std::uint8_t>(value >> (index * 8U));
    }
    snapshot_write_bytes(output, bytes, sizeof(bytes));
}

void snapshot_write_u16(SnapshotWriterImpl &output, std::uint16_t value)
{
    const std::uint8_t bytes[2] = {
        static_cast<std::uint8_t>(value),
        static_cast<std::uint8_t>(value >> 8U),
    };
    snapshot_write_bytes(output, bytes, sizeof(bytes));
}

std::uint16_t probability_to_unorm16(float value)
{
    if (!std::isfinite(value) || value < 0.0F || value > 1.0F) {
        throw SnapshotError("MethDB probability is outside [0, 1]");
    }
    constexpr double scale = static_cast<double>(UINT16_MAX);
    const auto encoded = static_cast<std::uint32_t>(
        static_cast<double>(value) * scale + 0.5);
    return static_cast<std::uint16_t>(encoded);
}

float unorm16_to_probability(std::uint16_t value)
{
    constexpr float scale = static_cast<float>(UINT16_MAX);
    return static_cast<float>(value) / scale;
}

void snapshot_write_varuint(SnapshotWriterImpl &output, std::uint64_t value)
{
    do {
        std::uint8_t byte = static_cast<std::uint8_t>(value & UINT64_C(0x7f));
        value >>= 7U;
        if (value != 0U) {byte |= UINT8_C(0x80);}
        snapshot_write_u8(output, byte);
    } while (value != 0U);
}

std::uint8_t snapshot_context_code(model::MethylationContext context)
{
    switch (context) {
    case model::MethylationContext::cg_c: return 0U;
    case model::MethylationContext::chg_c: return 1U;
    case model::MethylationContext::chh_c: return 2U;
    case model::MethylationContext::cg_g: return 3U;
    case model::MethylationContext::chg_g: return 4U;
    case model::MethylationContext::chh_g: return 5U;
    }
    throw SnapshotError("MethDB row has an invalid methylation context");
}

std::uint8_t snapshot_pack_site(
    model::MethylationContext context,
    model::MethylationSource source,
    model::MethylationAllele allele)
{
    const std::uint8_t source_value = static_cast<std::uint8_t>(source);
    const std::uint8_t allele_value = static_cast<std::uint8_t>(allele);
    if (source_value < 1U || source_value > 4U || allele_value > 2U) {
        throw SnapshotError("MethDB row has invalid source or allele metadata");
    }
    return static_cast<std::uint8_t>(
        snapshot_context_code(context)
        | static_cast<std::uint8_t>((source_value - 1U) << 3U)
        | static_cast<std::uint8_t>(allele_value << 5U));
}

void snapshot_write_metadata(
    SnapshotWriterImpl &output,
    const reference::ContigMetadata &metadata,
    std::uint8_t mode,
    std::uint32_t first_count,
    std::uint32_t second_count,
    std::uint32_t third_count)
{
    if (metadata.name.empty()
        || metadata.name.size() > std::numeric_limits<std::uint32_t>::max()
        || metadata.length > std::numeric_limits<std::uint32_t>::max()) {
        throw SnapshotError("MethDB contig metadata is outside the format boundary");
    }
    snapshot_write_u32(output, static_cast<std::uint32_t>(metadata.name.size()));
    snapshot_write_bytes(output, metadata.name.data(), metadata.name.size());
    snapshot_write_u32(output, static_cast<std::uint32_t>(metadata.length));
    snapshot_write_bytes(
        output, metadata.reference_sha256.data(), metadata.reference_sha256.size());
    snapshot_write_u8(output, mode);
    snapshot_write_u8(output, 0U);
    snapshot_write_u8(output, 0U);
    snapshot_write_u8(output, 0U);
    snapshot_write_u32(output, first_count);
    snapshot_write_u32(output, second_count);
    snapshot_write_u32(output, third_count);
}

void snapshot_write_catalog_sites(
    SnapshotWriterImpl &output,
    const std::vector<CatalogSite> &sites)
{
    std::uint64_t previous = 0U;
    bool first = true;
    for (const CatalogSite &site : sites) {
        if (!first && site.reference_position <= previous) {
            throw SnapshotError("MethDB reference sites are not strictly ordered");
        }
        const std::uint64_t delta = first
            ? site.reference_position
            : static_cast<std::uint64_t>(site.reference_position) - previous;
        snapshot_write_varuint(output, delta);
        snapshot_write_u8(output, snapshot_pack_site(
            site.context,
            site.methylation_source,
            model::MethylationAllele::shared));
        snapshot_write_u16(
            output, probability_to_unorm16(site.methylation_probability));
        previous = site.reference_position;
        first = false;
    }
}

void snapshot_write_diploid_sites(
    SnapshotWriterImpl &output,
    const std::vector<DiploidSite> &sites)
{
    std::uint64_t previous = 0U;
    bool first = true;
    for (const DiploidSite &site : sites) {
        if (!first && site.origin_id <= previous) {
            throw SnapshotError("MethDB diploid sites are not strictly ordered");
        }
        snapshot_write_varuint(
            output, first ? site.origin_id : site.origin_id - previous);
        snapshot_write_u8(output, snapshot_pack_site(
            site.context, site.methylation_source, site.allele));
        snapshot_write_u16(
            output, probability_to_unorm16(site.methylation_probability));
        previous = site.origin_id;
        first = false;
    }
}

class SnapshotReader {
public:
    explicit SnapshotReader(const std::string &path)
        : input_(path, std::ios::binary)
    {
        if (!input_) {throw SnapshotError("cannot open MethDB snapshot: " + path);}
        char magic[sizeof(methdb_magic) - 1U];
        read_raw(magic, sizeof(magic));
        if (std::memcmp(magic, methdb_magic, sizeof(magic)) != 0) {
            throw SnapshotError("MethDB snapshot magic is invalid");
        }
        std::uint8_t version = 0U;
        read_raw(&version, sizeof(version));
        if (version != methdb_version) {
            throw SnapshotError("MethDB snapshot version is unsupported");
        }
        if (inflateInit(&stream_) != Z_OK) {
            throw SnapshotError("cannot initialize MethDB decompression");
        }
        inflate_active_ = true;
    }

    ~SnapshotReader()
    {
        if (inflate_active_) {(void)inflateEnd(&stream_);}
    }

    void read(void *destination, std::size_t size)
    {
        auto *cursor = static_cast<std::uint8_t *>(destination);
        while (size != 0U) {
            if (inflate_finished_) {
                throw SnapshotError("MethDB snapshot is truncated");
            }
            const std::size_t chunk = std::min<std::size_t>(
                size, std::numeric_limits<uInt>::max());
            stream_.next_out = cursor;
            stream_.avail_out = static_cast<uInt>(chunk);
            while (stream_.avail_out != 0U && !inflate_finished_) {
                inflate_once();
            }
            const std::size_t produced = chunk - stream_.avail_out;
            if (produced == 0U) {
                throw SnapshotError("MethDB snapshot is truncated");
            }
            cursor += produced;
            size -= produced;
        }
    }

    std::uint8_t u8()
    {
        std::uint8_t value = 0U;
        read(&value, 1U);
        return value;
    }

    std::uint32_t u32()
    {
        std::uint8_t bytes[4];
        read(bytes, sizeof(bytes));
        std::uint32_t value = 0U;
        for (unsigned index = 0U; index < 4U; ++index) {
            value |= static_cast<std::uint32_t>(bytes[index]) << (index * 8U);
        }
        return value;
    }

    std::uint16_t u16()
    {
        std::uint8_t bytes[2];
        read(bytes, sizeof(bytes));
        return static_cast<std::uint16_t>(bytes[0])
            | static_cast<std::uint16_t>(
                static_cast<std::uint16_t>(bytes[1]) << 8U);
    }

    std::uint64_t varuint()
    {
        std::uint64_t value = 0U;
        for (unsigned index = 0U; index < 10U; ++index) {
            const std::uint8_t byte = u8();
            if (index == 9U && (byte & UINT8_C(0xfe)) != 0U) {
                throw SnapshotError("MethDB variable integer exceeds uint64");
            }
            value |= static_cast<std::uint64_t>(byte & UINT8_C(0x7f))
                << (index * 7U);
            if ((byte & UINT8_C(0x80)) == 0U) {
                if (index != 0U && byte == 0U) {
                    throw SnapshotError("MethDB variable integer is not canonical");
                }
                return value;
            }
        }
        throw SnapshotError("MethDB variable integer is unterminated");
    }

    crypto::Sha256Digest finish()
    {
        std::uint8_t extra = 0U;
        while (!inflate_finished_) {
            stream_.next_out = &extra;
            stream_.avail_out = 1U;
            inflate_once();
            if (stream_.avail_out == 0U) {
                throw SnapshotError("MethDB snapshot has trailing data");
            }
        }
        if (stream_.avail_in != 0U) {
            throw SnapshotError("MethDB snapshot has trailing bytes");
        }
        if (!input_exhausted_) {
            char raw_extra = 0;
            if (input_.get(raw_extra)) {
                throw SnapshotError("MethDB snapshot has trailing bytes");
            }
            if (!input_.eof()) {
                throw SnapshotError("MethDB snapshot read failed");
            }
        }
        if (inflateEnd(&stream_) != Z_OK) {
            inflate_active_ = false;
            throw SnapshotError("cannot finalize MethDB decompression");
        }
        inflate_active_ = false;
        return hash_.digest();
    }

private:
    void read_raw(void *destination, std::size_t size)
    {
        if (size == 0U) {return;}
        input_.read(
            static_cast<char *>(destination), static_cast<std::streamsize>(size));
        if (input_.gcount() != static_cast<std::streamsize>(size)) {
            throw SnapshotError("MethDB snapshot is truncated");
        }
        hash_.update(static_cast<const std::uint8_t *>(destination), size);
    }

    void fill_input()
    {
        if (input_exhausted_) {
            throw SnapshotError("MethDB compressed payload is truncated");
        }
        input_.read(
            reinterpret_cast<char *>(input_buffer_.data()),
            static_cast<std::streamsize>(input_buffer_.size()));
        const std::streamsize observed = input_.gcount();
        if (observed <= 0) {
            if (!input_.eof()) {
                throw SnapshotError("MethDB snapshot read failed");
            }
            input_exhausted_ = true;
            throw SnapshotError("MethDB compressed payload is truncated");
        }
        const std::size_t size = static_cast<std::size_t>(observed);
        hash_.update(input_buffer_.data(), size);
        stream_.next_in = input_buffer_.data();
        stream_.avail_in = static_cast<uInt>(size);
        if (size < input_buffer_.size()) {input_exhausted_ = true;}
    }

    void inflate_once()
    {
        if (stream_.avail_in == 0U) {fill_input();}
        const uInt before_in = stream_.avail_in;
        const uInt before_out = stream_.avail_out;
        const int status = inflate(&stream_, Z_NO_FLUSH);
        if (status == Z_STREAM_END) {
            inflate_finished_ = true;
        } else if (status != Z_OK) {
            throw SnapshotError("MethDB compressed payload is invalid");
        }
        if (before_in == stream_.avail_in && before_out == stream_.avail_out
            && !inflate_finished_) {
            throw SnapshotError("MethDB decompression made no progress");
        }
    }

    std::ifstream input_;
    crypto::Sha256 hash_;
    z_stream stream_ = {};
    std::array<std::uint8_t, 65536> input_buffer_ = {};
    bool inflate_active_ = false;
    bool inflate_finished_ = false;
    bool input_exhausted_ = false;
};

std::uint32_t checked_snapshot_count(std::size_t size)
{
    if (size > std::numeric_limits<std::uint32_t>::max()) {
        throw SnapshotError("MethDB row count exceeds uint32");
    }
    return static_cast<std::uint32_t>(size);
}

void snapshot_unpack_site(
    std::uint8_t packed,
    model::MethylationContext &context,
    model::MethylationSource &source,
    model::MethylationAllele &allele)
{
    if ((packed & UINT8_C(0x80)) != 0U) {
        throw SnapshotError("MethDB packed row reserved bit is nonzero");
    }
    switch (packed & UINT8_C(0x07)) {
    case 0U: context = model::MethylationContext::cg_c; break;
    case 1U: context = model::MethylationContext::chg_c; break;
    case 2U: context = model::MethylationContext::chh_c; break;
    case 3U: context = model::MethylationContext::cg_g; break;
    case 4U: context = model::MethylationContext::chg_g; break;
    case 5U: context = model::MethylationContext::chh_g; break;
    default: throw SnapshotError("MethDB packed context is invalid");
    }
    source = static_cast<model::MethylationSource>(
        ((packed >> 3U) & UINT8_C(0x03)) + 1U);
    const std::uint8_t allele_value =
        (packed >> 5U) & UINT8_C(0x03);
    if (allele_value > 2U) {
        throw SnapshotError("MethDB packed allele is invalid");
    }
    allele = static_cast<model::MethylationAllele>(allele_value);
}

std::uint64_t snapshot_read_delta(
    SnapshotReader &reader, std::uint64_t previous, bool first)
{
    const std::uint64_t delta = reader.varuint();
    if (!first && delta == 0U) {
        throw SnapshotError("MethDB delta is not strictly positive");
    }
    if (delta > std::numeric_limits<std::uint64_t>::max() - previous) {
        throw SnapshotError("MethDB delta overflows uint64");
    }
    return first ? delta : previous + delta;
}

CatalogSite snapshot_read_catalog_site(
    SnapshotReader &reader, std::uint64_t previous, bool first)
{
    const std::uint64_t position = snapshot_read_delta(reader, previous, first);
    if (position > std::numeric_limits<std::uint32_t>::max()) {
        throw SnapshotError("MethDB reference position exceeds uint32");
    }
    CatalogSite site;
    site.reference_position = static_cast<std::uint32_t>(position);
    model::MethylationAllele allele = model::MethylationAllele::shared;
    snapshot_unpack_site(
        reader.u8(), site.context, site.methylation_source, allele);
    if (allele != model::MethylationAllele::shared) {
        throw SnapshotError("MethDB reference row has a diploid allele");
    }
    site.methylation_probability = unorm16_to_probability(reader.u16());
    return site;
}

DiploidSite snapshot_read_diploid_site(
    SnapshotReader &reader, std::uint64_t previous, bool first)
{
    DiploidSite site;
    site.origin_id = snapshot_read_delta(reader, previous, first);
    snapshot_unpack_site(
        reader.u8(), site.context, site.methylation_source, site.allele);
    site.methylation_probability = unorm16_to_probability(reader.u16());
    return site;
}

std::string_view snapshot_context_name(model::MethylationContext context)
{
    switch (context) {
    case model::MethylationContext::cg_c: return "CG-C";
    case model::MethylationContext::chg_c: return "CHG-C";
    case model::MethylationContext::chh_c: return "CHH-C";
    case model::MethylationContext::cg_g: return "CG-G";
    case model::MethylationContext::chg_g: return "CHG-G";
    case model::MethylationContext::chh_g: return "CHH-G";
    }
    throw SnapshotError("MethDB BED row has an invalid context");
}

std::string_view snapshot_source_name(model::MethylationSource source)
{
    switch (source) {
    case model::MethylationSource::cgmap: return "cgmap";
    case model::MethylationSource::asm_source: return "asm";
    case model::MethylationSource::beta: return "beta";
    case model::MethylationSource::pooled_cgmap: return "pooled-cgmap";
    }
    throw SnapshotError("MethDB BED row has an invalid source");
}

std::string_view snapshot_allele_name(model::MethylationAllele allele)
{
    switch (allele) {
    case model::MethylationAllele::shared: return "shared";
    case model::MethylationAllele::reference_haplotype: return "reference";
    case model::MethylationAllele::alternate_haplotype: return "alternate";
    }
    throw SnapshotError("MethDB BED row has an invalid allele");
}

void snapshot_write_hex(
    std::ostream &sink, const crypto::Sha256Digest &digest)
{
    static constexpr char digits[] = "0123456789abcdef";
    for (const std::uint8_t byte : digest) {
        sink.put(digits[byte >> 4U]);
        sink.put(digits[byte & UINT8_C(0x0f)]);
    }
}

void snapshot_write_bed_record(
    std::ostream &sink,
    std::string_view contig_name,
    std::string_view set_name,
    std::uint64_t origin_id,
    model::MethylationContext context,
    model::MethylationSource source,
    model::MethylationAllele allele,
    float probability)
{
    constexpr std::uint64_t insertion_flag = UINT64_C(1) << 63U;
    const bool insertion = (origin_id & insertion_flag) != 0U;
    const std::uint16_t encoded_probability = probability_to_unorm16(probability);
    if (insertion) {
        const std::uint64_t payload = origin_id & ~insertion_flag;
        sink << "#insertion\t" << contig_name << '\t' << set_name << '\t'
             << origin_id << '\t' << (payload >> 2U) << '\t'
             << (payload & 3U);
    } else {
        const std::uint32_t score = (
            static_cast<std::uint32_t>(encoded_probability) * 1000U
            + 32767U) / 65535U;
        const char strand = static_cast<std::uint8_t>(context) < 8U
            ? '+'
            : '-';
        sink << contig_name << '\t' << origin_id << '\t' << (origin_id + 1U)
             << "\tmethdb:" << set_name << ':' << origin_id << '\t'
             << score << '\t' << strand << '\t' << set_name << '\t'
             << origin_id << "\treference\t.\t.";
    }
    sink << '\t' << snapshot_context_name(context)
         << '\t' << snapshot_source_name(source)
         << '\t' << snapshot_allele_name(allele)
         << '\t' << encoded_probability
         << '\t' << std::setprecision(std::numeric_limits<float>::max_digits10)
         << probability << '\n';
}

} // namespace

SnapshotWriter::SnapshotWriter(
    std::ostream &output,
    const crypto::Sha256Digest &binding,
    std::uint32_t contig_count)
    : output_(output),
      impl_(std::make_unique<SnapshotWriterImpl>(output)),
      contig_count_(contig_count)
{
    snapshot_write_raw(
        output_, methdb_magic, sizeof(methdb_magic) - 1U);
    snapshot_write_raw(output_, &methdb_version, sizeof(methdb_version));
    snapshot_write_bytes(*impl_, binding.data(), binding.size());
    snapshot_write_u32(*impl_, contig_count_);
}

SnapshotWriter::~SnapshotWriter() = default;

void SnapshotWriter::write_reference(
    const reference::ContigMetadata &metadata,
    const MethylationCatalog &catalog)
{
    if (finished_ || written_ >= contig_count_) {
        throw SnapshotError("MethDB snapshot received too many contigs");
    }
    const auto &sites = catalog.sites();
    snapshot_write_metadata(
        *impl_, metadata, 0U, checked_snapshot_count(sites.size()), 0U, 0U);
    snapshot_write_catalog_sites(*impl_, sites);
    ++written_;
}

void SnapshotWriter::write_diploid(
    const reference::ContigMetadata &metadata,
    const DiploidMethylationCatalog &catalog)
{
    if (finished_ || written_ >= contig_count_) {
        throw SnapshotError("MethDB snapshot received too many contigs");
    }
    const auto &shared = catalog.shared_sites();
    const auto &haplotype_0 = catalog.haplotype_sites(0U);
    const auto &haplotype_1 = catalog.haplotype_sites(1U);
    snapshot_write_metadata(
        *impl_,
        metadata,
        1U,
        checked_snapshot_count(shared.size()),
        checked_snapshot_count(haplotype_0.size()),
        checked_snapshot_count(haplotype_1.size()));
    snapshot_write_diploid_sites(*impl_, shared);
    snapshot_write_diploid_sites(*impl_, haplotype_0);
    snapshot_write_diploid_sites(*impl_, haplotype_1);
    ++written_;
}

void SnapshotWriter::finish()
{
    if (finished_) {throw SnapshotError("MethDB snapshot was finished twice");}
    if (written_ != contig_count_) {
        throw SnapshotError("MethDB snapshot contig count is incomplete");
    }
    impl_->finish();
    if (!output_) {throw SnapshotError("MethDB snapshot output failed");}
    finished_ = true;
}

void export_snapshot_bed(const std::string &path, std::ostream &sink)
{
    SnapshotReader reader(path);
    crypto::Sha256Digest binding = {};
    reader.read(binding.data(), binding.size());
    const std::uint32_t contig_count = reader.u32();

    sink << "#format\t" << methdb_bed_format << '\n'
         << "#binding_sha256\t";
    snapshot_write_hex(sink, binding);
    sink << '\n'
         << "#columns\tchrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
            "set\torigin_id\torigin_kind\tvariant_event\tinsertion_offset\t"
            "context\tsource\tallele\tprobability_u16\tprobability\n"
         << "#insertion_columns\tcontig\tset\torigin_id\tvariant_event\t"
            "insertion_offset\tcontext\tsource\tallele\tprobability_u16\t"
            "probability\n";

    for (std::uint32_t index = 0U; index < contig_count; ++index) {
        const std::uint32_t name_size = reader.u32();
        if (name_size == 0U || name_size > UINT32_C(1048576)) {
            throw SnapshotError("MethDB contig name length is invalid");
        }
        std::string name(name_size, '\0');
        reader.read(name.data(), name.size());
        if (name.find_first_of("\t\r\n") != std::string::npos) {
            throw SnapshotError("MethDB contig name cannot be represented in BED");
        }
        const std::uint32_t reference_length = reader.u32();
        crypto::Sha256Digest reference_sha256 = {};
        reader.read(reference_sha256.data(), reference_sha256.size());
        const std::uint8_t mode = reader.u8();
        if (reader.u8() != 0U || reader.u8() != 0U || reader.u8() != 0U) {
            throw SnapshotError("MethDB contig reserved bytes are nonzero");
        }
        const std::uint32_t first_count = reader.u32();
        const std::uint32_t second_count = reader.u32();
        const std::uint32_t third_count = reader.u32();
        const std::uint64_t generous_limit =
            static_cast<std::uint64_t>(reference_length) * 5U + 4U;
        if (first_count > generous_limit
            || second_count > generous_limit
            || third_count > generous_limit) {
            throw SnapshotError("MethDB row count is incompatible with contig length");
        }
        if (mode > 1U) {
            throw SnapshotError("MethDB contig mode is invalid");
        }
        if (mode == 0U && (second_count != 0U || third_count != 0U)) {
            throw SnapshotError("reference MethDB has diploid row counts");
        }

        sink << "#contig\t" << index << '\t' << name << '\t'
             << reference_length << '\t';
        snapshot_write_hex(sink, reference_sha256);
        sink << '\t' << (mode == 0U ? "reference" : "diploid") << '\n';

        if (mode == 0U) {
            std::uint64_t previous = 0U;
            for (std::uint32_t row = 0U; row < first_count; ++row) {
                const CatalogSite site = snapshot_read_catalog_site(
                    reader, previous, row == 0U);
                previous = site.reference_position;
                if (site.reference_position >= reference_length) {
                    throw SnapshotError(
                        "MethDB reference origin is outside the contig");
                }
                snapshot_write_bed_record(
                    sink,
                    name,
                    "reference",
                    site.reference_position,
                    site.context,
                    site.methylation_source,
                    model::MethylationAllele::shared,
                    site.methylation_probability);
            }
            continue;
        }

        const auto write_diploid_set = [&](
            std::uint32_t count, std::string_view set_name) {
            std::uint64_t previous = 0U;
            for (std::uint32_t row = 0U; row < count; ++row) {
                const DiploidSite site = snapshot_read_diploid_site(
                    reader, previous, row == 0U);
                previous = site.origin_id;
                if ((site.origin_id >> 63U) == 0U
                    && site.origin_id >= reference_length) {
                    throw SnapshotError(
                        "MethDB diploid reference origin is outside the contig");
                }
                snapshot_write_bed_record(
                    sink,
                    name,
                    set_name,
                    site.origin_id,
                    site.context,
                    site.methylation_source,
                    site.allele,
                    site.methylation_probability);
            }
        };
        write_diploid_set(first_count, "shared");
        write_diploid_set(second_count, "haplotype-1");
        write_diploid_set(third_count, "haplotype-2");
    }

    const crypto::Sha256Digest file_sha256 = reader.finish();
    sink << "#file_sha256\t";
    snapshot_write_hex(sink, file_sha256);
    sink << '\n';
    if (!sink) {throw SnapshotError("failed while writing MethDB BED");}
}

Snapshot::Snapshot(
    const std::string &path,
    const crypto::Sha256Digest &expected_binding,
    const std::vector<reference::ContigMetadata> &reference_catalog)
{
    SnapshotReader reader(path);
    crypto::Sha256Digest binding = {};
    reader.read(binding.data(), binding.size());
    if (binding != expected_binding) {
        throw SnapshotError("MethDB profile binding is incompatible");
    }
    const std::uint32_t contig_count = reader.u32();
    if (contig_count != reference_catalog.size()) {
        throw SnapshotError("MethDB snapshot contig count disagrees with reference");
    }
    contigs_.reserve(contig_count);
    for (std::uint32_t index = 0U; index < contig_count; ++index) {
        const std::uint32_t name_size = reader.u32();
        if (name_size == 0U || name_size > UINT32_C(1048576)) {
            throw SnapshotError("MethDB contig name length is invalid");
        }
        SnapshotContig entry;
        entry.name.resize(name_size);
        reader.read(entry.name.data(), entry.name.size());
        entry.reference_length = reader.u32();
        reader.read(entry.reference_sha256.data(), entry.reference_sha256.size());
        const auto &expected = reference_catalog[index];
        if (entry.name != expected.name
            || entry.reference_length != expected.length
            || entry.reference_sha256 != expected.reference_sha256) {
            throw SnapshotError("MethDB contig identity disagrees with reference");
        }
        const std::uint8_t mode = reader.u8();
        if (reader.u8() != 0U || reader.u8() != 0U || reader.u8() != 0U) {
            throw SnapshotError("MethDB contig reserved bytes are nonzero");
        }
        const std::uint32_t first_count = reader.u32();
        const std::uint32_t second_count = reader.u32();
        const std::uint32_t third_count = reader.u32();
        const std::uint64_t generous_limit =
            static_cast<std::uint64_t>(entry.reference_length) * 5U + 4U;
        if (first_count > generous_limit
            || second_count > generous_limit
            || third_count > generous_limit) {
            throw SnapshotError("MethDB row count is incompatible with contig length");
        }
        if (mode == 0U) {
            if (second_count != 0U || third_count != 0U) {
                throw SnapshotError("reference MethDB has diploid row counts");
            }
            entry.reference_sites.reserve(first_count);
            std::uint64_t previous = 0U;
            for (std::uint32_t row = 0U; row < first_count; ++row) {
                CatalogSite site = snapshot_read_catalog_site(
                    reader, previous, row == 0U);
                previous = site.reference_position;
                entry.reference_sites.push_back(std::move(site));
            }
            (void)MethylationCatalog(
                entry.reference_length, entry.reference_sites);
        } else if (mode == 1U) {
            entry.diploid = true;
            entry.shared_sites.reserve(first_count);
            entry.haplotype_sites[0].reserve(second_count);
            entry.haplotype_sites[1].reserve(third_count);
            std::uint64_t previous = 0U;
            for (std::uint32_t row = 0U; row < first_count; ++row) {
                DiploidSite site = snapshot_read_diploid_site(
                    reader, previous, row == 0U);
                previous = site.origin_id;
                entry.shared_sites.push_back(site);
            }
            previous = 0U;
            for (std::uint32_t row = 0U; row < second_count; ++row) {
                DiploidSite site = snapshot_read_diploid_site(
                    reader, previous, row == 0U);
                previous = site.origin_id;
                entry.haplotype_sites[0].push_back(site);
            }
            previous = 0U;
            for (std::uint32_t row = 0U; row < third_count; ++row) {
                DiploidSite site = snapshot_read_diploid_site(
                    reader, previous, row == 0U);
                previous = site.origin_id;
                entry.haplotype_sites[1].push_back(site);
            }
            (void)DiploidMethylationCatalog(
                index,
                entry.reference_length,
                entry.shared_sites,
                entry.haplotype_sites);
        } else {
            throw SnapshotError("MethDB contig mode is invalid");
        }
        contigs_.push_back(std::move(entry));
    }
    file_sha256_ = reader.finish();
}

const SnapshotContig &Snapshot::contig(std::uint32_t contig_index) const
{
    if (contig_index >= contigs_.size()) {
        throw SnapshotError("MethDB contig index is out of range");
    }
    return contigs_[contig_index];
}

} // namespace htsim::methdb
