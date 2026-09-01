#include "protocol.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <limits>
#include <unordered_set>
#include <utility>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>
#include <condition_variable>
#include <deque>
#include <future>
#include <iterator>
#include <memory>
#include <mutex>
#include <queue>
#include <string_view>
#include <thread>

#include "types.h"

// ---- wire --------------------------------------------------------

namespace htsim::protocol {
namespace {

constexpr std::array<std::uint8_t, 8> magic = {{
    'B', 'S', 'R', 'S', 'T', 'R', 'M', 0,
}};

constexpr std::array<std::uint32_t, 256> make_crc32c_table() noexcept
{
    std::array<std::uint32_t, 256> table = {};
    for (std::size_t value = 0; value < table.size(); ++value) {
        std::uint32_t crc = static_cast<std::uint32_t>(value);
        for (unsigned int bit = 0; bit < 8; ++bit) {
            crc = (crc >> 1)
                ^ ((crc & 1U) != 0 ? UINT32_C(0x82f63b78) : 0U);
        }
        table[value] = crc;
    }
    return table;
}

constexpr auto crc32c_table = make_crc32c_table();

bool valid_utf8(std::string_view text) noexcept
{
    std::size_t index = 0;
    while (index < text.size()) {
        const auto lead = static_cast<unsigned char>(text[index]);
        std::size_t length = 0;
        std::uint32_t codepoint = 0;
        if (lead <= 0x7fU) {
            length = 1;
            codepoint = lead;
        } else if (lead >= 0xc2U && lead <= 0xdfU) {
            length = 2;
            codepoint = lead & 0x1fU;
        } else if (lead >= 0xe0U && lead <= 0xefU) {
            length = 3;
            codepoint = lead & 0x0fU;
        } else if (lead >= 0xf0U && lead <= 0xf4U) {
            length = 4;
            codepoint = lead & 0x07U;
        } else {
            return false;
        }
        if (length > text.size() - index) {return false;}
        for (std::size_t offset = 1; offset < length; ++offset) {
            const auto continuation =
                static_cast<unsigned char>(text[index + offset]);
            if ((continuation & 0xc0U) != 0x80U) {return false;}
            codepoint = (codepoint << 6U) | (continuation & 0x3fU);
        }
        if ((length == 3 && codepoint < 0x800U)
            || (length == 4 && codepoint < 0x10000U)
            || (codepoint >= 0xd800U && codepoint <= 0xdfffU)
            || codepoint > 0x10ffffU) {
            return false;
        }
        index += length;
    }
    return true;
}

void validate_string(
    std::string_view name,
    std::string_view value,
    bool nonempty = false)
{
    if (value.size() > maximum_string_bytes) {
        throw ProtocolError(std::string(name) + " exceeds the 1 MiB string limit");
    }
    if (nonempty && value.empty()) {
        throw ProtocolError(std::string(name) + " must not be empty");
    }
    if (value.find('\0') != std::string_view::npos) {
        throw ProtocolError(std::string(name) + " must not contain NUL");
    }
    if (!valid_utf8(value)) {
        throw ProtocolError(std::string(name) + " must be valid UTF-8");
    }
}

bool canonical_uuid(std::string_view value) noexcept
{
    if (value.size() != 36) {return false;}
    for (std::size_t index = 0; index < value.size(); ++index) {
        if (index == 8 || index == 13 || index == 18 || index == 23) {
            if (value[index] != '-') {return false;}
            continue;
        }
        const char character = value[index];
        if (!((character >= '0' && character <= '9')
              || (character >= 'a' && character <= 'f'))) {
            return false;
        }
    }
    return true;
}

bool ascii_digit(char character) noexcept
{
    return character >= '0' && character <= '9';
}

bool identifier_character(char character) noexcept
{
    return ascii_digit(character)
        || (character >= 'A' && character <= 'Z')
        || (character >= 'a' && character <= 'z') || character == '-';
}

bool numeric_identifier(std::string_view value) noexcept
{
    return !value.empty() && (value.size() == 1 || value.front() != '0')
        && std::all_of(value.begin(), value.end(), ascii_digit);
}

bool identifier_list(
    std::string_view value,
    bool reject_numeric_leading_zero) noexcept
{
    if (value.empty()) {return false;}
    std::size_t start = 0;
    while (start <= value.size()) {
        const std::size_t end = value.find('.', start);
        const std::size_t stop =
            end == std::string_view::npos ? value.size() : end;
        const std::string_view item = value.substr(start, stop - start);
        if (item.empty()
            || !std::all_of(item.begin(), item.end(), identifier_character)) {
            return false;
        }
        const bool numeric =
            std::all_of(item.begin(), item.end(), ascii_digit);
        if (reject_numeric_leading_zero && numeric && item.size() > 1
            && item.front() == '0') {
            return false;
        }
        if (end == std::string_view::npos) {break;}
        start = end + 1;
    }
    return true;
}

bool semantic_version(std::string_view value) noexcept
{
    const std::size_t plus = value.find('+');
    if (plus != std::string_view::npos) {
        if (value.find('+', plus + 1) != std::string_view::npos
            || !identifier_list(value.substr(plus + 1), false)) {
            return false;
        }
        value = value.substr(0, plus);
    }
    const std::size_t hyphen = value.find('-');
    if (hyphen != std::string_view::npos) {
        if (!identifier_list(value.substr(hyphen + 1), true)) {return false;}
        value = value.substr(0, hyphen);
    }
    const std::size_t first = value.find('.');
    if (first == std::string_view::npos) {return false;}
    const std::size_t second = value.find('.', first + 1);
    if (second == std::string_view::npos
        || value.find('.', second + 1) != std::string_view::npos) {
        return false;
    }
    return numeric_identifier(value.substr(0, first))
        && numeric_identifier(value.substr(first + 1, second - first - 1))
        && numeric_identifier(value.substr(second + 1));
}

template <typename Value>
std::uint32_t checked_size(
    std::string_view name,
    const std::vector<Value> &values)
{
    if (values.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw ProtocolError(std::string(name) + " count exceeds u32");
    }
    return static_cast<std::uint32_t>(values.size());
}

void require_size(
    std::string_view name,
    std::size_t observed,
    std::size_t expected)
{
    if (observed != expected) {
        throw ProtocolError(std::string(name) + " length is inconsistent");
    }
}

void validate_prefix(
    std::string_view name,
    const std::vector<std::uint32_t> &values,
    std::uint32_t rows,
    std::uint32_t total)
{
    require_size(name, values.size(), static_cast<std::size_t>(rows) + 1U);
    if (values.front() != 0U || values.back() != total) {
        throw ProtocolError(
            std::string(name) + " must begin at zero and end at its flat count");
    }
    if (!std::is_sorted(values.begin(), values.end())) {
        throw ProtocolError(std::string(name) + " must be monotone");
    }
}

bool valid_capture_strand(std::uint8_t value) noexcept
{
    return value <= static_cast<std::uint8_t>(CaptureStrand::reverse);
}

bool valid_context(std::uint8_t value) noexcept
{
    switch (static_cast<MethylationContext>(value)) {
    case MethylationContext::cg_c:
    case MethylationContext::chg_c:
    case MethylationContext::chh_c:
    case MethylationContext::cg_g:
    case MethylationContext::chg_g:
    case MethylationContext::chh_g:
        return true;
    }
    return false;
}

bool cytosine_context(std::uint8_t value) noexcept
{
    return value == static_cast<std::uint8_t>(MethylationContext::cg_c)
        || value == static_cast<std::uint8_t>(MethylationContext::chg_c)
        || value == static_cast<std::uint8_t>(MethylationContext::chh_c);
}

bool guanine_context(std::uint8_t value) noexcept
{
    return value == static_cast<std::uint8_t>(MethylationContext::cg_g)
        || value == static_cast<std::uint8_t>(MethylationContext::chg_g)
        || value == static_cast<std::uint8_t>(MethylationContext::chh_g);
}

bool valid_source(std::uint8_t value) noexcept
{
    return value >= static_cast<std::uint8_t>(MethylationSource::cgmap)
        && value <= static_cast<std::uint8_t>(MethylationSource::pooled_cgmap);
}

bool valid_allele(std::uint8_t value) noexcept
{
    return value <=
        static_cast<std::uint8_t>(MethylationAllele::alternate_haplotype);
}

bool valid_variant_kind(std::uint8_t value) noexcept
{
    return value >= static_cast<std::uint8_t>(VariantKind::snv)
        && value <= static_cast<std::uint8_t>(VariantKind::deletion);
}

void validate_bases(
    std::string_view name,
    const Bases &values,
    std::uint8_t maximum)
{
    if (!std::all_of(values.begin(), values.end(), [maximum](std::uint8_t value) {
            return value <= maximum;
        })) {
        throw ProtocolError(std::string(name) + " contains an invalid base code");
    }
}

void validate_header(const Header &header)
{
    validate_string("header.run_id", header.run_id, true);
    if (!canonical_uuid(header.run_id)) {
        throw ProtocolError("header.run_id must be canonical lowercase UUID text");
    }
    validate_string("header.core_version", header.core_version, true);
    if (!semantic_version(header.core_version)) {
        throw ProtocolError("header.core_version must be a semantic version");
    }
    if (header.rng_contract != protocol::rng_contract) {
        throw ProtocolError("unsupported RNG contract");
    }
    switch (header.technology) {
    case Technology::wgbs:
    case Technology::rrbs:
    case Technology::tbs:
    case Technology::wgs:
    case Technology::wes:
    case Technology::ts:
        break;
    default:
        throw ProtocolError("header technology is invalid");
    }
    if (header.mates_per_fragment != 1U
        && header.mates_per_fragment != 2U) {
        throw ProtocolError("header mates_per_fragment must be one or two");
    }
    if (header.base_encoding != BaseEncoding::acgtn_u8) {
        throw ProtocolError("unsupported base encoding");
    }
    if (header.ambiguity_policy != AmbiguityPolicy::preserve_n) {
        throw ProtocolError("unsupported ambiguity policy");
    }
    if (header.read_length_r1 == 0U
        || (header.mates_per_fragment == 1U && header.read_length_r2 != 0U)
        || (header.mates_per_fragment == 2U && header.read_length_r2 == 0U)) {
        throw ProtocolError("header read lengths disagree with SE/PE mode");
    }
    checked_size("header.contigs", header.contigs);
    if (header.contigs.empty()) {
        throw ProtocolError("header.contigs must not be empty");
    }
    std::unordered_set<std::string> names;
    for (const Contig &contig : header.contigs) {
        validate_string("contig.name", contig.name, true);
        if (contig.length == 0U) {
            throw ProtocolError("contig length must be positive");
        }
        if (!names.insert(contig.name).second) {
            throw ProtocolError("contig names must be unique");
        }
    }
}

class Encoder {
public:
    explicit Encoder(std::vector<std::uint8_t> &bytes) noexcept : bytes_(bytes) {}

    void u8(std::uint8_t value) {bytes_.push_back(value);}

    void u16(std::uint16_t value)
    {
        for (unsigned int byte = 0; byte < 2; ++byte) {
            bytes_.push_back(static_cast<std::uint8_t>(value >> (8U * byte)));
        }
    }

    void u32(std::uint32_t value)
    {
        for (unsigned int byte = 0; byte < 4; ++byte) {
            bytes_.push_back(static_cast<std::uint8_t>(value >> (8U * byte)));
        }
    }

    void u64(std::uint64_t value)
    {
        for (unsigned int byte = 0; byte < 8; ++byte) {
            bytes_.push_back(static_cast<std::uint8_t>(value >> (8U * byte)));
        }
    }

    void f32(float value)
    {
        static_assert(sizeof(float) == sizeof(std::uint32_t),
                      "protocol requires IEEE-754 binary32");
        static_assert(std::numeric_limits<float>::is_iec559,
                      "protocol requires IEEE-754 binary32");
        std::uint32_t bits = 0;
        std::memcpy(&bits, &value, sizeof(bits));
        u32(bits);
    }

    void raw(const std::uint8_t *data, std::size_t size)
    {
        if (size != 0U) {bytes_.insert(bytes_.end(), data, data + size);}
    }

    template <std::size_t Size>
    void raw(const std::array<std::uint8_t, Size> &values)
    {
        raw(values.data(), values.size());
    }

    void raw(const std::vector<std::uint8_t> &values)
    {
        raw(values.data(), values.size());
    }

    void string(std::string_view value)
    {
        u32(static_cast<std::uint32_t>(value.size()));
        raw(reinterpret_cast<const std::uint8_t *>(value.data()), value.size());
    }

    void u32s(const std::vector<std::uint32_t> &values)
    {
        for (std::uint32_t value : values) {u32(value);}
    }

    void u64s(const std::vector<std::uint64_t> &values)
    {
        for (std::uint64_t value : values) {u64(value);}
    }

    void f32s(const std::vector<float> &values)
    {
        for (float value : values) {f32(value);}
    }

    void align4()
    {
        while (bytes_.size() % 4U != 0U) {bytes_.push_back(0U);}
    }

private:
    std::vector<std::uint8_t> &bytes_;
};

std::vector<std::uint8_t> encode_header_payload(const Header &header)
{
    validate_header(header);
    std::vector<std::uint8_t> payload;
    Encoder encoder(payload);
    encoder.string(header.run_id);
    encoder.string(header.core_version);
    encoder.string(header.rng_contract);
    encoder.u64(header.master_seed);
    encoder.raw(header.normalized_config_sha256);
    encoder.u8(static_cast<std::uint8_t>(header.technology));
    encoder.u8(static_cast<std::uint8_t>(header.has_details));
    encoder.u8(header.mates_per_fragment);
    encoder.u8(static_cast<std::uint8_t>(header.base_encoding));
    encoder.u8(static_cast<std::uint8_t>(header.ambiguity_policy));
    encoder.u8(0U);
    encoder.u8(0U);
    encoder.u8(0U);
    encoder.u32(header.read_length_r1);
    encoder.u32(header.read_length_r2);
    encoder.u32(static_cast<std::uint32_t>(header.contigs.size()));
    for (const Contig &contig : header.contigs) {
        encoder.string(contig.name);
        encoder.u32(contig.length);
        encoder.raw(contig.reference_sha256);
    }
    encoder.align4();
    return payload;
}

void validate_annotations(const FragmentBatch &, const Header &, const FragmentDetails &);

void validate_batch(const FragmentBatch &batch, const Header &header)
{
    validate_header(header);
    const std::uint32_t fragments =
        checked_size("batch.contig_indices", batch.contig_indices);
    const std::uint32_t bases =
        checked_size("batch.template_bases", batch.template_bases);
    const std::uint32_t mates =
        checked_size("batch.mate_indices", batch.mate_indices);
    const std::uint32_t sites =
        checked_size("batch.site_template_offsets", batch.site_template_offsets);
    if (fragments == 0U) {
        throw ProtocolError("fragment batch must not be empty");
    }
    if (static_cast<std::uint64_t>(batch.first_fragment_ordinal) + fragments
        > std::numeric_limits<std::uint32_t>::max()) {
        throw ProtocolError("batch fragment ordinal range exceeds u32");
    }
    require_size("batch.reference_starts", batch.reference_starts.size(), fragments);
    require_size("batch.reference_ends", batch.reference_ends.size(), fragments);
    require_size("batch.haplotypes", batch.haplotypes.size(), fragments);
    require_size("batch.capture_strands", batch.capture_strands.size(), fragments);
    require_size("batch.mate_template_starts", batch.mate_template_starts.size(), mates);
    require_size("batch.mate_template_ends", batch.mate_template_ends.size(), mates);
    require_size("batch.mate_reverse_complements",
                 batch.mate_reverse_complements.size(), mates);
    require_size("batch.site_probabilities", batch.site_probabilities.size(), sites);
    require_size("batch.site_contexts", batch.site_contexts.size(), sites);
    require_size("batch.methylation_sources", batch.methylation_sources.size(), sites);
    require_size("batch.site_alleles", batch.site_alleles.size(), sites);
    validate_prefix("batch.template_offsets", batch.template_offsets, fragments, bases);
    validate_prefix("batch.mate_offsets", batch.mate_offsets, fragments, mates);
    validate_prefix("batch.site_offsets", batch.site_offsets, fragments, sites);
    validate_bases("batch.template_bases", batch.template_bases, 4U);

    for (std::uint32_t row = 0; row < fragments; ++row) {
        const std::uint32_t contig_index = batch.contig_indices[row];
        if (contig_index >= header.contigs.size()) {
            throw ProtocolError("fragment contig index is outside the header");
        }
        if (batch.reference_starts[row] >= batch.reference_ends[row]
            || batch.reference_ends[row] > header.contigs[contig_index].length) {
            throw ProtocolError("fragment reference envelope is invalid");
        }
        if (batch.haplotypes[row] > 1U) {
            throw ProtocolError("fragment haplotype must be zero or one");
        }
        if (!valid_capture_strand(batch.capture_strands[row])) {
            throw ProtocolError("fragment capture strand is invalid");
        }
        const std::uint32_t template_length =
            batch.template_offsets[row + 1U] - batch.template_offsets[row];
        if (template_length == 0U) {
            throw ProtocolError("fragment template must not be empty");
        }
        const std::uint32_t mate_begin = batch.mate_offsets[row];
        const std::uint32_t mate_end = batch.mate_offsets[row + 1U];
        if (mate_end - mate_begin != header.mates_per_fragment) {
            throw ProtocolError("fragment mate count disagrees with header");
        }
        for (std::uint32_t local = 0; local < header.mates_per_fragment; ++local) {
            const std::uint32_t index = mate_begin + local;
            const std::uint32_t begin = batch.mate_template_starts[index];
            const std::uint32_t end = batch.mate_template_ends[index];
            const std::uint32_t expected =
                local == 0U ? header.read_length_r1 : header.read_length_r2;
            if (batch.mate_indices[index] != local
                || batch.mate_reverse_complements[index] > 1U
                || begin >= end || end > template_length
                || end - begin != expected) {
                throw ProtocolError("mate row disagrees with header or template");
            }
        }
        const std::uint32_t site_begin = batch.site_offsets[row];
        const std::uint32_t site_end = batch.site_offsets[row + 1U];
        std::optional<std::uint32_t> previous_site;
        for (std::uint32_t index = site_begin; index < site_end; ++index) {
            const std::uint32_t offset = batch.site_template_offsets[index];
            if (offset >= template_length
                || (previous_site && offset <= *previous_site)) {
                throw ProtocolError("site offsets must be in range and increasing");
            }
            previous_site = offset;
            const float probability = batch.site_probabilities[index];
            if (!std::isfinite(probability) || probability < 0.0F
                || probability > 1.0F) {
                throw ProtocolError("site probability must be finite and in [0,1]");
            }
            if (!valid_context(batch.site_contexts[index])
                || !valid_source(batch.methylation_sources[index])
                || !valid_allele(batch.site_alleles[index])) {
                throw ProtocolError("site enum value is invalid");
            }
            const std::uint8_t base = batch.template_bases[
                batch.template_offsets[row] + offset];
            if ((cytosine_context(batch.site_contexts[index]) && base != 1U)
                || (guanine_context(batch.site_contexts[index]) && base != 2U)) {
                throw ProtocolError("site context is incompatible with template base");
            }
        }
    }
    if (header.has_details) {
        if (!batch.details) {
            throw ProtocolError("header requires Full-Details batch columns");
        }
        validate_annotations(batch, header, *batch.details);
    } else if (batch.details) {
        throw ProtocolError("header forbids Full-Details batch columns");
    }
}

void validate_annotations(
    const FragmentBatch &batch,
    const Header &,
    const FragmentDetails &details)
{
    const std::uint32_t fragments = batch.fragment_count();
    const std::uint32_t sites = batch.methylation_site_count();
    const std::uint32_t projections = checked_size(
        "details.projection_template_starts",
        details.projection_template_starts);
    const std::uint32_t variants =
        checked_size("details.variant_indices", details.variant_indices);
    const std::uint32_t variant_id_bytes =
        checked_size("details.variant_ids", details.variant_ids);
    const std::uint32_t ref_bases =
        checked_size("details.variant_ref_bases", details.variant_ref_bases);
    const std::uint32_t alt_bases =
        checked_size("details.variant_alt_bases", details.variant_alt_bases);
    const std::uint32_t original_ns = checked_size(
        "details.original_n_template_offsets",
        details.original_n_template_offsets);

    validate_prefix(
        "details.projection_offsets", details.projection_offsets, fragments, projections);
    validate_prefix("details.variant_offsets", details.variant_offsets, fragments, variants);
    validate_prefix(
        "details.original_n_offsets", details.original_n_offsets, fragments, original_ns);
    validate_prefix(
        "details.variant_id_offsets",
        details.variant_id_offsets,
        variants,
        variant_id_bytes);
    validate_prefix(
        "details.variant_ref_offsets", details.variant_ref_offsets, variants, ref_bases);
    validate_prefix(
        "details.variant_alt_offsets", details.variant_alt_offsets, variants, alt_bases);
    require_size("details.projection_template_ends",
                 details.projection_template_ends.size(), projections);
    require_size("details.projection_reference_starts",
                 details.projection_reference_starts.size(), projections);
    require_size(
        "details.variant_reference_starts", details.variant_reference_starts.size(), variants);
    require_size(
        "details.variant_reference_ends", details.variant_reference_ends.size(), variants);
    require_size(
        "details.variant_template_starts", details.variant_template_starts.size(), variants);
    require_size(
        "details.variant_template_ends", details.variant_template_ends.size(), variants);
    require_size("details.site_reference_positions",
                 details.site_reference_positions.size(), sites);
    require_size("details.variant_kinds", details.variant_kinds.size(), variants);
    require_size("details.variant_sources", details.variant_sources.size(), variants);
    require_size("details.variant_phased_haplotypes",
                 details.variant_phased_haplotypes.size(), variants);
    validate_bases("details.variant_ref_bases", details.variant_ref_bases, 3U);
    validate_bases("details.variant_alt_bases", details.variant_alt_bases, 3U);

    for (std::uint32_t row = 0; row < fragments; ++row) {
        const std::uint32_t flat_template_begin = batch.template_offsets[row];
        const std::uint32_t template_length =
            batch.template_offsets[row + 1U] - flat_template_begin;
        const std::uint32_t reference_begin = batch.reference_starts[row];
        const std::uint32_t reference_end = batch.reference_ends[row];
        const std::uint8_t haplotype = batch.haplotypes[row];
        std::vector<std::uint8_t> projection_cover(template_length, 0U);
        std::vector<std::uint8_t> insertion_cover(template_length, 0U);
        std::vector<std::uint8_t> event_cover(template_length, 0U);
        std::vector<std::uint32_t> mapped_positions(
            template_length, no_reference_position);

        std::optional<std::uint32_t> previous_template_end;
        std::optional<std::uint32_t> previous_reference_begin;
        std::optional<std::uint32_t> previous_reference_end;
        for (std::uint32_t index = details.projection_offsets[row];
             index < details.projection_offsets[row + 1U];
             ++index) {
            const std::uint32_t template_begin =
                details.projection_template_starts[index];
            const std::uint32_t template_end =
                details.projection_template_ends[index];
            const std::uint32_t mapped_begin =
                details.projection_reference_starts[index];
            if (template_begin >= template_end || template_end > template_length) {
                throw ProtocolError("projection run has an invalid template interval");
            }
            const std::uint64_t mapped_end_wide =
                static_cast<std::uint64_t>(mapped_begin)
                + (template_end - template_begin);
            if (mapped_end_wide > std::numeric_limits<std::uint32_t>::max()
                || mapped_begin < reference_begin
                || mapped_end_wide > reference_end) {
                throw ProtocolError("projection run exceeds its reference envelope");
            }
            const auto mapped_end = static_cast<std::uint32_t>(mapped_end_wide);
            if (previous_template_end) {
                if (template_begin < *previous_template_end
                    || mapped_begin < *previous_reference_end
                    || mapped_begin <= *previous_reference_begin) {
                    throw ProtocolError("projection runs are not ordered");
                }
                if (template_begin == *previous_template_end
                    && mapped_begin == *previous_reference_end) {
                    throw ProtocolError("projection runs are not maximal");
                }
            }
            for (std::uint32_t offset = template_begin; offset < template_end;
                 ++offset) {
                if (projection_cover[offset] != 0U) {
                    throw ProtocolError("projection runs overlap");
                }
                projection_cover[offset] = 1U;
                mapped_positions[offset] = mapped_begin + (offset - template_begin);
            }
            previous_template_end = template_end;
            previous_reference_begin = mapped_begin;
            previous_reference_end = mapped_end;
        }

        std::optional<std::uint32_t> previous_variant_index;
        for (std::uint32_t index = details.variant_offsets[row];
             index < details.variant_offsets[row + 1U];
             ++index) {
            const std::uint32_t variant_index = details.variant_indices[index];
            if (variant_index == no_reference_position
                || (previous_variant_index && variant_index <= *previous_variant_index)) {
                throw ProtocolError("event ids must be strictly increasing");
            }
            previous_variant_index = variant_index;
            const std::uint8_t raw_source = details.variant_sources[index];
            if (raw_source < static_cast<std::uint8_t>(VariantSource::vcf)
                || raw_source
                    > static_cast<std::uint8_t>(
                        VariantSource::asm_profile)) {
                throw ProtocolError("variant source is invalid");
            }
            const std::uint32_t id_start = details.variant_id_offsets[index];
            const std::uint32_t id_end = details.variant_id_offsets[index + 1U];
            if (id_start == id_end) {
                throw ProtocolError("variant ID must not be empty");
            }
            const std::string_view id(
                reinterpret_cast<const char *>(details.variant_ids.data() + id_start),
                id_end - id_start);
            if (!valid_utf8(id)) {
                throw ProtocolError("variant ID must be valid UTF-8");
            }
            const std::uint8_t raw_kind = details.variant_kinds[index];
            if (!valid_variant_kind(raw_kind)) {
                throw ProtocolError("event kind is invalid");
            }
            const VariantKind kind = static_cast<VariantKind>(raw_kind);
            const std::uint8_t phased = details.variant_phased_haplotypes[index];
            if (phased != 255U && phased != haplotype) {
                throw ProtocolError("event phased haplotype disagrees with fragment");
            }
            const std::uint32_t event_reference_begin =
                details.variant_reference_starts[index];
            const std::uint32_t event_reference_end =
                details.variant_reference_ends[index];
            const std::uint32_t event_template_begin =
                details.variant_template_starts[index];
            const std::uint32_t event_template_end =
                details.variant_template_ends[index];
            if (event_reference_begin > event_reference_end
                || event_reference_begin < reference_begin
                || event_reference_end > reference_end) {
                throw ProtocolError("event reference interval exceeds its fragment");
            }
            if (event_template_begin > event_template_end
                || event_template_end > template_length) {
                throw ProtocolError("event template interval exceeds its fragment");
            }
            const std::uint32_t reference_span =
                event_reference_end - event_reference_begin;
            const std::uint32_t template_span =
                event_template_end - event_template_begin;
            const std::uint32_t ref_begin = details.variant_ref_offsets[index];
            const std::uint32_t ref_end = details.variant_ref_offsets[index + 1U];
            const std::uint32_t alt_begin = details.variant_alt_offsets[index];
            const std::uint32_t alt_end = details.variant_alt_offsets[index + 1U];
            if (ref_end - ref_begin != reference_span) {
                throw ProtocolError("event REF bases disagree with reference span");
            }
            if (alt_end - alt_begin != template_span) {
                throw ProtocolError("event ALT bases disagree with template span");
            }
            if ((kind == VariantKind::snv
                 && (reference_span == 0U || reference_span != template_span))
                || (kind == VariantKind::insertion
                    && (reference_span != 0U || template_span == 0U))
                || (kind == VariantKind::deletion
                    && (reference_span == 0U || template_span != 0U))) {
                throw ProtocolError("event kind and spans are inconsistent");
            }

            const auto previous_mapped = [&]() -> std::optional<std::uint32_t> {
                for (std::uint32_t offset = event_template_begin; offset != 0U;) {
                    --offset;
                    if (mapped_positions[offset] != no_reference_position) {
                        return mapped_positions[offset];
                    }
                }
                return std::nullopt;
            }();
            const auto next_mapped = [&]() -> std::optional<std::uint32_t> {
                const std::uint32_t begin = kind == VariantKind::deletion
                    ? event_template_begin
                    : event_template_end;
                for (std::uint32_t offset = begin; offset < template_length;
                     ++offset) {
                    if (mapped_positions[offset] != no_reference_position) {
                        return mapped_positions[offset];
                    }
                }
                return std::nullopt;
            }();
            if (kind == VariantKind::insertion) {
                if ((!previous_mapped && event_reference_begin != reference_begin)
                    || (previous_mapped
                        && *previous_mapped >= event_reference_begin)
                    || (!next_mapped && event_reference_begin != reference_end)
                    || (next_mapped && *next_mapped < event_reference_begin)) {
                    throw ProtocolError("insertion anchor disagrees with projection");
                }
            } else if (kind == VariantKind::deletion) {
                const bool maps_deleted_base = std::any_of(
                    mapped_positions.begin(),
                    mapped_positions.end(),
                    [event_reference_begin, event_reference_end](
                        std::uint32_t position) {
                        return position != no_reference_position
                            && position >= event_reference_begin
                            && position < event_reference_end;
                    });
                if ((previous_mapped
                     && *previous_mapped >= event_reference_begin)
                    || (next_mapped && *next_mapped < event_reference_end)
                    || maps_deleted_base) {
                    throw ProtocolError("deletion boundary disagrees with projection");
                }
            }

            for (std::uint32_t relative = 0; relative < template_span;
                 ++relative) {
                const std::uint32_t offset = event_template_begin + relative;
                if (event_cover[offset] != 0U) {
                    throw ProtocolError("event template spans overlap");
                }
                event_cover[offset] = 1U;
                if (batch.template_bases[flat_template_begin + offset]
                    != details.variant_alt_bases[alt_begin + relative]) {
                    throw ProtocolError("event ALT bases disagree with template");
                }
                if (kind == VariantKind::snv) {
                    if (mapped_positions[offset]
                        != event_reference_begin + relative) {
                        throw ProtocolError("SNV projection disagrees with event");
                    }
                } else if (kind == VariantKind::insertion) {
                    if (projection_cover[offset] != 0U) {
                        throw ProtocolError("insertion overlaps a projection run");
                    }
                    insertion_cover[offset] = 1U;
                }
            }
        }

        for (std::uint32_t offset = 0; offset < template_length; ++offset) {
            if (projection_cover[offset] + insertion_cover[offset] != 1U) {
                throw ProtocolError(
                    "projection and insertions do not cover template exactly");
            }
        }
        for (std::uint32_t index = batch.site_offsets[row];
             index < batch.site_offsets[row + 1U];
             ++index) {
            if (details.site_reference_positions[index]
                != mapped_positions[batch.site_template_offsets[index]]) {
                throw ProtocolError("site reference position disagrees with projection");
            }
        }

        std::vector<std::uint8_t> observed_n(template_length, 0U);
        std::optional<std::uint32_t> previous_n;
        for (std::uint32_t index = details.original_n_offsets[row];
             index < details.original_n_offsets[row + 1U];
             ++index) {
            const std::uint32_t offset = details.original_n_template_offsets[index];
            if (offset >= template_length || (previous_n && offset <= *previous_n)) {
                throw ProtocolError("original-N offsets are invalid");
            }
            previous_n = offset;
            observed_n[offset] = 1U;
            if (batch.template_bases[flat_template_begin + offset] != 4U) {
                throw ProtocolError("PRESERVE_N provenance does not point to N");
            }
        }
        for (std::uint32_t offset = 0; offset < template_length; ++offset) {
            const bool is_n = batch.template_bases[flat_template_begin + offset] == 4U;
            if ((observed_n[offset] != 0U) != is_n) {
                throw ProtocolError("PRESERVE_N provenance is incomplete");
            }
        }
    }
}

std::pair<std::uint8_t, std::vector<std::uint8_t>> encode_batch_payload(
    const FragmentBatch &batch,
    const Header &header)
{
    validate_batch(batch, header);
    std::vector<std::uint8_t> payload;
    Encoder encoder(payload);
    encoder.u32(batch.first_fragment_ordinal);
    encoder.u32(batch.fragment_count());
    encoder.u32(batch.template_base_count());
    encoder.u32(batch.mate_count());
    encoder.u32(batch.methylation_site_count());
    encoder.u32s(batch.contig_indices);
    encoder.u32s(batch.reference_starts);
    encoder.u32s(batch.reference_ends);
    encoder.u32s(batch.template_offsets);
    encoder.u32s(batch.mate_offsets);
    encoder.u32s(batch.site_offsets);
    encoder.u32s(batch.mate_template_starts);
    encoder.u32s(batch.mate_template_ends);
    encoder.u32s(batch.site_template_offsets);
    encoder.f32s(batch.site_probabilities);
    encoder.raw(batch.haplotypes);
    encoder.raw(batch.capture_strands);
    encoder.raw(batch.mate_indices);
    encoder.raw(batch.mate_reverse_complements);
    encoder.raw(batch.site_contexts);
    encoder.raw(batch.methylation_sources);
    encoder.raw(batch.site_alleles);
    encoder.raw(batch.template_bases);

    std::uint8_t flags = 0U;
    if (batch.details) {
        flags = details_present;
        const FragmentDetails &details = *batch.details;
        encoder.align4();
        encoder.u32(static_cast<std::uint32_t>(
            details.projection_template_starts.size()));
        encoder.u32(static_cast<std::uint32_t>(details.variant_indices.size()));
        encoder.u32(static_cast<std::uint32_t>(details.variant_ids.size()));
        encoder.u32(static_cast<std::uint32_t>(details.variant_ref_bases.size()));
        encoder.u32(static_cast<std::uint32_t>(details.variant_alt_bases.size()));
        encoder.u32(static_cast<std::uint32_t>(
            details.original_n_template_offsets.size()));
        encoder.u32s(details.projection_offsets);
        encoder.u32s(details.variant_offsets);
        encoder.u32s(details.original_n_offsets);
        encoder.u32s(details.projection_template_starts);
        encoder.u32s(details.projection_template_ends);
        encoder.u32s(details.projection_reference_starts);
        encoder.u32s(details.variant_indices);
        encoder.u32s(details.variant_id_offsets);
        encoder.u32s(details.variant_reference_starts);
        encoder.u32s(details.variant_reference_ends);
        encoder.u32s(details.variant_template_starts);
        encoder.u32s(details.variant_template_ends);
        encoder.u32s(details.variant_ref_offsets);
        encoder.u32s(details.variant_alt_offsets);
        encoder.u32s(details.site_reference_positions);
        encoder.u32s(details.original_n_template_offsets);
        encoder.raw(details.variant_sources);
        encoder.raw(details.variant_kinds);
        encoder.raw(details.variant_phased_haplotypes);
        encoder.raw(details.variant_ids);
        encoder.raw(details.variant_ref_bases);
        encoder.raw(details.variant_alt_bases);
    }
    encoder.align4();
    if (payload.size() > maximum_frame_payload) {
        throw ProtocolError("fragment batch exceeds the 64 MiB payload limit");
    }
    return {flags, std::move(payload)};
}

void checked_add(
    std::uint64_t current,
    std::uint64_t increment,
    std::string_view name)
{
    if (increment > std::numeric_limits<std::uint64_t>::max() - current) {
        throw ProtocolError(std::string(name) + " overflows u64");
    }
}

std::vector<std::uint8_t> encode_trailer_payload(const Trailer &trailer)
{
    checked_size(
        "trailer.per_contig_fragment_counts",
        trailer.per_contig_fragment_counts);
    std::uint64_t per_contig_sum = 0;
    for (std::uint64_t count : trailer.per_contig_fragment_counts) {
        checked_add(per_contig_sum, count, "trailer per-contig count sum");
        per_contig_sum += count;
    }
    if (per_contig_sum != trailer.fragment_count) {
        throw ProtocolError("trailer per-contig counts do not sum to fragments");
    }
    std::vector<std::uint8_t> payload;
    Encoder encoder(payload);
    encoder.u64(trailer.fragment_count);
    encoder.u64(trailer.fragment_batch_count);
    encoder.u64(trailer.mate_count);
    encoder.u64(trailer.template_base_count);
    encoder.u64(trailer.methylation_site_count);
    encoder.u64(trailer.skipped_fragment_count);
    encoder.u32(static_cast<std::uint32_t>(
        trailer.per_contig_fragment_counts.size()));
    encoder.u64s(trailer.per_contig_fragment_counts);
    encoder.raw(trailer.stream_sha256);
    encoder.align4();
    return payload;
}

std::vector<std::uint8_t> encode_error_payload(const ErrorFrame &error)
{
    validate_string("error.message", error.message);
    std::vector<std::uint8_t> payload;
    Encoder encoder(payload);
    encoder.u32(error.error_code);
    encoder.string(error.message);
    encoder.align4();
    return payload;
}

std::vector<std::uint8_t> encode_frame(
    FrameType type,
    std::uint8_t flags,
    std::uint64_t sequence,
    const std::vector<std::uint8_t> &payload)
{
    if (payload.size() > maximum_frame_payload) {
        throw ProtocolError("frame exceeds the 64 MiB payload limit");
    }
    std::vector<std::uint8_t> frame;
    frame.reserve(16U + payload.size() + 4U);
    Encoder encoder(frame);
    encoder.u32(static_cast<std::uint32_t>(payload.size()));
    encoder.u8(static_cast<std::uint8_t>(type));
    encoder.u8(flags);
    encoder.u16(0U);
    encoder.u64(sequence);
    encoder.raw(payload);
    encoder.u32(protocol::crc32c(frame));
    return frame;
}

std::vector<std::uint8_t> encode_preamble()
{
    std::vector<std::uint8_t> preamble;
    Encoder encoder(preamble);
    encoder.raw(magic);
    encoder.u16(protocol_major);
    encoder.u16(protocol_minor);
    encoder.u32(preamble_flags);
    return preamble;
}

} // namespace

std::uint32_t crc32c(const std::uint8_t *data, std::size_t size) noexcept
{
    std::uint32_t crc = UINT32_C(0xffffffff);
    for (std::size_t index = 0; index < size; ++index) {
        crc = crc32c_table[(crc ^ data[index]) & UINT32_C(0xff)] ^ (crc >> 8);
    }
    return crc ^ UINT32_C(0xffffffff);
}

std::uint32_t crc32c(const std::vector<std::uint8_t> &data) noexcept
{
    return crc32c(data.data(), data.size());
}

std::uint32_t FragmentBatch::fragment_count() const
{
    return checked_size("batch.contig_indices", contig_indices);
}

std::uint32_t FragmentBatch::template_base_count() const
{
    return checked_size("batch.template_bases", template_bases);
}

std::uint32_t FragmentBatch::mate_count() const
{
    return checked_size("batch.mate_indices", mate_indices);
}

std::uint32_t FragmentBatch::methylation_site_count() const
{
    return checked_size("batch.site_template_offsets", site_template_offsets);
}

PreparedFragmentBatch prepare_fragment_batch(
    const Header &header,
    FragmentBatch batch)
{
    auto encoded = encode_batch_payload(batch, header);
    PreparedFragmentBatch prepared;
    prepared.first_ordinal_ = batch.first_fragment_ordinal;
    prepared.fragment_count_ = batch.fragment_count();
    prepared.mate_count_ = batch.mate_count();
    prepared.template_base_count_ = batch.template_base_count();
    prepared.methylation_site_count_ = batch.methylation_site_count();
    prepared.frame_flags_ = encoded.first;
    prepared.header_payload_sha256_ = crypto::sha256(encode_header_payload(header));
    prepared.per_contig_fragment_counts_.assign(header.contigs.size(), 0U);
    for (std::uint32_t contig_index : batch.contig_indices) {
        ++prepared.per_contig_fragment_counts_[contig_index];
    }
    prepared.payload_ = std::move(encoded.second);
    return prepared;
}

Writer::Writer(std::ostream &sink) : sink_(&sink) {}

void Writer::ensure_open() const
{
    if (failed_) {
        throw ProtocolError("protocol writer is poisoned after a prior failure");
    }
    if (complete_) {
        throw ProtocolError("protocol writer is already complete");
    }
}

void Writer::write_raw(const std::uint8_t *bytes, std::size_t size)
{
    if (size == 0U) {return;}
    if (size > static_cast<std::size_t>(
                   std::numeric_limits<std::streamsize>::max())) {
        throw ProtocolError("protocol write exceeds streamsize");
    }
    sink_->write(
        reinterpret_cast<const char *>(bytes),
        static_cast<std::streamsize>(size));
    if (!*sink_) {throw ProtocolError("failed while writing protocol stream");}
}

void Writer::write_raw(const std::vector<std::uint8_t> &bytes)
{
    write_raw(bytes.data(), bytes.size());
}

void Writer::write_header(const Header &header)
{
    ensure_open();
    if (header_) {
        poison();
        throw ProtocolError("protocol header was written twice");
    }
    try {
        const std::vector<std::uint8_t> payload = encode_header_payload(header);
        const std::vector<std::uint8_t> preamble = encode_preamble();
        const std::vector<std::uint8_t> frame =
            encode_frame(FrameType::header, 0U, 0U, payload);
        write_raw(preamble);
        write_raw(frame);
        digest_state_.update(preamble);
        digest_state_.update(frame);
        header_ = header;
        header_payload_sha256_ = crypto::sha256(payload);
        per_contig_fragment_counts_.assign(header.contigs.size(), 0U);
        next_sequence_ = 1U;
    } catch (...) {
        poison();
        throw;
    }
}

void Writer::write_batch(FragmentBatch batch)
{
    ensure_open();
    if (!header_) {
        poison();
        throw ProtocolError("protocol header must precede batches");
    }
    try {
        write_prepared_batch(prepare_fragment_batch(*header_, std::move(batch)));
    } catch (...) {
        poison();
        throw;
    }
}

void Writer::write_prepared_batch(PreparedFragmentBatch &&batch)
{
    ensure_open();
    if (!header_) {
        poison();
        throw ProtocolError("protocol header must precede batches");
    }
    try {
        if (batch.first_ordinal_ != next_ordinal_) {
            throw ProtocolError("prepared batch ordinal range is not consecutive");
        }
        if (batch.header_payload_sha256_ != header_payload_sha256_) {
            throw ProtocolError("prepared batch header identity mismatch");
        }
        if (batch.per_contig_fragment_counts_.size()
            != per_contig_fragment_counts_.size()) {
            throw ProtocolError("prepared batch contig cardinality mismatch");
        }
        const bool annotation_present =
            (batch.frame_flags_ & details_present) != 0U;
        if ((batch.frame_flags_ & ~details_present) != 0U
            || annotation_present != header_->has_details) {
            throw ProtocolError("prepared batch details flags disagree with header");
        }
        checked_add(fragment_count_, batch.fragment_count_, "fragment count");
        checked_add(fragment_batch_count_, 1U, "fragment batch count");
        checked_add(mate_count_, batch.mate_count_, "mate count");
        checked_add(template_base_count_,
                    batch.template_base_count_,
                    "template base count");
        checked_add(methylation_site_count_,
                    batch.methylation_site_count_,
                    "methylation site count");
        if (next_ordinal_ + batch.fragment_count_
            > std::numeric_limits<std::uint32_t>::max()) {
            throw ProtocolError("stream fragment count exceeds u32");
        }
        for (std::size_t index = 0;
             index < per_contig_fragment_counts_.size();
             ++index) {
            checked_add(per_contig_fragment_counts_[index],
                        batch.per_contig_fragment_counts_[index],
                        "per-contig fragment count");
        }
        const std::vector<std::uint8_t> frame = encode_frame(
            FrameType::fragment_batch,
            batch.frame_flags_,
            next_sequence_,
            batch.payload_);
        write_raw(frame);
        digest_state_.update(frame);
        ++next_sequence_;
        next_ordinal_ += batch.fragment_count_;
        fragment_count_ += batch.fragment_count_;
        ++fragment_batch_count_;
        mate_count_ += batch.mate_count_;
        template_base_count_ += batch.template_base_count_;
        methylation_site_count_ += batch.methylation_site_count_;
        for (std::size_t index = 0;
             index < per_contig_fragment_counts_.size();
             ++index) {
            per_contig_fragment_counts_[index] +=
                batch.per_contig_fragment_counts_[index];
        }
    } catch (...) {
        poison();
        throw;
    }
}

Trailer Writer::finish(std::uint64_t skipped_fragment_count)
{
    ensure_open();
    if (!header_) {
        poison();
        throw ProtocolError("protocol header must precede trailer");
    }
    try {
        Trailer trailer;
        trailer.fragment_count = fragment_count_;
        trailer.fragment_batch_count = fragment_batch_count_;
        trailer.mate_count = mate_count_;
        trailer.template_base_count = template_base_count_;
        trailer.methylation_site_count = methylation_site_count_;
        trailer.skipped_fragment_count = skipped_fragment_count;
        trailer.per_contig_fragment_counts = per_contig_fragment_counts_;
        trailer.stream_sha256 = digest_state_.digest();
        const std::vector<std::uint8_t> frame = encode_frame(
            FrameType::trailer,
            0U,
            next_sequence_,
            encode_trailer_payload(trailer));
        write_raw(frame);
        ++next_sequence_;
        complete_ = true;
        return trailer;
    } catch (...) {
        poison();
        throw;
    }
}

void Writer::write_error(const ErrorFrame &error)
{
    ensure_open();
    if (!header_) {
        poison();
        throw ProtocolError("protocol header must precede error frame");
    }
    try {
        const std::vector<std::uint8_t> frame = encode_frame(
            FrameType::error,
            0U,
            next_sequence_,
            encode_error_payload(error));
        write_raw(frame);
        ++next_sequence_;
        complete_ = true;
    } catch (...) {
        poison();
        throw;
    }
}

} // namespace htsim::protocol

// ---- adapter --------------------------------------------------------

namespace htsim::protocol {
namespace {

[[noreturn]] void fail(const std::string &message)
{
    throw ProtocolError("protocol adapter: " + message);
}

std::uint32_t checked_u32(std::uint64_t value, const char *field)
{
    if (value > std::numeric_limits<std::uint32_t>::max()) {
        fail(std::string(field) + " exceeds u32");
    }
    return static_cast<std::uint32_t>(value);
}

std::uint32_t checked_size(std::size_t value, const char *field)
{
    return checked_u32(static_cast<std::uint64_t>(value), field);
}

void append_prefix(
    std::vector<std::uint32_t> &prefix,
    std::size_t size,
    const char *field)
{
    prefix.push_back(checked_size(size, field));
}

std::uint32_t wire_reference_position(std::int64_t position)
{
    if (position == -1) {return no_reference_position;}
    if (position < 0) {fail("reference position is below -1");}
    return checked_u32(static_cast<std::uint64_t>(position), "reference position");
}

std::pair<std::uint32_t, std::uint32_t> event_template_interval(
    const model::Fragment &fragment,
    const model::Variant &event)
{
    if (event.kind == model::VariantKind::deletion) {
        for (std::size_t offset = 0; offset < fragment.reference_positions.size();
             ++offset) {
            const std::int64_t position = fragment.reference_positions[offset];
            if (position >= 0
                && static_cast<std::uint64_t>(position) >= event.reference_end) {
                const auto boundary = checked_size(offset, "deletion boundary");
                return {boundary, boundary};
            }
        }
        const auto boundary = checked_size(
            fragment.template_bases.size(), "deletion boundary");
        return {boundary, boundary};
    }

    std::optional<std::size_t> first;
    std::size_t previous = 0U;
    std::size_t count = 0U;
    for (std::size_t offset = 0; offset < fragment.base_variant_indices.size();
         ++offset) {
        if (fragment.base_variant_indices[offset] != event.index) {continue;}
        if (!first) {
            first = offset;
        } else if (offset != previous + 1U) {
            fail("one event occupies non-contiguous template bases");
        }
        previous = offset;
        ++count;
    }
    if (!first || count != event.alt_bases.size()) {
        fail("event ALT span is incomplete in base_variant_indices");
    }
    return {
        checked_size(*first, "event template begin"),
        checked_size(*first + count, "event template end"),
    };
}

void append_projection_runs(
    FragmentDetails &details,
    const model::Fragment &fragment)
{
    bool run_open = false;
    std::size_t run_template_begin = 0U;
    std::uint32_t run_reference_begin = 0U;
    std::uint32_t previous_reference = 0U;

    const auto close_run = [&](std::size_t template_end) {
        if (!run_open) {return;}
        details.projection_template_starts.push_back(
            checked_size(run_template_begin, "projection template begin"));
        details.projection_template_ends.push_back(
            checked_size(template_end, "projection template end"));
        details.projection_reference_starts.push_back(run_reference_begin);
        run_open = false;
    };

    for (std::size_t offset = 0; offset < fragment.reference_positions.size();
         ++offset) {
        const std::int64_t raw_position = fragment.reference_positions[offset];
        if (raw_position == -1) {
            close_run(offset);
            continue;
        }
        const std::uint32_t position = wire_reference_position(raw_position);
        if (!run_open) {
            run_template_begin = offset;
            run_reference_begin = position;
            run_open = true;
        } else if (position != previous_reference + 1U) {
            close_run(offset);
            run_template_begin = offset;
            run_reference_begin = position;
            run_open = true;
        }
        previous_reference = position;
    }
    close_run(fragment.reference_positions.size());
}

void append_variants(
    FragmentDetails &details,
    const model::Fragment &fragment)
{
    std::vector<const model::Variant *> variants;
    variants.reserve(fragment.variants.size());
    for (const model::Variant &variant : fragment.variants) {
        variants.push_back(&variant);
    }
    std::sort(
        variants.begin(), variants.end(),
        [](const model::Variant *left,
           const model::Variant *right) {
            return left->index < right->index;
        });

    std::optional<std::uint32_t> previous_id;
    for (const model::Variant *variant : variants) {
        if (variant->index == model::no_variant_index
            || (previous_id && variant->index <= *previous_id)) {
            fail("variant indices are not unique u32 values");
        }
        if (variant->id.empty()) {fail("variant ID is empty");}
        previous_id = variant->index;
        const auto interval = event_template_interval(fragment, *variant);
        details.variant_indices.push_back(variant->index);
        details.variant_ids.insert(
            details.variant_ids.end(), variant->id.begin(), variant->id.end());
        append_prefix(
            details.variant_id_offsets,
            details.variant_ids.size(), "variant IDs");
        details.variant_sources.push_back(
            static_cast<std::uint8_t>(variant->source));
        details.variant_kinds.push_back(static_cast<std::uint8_t>(variant->kind));
        details.variant_phased_haplotypes.push_back(variant->phased_haplotype);
        details.variant_reference_starts.push_back(
            checked_u32(variant->reference_start, "variant reference start"));
        details.variant_reference_ends.push_back(
            checked_u32(variant->reference_end, "variant reference end"));
        details.variant_template_starts.push_back(interval.first);
        details.variant_template_ends.push_back(interval.second);
        details.variant_ref_bases.insert(
            details.variant_ref_bases.end(),
            variant->ref_bases.begin(), variant->ref_bases.end());
        details.variant_alt_bases.insert(
            details.variant_alt_bases.end(),
            variant->alt_bases.begin(), variant->alt_bases.end());
        append_prefix(
            details.variant_ref_offsets,
            details.variant_ref_bases.size(), "variant REF bases");
        append_prefix(
            details.variant_alt_offsets,
            details.variant_alt_bases.size(), "variant ALT bases");
    }
}

void validate_fragment_shape(
    const model::Fragment &fragment,
    bool include_annotations)
{
    if (fragment.template_bases.empty()) {
        fail("typed fragment template is empty");
    }
    const bool has_reference_positions = !fragment.reference_positions.empty();
    const bool has_base_variant_indices = !fragment.base_variant_indices.empty();
    if (has_reference_positions != has_base_variant_indices
        || (has_reference_positions
            && (fragment.reference_positions.size()
                    != fragment.template_bases.size()
                || fragment.base_variant_indices.size()
                    != fragment.template_bases.size()))) {
        fail("typed fragment arrays have inconsistent lengths");
    }
    if (!has_reference_positions) {
        if (include_annotations) {
            fail("Full Details fragment omitted typed projection arrays");
        }
        if (!fragment.variants.empty()) {
            fail("compact fragment retained unprojectable variant variants");
        }
    } else {
        std::optional<std::int64_t> previous_position;
        for (const std::int64_t position : fragment.reference_positions) {
            if (position < -1) {fail("reference position is below -1");}
            if (position >= 0) {
                if (previous_position && position <= *previous_position) {
                    fail("mapped reference positions are not increasing");
                }
                previous_position = position;
            }
        }
    }
    std::optional<std::uint32_t> previous_site_offset;
    for (std::size_t index = 0; index < fragment.methylation_sites.size(); ++index) {
        const model::MethylationSite &site = fragment.methylation_sites[index];
        if (site.site_index != index
            || site.template_offset >= fragment.template_bases.size()
            || (previous_site_offset
                && site.template_offset <= *previous_site_offset)
            || (has_reference_positions
                && site.reference_pos
                    != fragment.reference_positions[site.template_offset])) {
            fail("methylation sites are not canonical typed rows");
        }
        previous_site_offset = site.template_offset;
    }
}

} // namespace

FragmentBatch make_fragment_batch(
    const Header &header,
    const std::vector<model::Fragment> &fragments)
{
    if (fragments.empty()) {fail("fragment batch is empty");}

    FragmentBatch batch;
    batch.first_fragment_ordinal = checked_u32(
        fragments.front().fragment_ordinal, "first fragment ordinal");
    batch.template_offsets.push_back(0U);
    batch.mate_offsets.push_back(0U);
    batch.site_offsets.push_back(0U);

    const bool include_details = header.has_details;
    if (include_details) {
        batch.details.emplace();
        batch.details->projection_offsets.push_back(0U);
        batch.details->variant_offsets.push_back(0U);
        batch.details->original_n_offsets.push_back(0U);
        batch.details->variant_id_offsets.push_back(0U);
        batch.details->variant_ref_offsets.push_back(0U);
        batch.details->variant_alt_offsets.push_back(0U);
    }

    for (std::size_t row = 0; row < fragments.size(); ++row) {
        const model::Fragment &fragment = fragments[row];
        validate_fragment_shape(fragment, include_details);
        const std::uint64_t expected_ordinal =
            static_cast<std::uint64_t>(batch.first_fragment_ordinal) + row;
        if (fragment.fragment_ordinal != expected_ordinal) {
            fail("fragment ordinals are not consecutive");
        }

        batch.contig_indices.push_back(fragment.contig_index);
        batch.reference_starts.push_back(
            checked_u32(fragment.reference_start, "fragment reference begin"));
        batch.reference_ends.push_back(
            checked_u32(fragment.reference_end, "fragment reference end"));
        batch.haplotypes.push_back(fragment.haplotype);
        batch.capture_strands.push_back(
            static_cast<std::uint8_t>(fragment.capture_strand));
        batch.template_bases.insert(
            batch.template_bases.end(),
            fragment.template_bases.begin(), fragment.template_bases.end());
        append_prefix(
            batch.template_offsets,
            batch.template_bases.size(), "template bases");

        for (const model::Mate &mate : fragment.mates) {
            batch.mate_indices.push_back(mate.mate_index);
            batch.mate_reverse_complements.push_back(
                static_cast<std::uint8_t>(mate.reverse_complement));
            batch.mate_template_starts.push_back(mate.template_start);
            batch.mate_template_ends.push_back(mate.template_end);
        }
        append_prefix(batch.mate_offsets, batch.mate_indices.size(), "mates");

        for (const model::MethylationSite &site
             : fragment.methylation_sites) {
            batch.site_template_offsets.push_back(site.template_offset);
            batch.site_probabilities.push_back(site.methylation_probability);
            batch.site_contexts.push_back(
                static_cast<std::uint8_t>(site.context));
            batch.methylation_sources.push_back(
                static_cast<std::uint8_t>(site.methylation_source));
            batch.site_alleles.push_back(
                static_cast<std::uint8_t>(site.allele));
            if (include_details) {
                batch.details->site_reference_positions.push_back(
                    wire_reference_position(site.reference_pos));
            }
        }
        append_prefix(
            batch.site_offsets,
            batch.site_template_offsets.size(), "methylation sites");

        if (include_details) {
            append_projection_runs(*batch.details, fragment);
            append_prefix(
                batch.details->projection_offsets,
                batch.details->projection_template_starts.size(),
                "projection runs");
            append_variants(*batch.details, fragment);
            append_prefix(
                batch.details->variant_offsets,
                batch.details->variant_indices.size(), "variants");
            for (std::size_t offset = 0;
                 offset < fragment.template_bases.size(); ++offset) {
                if (fragment.template_bases[offset] == 4U) {
                    batch.details->original_n_template_offsets.push_back(
                        checked_size(offset, "original-N template offset"));
                }
            }
            append_prefix(
                batch.details->original_n_offsets,
                batch.details->original_n_template_offsets.size(),
                "original-N offsets");
        }
    }
    return batch;
}

} // namespace htsim::protocol

// ---- emitter --------------------------------------------------------

namespace htsim::protocol {
namespace {

constexpr std::size_t maximum_fragments_per_batch = 4096U;
using PreparedBatchGroup = std::vector<PreparedFragmentBatch>;

void append_prepared_batches(
    const Header &header,
    std::vector<model::Fragment> fragments,
    PreparedBatchGroup &prepared)
{
    try {
        prepared.push_back(prepare_fragment_batch(
            header, make_fragment_batch(header, fragments)));
        return;
    } catch (const ProtocolError &error) {
        if (fragments.size() <= 1U
            || std::string_view(error.what())
                != "fragment batch exceeds the 64 MiB payload limit") {
            throw;
        }
    }

    const std::size_t middle = fragments.size() / 2U;
    std::vector<model::Fragment> left;
    std::vector<model::Fragment> right;
    left.reserve(middle);
    right.reserve(fragments.size() - middle);
    std::move(
        fragments.begin(), fragments.begin() + static_cast<std::ptrdiff_t>(middle),
        std::back_inserter(left));
    std::move(
        fragments.begin() + static_cast<std::ptrdiff_t>(middle), fragments.end(),
        std::back_inserter(right));
    append_prepared_batches(header, std::move(left), prepared);
    append_prepared_batches(header, std::move(right), prepared);
}

PreparedBatchGroup prepare_batches(
    const Header &header,
    std::vector<model::Fragment> fragments)
{
    PreparedBatchGroup prepared;
    append_prepared_batches(header, std::move(fragments), prepared);
    return prepared;
}

class EncodePool {
public:
    explicit EncodePool(std::uint32_t worker_count)
    {
        if (worker_count == 0U) {
            throw ProtocolError("encode pool requires at least one worker");
        }
        threads_.reserve(worker_count);
        try {
            for (std::uint32_t index = 0U; index < worker_count; ++index) {
                threads_.emplace_back([this] {worker_loop();});
            }
        } catch (...) {
            shutdown();
            throw;
        }
    }

    EncodePool(const EncodePool &) = delete;
    EncodePool &operator=(const EncodePool &) = delete;

    ~EncodePool() {shutdown();}

    std::future<PreparedBatchGroup> submit(
        const Header &header,
        std::vector<model::Fragment> fragments)
    {
        const Header *header_pointer = &header;
        std::packaged_task<PreparedBatchGroup()> task(
            [header_pointer, fragments = std::move(fragments)]() mutable {
                return prepare_batches(*header_pointer, std::move(fragments));
            });
        auto result = task.get_future();
        {
            std::lock_guard<std::mutex> lock(mutex_);
            if (stopping_) {
                throw ProtocolError("encode pool is stopping");
            }
            tasks_.push(std::move(task));
        }
        ready_.notify_one();
        return result;
    }

private:
    using Task = std::packaged_task<PreparedBatchGroup()>;

    void shutdown() noexcept
    {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stopping_ = true;
        }
        ready_.notify_all();
        for (std::thread &thread : threads_) {
            if (thread.joinable()) {thread.join();}
        }
        threads_.clear();
    }

    void worker_loop()
    {
        while (true) {
            Task task;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                ready_.wait(lock, [this] {
                    return stopping_ || !tasks_.empty();
                });
                if (tasks_.empty()) {
                    if (stopping_) {return;}
                    continue;
                }
                task = std::move(tasks_.front());
                tasks_.pop();
            }
            task();
        }
    }

    std::mutex mutex_;
    std::condition_variable ready_;
    std::queue<Task> tasks_;
    std::vector<std::thread> threads_;
    bool stopping_ = false;
};

} // namespace

class BatchEmitter::Implementation {
public:
    Implementation(
        Writer &writer,
        const Header &header,
        std::uint32_t worker_count,
        std::uint32_t fragments_per_batch)
        : writer_(&writer), header_(&header),
          fragments_per_batch_(fragments_per_batch),
          maximum_pending_batches_(
              static_cast<std::size_t>(worker_count) * 2U)
    {
        if (worker_count == 0U || fragments_per_batch_ == 0U
            || fragments_per_batch_ > maximum_fragments_per_batch) {
            throw ProtocolError(
                "worker and batch fragment counts are outside bounds");
        }
        current_.reserve(fragments_per_batch_);
        if (worker_count > 1U) {
            pool_ = std::make_unique<EncodePool>(worker_count);
        }
    }

    void write(model::Fragment fragment)
    {
        if (finished_) {
            throw ProtocolError("batch emitter is already finished");
        }
        current_.push_back(std::move(fragment));
        if (current_.size() == fragments_per_batch_) {submit_current();}
    }

    void finish()
    {
        if (finished_) {
            throw ProtocolError("batch emitter was finished twice");
        }
        submit_current();
        while (!pending_.empty()) {write_oldest();}
        finished_ = true;
    }

private:
    void write_group(PreparedBatchGroup prepared)
    {
        for (PreparedFragmentBatch &batch : prepared) {
            writer_->write_prepared_batch(std::move(batch));
        }
    }

    void submit_current()
    {
        if (current_.empty()) {return;}
        std::vector<model::Fragment> fragments;
        fragments.swap(current_);
        current_.reserve(fragments_per_batch_);
        if (!pool_) {
            write_group(prepare_batches(*header_, std::move(fragments)));
            return;
        }
        while (pending_.size() >= maximum_pending_batches_) {
            write_oldest();
        }
        pending_.push_back(pool_->submit(*header_, std::move(fragments)));
    }

    void write_oldest()
    {
        auto result = std::move(pending_.front());
        pending_.pop_front();
        write_group(result.get());
    }

    Writer *writer_;
    const Header *header_;
    std::size_t fragments_per_batch_;
    std::size_t maximum_pending_batches_;
    std::unique_ptr<EncodePool> pool_;
    std::vector<model::Fragment> current_;
    std::deque<std::future<PreparedBatchGroup>> pending_;
    bool finished_ = false;
};

BatchEmitter::BatchEmitter(
    Writer &writer,
    const Header &header,
    std::uint32_t worker_count,
    std::uint32_t fragments_per_batch)
    : implementation_(std::make_unique<Implementation>(
          writer, header, worker_count, fragments_per_batch))
{}

BatchEmitter::~BatchEmitter() = default;

void BatchEmitter::write(model::Fragment fragment)
{
    implementation_->write(std::move(fragment));
}

void BatchEmitter::finish()
{
    implementation_->finish();
}

} // namespace htsim::protocol
