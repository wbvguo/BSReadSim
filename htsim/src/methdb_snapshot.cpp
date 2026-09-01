#include "methdb.h"

#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <tuple>
#include <unordered_map>
#include <utility>

#include <sys/types.h>

#include <zlib.h>

namespace htsim::methdb {
namespace {

constexpr std::array<char, 8> footer_magic = {{
    'm', 'd', 'b', '2', 'e', 'n', 'd', '\0',
}};
constexpr std::uint8_t little_endian_marker = 1U;
constexpr std::uint64_t maximum_section_bytes = UINT64_C(1) << 34U;

enum class SectionType : std::uint8_t {
    baseline = 1,
    events = 2,
    reference_overlay = 3,
    insertion_overlay = 4,
    asm_layer = 5,
};

enum class OverlayState : std::uint8_t {
    absent = 0,
    inherit = 1,
    explicit_site = 2,
};

struct OverlayRow {
    std::uint32_t key = 0U;
    std::array<OverlayState, 2> states = {{
        OverlayState::absent,
        OverlayState::absent,
    }};
    std::array<RuntimeSite, 2> payloads = {{0U, 0U}};
    std::vector<std::uint32_t> causes;
};

class ByteWriter {
public:
    void u8(std::uint8_t value) {bytes_.push_back(value);}

    void u16(std::uint16_t value)
    {
        u8(static_cast<std::uint8_t>(value));
        u8(static_cast<std::uint8_t>(value >> 8U));
    }

    void u32(std::uint32_t value)
    {
        for (unsigned shift = 0U; shift < 32U; shift += 8U) {
            u8(static_cast<std::uint8_t>(value >> shift));
        }
    }

    void u64(std::uint64_t value)
    {
        for (unsigned shift = 0U; shift < 64U; shift += 8U) {
            u8(static_cast<std::uint8_t>(value >> shift));
        }
    }

    void varuint(std::uint64_t value)
    {
        do {
            std::uint8_t byte = static_cast<std::uint8_t>(
                value & UINT64_C(0x7f));
            value >>= 7U;
            if (value != 0U) {byte |= UINT8_C(0x80);}
            u8(byte);
        } while (value != 0U);
    }

    void raw(const void *data, std::size_t size)
    {
        const auto *first = static_cast<const std::uint8_t *>(data);
        bytes_.insert(bytes_.end(), first, first + size);
    }

    void string(std::string_view value)
    {
        varuint(value.size());
        raw(value.data(), value.size());
    }

    const std::vector<std::uint8_t> &bytes() const noexcept {return bytes_;}
    std::vector<std::uint8_t> take() {return std::move(bytes_);}

private:
    std::vector<std::uint8_t> bytes_;
};

class ByteReader {
public:
    explicit ByteReader(const std::vector<std::uint8_t> &bytes)
        : bytes_(bytes)
    {}

    std::uint8_t u8()
    {
        require(1U);
        return bytes_[offset_++];
    }

    std::uint16_t u16()
    {
        const std::uint16_t first = u8();
        return static_cast<std::uint16_t>(first | (u8() << 8U));
    }

    std::uint32_t u32()
    {
        std::uint32_t value = 0U;
        for (unsigned shift = 0U; shift < 32U; shift += 8U) {
            value |= static_cast<std::uint32_t>(u8()) << shift;
        }
        return value;
    }

    std::uint64_t u64()
    {
        std::uint64_t value = 0U;
        for (unsigned shift = 0U; shift < 64U; shift += 8U) {
            value |= static_cast<std::uint64_t>(u8()) << shift;
        }
        return value;
    }

    std::uint64_t varuint()
    {
        std::uint64_t value = 0U;
        for (unsigned index = 0U; index < 10U; ++index) {
            const std::uint8_t byte = u8();
            if (index == 9U && (byte & UINT8_C(0xfe)) != 0U) {
                throw SnapshotError("MethDB varuint exceeds uint64");
            }
            value |= static_cast<std::uint64_t>(byte & UINT8_C(0x7f))
                << (index * 7U);
            if ((byte & UINT8_C(0x80)) == 0U) {
                if (index != 0U && byte == 0U) {
                    throw SnapshotError("MethDB varuint is not canonical");
                }
                return value;
            }
        }
        throw SnapshotError("MethDB varuint is unterminated");
    }

    std::string string()
    {
        const std::uint64_t size = varuint();
        if (size > UINT64_C(1048576)
            || size > std::numeric_limits<std::size_t>::max()) {
            throw SnapshotError("MethDB string length is invalid");
        }
        require(static_cast<std::size_t>(size));
        const auto *first = reinterpret_cast<const char *>(
            bytes_.data() + offset_);
        std::string result(first, static_cast<std::size_t>(size));
        offset_ += static_cast<std::size_t>(size);
        return result;
    }

    std::vector<std::uint8_t> bytes(std::size_t size)
    {
        require(size);
        std::vector<std::uint8_t> result(
            bytes_.begin() + static_cast<std::ptrdiff_t>(offset_),
            bytes_.begin() + static_cast<std::ptrdiff_t>(offset_ + size));
        offset_ += size;
        return result;
    }

    void finish() const
    {
        if (offset_ != bytes_.size()) {
            throw SnapshotError("MethDB section has trailing canonical bytes");
        }
    }

private:
    void require(std::size_t size) const
    {
        if (size > bytes_.size() - offset_) {
            throw SnapshotError("MethDB canonical section is truncated");
        }
    }

    const std::vector<std::uint8_t> &bytes_;
    std::size_t offset_ = 0U;
};

void write_raw(std::ostream &output, const void *data, std::size_t size)
{
    output.write(
        static_cast<const char *>(data),
        static_cast<std::streamsize>(size));
    if (!output) {throw SnapshotError("failed while writing MethDB v2");}
}

void write_u8(std::ostream &output, std::uint8_t value)
{
    write_raw(output, &value, sizeof(value));
}

void write_u16(std::ostream &output, std::uint16_t value)
{
    std::array<std::uint8_t, 2> bytes = {{
        static_cast<std::uint8_t>(value),
        static_cast<std::uint8_t>(value >> 8U),
    }};
    write_raw(output, bytes.data(), bytes.size());
}

void write_u32(std::ostream &output, std::uint32_t value)
{
    std::array<std::uint8_t, 4> bytes = {};
    for (unsigned shift = 0U; shift < 32U; shift += 8U) {
        bytes[shift / 8U] = static_cast<std::uint8_t>(value >> shift);
    }
    write_raw(output, bytes.data(), bytes.size());
}

void write_u64(std::ostream &output, std::uint64_t value)
{
    std::array<std::uint8_t, 8> bytes = {};
    for (unsigned shift = 0U; shift < 64U; shift += 8U) {
        bytes[shift / 8U] = static_cast<std::uint8_t>(value >> shift);
    }
    write_raw(output, bytes.data(), bytes.size());
}

void read_raw(std::istream &input, void *data, std::size_t size)
{
    input.read(static_cast<char *>(data), static_cast<std::streamsize>(size));
    if (input.gcount() != static_cast<std::streamsize>(size)) {
        throw SnapshotError("MethDB v2 is truncated");
    }
}

std::uint8_t read_u8(std::istream &input)
{
    std::uint8_t value = 0U;
    read_raw(input, &value, sizeof(value));
    return value;
}

std::uint16_t read_u16(std::istream &input)
{
    const std::uint16_t first = read_u8(input);
    return static_cast<std::uint16_t>(first | (read_u8(input) << 8U));
}

std::uint32_t read_u32(std::istream &input)
{
    std::uint32_t value = 0U;
    for (unsigned shift = 0U; shift < 32U; shift += 8U) {
        value |= static_cast<std::uint32_t>(read_u8(input)) << shift;
    }
    return value;
}

std::uint64_t read_u64(std::istream &input)
{
    std::uint64_t value = 0U;
    for (unsigned shift = 0U; shift < 64U; shift += 8U) {
        value |= static_cast<std::uint64_t>(read_u8(input)) << shift;
    }
    return value;
}

std::uint8_t metadata_byte(RuntimeSite site)
{
    const std::uint16_t metadata = static_cast<std::uint16_t>(site);
    if ((metadata & UINT16_C(0xff00)) != 0U) {
        throw SnapshotError("runtime site uses unassigned metadata bits");
    }
    return static_cast<std::uint8_t>(metadata);
}

RuntimeSite payload_site(
    std::uint32_t key,
    std::uint16_t probability,
    std::uint8_t metadata)
{
    const RuntimeSite site = (static_cast<RuntimeSite>(key) << 32U)
        | (static_cast<RuntimeSite>(probability) << 16U)
        | metadata;
    (void)runtime_site_context(site);
    (void)runtime_site_source(site);
    (void)runtime_site_allele(site);
    return site;
}

void write_payload(ByteWriter &writer, RuntimeSite site)
{
    writer.u8(metadata_byte(site));
    writer.u16(runtime_site_probability(site));
}

RuntimeSite read_payload(ByteReader &reader, std::uint32_t key)
{
    const std::uint8_t metadata = reader.u8();
    const std::uint16_t probability = reader.u16();
    return payload_site(key, probability, metadata);
}

std::vector<RuntimeSite> merge_site_sets(
    const std::vector<RuntimeSite> &shared,
    const std::vector<RuntimeSite> &specific)
{
    std::vector<RuntimeSite> result;
    result.reserve(shared.size() + specific.size());
    std::size_t left = 0U;
    std::size_t right = 0U;
    while (left < shared.size() || right < specific.size()) {
        if (right == specific.size()
            || (left < shared.size()
                && runtime_site_key(shared[left])
                    < runtime_site_key(specific[right]))) {
            result.push_back(shared[left++]);
        } else if (left == shared.size()
                   || runtime_site_key(specific[right])
                       < runtime_site_key(shared[left])) {
            result.push_back(specific[right++]);
        } else {
            throw SnapshotError(
                "shared and haplotype runtime sites overlap");
        }
    }
    return result;
}

bool baseline_payload_matches(const CatalogSite &baseline, RuntimeSite site)
{
    return baseline.probability_u16 == runtime_site_probability(site)
        && baseline.context == runtime_site_context(site)
        && baseline.methylation_source == runtime_site_source(site);
}

std::vector<std::uint32_t> overlay_causes(
    std::uint32_t position,
    const std::vector<variant::Variant> &variants)
{
    std::vector<std::uint32_t> causes;
    for (std::size_t index = 0U; index < variants.size(); ++index) {
        const variant::Variant &event = variants[index];
        const std::uint64_t position_with_halo =
            static_cast<std::uint64_t>(position) + 2U;
        const std::uint64_t event_end_with_halo =
            static_cast<std::uint64_t>(event.reference_end) + 2U;
        if (position_with_halo >= event.reference_start
            && position <= event_end_with_halo) {
            causes.push_back(static_cast<std::uint32_t>(index));
        }
    }
    if (causes.empty()) {
        throw SnapshotError(
            "variant reference overlay has no causal prepared event");
    }
    return causes;
}

std::vector<OverlayRow> build_reference_overlays(
    const std::vector<CatalogSite> &baseline,
    const DiploidMethylationCatalog &catalog,
    const std::vector<variant::Variant> &variants)
{
    const DiploidRuntimeArrays &arrays = catalog.runtime_arrays();
    const std::array<std::vector<RuntimeSite>, 2> resolved = {{
        merge_site_sets(
            arrays.reference_shared,
            arrays.reference_haplotypes[0]),
        merge_site_sets(
            arrays.reference_shared,
            arrays.reference_haplotypes[1]),
    }};
    std::array<std::size_t, 2> hap_index = {{0U, 0U}};
    std::size_t baseline_index = 0U;
    std::vector<OverlayRow> result;

    while (baseline_index < baseline.size()
           || hap_index[0] < resolved[0].size()
           || hap_index[1] < resolved[1].size()) {
        std::uint32_t key = std::numeric_limits<std::uint32_t>::max();
        if (baseline_index < baseline.size()) {
            key = baseline[baseline_index].reference_position;
        }
        for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
            if (hap_index[haplotype] < resolved[haplotype].size()) {
                key = std::min(
                    key,
                    runtime_site_key(resolved[haplotype][hap_index[haplotype]]));
            }
        }
        const CatalogSite *baseline_site =
            baseline_index < baseline.size()
                && baseline[baseline_index].reference_position == key
            ? &baseline[baseline_index]
            : nullptr;
        OverlayRow row;
        row.key = key;
        for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
            const RuntimeSite *site =
                hap_index[haplotype] < resolved[haplotype].size()
                    && runtime_site_key(
                        resolved[haplotype][hap_index[haplotype]]) == key
                ? &resolved[haplotype][hap_index[haplotype]]
                : nullptr;
            if (site == nullptr) {
                row.states[haplotype] = OverlayState::absent;
            } else if (runtime_site_reference_equivalent(*site)) {
                if (baseline_site == nullptr
                    || !baseline_payload_matches(*baseline_site, *site)) {
                    throw SnapshotError(
                        "reference-equivalent site disagrees with baseline");
                }
                row.states[haplotype] = OverlayState::inherit;
            } else {
                row.states[haplotype] = OverlayState::explicit_site;
                row.payloads[haplotype] = *site;
            }
            if (site != nullptr) {++hap_index[haplotype];}
        }
        if (baseline_site != nullptr) {++baseline_index;}

        const bool unchanged = baseline_site != nullptr
            && row.states[0] == OverlayState::inherit
            && row.states[1] == OverlayState::inherit;
        const bool empty = baseline_site == nullptr
            && row.states[0] == OverlayState::absent
            && row.states[1] == OverlayState::absent;
        if (!unchanged && !empty) {
            row.causes = overlay_causes(key, variants);
            result.push_back(std::move(row));
        }
    }
    return result;
}

std::vector<OverlayRow> build_insertion_overlays(
    const DiploidMethylationCatalog &catalog,
    const std::vector<variant::Variant> &variants)
{
    const DiploidRuntimeArrays &arrays = catalog.runtime_arrays();
    const std::array<std::vector<RuntimeSite>, 2> resolved = {{
        merge_site_sets(
            arrays.insertion_shared,
            arrays.insertion_haplotypes[0]),
        merge_site_sets(
            arrays.insertion_shared,
            arrays.insertion_haplotypes[1]),
    }};
    std::array<std::size_t, 2> indices = {{0U, 0U}};
    std::vector<OverlayRow> result;
    while (indices[0] < resolved[0].size()
           || indices[1] < resolved[1].size()) {
        std::uint32_t key = std::numeric_limits<std::uint32_t>::max();
        for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
            if (indices[haplotype] < resolved[haplotype].size()) {
                key = std::min(
                    key,
                    runtime_site_key(resolved[haplotype][indices[haplotype]]));
            }
        }
        const std::uint32_t event_ordinal = key >> 2U;
        const std::uint8_t insertion_offset = static_cast<std::uint8_t>(key & 3U);
        if (event_ordinal >= variants.size()
            || variants[event_ordinal].kind != model::VariantKind::insertion
            || insertion_offset >= variants[event_ordinal].alt_bases.size()) {
            throw SnapshotError(
                "insertion overlay does not identify an embedded event base");
        }
        OverlayRow row;
        row.key = key;
        row.causes.push_back(event_ordinal);
        for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
            if (indices[haplotype] < resolved[haplotype].size()
                && runtime_site_key(
                    resolved[haplotype][indices[haplotype]]) == key) {
                row.states[haplotype] = OverlayState::explicit_site;
                row.payloads[haplotype] =
                    resolved[haplotype][indices[haplotype]++];
            } else {
                row.states[haplotype] = OverlayState::absent;
            }
        }
        result.push_back(std::move(row));
    }
    return result;
}

std::uint8_t pack_states(const OverlayRow &row)
{
    return static_cast<std::uint8_t>(row.states[0])
        | static_cast<std::uint8_t>(
            static_cast<std::uint8_t>(row.states[1]) << 2U);
}

std::array<OverlayState, 2> unpack_states(std::uint8_t packed)
{
    if ((packed & UINT8_C(0xf0)) != 0U) {
        throw SnapshotError("MethDB overlay state reserved bits are nonzero");
    }
    std::array<OverlayState, 2> result = {{
        static_cast<OverlayState>(packed & 3U),
        static_cast<OverlayState>((packed >> 2U) & 3U),
    }};
    for (const OverlayState state : result) {
        if (state != OverlayState::absent
            && state != OverlayState::inherit
            && state != OverlayState::explicit_site) {
            throw SnapshotError("MethDB overlay state is invalid");
        }
    }
    return result;
}

std::vector<std::uint8_t> encode_baseline(
    const std::vector<CatalogSite> &sites)
{
    ByteWriter writer;
    writer.u8(1U);
    std::uint32_t previous = 0U;
    bool first = true;
    for (const CatalogSite &site : sites) {
        if (!first && site.reference_position <= previous) {
            throw SnapshotError("MethDB baseline is not strictly ordered");
        }
        writer.varuint(first
            ? site.reference_position
            : site.reference_position - previous);
        writer.u8(metadata_byte(pack_runtime_site(
            site.reference_position,
            site.probability_u16,
            site.context,
            site.methylation_source,
            model::MethylationAllele::shared,
            true)));
        writer.u16(site.probability_u16);
        previous = site.reference_position;
        first = false;
    }
    return writer.take();
}

std::vector<std::uint8_t> encode_events(
    const std::vector<variant::Variant> &variants,
    std::uint32_t contig_index)
{
    ByteWriter writer;
    writer.u8(1U);
    std::uint32_t previous_start = 0U;
    std::uint32_t previous_end = 0U;
    bool first = true;
    for (const variant::Variant &event : variants) {
        if (event.contig_index != contig_index
            || (!first && event.reference_start < previous_start)
            || (!first && event.reference_start < previous_end)
            || event.reference_end < event.reference_start
            || event.id.size() > UINT32_C(1048576)
            || event.id.find_first_of("\t\r\n") != std::string::npos
            || event.id.find('\0') != std::string::npos
            || event.ref_bases.size() > variant::maximum_indel_bases
            || event.alt_bases.size() > variant::maximum_indel_bases
            || !model::is_haplotype_mask(
                static_cast<std::uint8_t>(event.alt_haplotypes))
            || (event.source != model::VariantSource::vcf
                && event.source != model::VariantSource::de_novo
                && event.source != model::VariantSource::asm_profile)
            || std::any_of(
                event.ref_bases.begin(),
                event.ref_bases.end(),
                [](std::uint8_t base) {return base > 3U;})
            || std::any_of(
                event.alt_bases.begin(),
                event.alt_bases.end(),
                [](std::uint8_t base) {return base > 3U;})) {
            throw SnapshotError("prepared variant event is not canonical");
        }
        switch (event.kind) {
        case model::VariantKind::snv:
            if (event.reference_end != event.reference_start + 1U
                || event.ref_bases.size() != 1U
                || event.alt_bases.size() != 1U) {
                throw SnapshotError("prepared MethDB SNV shape is invalid");
            }
            break;
        case model::VariantKind::insertion:
            if (event.reference_end != event.reference_start
                || !event.ref_bases.empty()
                || event.alt_bases.empty()) {
                throw SnapshotError(
                    "prepared MethDB insertion shape is invalid");
            }
            break;
        case model::VariantKind::deletion:
            if (event.reference_end <= event.reference_start
                || event.ref_bases.empty()
                || !event.alt_bases.empty()
                || event.ref_bases.size()
                    != event.reference_end - event.reference_start) {
                throw SnapshotError(
                    "prepared MethDB deletion shape is invalid");
            }
            break;
        default:
            throw SnapshotError("prepared MethDB event kind is invalid");
        }
        writer.varuint(first
            ? event.reference_start
            : event.reference_start - previous_start);
        writer.varuint(event.reference_end - event.reference_start);
        writer.u8(static_cast<std::uint8_t>(event.kind));
        writer.u8(static_cast<std::uint8_t>(event.alt_haplotypes));
        writer.u8(static_cast<std::uint8_t>(event.source));
        writer.string(event.id);
        writer.varuint(event.ref_bases.size());
        if (!event.ref_bases.empty()) {
            writer.raw(event.ref_bases.data(), event.ref_bases.size());
        }
        writer.varuint(event.alt_bases.size());
        if (!event.alt_bases.empty()) {
            writer.raw(event.alt_bases.data(), event.alt_bases.size());
        }
        previous_start = event.reference_start;
        previous_end = event.reference_end;
        first = false;
    }
    return writer.take();
}

std::vector<std::uint8_t> encode_overlays(
    const std::vector<OverlayRow> &rows,
    bool insertion)
{
    ByteWriter writer;
    writer.u8(1U);
    std::uint32_t previous = 0U;
    bool first = true;
    for (const OverlayRow &row : rows) {
        if (!first && row.key <= previous) {
            throw SnapshotError("MethDB overlay keys are not strictly ordered");
        }
        writer.varuint(first ? row.key : row.key - previous);
        writer.u8(pack_states(row));
        if (!insertion) {
            writer.varuint(row.causes.size());
            std::uint32_t previous_cause = 0U;
            bool first_cause = true;
            for (const std::uint32_t cause : row.causes) {
                if (!first_cause && cause <= previous_cause) {
                    throw SnapshotError("MethDB cause set is not canonical");
                }
                writer.varuint(first_cause ? cause : cause - previous_cause);
                previous_cause = cause;
                first_cause = false;
            }
        }
        for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
            if (row.states[haplotype] == OverlayState::explicit_site) {
                write_payload(writer, row.payloads[haplotype]);
            }
        }
        previous = row.key;
        first = false;
    }
    return writer.take();
}

std::vector<std::uint8_t> encode_asm(
    const std::vector<AsmRecord> &records)
{
    ByteWriter writer;
    writer.u8(1U);
    std::uint32_t previous = 0U;
    bool first = true;
    for (const AsmRecord &record : records) {
        if (!first && record.target_reference_position <= previous) {
            throw SnapshotError("MethDB ASM records are not strictly ordered");
        }
        writer.varuint(first
            ? record.target_reference_position
            : record.target_reference_position - previous);
        writer.varuint(record.linked_variant_position);
        writer.u16(record.reference_probability_u16);
        writer.u16(record.alternate_probability_u16);
        writer.u8(static_cast<std::uint8_t>(record.context));
        writer.u8(record.dinucleotide_second);
        writer.u8(record.linked_reference_base);
        writer.u8(record.linked_alternate_base);
        previous = record.target_reference_position;
        first = false;
    }
    return writer.take();
}

std::vector<std::uint8_t> compress_section(
    const std::vector<std::uint8_t> &raw)
{
    if (raw.size() > std::numeric_limits<uLong>::max()) {
        throw SnapshotError("MethDB section exceeds zlib input boundary");
    }
    uLongf compressed_size = compressBound(static_cast<uLong>(raw.size()));
    std::vector<std::uint8_t> compressed(compressed_size);
    const std::uint8_t empty = 0U;
    const Bytef *input = raw.empty() ? &empty : raw.data();
    const int status = compress2(
        compressed.data(),
        &compressed_size,
        input,
        static_cast<uLong>(raw.size()),
        Z_DEFAULT_COMPRESSION);
    if (status != Z_OK) {
        throw SnapshotError("cannot compress MethDB v2 section");
    }
    compressed.resize(compressed_size);
    return compressed;
}

struct SectionEntry {
    SectionType type = SectionType::baseline;
    std::uint32_t record_count = 0U;
    std::uint64_t raw_size = 0U;
    std::uint64_t compressed_size = 0U;
    crypto::Sha256Digest raw_sha256 = {};
    std::uint64_t payload_offset = 0U;
};

struct ContigDirectory {
    reference::ContigMetadata metadata;
    bool diploid = false;
    std::vector<SectionEntry> sections;
};

struct ArchiveDirectory {
    crypto::Sha256Digest binding = {};
    crypto::Sha256Digest content_sha256 = {};
    crypto::Sha256Digest file_sha256 = {};
    std::vector<ContigDirectory> contigs;
};

std::uint64_t stream_position(std::istream &input)
{
    const std::streampos position = input.tellg();
    if (position < 0) {throw SnapshotError("cannot locate MethDB section");}
    return static_cast<std::uint64_t>(position);
}

void skip_bytes(std::istream &input, std::uint64_t size)
{
    if (size > static_cast<std::uint64_t>(
                   std::numeric_limits<std::streamoff>::max())) {
        throw SnapshotError("MethDB section offset exceeds stream boundary");
    }
    input.seekg(static_cast<std::streamoff>(size), std::ios::cur);
    if (!input) {throw SnapshotError("MethDB section payload is truncated");}
}

crypto::Sha256Digest hash_file(const std::string &path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {throw SnapshotError("cannot open MethDB snapshot: " + path);}
    crypto::Sha256 hash;
    std::array<std::uint8_t, 65536> buffer = {};
    while (input) {
        input.read(
            reinterpret_cast<char *>(buffer.data()),
            static_cast<std::streamsize>(buffer.size()));
        const std::streamsize observed = input.gcount();
        if (observed > 0) {
            hash.update(buffer.data(), static_cast<std::size_t>(observed));
        }
    }
    if (!input.eof()) {throw SnapshotError("failed while hashing MethDB v2");}
    return hash.digest();
}

bool valid_section_type(std::uint8_t value) noexcept
{
    return value >= static_cast<std::uint8_t>(SectionType::baseline)
        && value <= static_cast<std::uint8_t>(SectionType::asm_layer);
}

ArchiveDirectory parse_archive(
    const std::string &path,
    const crypto::Sha256Digest *expected_binding,
    const std::vector<reference::ContigMetadata> *expected_catalog)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {throw SnapshotError("cannot open MethDB snapshot: " + path);}
    std::array<char, sizeof(methdb_magic) - 1U> magic = {};
    read_raw(input, magic.data(), magic.size());
    if (std::memcmp(magic.data(), methdb_magic, magic.size()) != 0) {
        throw SnapshotError("MethDB snapshot magic is invalid");
    }
    if (read_u8(input) != methdb_version) {
        throw SnapshotError("MethDB snapshot version is unsupported");
    }
    if (read_u8(input) != little_endian_marker
        || read_u8(input) != 0U
        || read_u16(input) != 0U) {
        throw SnapshotError("MethDB v2 header flags are unsupported");
    }

    ArchiveDirectory archive;
    read_raw(input, archive.binding.data(), archive.binding.size());
    if (expected_binding != nullptr && archive.binding != *expected_binding) {
        throw SnapshotError("MethDB profile binding is incompatible");
    }
    const std::uint32_t contig_count = read_u32(input);
    if (expected_catalog != nullptr
        && contig_count != expected_catalog->size()) {
        throw SnapshotError(
            "MethDB snapshot contig count disagrees with reference");
    }
    crypto::Sha256 content_hash;
    static constexpr char content_domain[] =
        "BSReadSim/MethDB/canonical-v2";
    content_hash.update(
        reinterpret_cast<const std::uint8_t *>(content_domain),
        sizeof(content_domain) - 1U);
    content_hash.update(archive.binding.data(), archive.binding.size());
    ByteWriter encoded_count;
    encoded_count.u32(contig_count);
    content_hash.update(encoded_count.bytes());
    archive.contigs.reserve(contig_count);
    for (std::uint32_t contig_index = 0U;
         contig_index < contig_count;
         ++contig_index) {
        const std::uint32_t name_size = read_u32(input);
        if (name_size == 0U || name_size > UINT32_C(1048576)) {
            throw SnapshotError("MethDB contig name length is invalid");
        }
        ContigDirectory contig;
        contig.metadata.name.resize(name_size);
        read_raw(input, contig.metadata.name.data(), name_size);
        if (contig.metadata.name.find_first_of("\t\r\n")
                != std::string::npos
            || contig.metadata.name.find('\0') != std::string::npos) {
            throw SnapshotError("MethDB contig name is invalid");
        }
        contig.metadata.length = read_u32(input);
        read_raw(
            input,
            contig.metadata.reference_sha256.data(),
            contig.metadata.reference_sha256.size());
        const std::uint8_t mode = read_u8(input);
        const std::uint8_t section_count = read_u8(input);
        if (read_u16(input) != 0U || mode > 1U) {
            throw SnapshotError("MethDB contig header is invalid");
        }
        contig.diploid = mode == 1U;
        const std::uint8_t expected_sections = contig.diploid ? 5U : 1U;
        if (section_count != expected_sections) {
            throw SnapshotError("MethDB contig section set is incomplete");
        }
        if (expected_catalog != nullptr) {
            const reference::ContigMetadata &expected =
                expected_catalog->at(contig_index);
            if (contig.metadata.name != expected.name
                || contig.metadata.length != expected.length
                || contig.metadata.reference_sha256
                    != expected.reference_sha256) {
                throw SnapshotError(
                    "MethDB contig identity disagrees with reference");
            }
        }
        ByteWriter contig_identity;
        contig_identity.u32(contig_index);
        contig_identity.string(contig.metadata.name);
        contig_identity.u32(
            static_cast<std::uint32_t>(contig.metadata.length));
        contig_identity.raw(
            contig.metadata.reference_sha256.data(),
            contig.metadata.reference_sha256.size());
        contig_identity.u8(contig.diploid ? 1U : 0U);
        content_hash.update(contig_identity.bytes());
        std::array<bool, 6> observed_types = {{
            false, false, false, false, false, false,
        }};
        contig.sections.reserve(section_count);
        for (std::uint8_t section_index = 0U;
             section_index < section_count;
             ++section_index) {
            const std::uint8_t type_value = read_u8(input);
            static constexpr std::array<SectionType, 5> canonical_order = {{
                SectionType::baseline,
                SectionType::events,
                SectionType::reference_overlay,
                SectionType::insertion_overlay,
                SectionType::asm_layer,
            }};
            if (!valid_section_type(type_value)
                || type_value != static_cast<std::uint8_t>(
                    canonical_order[section_index])
                || read_u8(input) != 0U
                || read_u16(input) != 0U
                || read_u32(input) != contig_index) {
                throw SnapshotError("MethDB section header is invalid");
            }
            if (observed_types[type_value]) {
                throw SnapshotError("MethDB contig has duplicate sections");
            }
            observed_types[type_value] = true;
            SectionEntry section;
            section.type = static_cast<SectionType>(type_value);
            section.record_count = read_u32(input);
            section.raw_size = read_u64(input);
            section.compressed_size = read_u64(input);
            read_raw(
                input,
                section.raw_sha256.data(),
                section.raw_sha256.size());
            if (section.raw_size == 0U
                || section.raw_size > maximum_section_bytes
                || section.compressed_size == 0U
                || section.compressed_size > maximum_section_bytes) {
                throw SnapshotError("MethDB section size is invalid");
            }
            ByteWriter section_identity;
            section_identity.u8(type_value);
            section_identity.u32(contig_index);
            section_identity.u32(section.record_count);
            section_identity.u64(section.raw_size);
            content_hash.update(section_identity.bytes());
            content_hash.update(
                section.raw_sha256.data(), section.raw_sha256.size());
            section.payload_offset = stream_position(input);
            skip_bytes(input, section.compressed_size);
            contig.sections.push_back(section);
        }
        if (!observed_types[static_cast<std::uint8_t>(SectionType::baseline)]
            || (contig.diploid
                && (!observed_types[static_cast<std::uint8_t>(SectionType::events)]
                    || !observed_types[static_cast<std::uint8_t>(
                        SectionType::reference_overlay)]
                    || !observed_types[static_cast<std::uint8_t>(
                        SectionType::insertion_overlay)]
                    || !observed_types[static_cast<std::uint8_t>(
                        SectionType::asm_layer)]))) {
            throw SnapshotError("MethDB contig section types are incomplete");
        }
        archive.contigs.push_back(std::move(contig));
    }
    std::array<char, footer_magic.size()> observed_footer = {};
    read_raw(input, observed_footer.data(), observed_footer.size());
    if (observed_footer != footer_magic) {
        throw SnapshotError("MethDB v2 footer magic is invalid");
    }
    read_raw(
        input,
        archive.content_sha256.data(),
        archive.content_sha256.size());
    if (archive.content_sha256 != content_hash.digest()) {
        throw SnapshotError("MethDB canonical content root is invalid");
    }
    char trailing = 0;
    if (input.get(trailing)) {
        throw SnapshotError("MethDB v2 has trailing bytes");
    }
    if (!input.eof()) {throw SnapshotError("failed while closing MethDB v2");}
    archive.file_sha256 = hash_file(path);
    return archive;
}

const SectionEntry &find_section(
    const ContigDirectory &contig,
    SectionType type)
{
    const auto found = std::find_if(
        contig.sections.begin(),
        contig.sections.end(),
        [type](const SectionEntry &entry) {return entry.type == type;});
    if (found == contig.sections.end()) {
        throw SnapshotError("MethDB requested section is absent");
    }
    return *found;
}

std::vector<std::uint8_t> load_section(
    const std::string &path,
    const SectionEntry &section)
{
    if (section.compressed_size > std::numeric_limits<std::size_t>::max()
        || section.raw_size > std::numeric_limits<std::size_t>::max()
        || section.compressed_size > std::numeric_limits<uLong>::max()
        || section.raw_size > std::numeric_limits<uLongf>::max()) {
        throw SnapshotError("MethDB section exceeds this process boundary");
    }
    std::ifstream input(path, std::ios::binary);
    if (!input) {throw SnapshotError("cannot reopen MethDB snapshot: " + path);}
    if (section.payload_offset > static_cast<std::uint64_t>(
                                     std::numeric_limits<std::streamoff>::max())) {
        throw SnapshotError("MethDB payload offset exceeds stream boundary");
    }
    input.seekg(static_cast<std::streamoff>(section.payload_offset));
    if (!input) {throw SnapshotError("cannot seek to MethDB section");}
    std::vector<std::uint8_t> compressed(
        static_cast<std::size_t>(section.compressed_size));
    read_raw(input, compressed.data(), compressed.size());
    std::vector<std::uint8_t> raw(
        static_cast<std::size_t>(section.raw_size));
    uLongf raw_size = static_cast<uLongf>(raw.size());
    const int status = uncompress(
        raw.data(),
        &raw_size,
        compressed.data(),
        static_cast<uLong>(compressed.size()));
    if (status != Z_OK || raw_size != raw.size()) {
        throw SnapshotError("MethDB section decompression failed");
    }
    if (crypto::sha256(raw) != section.raw_sha256) {
        throw SnapshotError("MethDB section canonical digest is invalid");
    }
    return raw;
}

std::uint32_t checked_u32(std::uint64_t value, const char *label)
{
    if (value > std::numeric_limits<std::uint32_t>::max()) {
        throw SnapshotError(std::string("MethDB ") + label + " exceeds uint32");
    }
    return static_cast<std::uint32_t>(value);
}

std::uint32_t read_delta(
    ByteReader &reader,
    std::uint32_t previous,
    bool first,
    const char *label)
{
    const std::uint64_t delta = reader.varuint();
    if (!first && delta == 0U) {
        throw SnapshotError(std::string("MethDB ") + label
                            + " delta is not positive");
    }
    const std::uint64_t value = first ? delta : previous + delta;
    if (!first && value < previous) {
        throw SnapshotError(std::string("MethDB ") + label
                            + " delta overflows");
    }
    return checked_u32(value, label);
}

std::vector<CatalogSite> decode_baseline(
    const std::vector<std::uint8_t> &raw,
    std::uint32_t count,
    std::uint32_t reference_length)
{
    ByteReader reader(raw);
    if (reader.u8() != 1U) {
        throw SnapshotError("MethDB baseline schema is unsupported");
    }
    std::vector<CatalogSite> sites;
    sites.reserve(count);
    std::uint32_t previous = 0U;
    for (std::uint32_t index = 0U; index < count; ++index) {
        const std::uint32_t position = read_delta(
            reader, previous, index == 0U, "baseline position");
        const RuntimeSite payload = read_payload(reader, position);
        if (position >= reference_length
            || runtime_site_allele(payload)
                != model::MethylationAllele::shared) {
            throw SnapshotError("MethDB baseline site is invalid");
        }
        sites.push_back(CatalogSite{
            position,
            runtime_site_probability(payload),
            runtime_site_context(payload),
            runtime_site_source(payload),
        });
        previous = position;
    }
    reader.finish();
    (void)MethylationCatalog(reference_length, sites);
    return sites;
}

model::Bases read_bases(ByteReader &reader, const char *label)
{
    const std::uint64_t size = reader.varuint();
    if (size > variant::maximum_indel_bases
        || size > std::numeric_limits<std::size_t>::max()) {
        throw SnapshotError(std::string("MethDB ") + label
                            + " length is invalid");
    }
    model::Bases result;
    result.reserve(static_cast<std::size_t>(size));
    for (std::uint64_t index = 0U; index < size; ++index) {
        const std::uint8_t base = reader.u8();
        if (base > 3U) {
            throw SnapshotError(std::string("MethDB ") + label
                                + " contains an invalid base");
        }
        result.push_back(base);
    }
    return result;
}

std::vector<variant::Variant> decode_events(
    const std::vector<std::uint8_t> &raw,
    std::uint32_t count,
    std::uint32_t contig_index,
    std::uint32_t reference_length)
{
    ByteReader reader(raw);
    if (reader.u8() != 1U) {
        throw SnapshotError("MethDB event schema is unsupported");
    }
    std::vector<variant::Variant> variants;
    variants.reserve(count);
    std::uint32_t previous_start = 0U;
    std::uint32_t previous_end = 0U;
    for (std::uint32_t index = 0U; index < count; ++index) {
        variant::Variant event;
        event.contig_index = contig_index;
        event.reference_start = read_delta(
            reader, previous_start, index == 0U, "event position");
        event.reference_end = checked_u32(
            static_cast<std::uint64_t>(event.reference_start)
                + reader.varuint(),
            "event end");
        event.kind = static_cast<model::VariantKind>(reader.u8());
        event.alt_haplotypes = static_cast<model::HaplotypeMask>(reader.u8());
        event.source = static_cast<model::VariantSource>(reader.u8());
        event.id = reader.string();
        event.ref_bases = read_bases(reader, "REF");
        event.alt_bases = read_bases(reader, "ALT");
        if (event.reference_end > reference_length
            || (index != 0U && event.reference_start < previous_end)
            || event.id.find_first_of("\t\r\n") != std::string::npos
            || event.id.find('\0') != std::string::npos
            || !model::is_haplotype_mask(
                static_cast<std::uint8_t>(event.alt_haplotypes))
            || (event.source != model::VariantSource::vcf
                && event.source != model::VariantSource::de_novo
                && event.source != model::VariantSource::asm_profile)) {
            throw SnapshotError("MethDB prepared event is invalid");
        }
        switch (event.kind) {
        case model::VariantKind::snv:
            if (event.reference_end != event.reference_start + 1U
                || event.ref_bases.size() != 1U
                || event.alt_bases.size() != 1U) {
                throw SnapshotError("MethDB SNV shape is invalid");
            }
            break;
        case model::VariantKind::insertion:
            if (event.reference_end != event.reference_start
                || !event.ref_bases.empty()
                || event.alt_bases.empty()) {
                throw SnapshotError("MethDB insertion shape is invalid");
            }
            break;
        case model::VariantKind::deletion:
            if (event.reference_end <= event.reference_start
                || event.ref_bases.empty()
                || !event.alt_bases.empty()
                || event.ref_bases.size()
                    != event.reference_end - event.reference_start) {
                throw SnapshotError("MethDB deletion shape is invalid");
            }
            break;
        default:
            throw SnapshotError("MethDB event kind is invalid");
        }
        previous_start = event.reference_start;
        previous_end = event.reference_end;
        variants.push_back(std::move(event));
    }
    reader.finish();
    return variants;
}

std::vector<OverlayRow> decode_overlays(
    const std::vector<std::uint8_t> &raw,
    std::uint32_t count,
    bool insertion,
    std::uint32_t reference_length,
    const std::vector<variant::Variant> &variants)
{
    ByteReader reader(raw);
    if (reader.u8() != 1U) {
        throw SnapshotError("MethDB overlay schema is unsupported");
    }
    std::vector<OverlayRow> rows;
    rows.reserve(count);
    std::uint32_t previous = 0U;
    for (std::uint32_t index = 0U; index < count; ++index) {
        OverlayRow row;
        row.key = read_delta(
            reader, previous, index == 0U, "overlay key");
        row.states = unpack_states(reader.u8());
        if (!insertion) {
            if (row.key >= reference_length) {
                throw SnapshotError("MethDB reference overlay is out of range");
            }
            const std::uint32_t cause_count = checked_u32(
                reader.varuint(), "cause count");
            if (cause_count == 0U || cause_count > variants.size()) {
                throw SnapshotError("MethDB cause count is invalid");
            }
            row.causes.reserve(cause_count);
            std::uint32_t previous_cause = 0U;
            for (std::uint32_t cause_index = 0U;
                 cause_index < cause_count;
                 ++cause_index) {
                const std::uint32_t cause = read_delta(
                    reader,
                    previous_cause,
                    cause_index == 0U,
                    "cause ordinal");
                if (cause >= variants.size()) {
                    throw SnapshotError("MethDB cause ordinal is invalid");
                }
                row.causes.push_back(cause);
                previous_cause = cause;
            }
            if (row.causes != overlay_causes(row.key, variants)) {
                throw SnapshotError(
                    "MethDB reference overlay causes are not canonical");
            }
        } else {
            const std::uint32_t event_ordinal = row.key >> 2U;
            const std::uint8_t insertion_offset =
                static_cast<std::uint8_t>(row.key & 3U);
            if (event_ordinal >= variants.size()
                || variants[event_ordinal].kind
                    != model::VariantKind::insertion
                || insertion_offset
                    >= variants[event_ordinal].alt_bases.size()
                || row.states[0] == OverlayState::inherit
                || row.states[1] == OverlayState::inherit
                || (row.states[0] == OverlayState::absent
                    && row.states[1] == OverlayState::absent)) {
                throw SnapshotError("MethDB insertion overlay is invalid");
            }
            row.causes.push_back(event_ordinal);
        }
        for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
            if (row.states[haplotype] == OverlayState::explicit_site) {
                row.payloads[haplotype] = read_payload(reader, row.key);
                if (runtime_site_reference_equivalent(
                        row.payloads[haplotype])) {
                    throw SnapshotError(
                        "MethDB explicit overlay cannot inherit baseline");
                }
            }
        }
        rows.push_back(std::move(row));
        previous = rows.back().key;
    }
    reader.finish();
    return rows;
}

std::vector<AsmRecord> decode_asm(
    const std::vector<std::uint8_t> &raw,
    std::uint32_t count,
    std::uint32_t reference_length)
{
    ByteReader reader(raw);
    if (reader.u8() != 1U) {
        throw SnapshotError("MethDB ASM schema is unsupported");
    }
    std::vector<AsmRecord> records;
    records.reserve(count);
    std::uint32_t previous = 0U;
    for (std::uint32_t index = 0U; index < count; ++index) {
        AsmRecord record;
        record.target_reference_position = read_delta(
            reader, previous, index == 0U, "ASM target");
        record.linked_variant_position = checked_u32(
            reader.varuint(), "ASM linked position");
        record.reference_probability_u16 = reader.u16();
        record.alternate_probability_u16 = reader.u16();
        record.context = static_cast<model::MethylationContext>(reader.u8());
        record.dinucleotide_second = reader.u8();
        record.linked_reference_base = reader.u8();
        record.linked_alternate_base = reader.u8();
        if (record.target_reference_position >= reference_length
            || record.linked_variant_position >= reference_length) {
            throw SnapshotError("MethDB ASM record is outside the contig");
        }
        records.push_back(record);
        previous = record.target_reference_position;
    }
    reader.finish();
    return records;
}

RuntimeSite inherited_site(
    const CatalogSite &baseline,
    model::MethylationAllele allele)
{
    return pack_runtime_site(
        baseline.reference_position,
        baseline.probability_u16,
        baseline.context,
        baseline.methylation_source,
        allele,
        true);
}

RuntimeSite with_allele(
    RuntimeSite site,
    model::MethylationAllele allele)
{
    return pack_runtime_site(
        runtime_site_key(site),
        runtime_site_probability(site),
        runtime_site_context(site),
        runtime_site_source(site),
        allele,
        runtime_site_reference_equivalent(site));
}

std::optional<RuntimeSite> resolve_reference_state(
    OverlayState state,
    RuntimeSite payload,
    const CatalogSite *baseline)
{
    switch (state) {
    case OverlayState::absent:
        return std::nullopt;
    case OverlayState::inherit:
        if (baseline == nullptr) {
            throw SnapshotError(
                "MethDB overlay inherits a nonexistent baseline site");
        }
        return inherited_site(
            *baseline, model::MethylationAllele::reference_haplotype);
    case OverlayState::explicit_site:
        return payload;
    }
    throw SnapshotError("MethDB reference overlay state is invalid");
}

void append_resolved_pair(
    std::array<std::optional<RuntimeSite>, 2> pair,
    std::vector<RuntimeSite> &shared,
    std::array<std::vector<RuntimeSite>, 2> &haplotypes)
{
    if (pair[0] && pair[1] && *pair[0] == *pair[1]) {
        shared.push_back(with_allele(
            *pair[0], model::MethylationAllele::shared));
        return;
    }
    for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
        if (pair[haplotype]) {
            haplotypes[haplotype].push_back(*pair[haplotype]);
        }
    }
}

void compile_reference_layer(
    const std::vector<CatalogSite> &baseline,
    const std::vector<OverlayRow> &overlays,
    DiploidRuntimeArrays &arrays)
{
    std::size_t baseline_index = 0U;
    std::size_t overlay_index = 0U;
    while (baseline_index < baseline.size()
           || overlay_index < overlays.size()) {
        const std::uint32_t baseline_key = baseline_index < baseline.size()
            ? baseline[baseline_index].reference_position
            : std::numeric_limits<std::uint32_t>::max();
        const std::uint32_t overlay_key = overlay_index < overlays.size()
            ? overlays[overlay_index].key
            : std::numeric_limits<std::uint32_t>::max();
        if (baseline_key < overlay_key) {
            const CatalogSite &site = baseline[baseline_index++];
            arrays.reference_shared.push_back(pack_runtime_site(
                site.reference_position,
                site.probability_u16,
                site.context,
                site.methylation_source,
                model::MethylationAllele::shared,
                true));
            continue;
        }
        const OverlayRow &overlay = overlays[overlay_index++];
        const CatalogSite *baseline_site = nullptr;
        if (baseline_key == overlay_key) {
            baseline_site = &baseline[baseline_index++];
        }
        if ((baseline_site != nullptr
             && overlay.states[0] == OverlayState::inherit
             && overlay.states[1] == OverlayState::inherit)
            || (baseline_site == nullptr
                && overlay.states[0] == OverlayState::absent
                && overlay.states[1] == OverlayState::absent)) {
            throw SnapshotError("MethDB reference overlay is redundant");
        }
        std::array<std::optional<RuntimeSite>, 2> pair;
        for (std::uint8_t haplotype = 0U; haplotype < 2U; ++haplotype) {
            pair[haplotype] = resolve_reference_state(
                overlay.states[haplotype],
                overlay.payloads[haplotype],
                baseline_site);
        }
        append_resolved_pair(
            std::move(pair),
            arrays.reference_shared,
            arrays.reference_haplotypes);
    }
}

void compile_insertion_layer(
    const std::vector<OverlayRow> &overlays,
    DiploidRuntimeArrays &arrays)
{
    for (const OverlayRow &overlay : overlays) {
        std::array<std::optional<RuntimeSite>, 2> pair;
        for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
            if (overlay.states[haplotype] == OverlayState::inherit) {
                throw SnapshotError("insertion overlay cannot inherit");
            }
            if (overlay.states[haplotype] == OverlayState::explicit_site) {
                pair[haplotype] = overlay.payloads[haplotype];
            }
        }
        append_resolved_pair(
            std::move(pair),
            arrays.insertion_shared,
            arrays.insertion_haplotypes);
    }
}

SnapshotContig load_contig(
    const std::string &path,
    const ArchiveDirectory &archive,
    std::uint32_t contig_index)
{
    if (contig_index >= archive.contigs.size()) {
        throw SnapshotError("MethDB contig index is out of range");
    }
    const ContigDirectory &directory = archive.contigs[contig_index];
    SnapshotContig result;
    result.name = directory.metadata.name;
    result.reference_length = static_cast<std::uint32_t>(
        directory.metadata.length);
    result.reference_sha256 = directory.metadata.reference_sha256;
    result.diploid = directory.diploid;

    const SectionEntry &baseline_entry = find_section(
        directory, SectionType::baseline);
    result.reference_sites = decode_baseline(
        load_section(path, baseline_entry),
        baseline_entry.record_count,
        result.reference_length);
    if (!result.diploid) {return result;}

    const SectionEntry &event_entry = find_section(
        directory, SectionType::events);
    result.variants = decode_events(
        load_section(path, event_entry),
        event_entry.record_count,
        contig_index,
        result.reference_length);
    const SectionEntry &reference_overlay_entry = find_section(
        directory, SectionType::reference_overlay);
    const std::vector<OverlayRow> reference_overlays = decode_overlays(
        load_section(path, reference_overlay_entry),
        reference_overlay_entry.record_count,
        false,
        result.reference_length,
        result.variants);
    const SectionEntry &insertion_overlay_entry = find_section(
        directory, SectionType::insertion_overlay);
    const std::vector<OverlayRow> insertion_overlays = decode_overlays(
        load_section(path, insertion_overlay_entry),
        insertion_overlay_entry.record_count,
        true,
        result.reference_length,
        result.variants);
    compile_reference_layer(
        result.reference_sites,
        reference_overlays,
        result.diploid_sites);
    compile_insertion_layer(
        insertion_overlays,
        result.diploid_sites);
    const SectionEntry &asm_entry = find_section(
        directory, SectionType::asm_layer);
    const std::vector<AsmRecord> asm_records = decode_asm(
        load_section(path, asm_entry),
        asm_entry.record_count,
        result.reference_length);
    try {
        DiploidMethylationCatalog catalog(
            contig_index,
            result.reference_length,
            std::move(result.diploid_sites));
        catalog.apply_asm_layer(result.variants, asm_records);
        result.diploid_sites = std::move(catalog).take_runtime_arrays();
    } catch (const std::exception &error) {
        throw SnapshotError(
            std::string("MethDB runtime compilation failed: ") + error.what());
    }
    return result;
}

} // namespace

class SnapshotWriterImpl {
public:
    SnapshotWriterImpl(
        std::ostream &output,
        const crypto::Sha256Digest &binding,
        std::uint32_t contig_count)
        : output_(output), contig_count_(contig_count)
    {
        write_raw(output_, methdb_magic, sizeof(methdb_magic) - 1U);
        write_u8(output_, methdb_version);
        write_u8(output_, little_endian_marker);
        write_u8(output_, 0U);
        write_u16(output_, 0U);
        write_raw(output_, binding.data(), binding.size());
        write_u32(output_, contig_count_);
        static constexpr char domain[] = "BSReadSim/MethDB/canonical-v2";
        content_hash_.update(
            reinterpret_cast<const std::uint8_t *>(domain),
            sizeof(domain) - 1U);
        content_hash_.update(binding.data(), binding.size());
        ByteWriter count;
        count.u32(contig_count_);
        content_hash_.update(count.bytes());
    }

    void begin_contig(
        const reference::ContigMetadata &metadata,
        bool diploid,
        std::uint8_t section_count)
    {
        if (finished_ || written_ >= contig_count_
            || metadata.name.empty()
            || metadata.name.size() > UINT32_C(1048576)
            || metadata.name.find_first_of("\t\r\n") != std::string::npos
            || metadata.name.find('\0') != std::string::npos
            || metadata.length > std::numeric_limits<std::uint32_t>::max()) {
            throw SnapshotError("MethDB writer received invalid contig metadata");
        }
        write_u32(output_, static_cast<std::uint32_t>(metadata.name.size()));
        write_raw(output_, metadata.name.data(), metadata.name.size());
        write_u32(output_, static_cast<std::uint32_t>(metadata.length));
        write_raw(
            output_, metadata.reference_sha256.data(),
            metadata.reference_sha256.size());
        write_u8(output_, diploid ? 1U : 0U);
        write_u8(output_, section_count);
        write_u16(output_, 0U);

        ByteWriter identity;
        identity.u32(written_);
        identity.string(metadata.name);
        identity.u32(static_cast<std::uint32_t>(metadata.length));
        identity.raw(
            metadata.reference_sha256.data(),
            metadata.reference_sha256.size());
        identity.u8(diploid ? 1U : 0U);
        content_hash_.update(identity.bytes());
    }

    void section(
        SectionType type,
        std::uint32_t record_count,
        std::vector<std::uint8_t> raw)
    {
        if (raw.empty()) {throw SnapshotError("MethDB canonical section is empty");}
        const crypto::Sha256Digest raw_sha256 = crypto::sha256(raw);
        const std::vector<std::uint8_t> compressed = compress_section(raw);
        write_u8(output_, static_cast<std::uint8_t>(type));
        write_u8(output_, 0U);
        write_u16(output_, 0U);
        write_u32(output_, written_);
        write_u32(output_, record_count);
        write_u64(output_, raw.size());
        write_u64(output_, compressed.size());
        write_raw(output_, raw_sha256.data(), raw_sha256.size());
        write_raw(output_, compressed.data(), compressed.size());

        ByteWriter identity;
        identity.u8(static_cast<std::uint8_t>(type));
        identity.u32(written_);
        identity.u32(record_count);
        identity.u64(raw.size());
        content_hash_.update(identity.bytes());
        content_hash_.update(raw_sha256.data(), raw_sha256.size());
    }

    void finish_contig() {++written_;}

    void finish()
    {
        if (finished_) {throw SnapshotError("MethDB writer was finished twice");}
        if (written_ != contig_count_) {
            throw SnapshotError("MethDB writer contig count is incomplete");
        }
        write_raw(output_, footer_magic.data(), footer_magic.size());
        const crypto::Sha256Digest digest = content_hash_.digest();
        write_raw(output_, digest.data(), digest.size());
        if (!output_) {throw SnapshotError("failed while finalizing MethDB v2");}
        finished_ = true;
    }

private:
    std::ostream &output_;
    crypto::Sha256 content_hash_;
    std::uint32_t contig_count_ = 0U;
    std::uint32_t written_ = 0U;
    bool finished_ = false;
};

SnapshotWriter::SnapshotWriter(
    std::ostream &output,
    const crypto::Sha256Digest &binding,
    std::uint32_t contig_count)
    : output_(output),
      impl_(std::make_unique<SnapshotWriterImpl>(
          output, binding, contig_count)),
      contig_count_(contig_count)
{}

SnapshotWriter::~SnapshotWriter() = default;

void SnapshotWriter::write_reference(
    const reference::ContigMetadata &metadata,
    const MethylationCatalog &catalog)
{
    if (finished_ || written_ >= contig_count_) {
        throw SnapshotError("MethDB writer received too many contigs");
    }
    impl_->begin_contig(metadata, false, 1U);
    impl_->section(
        SectionType::baseline,
        static_cast<std::uint32_t>(catalog.sites().size()),
        encode_baseline(catalog.sites()));
    impl_->finish_contig();
    ++written_;
}

void SnapshotWriter::write_diploid(
    const reference::ContigMetadata &metadata,
    const MethylationCatalog &baseline,
    const DiploidMethylationCatalog &pre_asm_catalog,
    const std::vector<variant::Variant> &variants,
    const std::vector<AsmRecord> &asm_records)
{
    if (finished_ || written_ >= contig_count_) {
        throw SnapshotError("MethDB writer received too many contigs");
    }
    if (variants.size() > std::numeric_limits<std::uint32_t>::max()
        || baseline.sites().size() > std::numeric_limits<std::uint32_t>::max()
        || asm_records.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw SnapshotError("MethDB canonical record count exceeds uint32");
    }
    const std::vector<OverlayRow> reference_overlays =
        build_reference_overlays(
            baseline.sites(), pre_asm_catalog, variants);
    const std::vector<OverlayRow> insertion_overlays =
        build_insertion_overlays(pre_asm_catalog, variants);
    if (reference_overlays.size() > std::numeric_limits<std::uint32_t>::max()
        || insertion_overlays.size()
            > std::numeric_limits<std::uint32_t>::max()) {
        throw SnapshotError("MethDB overlay count exceeds uint32");
    }
    if (!asm_records.empty()) {
        pre_asm_catalog.validate_asm_layer(variants, asm_records);
    }
    impl_->begin_contig(metadata, true, 5U);
    impl_->section(
        SectionType::baseline,
        static_cast<std::uint32_t>(baseline.sites().size()),
        encode_baseline(baseline.sites()));
    impl_->section(
        SectionType::events,
        static_cast<std::uint32_t>(variants.size()),
        encode_events(variants, written_));
    impl_->section(
        SectionType::reference_overlay,
        static_cast<std::uint32_t>(reference_overlays.size()),
        encode_overlays(reference_overlays, false));
    impl_->section(
        SectionType::insertion_overlay,
        static_cast<std::uint32_t>(insertion_overlays.size()),
        encode_overlays(insertion_overlays, true));
    impl_->section(
        SectionType::asm_layer,
        static_cast<std::uint32_t>(asm_records.size()),
        encode_asm(asm_records));
    impl_->finish_contig();
    ++written_;
}

void SnapshotWriter::finish()
{
    if (finished_) {throw SnapshotError("MethDB writer was finished twice");}
    if (written_ != contig_count_) {
        throw SnapshotError("MethDB writer contig count is incomplete");
    }
    impl_->finish();
    finished_ = true;
}

struct Snapshot::Impl {
    Impl(
        std::string archive_path,
        const crypto::Sha256Digest &binding,
        const std::vector<reference::ContigMetadata> &reference_catalog)
        : path(std::move(archive_path)),
          archive(parse_archive(path, &binding, &reference_catalog))
    {}

    std::string path;
    ArchiveDirectory archive;
};

Snapshot::Snapshot(
    const std::string &path,
    const crypto::Sha256Digest &expected_binding,
    const std::vector<reference::ContigMetadata> &reference_catalog)
    : impl_(std::make_shared<Impl>(
          path, expected_binding, reference_catalog))
{
    file_sha256_ = impl_->archive.file_sha256;
    content_sha256_ = impl_->archive.content_sha256;
}

SnapshotContig Snapshot::contig(std::uint32_t contig_index) const
{
    return load_contig(impl_->path, impl_->archive, contig_index);
}

std::vector<variant::Variant> Snapshot::variants(
    std::uint32_t contig_index) const
{
    if (contig_index >= impl_->archive.contigs.size()) {
        throw SnapshotError("MethDB contig index is out of range");
    }
    const ContigDirectory &directory = impl_->archive.contigs[contig_index];
    if (!directory.diploid) {return {};}
    const SectionEntry &event_entry = find_section(
        directory, SectionType::events);
    return decode_events(
        load_section(impl_->path, event_entry),
        event_entry.record_count,
        contig_index,
        static_cast<std::uint32_t>(directory.metadata.length));
}

bool Snapshot::contig_is_diploid(std::uint32_t contig_index) const
{
    if (contig_index >= impl_->archive.contigs.size()) {
        throw SnapshotError("MethDB contig index is out of range");
    }
    return impl_->archive.contigs[contig_index].diploid;
}

bool Snapshot::has_diploid_contigs() const noexcept
{
    return std::any_of(
        impl_->archive.contigs.begin(),
        impl_->archive.contigs.end(),
        [](const ContigDirectory &contig) {return contig.diploid;});
}

namespace {

struct MethbedFileCloser {
    void operator()(std::FILE *file) const noexcept
    {
        if (file != nullptr) {(void)std::fclose(file);}
    }
};

using MethbedFilePointer = std::unique_ptr<std::FILE, MethbedFileCloser>;

struct MethbedSpoolEntry {
    std::uint64_t offset = 0U;
    std::uint64_t size = 0U;
    bool diploid = false;
};

std::vector<std::string_view> methbed_fields(std::string_view line)
{
    std::vector<std::string_view> fields;
    std::size_t begin = 0U;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(
            begin,
            end == std::string_view::npos ? line.size() - begin : end - begin));
        if (end == std::string_view::npos) {break;}
        begin = end + 1U;
    }
    return fields;
}

std::uint64_t methbed_u64(std::string_view text, const char *label)
{
    if (text.empty()) {
        throw SnapshotError(std::string("MethBED ") + label + " is empty");
    }
    std::uint64_t value = 0U;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (parsed.ec != std::errc{} || parsed.ptr != text.data() + text.size()) {
        throw SnapshotError(
            std::string("MethBED ") + label + " is not an unsigned decimal");
    }
    return value;
}

std::uint32_t methbed_u32(std::string_view text, const char *label)
{
    const std::uint64_t value = methbed_u64(text, label);
    if (value > std::numeric_limits<std::uint32_t>::max()) {
        throw SnapshotError(std::string("MethBED ") + label + " exceeds uint32");
    }
    return static_cast<std::uint32_t>(value);
}

std::uint16_t methbed_u16(std::string_view text, const char *label)
{
    const std::uint64_t value = methbed_u64(text, label);
    if (value > std::numeric_limits<std::uint16_t>::max()) {
        throw SnapshotError(std::string("MethBED ") + label + " exceeds uint16");
    }
    return static_cast<std::uint16_t>(value);
}

std::uint8_t methbed_hex_digit(char value)
{
    if (value >= '0' && value <= '9') {
        return static_cast<std::uint8_t>(value - '0');
    }
    if (value >= 'a' && value <= 'f') {
        return static_cast<std::uint8_t>(value - 'a' + 10);
    }
    throw SnapshotError("MethBED SHA-256 must use lowercase hexadecimal");
}

crypto::Sha256Digest methbed_digest(std::string_view text, const char *label)
{
    crypto::Sha256Digest result = {};
    if (text.size() != result.size() * 2U) {
        throw SnapshotError(
            std::string("MethBED ") + label + " must contain 64 hex digits");
    }
    for (std::size_t index = 0U; index < result.size(); ++index) {
        result[index] = static_cast<std::uint8_t>(
            (methbed_hex_digit(text[index * 2U]) << 4U)
            | methbed_hex_digit(text[index * 2U + 1U]));
    }
    return result;
}

float methbed_probability(
    std::string_view text,
    std::uint16_t probability_u16)
{
    double value = 0.0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value);
    if (text.empty() || parsed.ec != std::errc{}
        || parsed.ptr != text.data() + text.size()
        || !std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw SnapshotError(
            "MethBED probability must be finite and in [0, 1]");
    }
    const float result = value == 0.0 ? 0.0F : static_cast<float>(value);
    if (probability_to_u16(result) != probability_u16) {
        throw SnapshotError(
            "MethBED probability disagrees with probability_u16");
    }
    return result;
}

model::MethylationContext methbed_context(std::string_view text)
{
    if (text == "CG-C") {return model::MethylationContext::cg_c;}
    if (text == "CHG-C") {return model::MethylationContext::chg_c;}
    if (text == "CHH-C") {return model::MethylationContext::chh_c;}
    if (text == "CG-G") {return model::MethylationContext::cg_g;}
    if (text == "CHG-G") {return model::MethylationContext::chg_g;}
    if (text == "CHH-G") {return model::MethylationContext::chh_g;}
    throw SnapshotError("MethDB BED context is invalid");
}

model::MethylationSource methbed_source(std::string_view text)
{
    if (text == "input" || text == "cgmap") {
        return model::MethylationSource::cgmap;
    }
    if (text == "asm") {return model::MethylationSource::asm_source;}
    if (text == "beta") {return model::MethylationSource::beta;}
    if (text == "pooled-input" || text == "pooled-cgmap") {
        return model::MethylationSource::pooled_cgmap;
    }
    throw SnapshotError("MethBED methylation source is invalid");
}

model::MethylationAllele methbed_allele(std::string_view text)
{
    if (text == "shared") {return model::MethylationAllele::shared;}
    if (text == "reference") {
        return model::MethylationAllele::reference_haplotype;
    }
    if (text == "alternate") {
        return model::MethylationAllele::alternate_haplotype;
    }
    throw SnapshotError("MethBED methylation allele is invalid");
}

model::Bases methbed_bases(std::string_view text, const char *label)
{
    if (text == ".") {return {};}
    if (text.empty() || text.size() > variant::maximum_indel_bases) {
        throw SnapshotError(std::string("MethBED ") + label + " length is invalid");
    }
    model::Bases result;
    result.reserve(text.size());
    for (const char value : text) {
        switch (value) {
        case 'A': result.push_back(0U); break;
        case 'C': result.push_back(1U); break;
        case 'G': result.push_back(2U); break;
        case 'T': result.push_back(3U); break;
        default:
            throw SnapshotError(
                std::string("MethBED ") + label + " must contain only A/C/G/T");
        }
    }
    return result;
}

model::VariantKind methbed_variant_kind(std::string_view text)
{
    if (text == "SNV") {return model::VariantKind::snv;}
    if (text == "insertion") {return model::VariantKind::insertion;}
    if (text == "deletion") {return model::VariantKind::deletion;}
    throw SnapshotError("MethBED variant kind is invalid");
}

model::VariantSource methbed_variant_source(std::string_view text)
{
    if (text == "vcf") {return model::VariantSource::vcf;}
    if (text == "de-novo") {return model::VariantSource::de_novo;}
    if (text == "asm") {return model::VariantSource::asm_profile;}
    throw SnapshotError("MethBED variant source is invalid");
}

void methbed_append_runtime_site(
    SnapshotContig &contig,
    std::string_view set,
    bool insertion,
    RuntimeSite site)
{
    if (!contig.diploid) {
        if (insertion || set != "reference"
            || runtime_site_allele(site)
                != model::MethylationAllele::shared) {
            throw SnapshotError(
                "reference-mode MethBED row has an invalid set or allele");
        }
        contig.reference_sites.push_back(CatalogSite{
            runtime_site_key(site),
            runtime_site_probability(site),
            runtime_site_context(site),
            runtime_site_source(site),
        });
        return;
    }
    std::vector<RuntimeSite> *destination = nullptr;
    if (set == "shared") {
        destination = insertion
            ? &contig.diploid_sites.insertion_shared
            : &contig.diploid_sites.reference_shared;
    } else if (set == "haplotype-1") {
        destination = insertion
            ? &contig.diploid_sites.insertion_haplotypes[0]
            : &contig.diploid_sites.reference_haplotypes[0];
    } else if (set == "haplotype-2") {
        destination = insertion
            ? &contig.diploid_sites.insertion_haplotypes[1]
            : &contig.diploid_sites.reference_haplotypes[1];
    } else {
        throw SnapshotError("diploid MethBED row has an invalid set");
    }
    destination->push_back(site);
}

void methbed_write_runtime_vector(
    ByteWriter &writer,
    const std::vector<RuntimeSite> &sites)
{
    if (sites.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw SnapshotError("MethBED per-set site count exceeds uint32");
    }
    writer.u32(static_cast<std::uint32_t>(sites.size()));
    for (const RuntimeSite site : sites) {writer.u64(site);}
}

std::vector<RuntimeSite> methbed_read_runtime_vector(ByteReader &reader)
{
    const std::uint32_t count = reader.u32();
    std::vector<RuntimeSite> result;
    result.reserve(count);
    for (std::uint32_t index = 0U; index < count; ++index) {
        result.push_back(reader.u64());
    }
    return result;
}

std::string methbed_spool_error(const char *operation)
{
    return std::string("cannot ") + operation + " MethBED spool: "
        + std::strerror(errno);
}

} // namespace

class MethbedSnapshot::Impl {
public:
    Impl(
        const std::string &path,
        const crypto::Sha256Digest &expected_binding,
        const std::vector<reference::ContigMetadata> &reference_catalog)
        : reference_catalog_(reference_catalog)
    {
        try {
            if (reference_catalog_.empty()) {
                throw SnapshotError(
                    "MethBED requires a non-empty reference catalog");
            }
            spool_.reset(std::tmpfile());
            if (!spool_) {throw SnapshotError(methbed_spool_error("create"));}
            text::TextSnapshot snapshot(path);
            file_sha256_ = snapshot.file_sha256();
            bool saw_format = false;
            bool saw_binding = false;
            bool saw_content = false;
            bool saw_columns = false;
            bool saw_footer = false;
            bool current_has_sites = false;
            std::optional<SnapshotContig> current;

            const auto flush = [&]() {
                if (!current) {return;}
                if (!current->diploid && !current->variants.empty()) {
                    throw SnapshotError(
                        "reference-mode MethBED contig contains variants");
                }
                append(std::move(*current));
                current.reset();
                current_has_sites = false;
            };

            snapshot.visit_lines(
                [&](std::string_view line, std::uint64_t line_number) {
                    if (line.empty()) {return;}
                    try {
                        if (saw_footer) {
                            throw SnapshotError(
                                "MethBED contains data after #file_sha256");
                        }
                        const auto fields = methbed_fields(line);
                        if (fields[0] == "#format") {
                            if (saw_format || fields.size() != 2U
                                || (fields[1] != methdb_bed_format
                                    && fields[1]
                                        != legacy_methbed_snapshot_format)) {
                                throw SnapshotError(
                                    "MethBED format marker is invalid");
                            }
                            saw_format = true;
                            return;
                        }
                        if (fields[0] == "#binding_sha256") {
                            if (!saw_format || saw_binding
                                || fields.size() != 2U
                                || methbed_digest(fields[1], "binding SHA-256")
                                    != expected_binding) {
                                throw SnapshotError(
                                    "MethBED reference binding is incompatible");
                            }
                            saw_binding = true;
                            return;
                        }
                        if (fields[0] == "#content_sha256") {
                            if (!saw_binding || saw_content
                                || fields.size() != 2U) {
                                throw SnapshotError(
                                    "MethBED content SHA-256 record is invalid");
                            }
                            source_content_sha256_ = methbed_digest(
                                fields[1], "content SHA-256");
                            saw_content = true;
                            return;
                        }
                        if (fields[0] == "#columns") {
                            static constexpr std::array<std::string_view, 17U>
                                expected = {{
                                    "#columns", "chrom", "chromStart", "chromEnd",
                                    "name", "score", "strand", "set", "origin_id",
                                    "origin_kind", "variant_event",
                                    "insertion_offset", "context", "source", "allele",
                                    "probability_u16", "probability",
                                }};
                            if (!saw_content || saw_columns
                                || fields.size() != expected.size()
                                || !std::equal(
                                    fields.begin(), fields.end(), expected.begin())) {
                                throw SnapshotError(
                                    "MethBED column declaration is invalid");
                            }
                            saw_columns = true;
                            return;
                        }
                        if (fields[0] == "#contig") {
                            if (!saw_columns || fields.size() != 6U) {
                                throw SnapshotError(
                                    "MethBED contig record is invalid");
                            }
                            flush();
                            const std::uint32_t index = methbed_u32(
                                fields[1], "contig index");
                            if (index != entries_.size()
                                || index >= reference_catalog_.size()) {
                                throw SnapshotError(
                                    "MethBED contigs must follow reference order");
                            }
                            const reference::ContigMetadata &expected_contig =
                                reference_catalog_[index];
                            const std::uint32_t length = methbed_u32(
                                fields[3], "contig length");
                            const crypto::Sha256Digest digest = methbed_digest(
                                fields[4], "contig SHA-256");
                            if (fields[2] != expected_contig.name
                                || length != expected_contig.length
                                || digest != expected_contig.reference_sha256
                                || (fields[5] != "diploid"
                                    && fields[5] != "reference")) {
                                throw SnapshotError(
                                    "MethBED contig identity disagrees with reference");
                            }
                            current.emplace();
                            current->name = expected_contig.name;
                            current->reference_length = length;
                            current->reference_sha256 = digest;
                            current->diploid = fields[5] == "diploid";
                            return;
                        }
                        if (fields[0] == "#variant") {
                            if (!current || !current->diploid
                                || current_has_sites || fields.size() != 11U
                                || fields[1] != current->name) {
                                throw SnapshotError(
                                    "MethBED variant record is misplaced or invalid");
                            }
                            const std::uint32_t ordinal = methbed_u32(
                                fields[2], "variant ordinal");
                            if (ordinal != current->variants.size()) {
                                throw SnapshotError(
                                    "MethBED variant ordinals must be consecutive");
                            }
                            const std::uint64_t mask_value = methbed_u64(
                                fields[8], "haplotype mask");
                            if (mask_value > 3U
                                || !model::is_haplotype_mask(
                                    static_cast<std::uint8_t>(mask_value))) {
                                throw SnapshotError(
                                    "MethBED haplotype mask is invalid");
                            }
                            variant::Variant event;
                            event.contig_index = static_cast<std::uint32_t>(
                                entries_.size());
                            event.reference_start = methbed_u32(
                                fields[3], "variant start");
                            event.reference_end = methbed_u32(
                                fields[4], "variant end");
                            event.kind = methbed_variant_kind(fields[5]);
                            event.ref_bases = methbed_bases(fields[6], "REF");
                            event.alt_bases = methbed_bases(fields[7], "ALT");
                            event.alt_haplotypes =
                                static_cast<model::HaplotypeMask>(mask_value);
                            event.id = fields[9] == "."
                                ? std::string{}
                                : std::string(fields[9]);
                            event.source = methbed_variant_source(fields[10]);
                            current->variants.push_back(std::move(event));
                            return;
                        }
                        if (fields[0] == "#insertion") {
                            if (!current || !current->diploid
                                || fields.size() != 11U
                                || fields[1] != current->name) {
                                throw SnapshotError(
                                    "MethBED insertion record is invalid");
                            }
                            current_has_sites = true;
                            const std::uint32_t event_ordinal = methbed_u32(
                                fields[4], "insertion event ordinal");
                            const std::uint32_t insertion_offset = methbed_u32(
                                fields[5], "insertion offset");
                            if (event_ordinal >= current->variants.size()
                                || current->variants[event_ordinal].kind
                                    != model::VariantKind::insertion
                                || insertion_offset
                                    >= current->variants[event_ordinal]
                                           .alt_bases.size()
                                || insertion_offset >= 4U
                                || event_ordinal
                                    > (std::numeric_limits<std::uint32_t>::max()
                                       >> 2U)) {
                                throw SnapshotError(
                                    "MethBED insertion origin is invalid");
                            }
                            const std::uint32_t key =
                                (event_ordinal << 2U) | insertion_offset;
                            const std::uint64_t origin_id = methbed_u64(
                                fields[3], "insertion origin ID");
                            if (origin_id != ((UINT64_C(1) << 63U) | key)) {
                                throw SnapshotError(
                                    "MethBED insertion origin ID is invalid");
                            }
                            const std::uint16_t probability = methbed_u16(
                                fields[9], "probability_u16");
                            (void)methbed_probability(fields[10], probability);
                            const RuntimeSite site = pack_runtime_site(
                                key,
                                probability,
                                methbed_context(fields[6]),
                                methbed_source(fields[7]),
                                methbed_allele(fields[8]),
                                false);
                            methbed_append_runtime_site(
                                *current, fields[2], true, site);
                            return;
                        }
                        if (fields[0] == "#file_sha256") {
                            if (!current || fields.size() != 2U) {
                                throw SnapshotError(
                                    "MethBED source file SHA-256 record is invalid");
                            }
                            flush();
                            source_file_sha256_ = methbed_digest(
                                fields[1], "source file SHA-256");
                            saw_footer = true;
                            return;
                        }
                        if (!fields[0].empty() && fields[0].front() == '#') {
                            throw SnapshotError("MethBED contains an unknown record");
                        }
                        if (!current || fields.size() != 16U
                            || fields[0] != current->name) {
                            throw SnapshotError("MethBED site row is invalid");
                        }
                        current_has_sites = true;
                        const std::uint32_t start = methbed_u32(
                            fields[1], "chromStart");
                        const std::uint32_t end = methbed_u32(
                            fields[2], "chromEnd");
                        if (start >= current->reference_length
                            || end != static_cast<std::uint64_t>(start) + 1U
                            || fields[7] != fields[1]
                            || fields[8] != "reference"
                            || fields[9] != "." || fields[10] != ".") {
                            throw SnapshotError(
                                "MethBED reference origin columns are invalid");
                        }
                        const std::string expected_name =
                            "methdb:" + std::string(fields[6]) + ":"
                            + std::to_string(start);
                        const std::string alternate_name =
                            "methbed:" + std::string(fields[6]) + ":"
                            + std::to_string(start);
                        if (fields[3] != expected_name
                            && fields[3] != alternate_name) {
                            throw SnapshotError("MethBED site name is invalid");
                        }
                        const std::uint16_t probability = methbed_u16(
                            fields[14], "probability_u16");
                        (void)methbed_probability(fields[15], probability);
                        const std::uint32_t expected_score =
                            (static_cast<std::uint32_t>(probability) * 1000U
                             + 32767U)
                            / 65535U;
                        if (methbed_u32(fields[4], "score") != expected_score) {
                            throw SnapshotError(
                                "MethBED score disagrees with probability_u16");
                        }
                        const model::MethylationContext context =
                            methbed_context(fields[11]);
                        const bool reverse =
                            static_cast<std::uint8_t>(context) >= 8U;
                        if (fields[5] != (reverse ? "-" : "+")) {
                            throw SnapshotError(
                                "MethBED strand disagrees with context");
                        }
                        const RuntimeSite site = pack_runtime_site(
                            start,
                            probability,
                            context,
                            methbed_source(fields[12]),
                            methbed_allele(fields[13]),
                            false);
                        methbed_append_runtime_site(
                            *current, fields[6], false, site);
                    } catch (const SnapshotError &error) {
                        throw SnapshotError(
                            "MethBED line " + std::to_string(line_number)
                            + ": " + error.what());
                    }
                });
            if (!saw_format || !saw_binding || !saw_content || !saw_columns
                || !saw_footer || current
                || entries_.size() != reference_catalog_.size()) {
                throw SnapshotError(
                    "MethBED header, footer, or contig set is incomplete");
            }
            if (std::fflush(spool_.get()) != 0) {
                throw SnapshotError(methbed_spool_error("flush"));
            }
        } catch (const SnapshotError &) {
            throw;
        } catch (const std::exception &error) {
            throw SnapshotError(error.what());
        }
    }

    SnapshotContig load(std::uint32_t contig_index) const
    {
        if (contig_index >= entries_.size()) {
            throw SnapshotError("MethBED contig index is out of range");
        }
        const MethbedSpoolEntry &entry = entries_[contig_index];
        if (entry.offset > static_cast<std::uint64_t>(
                               std::numeric_limits<off_t>::max())
            || entry.size > std::numeric_limits<std::size_t>::max()) {
            throw SnapshotError("MethBED spool boundary exceeds this process");
        }
        if (::fseeko(
                spool_.get(), static_cast<off_t>(entry.offset), SEEK_SET) != 0) {
            throw SnapshotError(methbed_spool_error("seek"));
        }
        std::vector<std::uint8_t> bytes(static_cast<std::size_t>(entry.size));
        if (!bytes.empty()
            && std::fread(bytes.data(), bytes.size(), 1U, spool_.get()) != 1U) {
            throw SnapshotError("MethBED spool is truncated or unreadable");
        }
        ByteReader reader(bytes);
        if (reader.u8() != 1U || (reader.u8() != 0U) != entry.diploid) {
            throw SnapshotError("MethBED spool header is corrupt");
        }
        SnapshotContig result;
        const reference::ContigMetadata &metadata =
            reference_catalog_[contig_index];
        result.name = metadata.name;
        result.reference_length = static_cast<std::uint32_t>(metadata.length);
        result.reference_sha256 = metadata.reference_sha256;
        result.diploid = entry.diploid;

        const std::uint32_t event_count = reader.u32();
        const std::uint64_t event_size = reader.u64();
        if (event_size > std::numeric_limits<std::size_t>::max()) {
            throw SnapshotError("MethBED event spool exceeds this process");
        }
        result.variants = decode_events(
            reader.bytes(static_cast<std::size_t>(event_size)),
            event_count,
            contig_index,
            result.reference_length);
        if (!result.diploid) {
            const std::uint32_t site_count = reader.u32();
            const std::uint64_t site_size = reader.u64();
            if (site_size > std::numeric_limits<std::size_t>::max()) {
                throw SnapshotError("MethBED site spool exceeds this process");
            }
            result.reference_sites = decode_baseline(
                reader.bytes(static_cast<std::size_t>(site_size)),
                site_count,
                result.reference_length);
        } else {
            result.diploid_sites.reference_shared =
                methbed_read_runtime_vector(reader);
            result.diploid_sites.insertion_shared =
                methbed_read_runtime_vector(reader);
            for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
                result.diploid_sites.reference_haplotypes[haplotype] =
                    methbed_read_runtime_vector(reader);
                result.diploid_sites.insertion_haplotypes[haplotype] =
                    methbed_read_runtime_vector(reader);
            }
            DiploidMethylationCatalog validated(
                contig_index,
                result.reference_length,
                std::move(result.diploid_sites));
            result.diploid_sites =
                std::move(validated).take_runtime_arrays();
        }
        reader.finish();
        return result;
    }

    const crypto::Sha256Digest &file_sha256() const noexcept
    {
        return file_sha256_;
    }

    bool contig_is_diploid(std::uint32_t contig_index) const
    {
        if (contig_index >= entries_.size()) {
            throw SnapshotError("MethBED contig index is out of range");
        }
        return entries_[contig_index].diploid;
    }

    bool has_diploid_contigs() const noexcept
    {
        return std::any_of(
            entries_.begin(),
            entries_.end(),
            [](const MethbedSpoolEntry &entry) {return entry.diploid;});
    }

private:
    void append(SnapshotContig contig)
    {
        const std::uint32_t contig_index = static_cast<std::uint32_t>(
            entries_.size());
        if (contig_index >= reference_catalog_.size()) {
            throw SnapshotError("MethBED has too many contigs");
        }
        ByteWriter writer;
        writer.u8(1U);
        writer.u8(contig.diploid ? 1U : 0U);
        const std::vector<std::uint8_t> events = encode_events(
            contig.variants, contig_index);
        writer.u32(static_cast<std::uint32_t>(contig.variants.size()));
        writer.u64(events.size());
        writer.raw(events.data(), events.size());
        if (!contig.diploid) {
            if (!contig.variants.empty()) {
                throw SnapshotError(
                    "reference-mode MethBED contig contains variants");
            }
            const MethylationCatalog validated(
                contig.reference_length,
                std::move(contig.reference_sites));
            const std::vector<std::uint8_t> sites = encode_baseline(
                validated.sites());
            writer.u32(static_cast<std::uint32_t>(validated.sites().size()));
            writer.u64(sites.size());
            writer.raw(sites.data(), sites.size());
        } else {
            if (!contig.reference_sites.empty()) {
                throw SnapshotError(
                    "diploid MethBED contig contains reference-mode sites");
            }
            DiploidMethylationCatalog validated(
                contig_index,
                contig.reference_length,
                std::move(contig.diploid_sites));
            const DiploidRuntimeArrays arrays =
                std::move(validated).take_runtime_arrays();
            methbed_write_runtime_vector(writer, arrays.reference_shared);
            methbed_write_runtime_vector(writer, arrays.insertion_shared);
            for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
                methbed_write_runtime_vector(
                    writer, arrays.reference_haplotypes[haplotype]);
                methbed_write_runtime_vector(
                    writer, arrays.insertion_haplotypes[haplotype]);
            }
        }
        const std::vector<std::uint8_t> bytes = writer.take();
        const off_t position = ::ftello(spool_.get());
        if (position < 0) {throw SnapshotError(methbed_spool_error("locate"));}
        if (!bytes.empty()
            && std::fwrite(bytes.data(), bytes.size(), 1U, spool_.get()) != 1U) {
            throw SnapshotError(methbed_spool_error("write"));
        }
        entries_.push_back(MethbedSpoolEntry{
            static_cast<std::uint64_t>(position),
            static_cast<std::uint64_t>(bytes.size()),
            contig.diploid,
        });
    }

    std::vector<reference::ContigMetadata> reference_catalog_;
    std::vector<MethbedSpoolEntry> entries_;
    MethbedFilePointer spool_;
    crypto::Sha256Digest file_sha256_ = {};
    crypto::Sha256Digest source_content_sha256_ = {};
    crypto::Sha256Digest source_file_sha256_ = {};
};

MethbedSnapshot::MethbedSnapshot(
    const std::string &path,
    const crypto::Sha256Digest &expected_binding,
    const std::vector<reference::ContigMetadata> &reference_catalog)
    : impl_(std::make_shared<Impl>(
          path, expected_binding, reference_catalog))
{}

MethbedSnapshot::~MethbedSnapshot() = default;

SnapshotContig MethbedSnapshot::contig(
    const reference::Contig &contig) const
{
    SnapshotContig result = impl_->load(contig.index);
    normalize_methbed_contig(contig, result);
    return result;
}

std::vector<variant::Variant> MethbedSnapshot::variants(
    std::uint32_t contig_index) const
{
    return impl_->load(contig_index).variants;
}

bool MethbedSnapshot::contig_is_diploid(
    std::uint32_t contig_index) const
{
    return impl_->contig_is_diploid(contig_index);
}

bool MethbedSnapshot::has_diploid_contigs() const noexcept
{
    return impl_->has_diploid_contigs();
}

const crypto::Sha256Digest &MethbedSnapshot::file_sha256() const noexcept
{
    return impl_->file_sha256();
}

namespace {

std::string_view context_name(model::MethylationContext context)
{
    switch (context) {
    case model::MethylationContext::cg_c: return "CG-C";
    case model::MethylationContext::chg_c: return "CHG-C";
    case model::MethylationContext::chh_c: return "CHH-C";
    case model::MethylationContext::cg_g: return "CG-G";
    case model::MethylationContext::chg_g: return "CHG-G";
    case model::MethylationContext::chh_g: return "CHH-G";
    }
    throw SnapshotError("MethBED context is invalid");
}

std::string_view source_name(model::MethylationSource source)
{
    switch (source) {
    case model::MethylationSource::cgmap: return "input";
    case model::MethylationSource::asm_source: return "asm";
    case model::MethylationSource::beta: return "beta";
    case model::MethylationSource::pooled_cgmap: return "pooled-input";
    }
    throw SnapshotError("MethDB BED source is invalid");
}

std::string_view allele_name(model::MethylationAllele allele)
{
    switch (allele) {
    case model::MethylationAllele::shared: return "shared";
    case model::MethylationAllele::reference_haplotype: return "reference";
    case model::MethylationAllele::alternate_haplotype: return "alternate";
    }
    throw SnapshotError("MethDB BED allele is invalid");
}

std::string bases_text(const model::Bases &bases)
{
    if (bases.empty()) {return ".";}
    static constexpr char alphabet[] = "ACGT";
    std::string result;
    result.reserve(bases.size());
    for (const std::uint8_t base : bases) {
        if (base > 3U) {
            throw SnapshotError("MethDB event contains an invalid base");
        }
        result.push_back(alphabet[base]);
    }
    return result;
}

std::string_view variant_kind_name(model::VariantKind kind)
{
    switch (kind) {
    case model::VariantKind::snv: return "SNV";
    case model::VariantKind::insertion: return "insertion";
    case model::VariantKind::deletion: return "deletion";
    }
    throw SnapshotError("MethDB event kind is invalid");
}

std::string_view variant_source_name(model::VariantSource source)
{
    switch (source) {
    case model::VariantSource::vcf: return "vcf";
    case model::VariantSource::de_novo: return "de-novo";
    case model::VariantSource::asm_profile: return "asm";
    }
    throw SnapshotError("MethDB event source is invalid");
}

void write_hex(std::ostream &sink, const crypto::Sha256Digest &digest)
{
    static constexpr char digits[] = "0123456789abcdef";
    for (const std::uint8_t byte : digest) {
        sink.put(digits[byte >> 4U]);
        sink.put(digits[byte & UINT8_C(0x0f)]);
    }
}

void write_bed_site(
    std::ostream &sink,
    std::string_view contig,
    std::string_view set_name,
    RuntimeSite site,
    bool insertion)
{
    constexpr std::uint64_t insertion_flag = UINT64_C(1) << 63U;
    const std::uint32_t key = runtime_site_key(site);
    const std::uint64_t origin_id = insertion
        ? insertion_flag | key
        : key;
    const std::uint16_t probability = runtime_site_probability(site);
    if (insertion) {
        sink << "#insertion\t" << contig << '\t' << set_name << '\t'
             << origin_id << '\t' << (key >> 2U) << '\t' << (key & 3U);
    } else {
        const std::uint32_t score =
            (static_cast<std::uint32_t>(probability) * 1000U + 32767U)
            / 65535U;
        sink << contig << '\t' << key << '\t'
             << (key + 1U) << "\tmethdb:" << set_name << ':'
             << key << '\t' << score << '\t'
             << (static_cast<std::uint8_t>(runtime_site_context(site)) < 8U
                    ? '+' : '-')
             << '\t' << set_name << '\t' << key
             << "\t" << (insertion ? "insertion" : "reference") << "\t.\t.";
    }
    sink << '\t' << context_name(runtime_site_context(site))
         << '\t' << source_name(runtime_site_source(site))
         << '\t' << allele_name(runtime_site_allele(site))
         << '\t' << probability
         << '\t' << std::setprecision(std::numeric_limits<float>::max_digits10)
         << probability_from_u16(probability) << '\n';
}

} // namespace

void export_snapshot_bed(const std::string &path, std::ostream &sink)
{
    const ArchiveDirectory archive = parse_archive(path, nullptr, nullptr);
    sink << "#format\t" << methdb_bed_format << '\n'
         << "#binding_sha256\t";
    write_hex(sink, archive.binding);
    sink << "\n#content_sha256\t";
    write_hex(sink, archive.content_sha256);
    sink << "\n#columns\tchrom\tchromStart\tchromEnd\tname\tscore\tstrand\t"
            "set\torigin_id\torigin_kind\tvariant_event\tinsertion_offset\t"
            "context\tsource\tallele\tprobability_u16\tprobability\n";
    for (std::uint32_t index = 0U; index < archive.contigs.size(); ++index) {
        SnapshotContig contig = load_contig(path, archive, index);
        sink << "#contig\t" << index << '\t' << contig.name << '\t'
             << contig.reference_length << '\t';
        write_hex(sink, contig.reference_sha256);
        sink << '\t' << (contig.diploid ? "diploid" : "reference") << '\n';
        for (std::size_t ordinal = 0U;
             ordinal < contig.variants.size();
             ++ordinal) {
            const variant::Variant &event = contig.variants[ordinal];
            sink << "#variant\t" << contig.name << '\t' << ordinal << '\t'
                 << event.reference_start << '\t' << event.reference_end
                 << '\t' << variant_kind_name(event.kind) << '\t'
                 << bases_text(event.ref_bases) << '\t'
                 << bases_text(event.alt_bases) << '\t'
                 << static_cast<unsigned>(event.alt_haplotypes) << '\t'
                 << (event.id.empty() ? "." : event.id) << '\t'
                 << variant_source_name(event.source) << '\n';
        }
        if (!contig.diploid) {
            for (const CatalogSite &site : contig.reference_sites) {
                write_bed_site(
                    sink,
                    contig.name,
                    "reference",
                    pack_runtime_site(
                        site.reference_position,
                        site.probability_u16,
                        site.context,
                        site.methylation_source,
                        model::MethylationAllele::shared,
                        true),
                    false);
            }
        } else {
            for (const RuntimeSite site
                 : contig.diploid_sites.reference_shared) {
                write_bed_site(
                    sink, contig.name, "shared", site, false);
            }
            for (const RuntimeSite site
                 : contig.diploid_sites.insertion_shared) {
                write_bed_site(
                    sink, contig.name, "shared", site, true);
            }
            for (std::size_t haplotype = 0U; haplotype < 2U; ++haplotype) {
                const std::string set_name = "haplotype-"
                    + std::to_string(haplotype + 1U);
                for (const RuntimeSite site
                     : contig.diploid_sites
                           .reference_haplotypes[haplotype]) {
                    write_bed_site(
                        sink, contig.name, set_name, site, false);
                }
                for (const RuntimeSite site
                     : contig.diploid_sites
                           .insertion_haplotypes[haplotype]) {
                    write_bed_site(
                        sink, contig.name, set_name, site, true);
                }
            }
        }
    }
    sink << "#file_sha256\t";
    write_hex(sink, archive.file_sha256);
    sink << '\n';
    if (!sink) {throw SnapshotError("failed while exporting MethDB BED");}
}

} // namespace htsim::methdb
