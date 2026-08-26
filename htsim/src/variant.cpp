#include "variant.h"

#include <algorithm>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <memory>
#include <mutex>
#include <ostream>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <cmath>

#include <htslib/hts.h>
#include <htslib/hts_log.h>
#include <htslib/vcf.h>

#include "utilities.h"
#include "types.h"

// ---- variant_catalog --------------------------------------------------------

namespace htsim::variant {
namespace {

std::string hexadecimal(std::uint64_t value)
{
    constexpr char digits[] = "0123456789abcdef";
    std::string result;
    do {
        result.push_back(digits[value & UINT64_C(0xf)]);
        value >>= 4U;
    } while (value != 0U);
    std::reverse(result.begin(), result.end());
    return result;
}

std::uint64_t catalog_address(
    std::uint32_t contig_index,
    std::size_t local_index)
{
    if (local_index > std::numeric_limits<std::uint32_t>::max()) {
        throw VariantCatalogError("variant index exceeds uint32");
    }
    return (static_cast<std::uint64_t>(contig_index) << 32U)
        | static_cast<std::uint32_t>(local_index);
}

std::string unique_variant_id(
    std::string requested,
    std::string_view fallback_prefix,
    std::uint64_t address,
    std::unordered_set<std::string> &used)
{
    if (requested.empty() || requested == ".") {
        requested = std::string(fallback_prefix) + hexadecimal(address);
    }
    if (used.insert(requested).second) {return requested;}

    const std::string disambiguated = requested + "@" + hexadecimal(address);
    if (used.insert(disambiguated).second) {return disambiguated;}
    throw VariantCatalogError("variant ID cannot be made unique");
}

model::Bases parse_allele(std::string_view text, const char *field)
{
    if (text.empty()) {
        throw VariantCatalogError(std::string("VCF ") + field + " is empty");
    }
    if (text.find(',') != std::string_view::npos) {
        throw VariantCatalogError("VCF subset requires one ALT allele");
    }
    model::Bases bases;
    bases.reserve(text.size());
    for (const char base : text) {
        switch (base) {
        case 'A': bases.push_back(0U); break;
        case 'C': bases.push_back(1U); break;
        case 'G': bases.push_back(2U); break;
        case 'T': bases.push_back(3U); break;
        case 'N': bases.push_back(4U); break;
        default:
            throw VariantCatalogError(
                std::string("VCF ") + field
                + " must contain only uppercase A/C/G/T/N");
        }
    }
    return bases;
}

struct Genotype {
    std::uint8_t first = 0;
    std::uint8_t second = 0;
    bool phased = false;
};

struct GenotypeBuffer {
    ~GenotypeBuffer() {std::free(values);}

    GenotypeBuffer(const GenotypeBuffer &) = delete;
    GenotypeBuffer &operator=(const GenotypeBuffer &) = delete;
    GenotypeBuffer() = default;

    int32_t *values = nullptr;
    int capacity = 0;
};

class HtsLogGuard {
public:
    HtsLogGuard() : lock_(mutex()), previous_(hts_get_log_level())
    {
        // Undefined contig/unused FORMAT declarations are intentionally
        // tolerated by the supported subset. Keep HTSlib errors, but do not turn
        // those expected warnings into htsim-core stderr on a valid run.
        hts_set_log_level(HTS_LOG_ERROR);
    }

    ~HtsLogGuard() {hts_set_log_level(previous_);}

    HtsLogGuard(const HtsLogGuard &) = delete;
    HtsLogGuard &operator=(const HtsLogGuard &) = delete;

private:
    static std::mutex &mutex()
    {
        static std::mutex value;
        return value;
    }

    std::unique_lock<std::mutex> lock_;
    htsLogLevel previous_;
};

void validate_format_shape(std::string_view line)
{
    std::size_t format_begin = 0;
    for (std::size_t field = 0; field < 8U; ++field) {
        const std::size_t separator = line.find('\t', format_begin);
        if (separator == std::string_view::npos) {
            throw VariantCatalogError(
                "VCF data row must contain exactly ten fields");
        }
        format_begin = separator + 1U;
    }
    const std::size_t format_end = line.find('\t', format_begin);
    if (format_end == std::string_view::npos) {
        throw VariantCatalogError(
            "VCF data row must contain exactly ten fields");
    }
    const std::string_view format =
        line.substr(format_begin, format_end - format_begin);
    const std::string_view sample = line.substr(format_end + 1U);
    if (std::count(format.begin(), format.end(), ':')
        != std::count(sample.begin(), sample.end(), ':')) {
        throw VariantCatalogError(
            "VCF FORMAT and sample field counts must match exactly");
    }

    std::size_t gt_fields = 0;
    std::size_t token_begin = 0;
    while (true) {
        const std::size_t token_end = format.find(':', token_begin);
        const std::string_view token = format.substr(
            token_begin,
            token_end == std::string_view::npos
                ? format.size() - token_begin
                : token_end - token_begin);
        if (token == "GT") {++gt_fields;}
        if (token_end == std::string_view::npos) {break;}
        token_begin = token_end + 1U;
    }
    if (gt_fields == 0U) {
        throw VariantCatalogError("VCF FORMAT does not contain GT");
    }
    if (gt_fields != 1U) {
        throw VariantCatalogError("VCF FORMAT contains duplicate GT fields");
    }
}

Genotype parse_genotype(
    const bcf_hdr_t *header,
    bcf1_t *record,
    GenotypeBuffer &buffer)
{
    const int gt_id = bcf_hdr_id2int(header, BCF_DT_ID, "GT");
    std::size_t gt_fields = 0;
    for (std::size_t index = 0; index < record->n_fmt; ++index) {
        if (record->d.fmt[index].id == gt_id) {++gt_fields;}
    }
    if (gt_fields == 0U) {
        throw VariantCatalogError("VCF FORMAT does not contain GT");
    }
    if (gt_fields != 1U) {
        throw VariantCatalogError("VCF FORMAT contains duplicate GT fields");
    }
    const int count = bcf_get_genotypes(
        header, record, &buffer.values, &buffer.capacity);
    if (count != 2 || buffer.values == nullptr
        || bcf_gt_is_missing(buffer.values[0])
        || bcf_gt_is_missing(buffer.values[1])
        || buffer.values[0] == bcf_int32_vector_end
        || buffer.values[1] == bcf_int32_vector_end) {
        throw VariantCatalogError(
            "VCF GT must be diploid 0/0, 0/1, 1/0, or 1/1");
    }
    const int first = bcf_gt_allele(buffer.values[0]);
    const int second = bcf_gt_allele(buffer.values[1]);
    if ((first != 0 && first != 1) || (second != 0 && second != 1)) {
        throw VariantCatalogError(
            "VCF GT must be diploid 0/0, 0/1, 1/0, or 1/1");
    }
    return {
        static_cast<std::uint8_t>(first),
        static_cast<std::uint8_t>(second),
        bcf_gt_is_phased(buffer.values[1]) != 0,
    };
}

void validate_vcf_surface(text::TextSnapshot &snapshot)
{
    constexpr std::string_view header_prefix =
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t";
    bool saw_fileformat = false;
    bool saw_header = false;
    snapshot.visit_lines(
        [&](std::string_view line, std::uint64_t line_number) {
            try {
                if (!saw_header && line.substr(0, 2U) == "##") {
                    if (line.substr(0, 17U) == "##fileformat=VCFv") {
                        if (saw_fileformat
                            || (line != "##fileformat=VCFv4.2"
                                && line != "##fileformat=VCFv4.3")) {
                            throw VariantCatalogError(
                                "VCF requires one supported fileformat "
                                "declaration");
                        }
                        saw_fileformat = true;
                    }
                    return;
                }
                if (!saw_header && line.substr(0, header_prefix.size())
                        == header_prefix
                    && line.size() > header_prefix.size()
                    && line.find('\t', header_prefix.size())
                        == std::string_view::npos) {
                    if (!saw_fileformat) {
                        throw VariantCatalogError(
                            "VCF header must follow a supported fileformat "
                            "declaration");
                    }
                    saw_header = true;
                    return;
                }
                if (!saw_header) {
                    throw VariantCatalogError(
                        "VCF data appeared before its header");
                }
                if (line.empty() || line.front() == '#') {
                    throw VariantCatalogError(
                        "VCF data section contains an empty or header line");
                }
                if (static_cast<std::size_t>(
                        std::count(line.begin(), line.end(), '\t'))
                    != 9U) {
                    throw VariantCatalogError(
                        "VCF data row must contain exactly ten fields");
                }
                validate_format_shape(line);
            } catch (const VariantCatalogError &error) {
                throw VariantCatalogError(
                    "VCF line " + std::to_string(line_number) + ": "
                    + error.what());
            }
        });
    if (!saw_header) {throw VariantCatalogError("VCF header is missing");}
}

model::HaplotypeMask resolve_alt_haplotypes(
    const Genotype &genotype,
    std::uint64_t key,
    std::uint64_t event_ordinal)
{
    const std::uint8_t alt_count = static_cast<std::uint8_t>(
        genotype.first + genotype.second);
    if (alt_count == 0U) {
        throw VariantCatalogError("reference genotype has no alternate haplotype");
    }
    if (alt_count == 2U) {return model::HaplotypeMask::both;}
    if (genotype.phased) {
        return genotype.first == 1U
            ? model::HaplotypeMask::haplotype_1
            : model::HaplotypeMask::haplotype_2;
    }
    // Addressed random phasing places ALT on haplotype 1 for true and on
    // haplotype 2 for false.
    return rng::bernoulli(key, event_ordinal, UINT64_C(1), 0.5)
        ? model::HaplotypeMask::haplotype_1
        : model::HaplotypeMask::haplotype_2;
}

Variant normalize_event(
    std::uint32_t contig_index,
    std::uint32_t one_based_position,
    model::Bases reference,
    model::Bases alternate,
    model::HaplotypeMask alt_haplotypes,
    std::uint64_t contig_length)
{
    if (reference == alternate) {
        throw VariantCatalogError("VCF REF and ALT must differ");
    }
    const std::uint64_t original_start = one_based_position - 1U;
    const std::uint64_t original_end = original_start + reference.size();
    if (original_end > contig_length) {
        throw VariantCatalogError("VCF REF interval exceeds its reference contig");
    }

    std::size_t prefix = 0;
    while (prefix < reference.size() && prefix < alternate.size()
           && reference[prefix] == alternate[prefix]) {
        ++prefix;
    }
    std::size_t suffix = 0;
    while (suffix < reference.size() - prefix
           && suffix < alternate.size() - prefix
           && reference[reference.size() - 1U - suffix]
               == alternate[alternate.size() - 1U - suffix]) {
        ++suffix;
    }
    model::Bases normalized_ref(
        reference.begin() + static_cast<std::ptrdiff_t>(prefix),
        reference.end() - static_cast<std::ptrdiff_t>(suffix));
    model::Bases normalized_alt(
        alternate.begin() + static_cast<std::ptrdiff_t>(prefix),
        alternate.end() - static_cast<std::ptrdiff_t>(suffix));
    if (normalized_ref.empty() && normalized_alt.empty()) {
        throw VariantCatalogError("VCF normalization produced no variant");
    }
    const auto resolved = [](const model::Bases &bases) {
        return std::all_of(
            bases.begin(), bases.end(),
            [](std::uint8_t base) {return base < 4U;});
    };
    if (!resolved(normalized_ref) || !resolved(normalized_alt)) {
        throw VariantCatalogError(
            "VCF N is allowed only in a common indel anchor");
    }

    model::VariantKind kind = model::VariantKind::snv;
    if (normalized_ref.empty()) {
        if (normalized_alt.size() > maximum_indel_bases) {
            throw VariantCatalogError("VCF insertion exceeds four bases");
        }
        kind = model::VariantKind::insertion;
    } else if (normalized_alt.empty()) {
        if (normalized_ref.size() > maximum_indel_bases) {
            throw VariantCatalogError("VCF deletion exceeds four bases");
        }
        kind = model::VariantKind::deletion;
    } else if (normalized_ref.size() != 1U || normalized_alt.size() != 1U) {
        throw VariantCatalogError(
            "VCF subset rejects MNP and complex replacement variants");
    }

    const std::uint64_t normalized_start = original_start + prefix;
    const std::uint64_t normalized_end =
        normalized_start + normalized_ref.size();
    if (normalized_start > std::numeric_limits<std::uint32_t>::max()
        || normalized_end > std::numeric_limits<std::uint32_t>::max()) {
        throw VariantCatalogError("VCF normalized coordinates exceed uint32");
    }
    return {
        contig_index,
        static_cast<std::uint32_t>(normalized_start),
        static_cast<std::uint32_t>(normalized_end),
        kind,
        std::move(normalized_ref),
        std::move(normalized_alt),
        alt_haplotypes,
    };
}

bool conflicts(const Variant &left, const Variant &right) noexcept
{
    const bool left_insertion = left.reference_start == left.reference_end;
    const bool right_insertion = right.reference_start == right.reference_end;
    if (left_insertion && right_insertion) {
        return left.reference_start == right.reference_start;
    }
    if (left_insertion) {
        return left.reference_start >= right.reference_start
            && left.reference_start < right.reference_end;
    }
    if (right_insertion) {
        return right.reference_start >= left.reference_start
            && right.reference_start < left.reference_end;
    }
    return left.reference_start < right.reference_end
        && right.reference_start < left.reference_end;
}

bool event_less(const Variant &left, const Variant &right) noexcept
{
    return std::tie(left.reference_start, left.reference_end, left.kind)
        < std::tie(right.reference_start, right.reference_end, right.kind);
}

} // namespace

VariantFile::VariantFile(
    const std::string &path,
    const std::vector<reference::ContigMetadata> &reference_catalog,
    std::uint64_t master_seed)
    : variants_by_contig_(reference_catalog.size())
{
    try {
        std::unordered_map<std::string, std::uint32_t> contig_indices;
        contig_indices.reserve(reference_catalog.size());
        for (std::size_t index = 0; index < reference_catalog.size(); ++index) {
            if (index > std::numeric_limits<std::uint32_t>::max()) {
                throw VariantCatalogError("reference contig index exceeds uint32");
            }
            if (reference_catalog[index].length
                > std::numeric_limits<std::uint32_t>::max()) {
                throw VariantCatalogError("VCF reference contig exceeds uint32");
            }
            if (!contig_indices.emplace(
                    reference_catalog[index].name,
                    static_cast<std::uint32_t>(index)).second) {
                throw VariantCatalogError("reference catalog contains duplicate names");
            }
        }

        text::TextSnapshot snapshot(path);
        file_sha256_ = snapshot.file_sha256();
        validate_vcf_surface(snapshot);

        std::uint32_t previous_contig = 0;
        std::uint32_t previous_position = 0;
        bool have_previous_record = false;
        std::unordered_set<std::string> used_variant_ids;
        snapshot.visit_hts([&](htsFile *file) {
            HtsLogGuard log_guard;
            const htsFormat *format = hts_get_format(file);
            if (format == nullptr || format->format != ::vcf) {
                throw VariantCatalogError(
                    "variant input must be textual VCF, not BCF or another "
                    "HTS format");
            }

            using HeaderPtr = std::unique_ptr<
                bcf_hdr_t, decltype(&bcf_hdr_destroy)>;
            HeaderPtr header(bcf_hdr_read(file), &bcf_hdr_destroy);
            if (!header) {
                throw VariantCatalogError("HTSlib could not parse the VCF header");
            }
            const char *version = bcf_hdr_get_version(header.get());
            if (version == nullptr
                || (std::string_view(version) != "VCFv4.2"
                    && std::string_view(version) != "VCFv4.3")) {
                throw VariantCatalogError(
                    "VCF requires version 4.2 or 4.3");
            }
            if (bcf_hdr_nsamples(header.get()) != 1) {
                throw VariantCatalogError(
                    "VCF header must declare exactly one named sample");
            }
            // The BSReadSim subset historically permits compact VCFs without
            // redundant ##FORMAT metadata.  Give HTSlib the standard GT type
            // before record decoding while leaving the source bytes untouched.
            if (bcf_hdr_id2int(header.get(), BCF_DT_ID, "GT") < 0) {
                if (bcf_hdr_append(
                        header.get(),
                        "##FORMAT=<ID=GT,Number=1,Type=String,"
                        "Description=\"Genotype\">") != 0
                    || bcf_hdr_sync(header.get()) != 0) {
                    throw VariantCatalogError(
                        "HTSlib could not install the standard GT definition");
                }
            }

            using RecordPtr = std::unique_ptr<bcf1_t, decltype(&bcf_destroy)>;
            RecordPtr record(bcf_init(), &bcf_destroy);
            if (!record) {
                throw VariantCatalogError(
                    "HTSlib could not allocate a VCF record");
            }
            GenotypeBuffer genotype_buffer;
            while (true) {
                const int status = bcf_read(file, header.get(), record.get());
                if (status == -1) {break;}
                if (status < -1) {
                    throw VariantCatalogError(
                        "VCF line " + std::to_string(file->lineno)
                        + ": HTSlib could not parse the record");
                }
                try {
                    constexpr int tolerated_missing_declarations =
                        BCF_ERR_CTG_UNDEF | BCF_ERR_TAG_UNDEF;
                    if ((record->errcode & ~tolerated_missing_declarations)
                        != 0) {
                        throw VariantCatalogError(
                            "HTSlib reported a malformed VCF record");
                    }
                    if (bcf_unpack(record.get(), BCF_UN_ALL) != 0) {
                        throw VariantCatalogError(
                            "HTSlib could not unpack the VCF record");
                    }
                    if (record->n_sample != 1U) {
                        throw VariantCatalogError(
                            "VCF record must contain exactly one sample");
                    }
                    if (record->n_allele != 2U) {
                        throw VariantCatalogError(
                            "VCF subset requires one ALT allele");
                    }
                    const char *contig_name =
                        bcf_hdr_id2name(header.get(), record->rid);
                    if (contig_name == nullptr) {
                        throw VariantCatalogError(
                            "VCF row has no valid contig name");
                    }
                    const auto found = contig_indices.find(contig_name);
                    if (found == contig_indices.end()) {
                        throw VariantCatalogError(
                            "VCF row names an unknown contig");
                    }
                    const std::uint32_t contig_index = found->second;
                    if (record->pos < 0
                        || static_cast<std::uint64_t>(record->pos)
                            >= std::numeric_limits<std::uint32_t>::max()) {
                        throw VariantCatalogError(
                            "VCF POS must be a positive uint32 decimal");
                    }
                    const std::uint32_t position =
                        static_cast<std::uint32_t>(record->pos) + 1U;
                    if (have_previous_record
                        && (contig_index < previous_contig
                            || (contig_index == previous_contig
                                && position < previous_position))) {
                        throw VariantCatalogError(
                            "VCF rows must be sorted in reference "
                            "order");
                    }
                    previous_contig = contig_index;
                    previous_position = position;
                    have_previous_record = true;

                    model::Bases reference = parse_allele(
                        record->d.allele[0], "REF");
                    model::Bases alternate = parse_allele(
                        record->d.allele[1], "ALT");
                    const Genotype genotype = parse_genotype(
                        header.get(), record.get(), genotype_buffer);
                    if (genotype.first == 0U && genotype.second == 0U) {
                        continue;
                    }

                    auto &contig_events = variants_by_contig_[contig_index];
                    if (contig_events.size()
                        == std::numeric_limits<std::uint32_t>::max()) {
                        throw VariantCatalogError(
                            "VCF event count for one contig exceeds uint32");
                    }
                    const std::uint64_t key = rng::derive_key(
                        master_seed,
                        rng::Stage::haplotype,
                        static_cast<std::uint32_t>(contig_index));
                    const auto mask = resolve_alt_haplotypes(
                        genotype, key, contig_events.size());
                    Variant event = normalize_event(
                        contig_index,
                        position,
                        std::move(reference),
                        std::move(alternate),
                        mask,
                        reference_catalog[contig_index].length);
                    const std::uint64_t address = catalog_address(
                        contig_index, contig_events.size());
                    event.id = unique_variant_id(
                        record->d.id == nullptr ? std::string() : record->d.id,
                        "aaracf_",
                        address,
                        used_variant_ids);
                    event.source = model::VariantSource::vcf;
                    if (!contig_events.empty()) {
                        const Variant &previous = contig_events.back();
                        if (event_less(event, previous)) {
                            throw VariantCatalogError(
                                "VCF normalized variants are not in canonical "
                                "order");
                        }
                        if (conflicts(previous, event)) {
                            throw VariantCatalogError(
                                "VCF normalized variants overlap or share an "
                                "insertion anchor");
                        }
                    }
                    contig_events.push_back(std::move(event));
                    if (variant_count_
                        == std::numeric_limits<std::uint64_t>::max()) {
                        throw VariantCatalogError(
                            "VCF event count exceeds uint64");
                    }
                    ++variant_count_;
                } catch (const VariantCatalogError &error) {
                    throw VariantCatalogError(
                        "VCF line " + std::to_string(file->lineno) + ": "
                        + error.what());
                }
            }
        });
        if (snapshot.file_sha256() != file_sha256_) {
            throw VariantCatalogError("VCF digest changed");
        }
    } catch (const VariantCatalogError &) {
        throw;
    } catch (const std::exception &error) {
        throw VariantCatalogError(error.what());
    }
}

const std::vector<Variant> &VariantFile::variants(
    std::uint32_t contig_index) const
{
    if (contig_index >= variants_by_contig_.size()) {
        throw VariantCatalogError("VCF contig index is out of range");
    }
    return variants_by_contig_[contig_index];
}

std::uint64_t VariantFile::variant_count() const noexcept
{
    return variant_count_;
}

const crypto::Sha256Digest &VariantFile::file_sha256() const noexcept
{
    return file_sha256_;
}

ContigVariants::ContigVariants(
    const model::Bases &reference_bases,
    const std::vector<Variant> &variants,
    std::uint32_t contig_index)
    : contig_index_(contig_index), variants_(variants)
{
    if (reference_bases.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw VariantCatalogError("VCF reference contig exceeds uint32");
    }
    reference_length_ = static_cast<std::uint32_t>(reference_bases.size());
    for (std::size_t index = 0; index < variants_.size(); ++index) {
        const Variant &event = variants_[index];
        if (event.contig_index != contig_index) {
            throw VariantCatalogError("VCF event belongs to another contig");
        }
        if (!model::is_haplotype_mask(
                static_cast<std::uint8_t>(event.alt_haplotypes))) {
            throw VariantCatalogError("VCF event has an invalid haplotype mask");
        }
        const auto valid_bases = [](const model::Bases &bases) {
            return std::all_of(
                bases.begin(), bases.end(),
                [](std::uint8_t base) {return base <= 3U;});
        };
        if (!valid_bases(event.ref_bases) || !valid_bases(event.alt_bases)) {
            throw VariantCatalogError("VCF event contains an invalid base code");
        }
        if (event.reference_start > event.reference_end
            || event.reference_end > reference_bases.size()) {
            throw VariantCatalogError("VCF event interval exceeds the contig");
        }
        if (event.ref_bases.size()
            != event.reference_end - event.reference_start) {
            throw VariantCatalogError("VCF REF length disagrees with its interval");
        }
        switch (event.kind) {
        case model::VariantKind::snv:
            if (event.ref_bases.size() != 1U
                || event.alt_bases.size() != 1U) {
                throw VariantCatalogError("VCF SNV is not one base to one base");
            }
            break;
        case model::VariantKind::insertion:
            if (!event.ref_bases.empty()
                || event.reference_start != event.reference_end
                || event.alt_bases.empty()
                || event.alt_bases.size() > maximum_indel_bases) {
                throw VariantCatalogError("VCF insertion shape is invalid");
            }
            break;
        case model::VariantKind::deletion:
            if (event.ref_bases.empty()
                || event.ref_bases.size() > maximum_indel_bases
                || !event.alt_bases.empty()) {
                throw VariantCatalogError("VCF deletion shape is invalid");
            }
            break;
        default:
            throw VariantCatalogError("VCF event kind is invalid");
        }
        if (!std::equal(
                event.ref_bases.begin(),
                event.ref_bases.end(),
                reference_bases.begin()
                    + static_cast<std::ptrdiff_t>(event.reference_start))) {
            throw VariantCatalogError("VCF REF does not match the reference snapshot");
        }
        if (index != 0U && conflicts(variants_[index - 1U], event)) {
            throw VariantCatalogError("VCF contig variants overlap");
        }
        if (index != 0U && event_less(event, variants_[index - 1U])) {
            throw VariantCatalogError("VCF contig variants are not canonical");
        }
    }
}

const std::vector<Variant> &ContigVariants::variants() const noexcept
{
    return variants_;
}

std::uint32_t ContigVariants::contig_index() const noexcept
{
    return contig_index_;
}

std::uint32_t ContigVariants::reference_length() const noexcept
{
    return reference_length_;
}

} // namespace htsim::variant

// ---- haplotype_projector --------------------------------------------------------

namespace htsim::haplotype {
namespace {

void require_append_capacity(std::size_t current, std::size_t additional)
{
    const std::size_t maximum = std::numeric_limits<std::uint32_t>::max();
    if (current > maximum || additional > maximum - current) {
        throw ProjectionError(
            ProjectionFailure::capacity,
            "projected haplotype interval exceeds uint32 bases");
    }
}

std::uint8_t phased_haplotype(model::HaplotypeMask mask)
{
    switch (mask) {
    case model::HaplotypeMask::haplotype_1: return 0U;
    case model::HaplotypeMask::haplotype_2: return 1U;
    case model::HaplotypeMask::both: return 255U;
    }
    throw ProjectionError(
        ProjectionFailure::invariant,
        "variant has an invalid haplotype mask");
}

model::Variant to_model_variant(
    const variant::Variant &event,
    std::uint32_t variant_index)
{
    return {
        variant_index,
        event.id,
        event.source,
        event.kind,
        phased_haplotype(event.alt_haplotypes),
        event.reference_start,
        event.reference_end,
        event.ref_bases,
        event.alt_bases,
    };
}

void append_reference(
    ProjectedInterval &projection,
    const model::Bases &reference,
    std::uint32_t begin,
    std::uint32_t end)
{
    if (begin > end) {
        throw ProjectionError(
            ProjectionFailure::invariant,
            "variant projection moved the reference cursor backwards");
    }
    const std::size_t count = static_cast<std::size_t>(end - begin);
    require_append_capacity(projection.template_bases.size(), count);
    const auto first = reference.begin() + static_cast<std::ptrdiff_t>(begin);
    const auto last = reference.begin() + static_cast<std::ptrdiff_t>(end);
    projection.template_bases.insert(
        projection.template_bases.end(), first, last);
    for (std::uint32_t position = begin; position < end; ++position) {
        projection.reference_positions.push_back(position);
        projection.base_variant_indices.push_back(model::no_variant_index);
    }
}

void append_alternate(
    ProjectedInterval &projection,
    const variant::Variant &event,
    std::uint32_t variant_index)
{
    require_append_capacity(
        projection.template_bases.size(), event.alt_bases.size());
    projection.template_bases.insert(
        projection.template_bases.end(),
        event.alt_bases.begin(),
        event.alt_bases.end());
    if (event.kind == model::VariantKind::insertion) {
        projection.reference_positions.insert(
            projection.reference_positions.end(), event.alt_bases.size(), -1);
    } else {
        for (std::uint32_t position = event.reference_start;
             position < event.reference_end;
             ++position) {
            projection.reference_positions.push_back(position);
        }
    }
    projection.base_variant_indices.insert(
        projection.base_variant_indices.end(), event.alt_bases.size(), variant_index);
}

bool insertion_belongs_to_interval(
    std::uint32_t anchor,
    std::uint32_t interval_start,
    std::uint32_t interval_end,
    ProjectionBoundaryPolicy boundary_policy) noexcept
{
    if (anchor < interval_start || anchor > interval_end) {return false;}
    if (anchor == interval_start && anchor == interval_end) {
        return boundary_policy.include_start_anchor_insertion
            || boundary_policy.include_end_anchor_insertion;
    }
    if (anchor == interval_start) {
        return boundary_policy.include_start_anchor_insertion;
    }
    if (anchor == interval_end) {
        return boundary_policy.include_end_anchor_insertion;
    }
    return true;
}

} // namespace

ProjectedInterval project_interval(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint8_t zero_based_haplotype,
    std::uint32_t reference_start,
    std::uint32_t reference_end)
{
    if (contig.length > std::numeric_limits<std::uint32_t>::max()) {
        throw ProjectionError(
            ProjectionFailure::invalid_input,
            "materialized contig length is inconsistent");
    }
    if (reference_start >= reference_end) {
        throw ProjectionError(
            ProjectionFailure::invalid_input,
            "projection interval is empty or reversed");
    }
    return project_interval(
        contig,
        variants,
        zero_based_haplotype,
        reference_start,
        reference_end,
        ProjectionBoundaryPolicy{
            true,
            reference_end == static_cast<std::uint32_t>(contig.length),
        });
}

ProjectedInterval project_interval(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint8_t zero_based_haplotype,
    std::uint32_t reference_start,
    std::uint32_t reference_end,
    ProjectionBoundaryPolicy boundary_policy)
{
    if (zero_based_haplotype > 1U) {
        throw ProjectionError(
            ProjectionFailure::invalid_input,
            "haplotype must be zero or one");
    }
    if (contig.length != contig.bases.size()
        || contig.length > std::numeric_limits<std::uint32_t>::max()) {
        throw ProjectionError(
            ProjectionFailure::invalid_input,
            "materialized contig length is inconsistent");
    }
    const auto contig_length = static_cast<std::uint32_t>(contig.length);
    if (variants.contig_index() != contig.index
        || variants.reference_length() != contig_length) {
        throw ProjectionError(
            ProjectionFailure::invalid_input,
            "variant set does not match the materialized contig");
    }
    if (reference_start > reference_end || reference_end > contig_length) {
        throw ProjectionError(
            ProjectionFailure::invalid_input,
            "projection interval is reversed or outside its contig");
    }

    ProjectedInterval projection;
    projection.contig_index = contig.index;
    projection.haplotype = zero_based_haplotype;
    projection.reference_start = reference_start;
    projection.reference_end = reference_end;
    projection.template_bases.reserve(reference_end - reference_start);
    projection.reference_positions.reserve(reference_end - reference_start);
    projection.base_variant_indices.reserve(reference_end - reference_start);

    std::uint32_t cursor = reference_start;
    const auto &variant_records = variants.variants();
    if (variant_records.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw ProjectionError(
            ProjectionFailure::capacity,
            "variant event count exceeds uint32");
    }
    auto event = std::lower_bound(
        variant_records.begin(),
        variant_records.end(),
        reference_start,
        [](const variant::Variant &candidate, std::uint32_t position) {
            return candidate.reference_start < position;
        });
    if (event != variant_records.begin()) {--event;}
    for (; event != variant_records.end(); ++event) {
        const std::size_t ordinal = static_cast<std::size_t>(
            std::distance(variant_records.begin(), event));
        const variant::Variant &current = *event;
        if (current.contig_index != contig.index
            || current.reference_end > contig.bases.size()
            || !std::equal(
                current.ref_bases.begin(),
                current.ref_bases.end(),
                contig.bases.begin()
                    + static_cast<std::ptrdiff_t>(current.reference_start))) {
            throw ProjectionError(
                ProjectionFailure::invariant,
                "variant set does not match the materialized contig");
        }
        const bool selected_end_insertion =
            current.kind == model::VariantKind::insertion
            && current.reference_start == reference_end
            && boundary_policy.include_end_anchor_insertion;
        if (current.reference_start >= reference_end
            && !selected_end_insertion) {
            break;
        }
        if (!model::mask_contains(
                current.alt_haplotypes, zero_based_haplotype)) {
            continue;
        }
        const std::uint32_t variant_index = static_cast<std::uint32_t>(ordinal);
        if (variant_index == model::no_variant_index) {
            throw ProjectionError(
                ProjectionFailure::capacity,
                "variant event ordinal uses the no-event sentinel");
        }

        if (current.kind == model::VariantKind::insertion) {
            if (!insertion_belongs_to_interval(
                    current.reference_start,
                    reference_start,
                    reference_end,
                    boundary_policy)) {
                continue;
            }
            append_reference(
                projection, contig.bases, cursor, current.reference_start);
            projection.variants.push_back(to_model_variant(current, variant_index));
            append_alternate(projection, current, variant_index);
            cursor = current.reference_start;
            continue;
        }

        if (current.reference_end <= reference_start) {continue;}
        if (current.reference_start >= reference_end) {break;}
        if (current.reference_start < reference_start
            || current.reference_end > reference_end) {
            throw ProjectionError(
                ProjectionFailure::boundary_cut,
                "projection interval cuts through an active variant");
        }
        append_reference(
            projection, contig.bases, cursor, current.reference_start);
        projection.variants.push_back(to_model_variant(current, variant_index));
        if (current.kind == model::VariantKind::snv) {
            append_alternate(projection, current, variant_index);
        }
        cursor = current.reference_end;
    }
    append_reference(projection, contig.bases, cursor, reference_end);

    if (projection.template_bases.empty()) {
        throw ProjectionError(
            ProjectionFailure::empty_projection,
            "active variants removed the complete interval");
    }
    if (projection.reference_positions.size() != projection.template_bases.size()
        || projection.base_variant_indices.size() != projection.template_bases.size()) {
        throw ProjectionError(
            ProjectionFailure::invariant,
            "projected per-base arrays are inconsistent");
    }
    return projection;
}

} // namespace htsim::haplotype

// ---- haplotype_layout --------------------------------------------------------

namespace htsim::haplotype {
namespace {

struct StartInterval {
    std::uint32_t lower = 0;
    std::uint32_t upper = 0;
    std::uint64_t payload_bytes = 0;
};

std::uint32_t population_count(std::uint64_t value) noexcept
{
    std::uint32_t count = 0U;
    while (value != 0U) {
        value &= value - 1U;
        ++count;
    }
    return count;
}

void require_identity(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint8_t haplotype)
{
    if (haplotype > 1U) {
        throw HaplotypeLayoutError("haplotype must be zero or one");
    }
    if (contig.length != contig.bases.size()
        || contig.length > std::numeric_limits<std::uint32_t>::max()
        || variants.contig_index() != contig.index
        || variants.reference_length() != contig.length) {
        throw HaplotypeLayoutError(
            "haplotype layout inputs do not identify one uint32 contig");
    }
}

std::uint64_t event_payload_bytes(const variant::Variant &event)
{
    const std::uint64_t bases = static_cast<std::uint64_t>(
        event.ref_bases.size() + event.alt_bases.size());
    if (bases > std::numeric_limits<std::uint64_t>::max() - UINT64_C(32)) {
        throw HaplotypeLayoutError("variant payload size exceeds uint64");
    }
    return UINT64_C(32) + bases;
}

} // namespace

HaplotypeLayout::HaplotypeLayout(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint8_t zero_based_haplotype,
    bool materialize_bases)
    : contig_index_(contig.index),
      haplotype_(zero_based_haplotype),
      reference_length_(static_cast<std::uint32_t>(contig.length)),
      materialized_(materialize_bases)
{
    require_identity(contig, variants, zero_based_haplotype);
    if (materialized_) {bases_.reserve(contig.bases.size());}

    std::uint32_t reference_cursor = 0U;
    std::int64_t cumulative_delta = 0;
    const auto &variant_records = variants.variants();
    for (std::size_t variant_index = 0U;
         variant_index < variant_records.size();
         ++variant_index) {
        if (variant_index >= model::no_variant_index) {
            throw HaplotypeLayoutError(
                "variant event ordinal uses the no-event sentinel");
        }
        const variant::Variant &event = variant_records[variant_index];
        if (!model::mask_contains(
                event.alt_haplotypes, zero_based_haplotype)) {
            continue;
        }
        append_reference_run(
            contig.bases, reference_cursor, event.reference_start);
        const std::uint32_t haplotype_before = length_;
        switch (event.kind) {
        case model::VariantKind::insertion: {
            FragmentBoundary &before = upsert_exception(
                haplotype_before, event.reference_start);
            before.right_reference_start = event.reference_start;
            before.include_insertion_in_right_fragment = true;
            append_alternate_run(event, event.reference_start);
            FragmentBoundary &after = upsert_exception(
                length_, event.reference_start);
            after.left_reference_end = event.reference_start;
            after.right_reference_start = event.reference_start;
            after.include_insertion_in_left_fragment = true;
            after.include_insertion_in_right_fragment = false;
            cumulative_delta += static_cast<std::int64_t>(
                event.alt_bases.size());
            reference_cursor = event.reference_start;
            break;
        }
        case model::VariantKind::snv:
            append_alternate_run(event, event.reference_start);
            reference_cursor = event.reference_end;
            break;
        case model::VariantKind::deletion: {
            FragmentBoundary &boundary = upsert_exception(
                haplotype_before, event.reference_start);
            // Consecutive deletions collapse at the same physical boundary.
            // Preserve the left edge of the first and extend the right edge.
            boundary.right_reference_start = event.reference_end;
            boundary.include_insertion_in_right_fragment = true;
            cumulative_delta -= static_cast<std::int64_t>(
                event.reference_end - event.reference_start);
            reference_cursor = event.reference_end;
            break;
        }
        default:
            throw HaplotypeLayoutError(
                "variant kind is outside the prepared variant set");
        }
        active_events_.push_back({
            event.reference_start,
            event.reference_end,
            event.kind,
            cumulative_delta,
            haplotype_before,
            length_,
            static_cast<std::uint32_t>(variant_index),
        });
    }
    append_reference_run(contig.bases, reference_cursor, reference_length_);

    ambiguous_prefix_.reserve(ambiguous_words_.size() + 1U);
    ambiguous_prefix_.push_back(0U);
    for (const std::uint64_t word : ambiguous_words_) {
        const std::uint32_t count = population_count(word);
        if (count > std::numeric_limits<std::uint32_t>::max()
                - ambiguous_prefix_.back()) {
            throw HaplotypeLayoutError("haplotype N count exceeds uint32");
        }
        ambiguous_prefix_.push_back(ambiguous_prefix_.back() + count);
    }
}

const model::Bases &HaplotypeLayout::bases() const
{
    if (!materialized_) {
        throw HaplotypeLayoutError("haplotype bases were not materialized");
    }
    return bases_;
}

void HaplotypeLayout::append_base(std::uint8_t base)
{
    if (base > 4U) {
        throw HaplotypeLayoutError(
            "haplotype contains a base outside protocol encoding");
    }
    if (length_ == std::numeric_limits<std::uint32_t>::max()) {
        throw HaplotypeLayoutError("haplotype length exceeds uint32");
    }
    const std::size_t word_index = static_cast<std::size_t>(length_ / 64U);
    if (word_index == ambiguous_words_.size()) {
        ambiguous_words_.push_back(0U);
    }
    if (base == 4U) {
        ambiguous_words_.back() |=
            UINT64_C(1) << static_cast<unsigned>(length_ % 64U);
    }
    if (materialized_) {bases_.push_back(base);}
    ++length_;
}

void HaplotypeLayout::append_reference_run(
    const model::Bases &reference,
    std::uint32_t begin,
    std::uint32_t end)
{
    if (begin > end || end > reference.size()) {
        throw HaplotypeLayoutError("reference cursor is outside the contig");
    }
    if (begin == end) {return;}
    const std::uint32_t haplotype_begin = length_;
    for (std::uint32_t position = begin; position < end; ++position) {
        append_base(reference[position]);
    }
    reference_runs_.push_back({haplotype_begin, length_, begin});
}

void HaplotypeLayout::append_alternate_run(
    const variant::Variant &event,
    std::uint32_t reference_begin)
{
    const std::uint32_t haplotype_begin = length_;
    for (const std::uint8_t base : event.alt_bases) {append_base(base);}
    if (event.kind == model::VariantKind::snv) {
        reference_runs_.push_back({haplotype_begin, length_, reference_begin});
    }
}

FragmentBoundary &HaplotypeLayout::upsert_exception(
    std::uint32_t haplotype_offset,
    std::uint32_t reference_offset)
{
    if (!boundary_exceptions_.empty()
        && boundary_exceptions_.back().haplotype_offset == haplotype_offset) {
        return boundary_exceptions_.back().boundary;
    }
    boundary_exceptions_.push_back({
        haplotype_offset,
        FragmentBoundary{reference_offset, reference_offset, false, true},
    });
    return boundary_exceptions_.back().boundary;
}

std::optional<FragmentBoundary> HaplotypeLayout::boundary(
    std::uint32_t haplotype_offset) const
{
    if (haplotype_offset > length_) {
        throw HaplotypeLayoutError("haplotype boundary exceeds its sequence");
    }
    const auto exceptional = std::lower_bound(
        boundary_exceptions_.begin(),
        boundary_exceptions_.end(),
        haplotype_offset,
        [](const BoundaryException &candidate, std::uint32_t offset) {
            return candidate.haplotype_offset < offset;
        });
    if (exceptional != boundary_exceptions_.end()
        && exceptional->haplotype_offset == haplotype_offset) {
        return exceptional->boundary;
    }

    const auto run = std::upper_bound(
        reference_runs_.begin(),
        reference_runs_.end(),
        haplotype_offset,
        [](std::uint32_t offset, const ReferenceRun &candidate) {
            return offset < candidate.haplotype_begin;
        });
    if (run == reference_runs_.begin()) {return std::nullopt;}
    const ReferenceRun &candidate = *(run - 1);
    if (haplotype_offset < candidate.haplotype_begin
        || haplotype_offset > candidate.haplotype_end) {
        return std::nullopt;
    }
    const std::uint32_t reference_offset = candidate.reference_begin
        + (haplotype_offset - candidate.haplotype_begin);
    return FragmentBoundary{
        reference_offset, reference_offset, false, true};
}

std::optional<std::uint32_t> HaplotypeLayout::boundary_before_reference(
    std::uint32_t reference_offset) const
{
    if (reference_offset > reference_length_) {
        throw HaplotypeLayoutError("reference boundary exceeds its contig");
    }
    const auto first_not_before = std::lower_bound(
        active_events_.begin(),
        active_events_.end(),
        reference_offset,
        [](const ActiveEventCoordinate &event, std::uint32_t offset) {
            return event.reference_start < offset;
        });
    std::int64_t delta = 0;
    if (first_not_before != active_events_.begin()) {
        const ActiveEventCoordinate &previous = *(first_not_before - 1);
        if (previous.kind == model::VariantKind::deletion
            && reference_offset > previous.reference_start
            && reference_offset < previous.reference_end) {
            return std::nullopt;
        }
        delta = previous.cumulative_delta_after;
    }
    const std::int64_t projected = static_cast<std::int64_t>(reference_offset)
        + delta;
    if (projected < 0
        || projected > static_cast<std::int64_t>(
            std::numeric_limits<std::uint32_t>::max())) {
        throw HaplotypeLayoutError(
            "reference boundary projection exceeds uint32");
    }
    return static_cast<std::uint32_t>(projected);
}

std::uint32_t HaplotypeLayout::n_rank(std::uint32_t haplotype_end) const
{
    if (haplotype_end > length_
        || ambiguous_prefix_.size() != ambiguous_words_.size() + 1U) {
        throw HaplotypeLayoutError("haplotype N-rank boundary is invalid");
    }
    const std::size_t word = static_cast<std::size_t>(haplotype_end / 64U);
    const unsigned bit = static_cast<unsigned>(haplotype_end % 64U);
    std::uint32_t result = ambiguous_prefix_.at(word);
    if (bit != 0U) {
        const std::uint64_t mask = (UINT64_C(1) << bit) - 1U;
        result += population_count(ambiguous_words_.at(word) & mask);
    }
    return result;
}

std::uint32_t HaplotypeLayout::ambiguous_count(
    std::uint32_t haplotype_begin,
    std::uint32_t haplotype_end) const
{
    if (haplotype_begin > haplotype_end || haplotype_end > length_) {
        throw HaplotypeLayoutError("haplotype ambiguity interval is invalid");
    }
    return n_rank(haplotype_end) - n_rank(haplotype_begin);
}

std::uint64_t HaplotypeLayout::maximum_variant_payload_bytes(
    const variant::ContigVariants &variants,
    std::uint32_t physical_span) const
{
    if (variants.contig_index() != contig_index_
        || variants.reference_length() != reference_length_
        || physical_span == 0U) {
        throw HaplotypeLayoutError(
            "variant payload window has the wrong identity or length");
    }
    if (physical_span > length_ || active_events_.empty()) {return 0U;}
    const std::uint32_t last_start = length_ - physical_span;
    const auto interval_at = [&](std::size_t index)
        -> std::optional<StartInterval> {
        const ActiveEventCoordinate &coordinate = active_events_.at(index);
        if (coordinate.event_ordinal >= variants.variants().size()
            || coordinate.haplotype_begin > coordinate.haplotype_end
            || coordinate.haplotype_end > length_) {
            throw HaplotypeLayoutError(
                "active variant coordinate is inconsistent");
        }
        std::uint32_t lower = 0U;
        std::uint32_t upper = 0U;
        if (coordinate.haplotype_begin == coordinate.haplotype_end) {
            // A deletion belongs only to a window that crosses its collapsed
            // boundary; a window starting or ending there excludes it.
            if (coordinate.haplotype_begin == 0U) {return std::nullopt;}
            lower = coordinate.haplotype_begin >= physical_span
                ? coordinate.haplotype_begin - physical_span + 1U
                : 0U;
            upper = std::min(
                coordinate.haplotype_begin - 1U, last_start);
        } else {
            const std::uint32_t event_span =
                coordinate.haplotype_end - coordinate.haplotype_begin;
            if (event_span > physical_span) {return std::nullopt;}
            lower = coordinate.haplotype_end > physical_span
                ? coordinate.haplotype_end - physical_span
                : 0U;
            upper = std::min(coordinate.haplotype_begin, last_start);
        }
        if (lower > upper) {return std::nullopt;}
        return StartInterval{
            lower,
            upper,
            event_payload_bytes(
                variants.variants()[coordinate.event_ordinal]),
        };
    };

    // Canonical non-overlapping variants yield monotone start intervals.  Check
    // that invariant explicitly so a future event kind cannot invalidate the
    // allocation-free two-cursor sweep.
    std::optional<StartInterval> previous;
    for (std::size_t index = 0U; index < active_events_.size(); ++index) {
        const auto current = interval_at(index);
        if (!current) {continue;}
        if (previous
            && (current->lower < previous->lower
                || current->upper < previous->upper)) {
            throw HaplotypeLayoutError(
                "variant payload start intervals are not monotone");
        }
        previous = current;
    }

    std::uint64_t current_bytes = 0U;
    std::uint64_t maximum_bytes = 0U;
    std::size_t add_index = 0U;
    std::size_t remove_index = 0U;
    while (add_index < active_events_.size()) {
        auto first = interval_at(add_index);
        if (!first) {
            ++add_index;
            continue;
        }
        const std::uint32_t position = first->lower;
        while (remove_index < add_index) {
            const auto expired = interval_at(remove_index);
            if (!expired) {
                ++remove_index;
                continue;
            }
            if (expired->upper >= position) {break;}
            if (current_bytes < expired->payload_bytes) {
                throw HaplotypeLayoutError(
                    "variant payload sweep underflowed");
            }
            current_bytes -= expired->payload_bytes;
            ++remove_index;
        }
        while (add_index < active_events_.size()) {
            const auto added = interval_at(add_index);
            if (!added) {
                ++add_index;
                continue;
            }
            if (added->lower != position) {break;}
            if (added->payload_bytes
                > std::numeric_limits<std::uint64_t>::max() - current_bytes) {
                throw HaplotypeLayoutError(
                    "variant payload sweep exceeds uint64");
            }
            current_bytes += added->payload_bytes;
            ++add_index;
        }
        maximum_bytes = std::max(maximum_bytes, current_bytes);
    }
    return maximum_bytes;
}

ProjectedInterval HaplotypeLayout::project(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint32_t haplotype_begin,
    std::uint32_t haplotype_end) const
{
    require_identity(contig, variants, haplotype_);
    if (contig.index != contig_index_
        || variants.reference_length() != reference_length_
        || haplotype_begin >= haplotype_end
        || haplotype_end > length_) {
        throw HaplotypeLayoutError(
            "haplotype projection slice is empty or has the wrong identity");
    }
    const auto left = boundary(haplotype_begin);
    const auto right = boundary(haplotype_end);
    if (!left || !right) {
        throw ProjectionError(
            ProjectionFailure::boundary_cut,
            "haplotype fragment boundary cuts through an insertion");
    }
    if (left->right_reference_start > right->left_reference_end) {
        throw HaplotypeLayoutError(
            "haplotype fragment reference envelope is reversed");
    }
    ProjectedInterval projection = project_interval(
        contig,
        variants,
        haplotype_,
        left->right_reference_start,
        right->left_reference_end,
        ProjectionBoundaryPolicy{
            left->include_insertion_in_right_fragment,
            right->include_insertion_in_left_fragment,
        });
    if (projection.template_bases.size()
        != static_cast<std::size_t>(haplotype_end - haplotype_begin)) {
        throw HaplotypeLayoutError(
            "haplotype slice length disagrees with its typed projection");
    }
    return projection;
}

} // namespace htsim::haplotype

// ---- mutation_catalog --------------------------------------------------------

namespace htsim::variant {
namespace {

inline constexpr std::uint64_t mutation_draw = 0U;
inline constexpr std::uint64_t indel_draw = 1U;
inline constexpr std::uint64_t substitution_draw = 2U;
inline constexpr std::uint64_t deletion_draw = 3U;
inline constexpr std::uint64_t extension_draw_begin = 4U;
inline constexpr std::uint64_t insertion_base_draw_begin = 8U;
inline constexpr std::uint64_t homozygous_draw = 16U;
inline constexpr std::uint64_t haplotype_draw = 17U;
inline constexpr double one_third = 0x1.5555555555555p-2;

void validate_probability(double value, std::string_view name)
{
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
        throw MutationCatalogError(
            std::string(name) + " must be finite and in [0, 1]");
    }
}

std::uint32_t validate_inputs(
    const reference::Contig &contig,
    const MutationParameters &parameters)
{
    if (contig.length != contig.bases.size()
        || contig.length > std::numeric_limits<std::uint32_t>::max()) {
        throw MutationCatalogError(
            "de novo mutation input is not one materialized uint32 contig");
    }
    if (!std::all_of(
            contig.bases.begin(),
            contig.bases.end(),
            [](std::uint8_t base) {return base <= 4U;})) {
        throw MutationCatalogError(
            "de novo mutation input contains an invalid reference base");
    }
    validate_probability(parameters.mutation_rate, "mutation rate");
    validate_probability(parameters.indel_fraction, "indel fraction");
    validate_probability(
        parameters.indel_extension_probability,
        "indel extension probability");
    return static_cast<std::uint32_t>(contig.length);
}

model::HaplotypeMask choose_haplotype_mask(
    std::uint64_t key,
    std::uint32_t reference_position,
    bool homozygous_only)
{
    if (homozygous_only
        || rng::bernoulli(
            key,
            reference_position,
            homozygous_draw,
            one_third)) {
        return model::HaplotypeMask::both;
    }
    return rng::bernoulli(
               key,
               reference_position,
               haplotype_draw,
               0.5)
        ? model::HaplotypeMask::haplotype_1
        : model::HaplotypeMask::haplotype_2;
}

std::uint32_t draw_indel_length(
    std::uint64_t key,
    std::uint32_t reference_position,
    double extension_probability)
{
    std::uint32_t length = 1U;
    while (length < maximum_indel_bases
           && rng::bernoulli(
               key,
               reference_position,
               extension_draw_begin + length - 1U,
               extension_probability)) {
        ++length;
    }
    return length;
}

model::Bases insertion_bases(
    std::uint64_t key,
    std::uint32_t reference_position,
    std::uint32_t length)
{
    model::Bases result;
    result.reserve(length);
    for (std::uint32_t offset = 0U; offset < length; ++offset) {
        result.push_back(static_cast<std::uint8_t>(rng::bounded_integer(
            key,
            reference_position,
            insertion_base_draw_begin + offset,
            4U)));
    }
    return result;
}

std::uint32_t resolved_deletion_length(
    const model::Bases &reference,
    std::uint32_t reference_position,
    std::uint32_t requested_length)
{
    std::uint32_t length = 0U;
    const std::size_t available =
        reference.size() - static_cast<std::size_t>(reference_position);
    while (length < requested_length
           && length < available
           && reference[reference_position + length] < 4U) {
        ++length;
    }
    if (length == 0U) {
        throw MutationCatalogError(
            "de novo deletion does not begin on a resolved reference base");
    }
    return length;
}

void require_event_capacity(const std::vector<Variant> &variants)
{
    if (variants.size()
        >= static_cast<std::size_t>(model::no_variant_index)) {
        throw MutationCatalogError(
            "de novo event count would use the no-event sentinel");
    }
}

} // namespace

std::vector<Variant> generate_de_novo_events(
    const reference::Contig &contig,
    std::uint64_t master_seed,
    const MutationParameters &parameters)
{
    try {
        const std::uint32_t reference_length =
            validate_inputs(contig, parameters);
        const std::uint64_t key = rng::derive_key(
            master_seed, rng::Stage::mutation, contig.index);
        std::vector<Variant> variants;
        std::uint32_t position = 0U;
        while (position < reference_length) {
            const std::uint8_t reference_base = contig.bases[position];
            if (reference_base == 4U
                || !rng::bernoulli(
                    key,
                    position,
                    mutation_draw,
                    parameters.mutation_rate)) {
                ++position;
                continue;
            }

            require_event_capacity(variants);
            const model::HaplotypeMask haplotypes =
                choose_haplotype_mask(
                    key, position, parameters.homozygous_only);
            if (!rng::bernoulli(
                    key,
                    position,
                    indel_draw,
                    parameters.indel_fraction)) {
                const std::uint8_t offset = static_cast<std::uint8_t>(
                    rng::bounded_integer(
                        key, position, substitution_draw, 3U));
                const std::uint8_t alternate_base =
                    static_cast<std::uint8_t>(
                        (reference_base + offset + 1U) & 3U);
                variants.push_back({
                    contig.index,
                    position,
                    position + 1U,
                    model::VariantKind::snv,
                    {reference_base},
                    {alternate_base},
                    haplotypes,
                    "varsim_" + hexadecimal(catalog_address(
                        contig.index, variants.size())),
                    model::VariantSource::de_novo,
                });
                ++position;
                continue;
            }

            const std::uint32_t requested_length = draw_indel_length(
                key,
                position,
                parameters.indel_extension_probability);
            if (rng::bernoulli(
                    key, position, deletion_draw, 0.5)) {
                const std::uint32_t deletion_length =
                    resolved_deletion_length(
                        contig.bases, position, requested_length);
                const std::uint32_t end = position + deletion_length;
                variants.push_back({
                    contig.index,
                    position,
                    end,
                    model::VariantKind::deletion,
                    model::Bases(
                        contig.bases.begin()
                            + static_cast<std::ptrdiff_t>(position),
                        contig.bases.begin()
                            + static_cast<std::ptrdiff_t>(end)),
                    {},
                    haplotypes,
                    "varsim_" + hexadecimal(catalog_address(
                        contig.index, variants.size())),
                    model::VariantSource::de_novo,
                });
                position = end;
                continue;
            }

            const std::uint32_t anchor = position + 1U;
            variants.push_back({
                contig.index,
                anchor,
                anchor,
                model::VariantKind::insertion,
                {},
                insertion_bases(key, position, requested_length),
                haplotypes,
                "varsim_" + hexadecimal(catalog_address(
                    contig.index, variants.size())),
                model::VariantSource::de_novo,
            });
            // The typed catalog forbids an insertion and another event at the
            // same anchor. The insertion follows `position`, so the immediately
            // following reference base is the one suppressed.
            position = anchor < reference_length ? anchor + 1U : anchor;
        }
        return variants;
    } catch (const MutationCatalogError &) {
        throw;
    } catch (const std::exception &error) {
        throw MutationCatalogError(error.what());
    }
}

namespace {

char vcf_base(std::uint8_t base)
{
    constexpr char alphabet[] = "ACGTN";
    if (base >= sizeof(alphabet) - 1U) {
        throw MutationCatalogError("VCF export received an invalid allele base");
    }
    return alphabet[base];
}

void write_vcf_bases(std::ostream &sink, const model::Bases &bases)
{
    for (const std::uint8_t base : bases) {sink.put(vcf_base(base));}
}

const char *vcf_genotype(model::HaplotypeMask mask)
{
    switch (mask) {
    case model::HaplotypeMask::haplotype_1: return "1|0";
    case model::HaplotypeMask::haplotype_2: return "0|1";
    case model::HaplotypeMask::both: return "1|1";
    }
    throw MutationCatalogError("VCF export received an invalid haplotype mask");
}

void require_vcf_token(std::string_view value, const char *field)
{
    if (value.empty()
        || value.find_first_of("\t\r\n") != std::string_view::npos) {
        throw MutationCatalogError(
            std::string("VCF export received an invalid ") + field);
    }
}

} // namespace

void write_vcf_header(std::ostream &sink)
{
    sink
        << "##fileformat=VCFv4.3\n"
        << "##source=BSReadSim\n"
        << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIMULATED\n";
    if (!sink) {throw MutationCatalogError("failed while writing the VCF header");}
}

void write_vcf_contig(
    std::ostream &sink,
    const reference::Contig &contig,
    const std::vector<Variant> &variants)
{
    // Reuse the typed catalog boundary before emitting any row for a contig.
    const ContigVariants checked(contig.bases, variants, contig.index);
    (void)checked;
    require_vcf_token(contig.name, "contig name");

    for (const Variant &event : variants) {
        require_vcf_token(event.id, "variant ID");
        std::uint32_t position = 0U;
        model::Bases reference;
        model::Bases alternate;
        switch (event.kind) {
        case model::VariantKind::snv:
            position = event.reference_start + 1U;
            reference = event.ref_bases;
            alternate = event.alt_bases;
            break;
        case model::VariantKind::insertion: {
            const std::uint32_t anchor = event.reference_start;
            if (anchor == 0U || anchor > contig.bases.size()) {
                throw MutationCatalogError(
                    "VCF export cannot anchor an insertion at this boundary");
            }
            const std::uint8_t anchor_base = contig.bases[anchor - 1U];
            (void)vcf_base(anchor_base);
            position = anchor;
            reference = {anchor_base};
            alternate = reference;
            alternate.insert(
                alternate.end(), event.alt_bases.begin(), event.alt_bases.end());
            break;
        }
        case model::VariantKind::deletion:
            if (event.reference_start > 0U) {
                const std::uint8_t anchor =
                    contig.bases[event.reference_start - 1U];
                position = event.reference_start;
                reference = {anchor};
                reference.insert(
                    reference.end(),
                    event.ref_bases.begin(),
                    event.ref_bases.end());
                alternate = {anchor};
            } else if (event.reference_end < contig.bases.size()) {
                const std::uint8_t anchor = contig.bases[event.reference_end];
                position = event.reference_start + 1U;
                reference = event.ref_bases;
                reference.push_back(anchor);
                alternate = {anchor};
            } else {
                throw MutationCatalogError(
                    "VCF export cannot find a resolved deletion anchor");
            }
            break;
        default:
            throw MutationCatalogError("VCF export received an invalid event kind");
        }

        sink << contig.name << '\t' << position << '\t' << event.id << '\t';
        write_vcf_bases(sink, reference);
        sink.put('\t');
        write_vcf_bases(sink, alternate);
        sink << "\t.\tPASS\t.\tGT\t" << vcf_genotype(event.alt_haplotypes)
             << '\n';
        if (!sink) {
            throw MutationCatalogError("failed while writing a VCF record");
        }
    }
}

} // namespace htsim::variant
