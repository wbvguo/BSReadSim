#ifndef HTSIM_TYPES_H
#define HTSIM_TYPES_H

#include <cstdint>
#include <string>
#include <vector>

namespace htsim::model {

using Bases = std::vector<std::uint8_t>;

inline constexpr std::uint8_t maximum_insertion_bases = 4;

enum class HaplotypeMask : std::uint8_t {
    haplotype_1 = 1,
    haplotype_2 = 2,
    both = 3,
};

inline constexpr bool is_haplotype_mask(std::uint8_t value) noexcept
{
    return value >= 1 && value <= 3;
}

inline constexpr bool mask_contains(
    HaplotypeMask mask,
    std::uint8_t zero_based_haplotype) noexcept
{
    return zero_based_haplotype < 2
        && (static_cast<std::uint8_t>(mask)
            & static_cast<std::uint8_t>(1U << zero_based_haplotype)) != 0;
}

inline constexpr std::uint32_t no_variant_index = UINT32_C(0xffffffff);

enum class CaptureStrand : std::uint8_t {
    unknown = 0,
    forward = 1,
    reverse = 2,
};

enum class VariantKind : std::uint8_t {
    snv = 1,
    insertion = 2,
    deletion = 3,
};

enum class VariantSource : std::uint8_t {
    vcf = 1,
    de_novo = 2,
    asm_profile = 3,
};

enum class MethylationContext : std::uint8_t {
    cg_c = 1,
    chg_c = 3,
    chh_c = 7,
    cg_g = 9,
    chg_g = 11,
    chh_g = 15,
};

enum class MethylationSource : std::uint8_t {
    cgmap = 1,
    asm_source = 2,
    beta = 3,
    pooled_cgmap = 4,
};

enum class MethylationAllele : std::uint8_t {
    shared = 0,
    reference_haplotype = 1,
    alternate_haplotype = 2,
};

struct Variant {
    std::uint32_t index = 0;
    std::string id{};
    VariantSource source = VariantSource::vcf;
    VariantKind kind = VariantKind::snv;
    std::uint8_t phased_haplotype = 255;
    std::uint64_t reference_start = 0;
    std::uint64_t reference_end = 0;
    Bases ref_bases;
    Bases alt_bases;
};

struct MethylationSite {
    std::uint32_t site_index = 0;
    std::uint32_t template_offset = 0;
    std::int64_t reference_pos = 0;
    MethylationContext context = MethylationContext::cg_c;
    MethylationSource methylation_source = MethylationSource::cgmap;
    MethylationAllele allele = MethylationAllele::shared;
    float methylation_probability = 0.0F;
};

struct SiteReference {
    std::uint32_t read_offset = 0;
    std::uint32_t site_index = 0;
};

struct Mate {
    std::uint8_t mate_index = 0;
    bool reverse_complement = false;
    std::uint32_t template_start = 0;
    std::uint32_t template_end = 0;
    std::uint64_t reference_start = 0;
    std::uint64_t reference_end = 0;
    std::vector<SiteReference> site_refs;
};

struct Fragment {
    std::uint64_t fragment_ordinal = 0;
    std::uint32_t contig_index = 0;
    std::uint8_t haplotype = 0;
    CaptureStrand capture_strand = CaptureStrand::unknown;
    std::uint64_t reference_start = 0;
    std::uint64_t reference_end = 0;
    Bases template_bases;
    std::vector<std::int64_t> reference_positions;
    std::vector<std::uint32_t> base_variant_indices;
    std::vector<Variant> variants;
    std::vector<MethylationSite> methylation_sites;
    std::vector<Mate> mates;
};

} // namespace htsim::model

#endif // HTSIM_TYPES_H
