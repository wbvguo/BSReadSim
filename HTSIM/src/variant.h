#ifndef HTSIM_VARIANT_H
#define HTSIM_VARIANT_H

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>
#include <optional>

#include "types.h"
#include "reference.h"
#include "utilities.h"

// ---- variant_catalog --------------------------------------------------------

namespace htsim::variant {

inline constexpr std::uint8_t maximum_indel_bases = 4;

class VariantCatalogError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// One normalized, biallelic event in original-reference coordinates. REF and
// ALT have already had their common prefix and suffix removed. Insertions use
// a zero-width interval and empty REF; deletions use empty ALT.
struct Variant {
    std::uint32_t contig_index = 0;
    std::uint32_t reference_start = 0;
    std::uint32_t reference_end = 0;
    model::VariantKind kind = model::VariantKind::snv;
    model::Bases ref_bases;
    model::Bases alt_bases;
    model::HaplotypeMask alt_haplotypes =
        model::HaplotypeMask::both;
    std::string id;
    model::VariantSource source = model::VariantSource::vcf;
};

// Verified text VCF snapshot projected into reference order. The frozen v1
// subset requires one diploid sample, biallelic A/C/G/T records, and GT alleles
// 0 or 1. Unphased heterozygotes receive deterministic phase from the
// haplotype RNG domain before they enter this catalog.
class VariantFile {
public:
    VariantFile(
        const std::string &path,
        const crypto::Sha256Digest &expected_file_sha256,
        const std::vector<reference::ContigMetadata> &reference_catalog,
        std::uint64_t master_seed);

    const std::vector<Variant> &variants(std::uint32_t contig_index) const;
    std::uint64_t variant_count() const noexcept;
    const crypto::Sha256Digest &file_sha256() const noexcept;

private:
    crypto::Sha256Digest file_sha256_ = {};
    std::vector<std::vector<Variant>> variants_by_contig_;
    std::uint64_t variant_count_ = 0;
};

// Per-contig reference validation boundary. Construction verifies canonical
// order, non-overlap, contig identity, and exact REF bases against the opened
// reference snapshot. This component retains only typed variants.
class ContigVariants {
public:
    ContigVariants(
        const model::Bases &reference_bases,
        const std::vector<Variant> &variants,
        std::uint32_t contig_index);

    const std::vector<Variant> &variants() const noexcept;
    std::uint32_t contig_index() const noexcept;
    std::uint32_t reference_length() const noexcept;

private:
    std::uint32_t contig_index_ = 0;
    std::uint32_t reference_length_ = 0;
    std::vector<Variant> variants_;
};

} // namespace htsim::variant

// ---- haplotype_projector --------------------------------------------------------

namespace htsim::haplotype {

enum class ProjectionFailure {
    invalid_input,
    boundary_cut,
    empty_projection,
    capacity,
    invariant,
};

class ProjectionError : public std::runtime_error {
public:
    ProjectionError(ProjectionFailure failure, const char *message)
        : std::runtime_error(message), failure_(failure)
    {}

    ProjectionFailure failure() const noexcept {return failure_;}

private:
    ProjectionFailure failure_;
};

// A selected haplotype projected over one non-empty, 0-based half-open
// original-reference interval. Coordinates remain uint32_t at this contig-local
// boundary. Inserted bases use reference position -1 and all event ids are the
// stable, per-contig ordinals from ContigVariants.
struct ProjectedInterval {
    std::uint32_t contig_index = 0;
    std::uint8_t haplotype = 0;
    std::uint32_t reference_start = 0;
    std::uint32_t reference_end = 0;
    model::Bases template_bases;
    std::vector<std::int64_t> reference_positions;
    std::vector<std::uint32_t> base_variant_indices;
    std::vector<model::Variant> variants;
};

// Insertions are anchored at a reference boundary.  A haplotype-coordinate
// fragment may begin immediately after, or end immediately after, all bases of
// an insertion.  Those two cases cannot be expressed by a bare half-open
// reference interval, so the boundary policy records which anchored insertion
// belongs to the projected template.  It never permits cutting through the
// inserted bases themselves.
struct ProjectionBoundaryPolicy {
    bool include_start_anchor_insertion = true;
    bool include_end_anchor_insertion = false;
};

// Apply only variants whose HaplotypeMask contains zero_based_haplotype. An
// interval that cuts through an active deletion is rejected; callers must not
// silently truncate a deletion at a fragment boundary.
ProjectedInterval project_interval(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint8_t zero_based_haplotype,
    std::uint32_t reference_start,
    std::uint32_t reference_end);

ProjectedInterval project_interval(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    std::uint8_t zero_based_haplotype,
    std::uint32_t reference_start,
    std::uint32_t reference_end,
    ProjectionBoundaryPolicy boundary_policy);

} // namespace htsim::haplotype

// ---- haplotype_layout --------------------------------------------------------

namespace htsim::haplotype {

class HaplotypeLayoutError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// One physical boundary on a selected haplotype.  A deletion can make the
// reference envelope on the left end before the envelope on the right begins.
// An insertion has two distinct physical boundaries at the same reference
// anchor; the flags say which adjacent fragment owns the complete insertion.
struct FragmentBoundary {
    std::uint32_t left_reference_end = 0;
    std::uint32_t right_reference_start = 0;
    bool include_insertion_in_left_fragment = false;
    bool include_insertion_in_right_fragment = true;
};

// Compact coordinate and ambiguity index for one selected haplotype.  All
// contig-local coordinates and counts remain uint32_t.  Only exceptional
// variant boundaries and an N bit-vector are retained; materialized bases are
// optional and used by motif discovery, not by ordinary fragment projection.
class HaplotypeLayout {
public:
    HaplotypeLayout(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        std::uint8_t zero_based_haplotype,
        bool materialize_bases = false);

    std::uint32_t contig_index() const noexcept {return contig_index_;}
    std::uint8_t haplotype() const noexcept {return haplotype_;}
    std::uint32_t reference_length() const noexcept {return reference_length_;}
    std::uint32_t length() const noexcept {return length_;}

    bool has_materialized_bases() const noexcept {return materialized_;}
    const model::Bases &bases() const;

    // Return nullopt only for a boundary strictly inside inserted ALT bases.
    std::optional<FragmentBoundary> boundary(
        std::uint32_t haplotype_offset) const;

    // BED/reference boundary projected immediately before an insertion at the
    // same anchor.  A boundary strictly inside a deletion has no haplotype
    // counterpart and returns nullopt.
    std::optional<std::uint32_t> boundary_before_reference(
        std::uint32_t reference_offset) const;

    std::uint32_t ambiguous_count(
        std::uint32_t haplotype_begin,
        std::uint32_t haplotype_end) const;

    // Exact upper bound for the encoded variant-event records contained by
    // any representable physical window of this length.  Per-base protocol
    // arrays are deliberately excluded because their cost is already bounded
    // by the physical template length.
    std::uint64_t maximum_variant_payload_bytes(
        const variant::ContigVariants &variants,
        std::uint32_t physical_span) const;

    // Project one non-empty physical haplotype slice back to typed reference
    // details.  A slice boundary inside an insertion fails closed.
    ProjectedInterval project(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        std::uint32_t haplotype_begin,
        std::uint32_t haplotype_end) const;

private:
    struct ReferenceRun {
        std::uint32_t haplotype_begin = 0;
        std::uint32_t haplotype_end = 0;
        std::uint32_t reference_begin = 0;
    };

    struct BoundaryException {
        std::uint32_t haplotype_offset = 0;
        FragmentBoundary boundary;
    };

    struct ActiveEventCoordinate {
        std::uint32_t reference_start = 0;
        std::uint32_t reference_end = 0;
        model::VariantKind kind = model::VariantKind::snv;
        std::int64_t cumulative_delta_after = 0;
        std::uint32_t haplotype_begin = 0;
        std::uint32_t haplotype_end = 0;
        std::uint32_t event_ordinal = 0;
    };

    void append_base(std::uint8_t base);
    void append_reference_run(
        const model::Bases &reference,
        std::uint32_t begin,
        std::uint32_t end);
    void append_alternate_run(
        const variant::Variant &event,
        std::uint32_t reference_begin);
    FragmentBoundary &upsert_exception(
        std::uint32_t haplotype_offset,
        std::uint32_t reference_offset);
    std::uint32_t n_rank(std::uint32_t haplotype_end) const;

    std::uint32_t contig_index_ = 0;
    std::uint8_t haplotype_ = 0;
    std::uint32_t reference_length_ = 0;
    std::uint32_t length_ = 0;
    bool materialized_ = false;
    model::Bases bases_;
    std::vector<std::uint64_t> ambiguous_words_;
    std::vector<std::uint32_t> ambiguous_prefix_;
    std::vector<ReferenceRun> reference_runs_;
    std::vector<BoundaryException> boundary_exceptions_;
    std::vector<ActiveEventCoordinate> active_events_;
};

} // namespace htsim::haplotype

// ---- mutation_catalog --------------------------------------------------------

namespace htsim::variant {

class MutationCatalogError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct MutationParameters {
    double mutation_rate = 0.0;
    double indel_fraction = 0.0;
    double indel_extension_probability = 0.0;
    bool homozygous_only = false;
};

// Generate one canonical, non-overlapping typed event stream for a materialized
// contig. Reference coordinates and per-contig event ordinals remain uint32;
// RNG entities widen the reference coordinate to uint64 without packing it
// together with another field.
std::vector<Variant> generate_de_novo_events(
    const reference::Contig &contig,
    std::uint64_t master_seed,
    const MutationParameters &parameters);

} // namespace htsim::variant

#endif // HTSIM_VARIANT_H
