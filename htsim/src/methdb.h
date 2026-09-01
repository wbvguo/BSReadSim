#ifndef HTSIM_METHDB_H
#define HTSIM_METHDB_H

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <memory>
#include <string>
#include <vector>
#include <string_view>
#include <type_traits>
#include <utility>
#include <array>

#include "types.h"
#include "reference.h"
#include "utilities.h"
#include "variant.h"

// ---- site --------------------------------------------------------

namespace htsim::methdb {

// MethDB probabilities are unsigned-normalized 16-bit values. This is the
// sole catalog/runtime authority; floating point exists only while parsing
// text/model input and while emitting protocol/text output.
using ProbabilityU16 = std::uint16_t;

ProbabilityU16 probability_to_u16(float probability);
float probability_from_u16(ProbabilityU16 probability) noexcept;

class ContextError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct ContextNeighborhood {
    std::optional<std::uint8_t> upstream_second;
    std::optional<std::uint8_t> upstream_first;
    std::uint8_t center = 0;
    std::optional<std::uint8_t> downstream_first;
    std::optional<std::uint8_t> downstream_second;
};

// Classify from an explicit haplotype-local five-base neighborhood. Missing
// values represent true contig boundaries, not an arbitrary fragment edge.
std::optional<model::MethylationContext> classify_context(
    const ContextNeighborhood &neighborhood,
    bool collect_non_cpg);

// Classify a C/G at reference_position using the entire normalized contig,
// never a fragment-local slice.  A missing flank or N in a required flank is
// unclassifiable and returns nullopt.  When collect_non_cpg is false, only the
// two CpG contexts are returned.
std::optional<model::MethylationContext> classify_context(
    const model::Bases &contig_bases,
    std::uint64_t reference_position,
    bool collect_non_cpg);

class EntityError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

enum class SiteEntityKind : std::uint8_t {
    reference_baseline = 0,
    variant_reference_shared = 1,
    variant_reference_haplotype_0 = 2,
    variant_reference_haplotype_1 = 3,
    insertion_shared = 4,
    insertion_haplotype_0 = 5,
    insertion_haplotype_1 = 6,
};

// A validated 64-bit RNG entity. Bits 63..60 are SiteEntityKind and bits
// 59..0 are a kind-specific payload. Construction is intentionally restricted
// to the checked factories below.
class SiteEntity {
public:
    std::uint64_t value() const noexcept {return value_;}

    friend constexpr bool operator==(SiteEntity left, SiteEntity right) noexcept
    {
        return left.value_ == right.value_;
    }
    friend constexpr bool operator!=(SiteEntity left, SiteEntity right) noexcept
    {
        return !(left == right);
    }

private:
    explicit constexpr SiteEntity(std::uint64_t value) noexcept : value_(value) {}

    std::uint64_t value_;

    friend SiteEntity reference_site_entity(std::uint32_t);
    friend SiteEntity variant_reference_site_entity(
        std::uint32_t, model::HaplotypeMask, std::uint8_t);
    friend SiteEntity insertion_site_entity(
        std::uint32_t, std::uint8_t, model::HaplotypeMask, std::uint8_t);
    friend SiteEntity decode_site_entity(std::uint64_t);
};

struct DecodedSiteEntity {
    SiteEntityKind kind = SiteEntityKind::reference_baseline;
    std::uint32_t reference_position = 0;
    std::uint32_t event_ordinal = 0;
    std::uint8_t insertion_offset = 0;
};

SiteEntity reference_site_entity(std::uint32_t reference_position);
SiteEntity variant_reference_site_entity(
    std::uint32_t reference_position,
    model::HaplotypeMask alt_haplotypes,
    std::uint8_t zero_based_haplotype);
SiteEntity insertion_site_entity(
    std::uint32_t event_ordinal,
    std::uint8_t insertion_offset,
    model::HaplotypeMask alt_haplotypes,
    std::uint8_t zero_based_haplotype);

DecodedSiteEntity decode_site_entity(SiteEntity entity) noexcept;
SiteEntity decode_site_entity(std::uint64_t encoded);
bool entity_is_insertion(SiteEntity entity) noexcept;

} // namespace htsim::methdb

// ---- cgmap_profile --------------------------------------------------------

namespace htsim::methdb {

class CgmapProfileError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// All supported text formats normalize to the same position-specific MethDB
// overlay. The
// explicit selector is part of the launch contract; paths and suffixes never
// influence parsing.
enum class MethylationProfileFormat : std::uint8_t {
    cgmap = 0,
    bed_methyl = 1,
    methbg = 2,
    methbed = 3,
};

struct CgmapRecord {
    std::uint32_t reference_position = 0;
    ProbabilityU16 probability_u16 = 0U;
    model::MethylationContext context = model::MethylationContext::cg_c;
    bool has_probability = false;
    // Protocol base encoding for the second base in the cytosine-oriented
    // DINUC field: 0=A, 1=C, 2=G, 3=T.
    std::uint8_t dinucleotide_second = 0U;
};

static_assert(sizeof(CgmapRecord) == 12U,
              "CGmap in-memory record must remain 12 bytes");

// Validate a sorted per-contig record set against the complete reference
// sequence. This is shared by the profile snapshot and both MethDB builders.
void validate_cgmap_records(
    const model::Bases &contig_bases,
    const std::vector<CgmapRecord> &records);

// Verified text methylation input with a bounded-RAM per-contig access
// boundary. Parsed rows are normalized into an unlinked fixed-record spool;
// records() only materializes and reference-validates the requested contig.
// Instances are single-owner and not thread-safe because the internal spool
// has one seek cursor.
class CgmapProfile {
public:
    CgmapProfile(
        const std::string &path,
        const std::vector<reference::ContigMetadata> &reference_catalog,
        MethylationProfileFormat format = MethylationProfileFormat::cgmap);
    ~CgmapProfile();

    CgmapProfile(const CgmapProfile &) = delete;
    CgmapProfile &operator=(const CgmapProfile &) = delete;
    CgmapProfile(CgmapProfile &&) = delete;
    CgmapProfile &operator=(CgmapProfile &&) = delete;

    const crypto::Sha256Digest &file_sha256() const noexcept;
    std::uint64_t row_count() const noexcept;
    std::uint64_t defined_probability_count() const noexcept;

    std::vector<CgmapRecord> records(const reference::Contig &contig) const;
    void validate_contig(const reference::Contig &contig) const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace htsim::methdb

// ---- asm_profile --------------------------------------------------------

namespace htsim::methdb {

class AsmProfileError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class AsmProfileFormat : std::uint8_t {
    cgmaptools_ass = 0,
    bed = 1,
};

// One normalized allele-specific methylation row. Both positions are zero-based
// contig-local coordinates. The linked variant is an SNV whose HaplotypeMask
// later determines which haplotype receives each probability.
struct AsmRecord {
    std::uint32_t target_reference_position = 0;
    std::uint32_t linked_variant_position = 0;
    ProbabilityU16 reference_probability_u16 = 0U;
    ProbabilityU16 alternate_probability_u16 = 0U;
    model::MethylationContext context = model::MethylationContext::cg_c;
    std::uint8_t dinucleotide_second = 0;
    std::uint8_t linked_reference_base = 0;
    std::uint8_t linked_alternate_base = 0;
};

static_assert(sizeof(AsmRecord) == 16U,
              "ASM profile record must remain a compact 16-byte record");

// Validate already-normalized rows against one materialized reference contig.
// VCF linkage and diploid-site availability are validated by the overlay
// boundary, which owns the typed variant catalog.
void validate_asm_records(
    const model::Bases &contig_bases,
    const std::vector<AsmRecord> &records);

// Build the minimal heterozygous SNV authority implied by an ASM profile when
// no VCF is supplied. Repeated links to one identical SNV collapse to one
// event; conflicting allele definitions fail. Phase is deterministic under
// phasing_seed and stable reference order.
std::vector<variant::Variant> variants_from_asm(
    const reference::Contig &contig,
    const std::vector<AsmRecord> &records,
    std::uint64_t phasing_seed);

// Verified plain/gzip CGmapTools ASS or BSReadSim ASM BED snapshot. Parsing writes
// fixed-width records to an unlinked temporary spool and materializes only one
// reference-validated contig at a time. Native ASS rows are normalized from
// SNP order into target order at that boundary. Calls are intentionally not
// concurrent because the spool has one seek cursor.
class AsmProfile {
public:
    AsmProfile(
        const std::string &path,
        const std::vector<reference::ContigMetadata> &reference_catalog,
        AsmProfileFormat format = AsmProfileFormat::cgmaptools_ass);
    ~AsmProfile();

    AsmProfile(const AsmProfile &) = delete;
    AsmProfile &operator=(const AsmProfile &) = delete;
    AsmProfile(AsmProfile &&) = delete;
    AsmProfile &operator=(AsmProfile &&) = delete;

    const crypto::Sha256Digest &file_sha256() const noexcept;
    std::uint64_t row_count() const noexcept;

    std::vector<AsmRecord> records(const reference::Contig &contig) const;
    void validate_contig(const reference::Contig &contig) const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace htsim::methdb

// ---- beta_sampler --------------------------------------------------------

namespace htsim::beta_sampler {

// This identifier freezes both the distribution algorithm and the counter
// layout documented below.  Any change to either requires a new identifier.
inline constexpr std::string_view algorithm_id =
    "marsaglia-tsang-box-muller-beta";
inline constexpr std::string_view site_entity_algorithm_id =
    "marsaglia-tsang-box-muller-beta-site-entity";

inline constexpr std::uint32_t default_max_gamma_attempts = 1024;

class SamplingError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The default is also the absolute hard maximum in this algorithm. A smaller
// value is useful for testing fail-closed exhaustion; changing the value cannot
// change a successful draw, because every attempt has an explicit address.
struct Options {
    std::uint32_t max_gamma_attempts = default_max_gamma_attempts;
};

// Draw one Beta(alpha, beta) methylation level.
//
// RNG address contract
// --------------------
// key              = rng::derive_key(
//                        master_seed,
//                        rng::Stage::methylation_level,
//                        contig_index)
// entity_ordinal   = reference_position (the zero-based reference/site scan
//                    position for this baseline)
// local_index bits = [63:62 role][61:0 role-local offset]
//
// Roles are 0=alpha Gamma candidates, 1=alpha shape<1 boost uniform,
// 2=beta Gamma candidates, and 3=beta shape<1 boost uniform.  For candidate
// attempt i, offset 2*i is a Philox block whose pair 0 supplies the Box-Muller
// radius uniform U=(top53(u64)+1)/2^53 in (0,1] and pair 1 supplies its
// angle uniform in [0,1) using rng::uniform01.
// Offset 2*i+1 is a block whose pair 0 supplies the Marsaglia-Tsang acceptance
// uniform with the same explicit (0,1] conversion; its pair 1 is reserved.  A
// boost role uses offset zero, pair zero and that same conversion; all other
// words in that role are reserved.  Attempts are
// zero-based.  This fixed layout makes results independent of call order,
// chunking, and threads, including when either Gamma draw rejects candidates.
//
// Algorithm contract
// ------------------
// Gamma(shape>=1) uses the Marsaglia-Tsang d=shape-1/3,
// c=(1/3)/sqrt(d) evaluation
// rejection algorithm, a non-cached Box-Muller normal
// sqrt(-2*log(U1))*cos(2*pi*U2), and pi rounded to IEEE binary64.  The quick
// acceptance constant is binary64 0.0331.  Gamma(shape<1) uses
// Gamma(shape+1) * U^(1/shape).  The two independent Gamma values are
// combined with an overflow-resistant ratio and converted to float.  A zero
// attempt cap, exhausted cap, non-finite intermediate, or indeterminate 0/0
// ratio fails closed with SamplingError.  Shapes must be finite and >0.
// IEEE-754 round-to-nearest must be active; other rounding modes fail closed.
// contig_index is the frozen zero-based reference-catalog ordinal and is part
// of the RNG domain identity.
//
// This reference-coordinate address means a reference-derived site with no
// variant or allele-specific identity. A future variant/allele-aware sampler MUST define
// a new address contract and algorithm_id; it MUST NOT reuse this address while
// silently changing the entity semantics.
//
// The counter mapping is bit-exact across implementations.  The numerical path
// calls the platform log, sqrt, cos, and pow functions, so bitwise equality is
// only guaranteed for a frozen compiler/libm environment.  Cross-platform
// implementations must use the released numerical tolerance unless a portable
// correctly-rounded math contract is added in a future algorithm change.
float sample_beta(
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    std::uint64_t reference_position,
    double alpha,
    double beta,
    Options options = {});

// Variant-aware entry point. The numerical algorithm and local-index roles are
// identical to the reference-coordinate entry point, but entity_ordinal is
// the validated 64-bit SiteEntity from
// methylation-site-entity rather than an untyped reference position. Tag-0
// reference entities deliberately preserve every reference counter address.
float sample_beta_for_site(
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    methdb::SiteEntity entity,
    double alpha,
    double beta,
    Options options = {});

// Prevent accidental use of C++ booleans as shape parameters.  Explicitly
// converting a boolean to double remains an equally explicit caller decision.
template <
    typename Alpha,
    typename Beta,
    std::enable_if_t<
        std::is_same_v<std::decay_t<Alpha>, bool>
            || std::is_same_v<std::decay_t<Beta>, bool>,
        int> = 0>
float sample_beta(
    std::uint64_t,
    std::uint32_t,
    std::uint64_t,
    Alpha,
    Beta,
    Options = {}) = delete;

template <
    typename Alpha,
    typename Beta,
    std::enable_if_t<
        std::is_same_v<std::decay_t<Alpha>, bool>
            || std::is_same_v<std::decay_t<Beta>, bool>,
        int> = 0>
float sample_beta_for_site(
    std::uint64_t,
    std::uint32_t,
    methdb::SiteEntity,
    Alpha,
    Beta,
    Options = {}) = delete;

} // namespace htsim::beta_sampler

// ---- catalog --------------------------------------------------------

namespace htsim::methdb {

class CatalogError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct ShapePair {
    double alpha = 0.0;
    double beta = 0.0;
};

struct ContextShapes {
    ShapePair cg;
    ShapePair chg;
    ShapePair chh;
};

void validate_context_shapes(const ContextShapes &shapes);
const ShapePair &shape_for_context(
    model::MethylationContext context,
    const ContextShapes &shapes);

// One cache-oriented runtime word. Numeric order is key order because the key
// occupies the high 32 bits:
//   [63:32 key][31:16 probability][15:0 metadata]
using RuntimeSite = std::uint64_t;

RuntimeSite pack_runtime_site(
    std::uint32_t key,
    ProbabilityU16 probability,
    model::MethylationContext context,
    model::MethylationSource source,
    model::MethylationAllele allele,
    bool reference_equivalent = false);
std::uint32_t runtime_site_key(RuntimeSite site) noexcept;
ProbabilityU16 runtime_site_probability(RuntimeSite site) noexcept;
model::MethylationContext runtime_site_context(RuntimeSite site);
model::MethylationSource runtime_site_source(RuntimeSite site);
model::MethylationAllele runtime_site_allele(RuntimeSite site);
bool runtime_site_reference_equivalent(RuntimeSite site) noexcept;

struct CatalogSite {
    std::uint32_t reference_position = 0;
    ProbabilityU16 probability_u16 = 0U;
    model::MethylationContext context = model::MethylationContext::cg_c;
    model::MethylationSource methylation_source = model::MethylationSource::beta;
};

static_assert(sizeof(CatalogSite) == 8U,
              "methylation catalog site must remain an 8-byte record");

// Per-contig methylation-level catalog. Sites are stored once in ascending
// reference order, so overlapping fragments reuse the same deterministic
// genomic probability. The catalog owns no copy of the reference bases.
class MethylationCatalog {
public:
    MethylationCatalog(
        std::uint32_t reference_length,
        std::vector<CatalogSite> sites);
    using const_iterator = std::vector<CatalogSite>::const_iterator;

    MethylationCatalog(
        const model::Bases &contig_bases,
        std::uint32_t contig_index,
        std::uint64_t master_seed,
        bool collect_non_cpg,
        const ContextShapes &shapes,
        const std::vector<CgmapRecord> *cgmap_records = nullptr,
        bool pool_cgmap = false);

    const std::vector<CatalogSite> &sites() const noexcept;
    std::pair<const_iterator, const_iterator> sites_in_range(
        std::uint32_t begin,
        std::uint32_t end) const;

private:
    std::vector<CatalogSite> sites_;
};

} // namespace htsim::methdb

// ---- cgmap_pool --------------------------------------------------------

namespace htsim::methdb {

inline constexpr std::string_view cgmap_pool_algorithm_id =
    "cgmap-context-pool";
inline constexpr std::uint64_t cgmap_pool_local_index = UINT64_MAX;

class CgmapPoolError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// One-contig, three-class empirical methylation-level pool. The constructor
// accepts normalized profile rows and retains only defined q16 values.
// C/G-oriented protocol contexts share a class through an explicit enum
// mapping.
class CgmapPool {
public:
    explicit CgmapPool(const std::vector<CgmapRecord> &records);

    std::uint32_t size(model::MethylationContext context) const;

    // Returns no value when this contig has no defined input for the requested
    // context class; callers then retain the addressed Beta fallback.
    std::optional<ProbabilityU16> sample(
        model::MethylationContext context,
        std::uint64_t master_seed,
        std::uint32_t contig_index,
        SiteEntity entity) const;

private:
    std::array<std::vector<ProbabilityU16>, 3> values_;
};

} // namespace htsim::methdb

// ---- diploid_catalog --------------------------------------------------------

namespace htsim::methdb {

class DiploidCatalogError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Fail-closed ASM preflight before protocol output. This validates exact SNV
// linkage and streams both haplotypes only far enough to prove that each ASM
// target is one shared reference-equivalent context. It does not sample
// probabilities or materialize a methylation catalog.
void validate_asm_targets(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    bool collect_non_cpg,
    const std::vector<AsmRecord> &asm_records);

// The complete active-contig diploid representation. Reference-backed and
// inserted origins have separate key domains, and shared sites are stored only
// once. Every row is the authoritative 64-bit RuntimeSite; no second decoded
// float-bearing site representation is retained.
struct DiploidRuntimeArrays {
    std::vector<RuntimeSite> reference_shared;
    std::array<std::vector<RuntimeSite>, 2> reference_haplotypes;
    std::vector<RuntimeSite> insertion_shared;
    std::array<std::vector<RuntimeSite>, 2> insertion_haplotypes;
};

// One per-contig diploid MethDB. Context is discovered on complete haplotypes
// through a streaming five-base window. Sites identical on both haplotypes are
// stored once unless an ASM profile deliberately splits their probabilities;
// other differences are stored in haplotype overlays. No whole-haplotype
// position/event arrays are materialized.
class DiploidMethylationCatalog {
public:
    DiploidMethylationCatalog(
        std::uint32_t contig_index,
        std::uint32_t reference_length,
        DiploidRuntimeArrays runtime_arrays);
    DiploidMethylationCatalog(
        const reference::Contig &contig,
        const variant::ContigVariants &variants,
        std::uint64_t master_seed,
        bool collect_non_cpg,
        const ContextShapes &shapes,
        const std::vector<CgmapRecord> *cgmap_records = nullptr,
        const std::vector<AsmRecord> *asm_records = nullptr,
        bool pool_cgmap = false);

    const DiploidRuntimeArrays &runtime_arrays() const noexcept;
    void validate_asm_layer(
        const std::vector<variant::Variant> &variants,
        const std::vector<AsmRecord> &asm_records) const;
    void apply_asm_layer(
        const std::vector<variant::Variant> &variants,
        const std::vector<AsmRecord> &asm_records);
    DiploidRuntimeArrays take_runtime_arrays() && noexcept;

    // Map catalog sites into one already validated haplotype projection. The
    // returned sites are consecutive and sorted by template_offset.
    std::vector<model::MethylationSite> sites_for_projection(
        const haplotype::ProjectedInterval &projection) const;

private:
    std::uint32_t contig_index_ = 0;
    std::uint32_t reference_length_ = 0;
    DiploidRuntimeArrays runtime_arrays_;
};

} // namespace htsim::methdb

// ---- fixed snapshot --------------------------------------------------------

namespace htsim::methdb {

inline constexpr char methdb_magic[] = "methdb";
inline constexpr std::uint8_t methdb_version = 2U;
inline constexpr std::string_view methdb_bed_format = "methdb-bed-v2";
inline constexpr std::string_view legacy_methbed_snapshot_format = "methbed-v1";

class SnapshotError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct SnapshotContig {
    std::string name;
    std::uint32_t reference_length = 0;
    crypto::Sha256Digest reference_sha256 = {};
    bool diploid = false;
    std::vector<variant::Variant> variants;
    std::vector<CatalogSite> reference_sites;
    DiploidRuntimeArrays diploid_sites;
};

class SnapshotWriterImpl;

class SnapshotWriter {
public:
    SnapshotWriter(
        std::ostream &output,
        const crypto::Sha256Digest &binding,
        std::uint32_t contig_count);
    ~SnapshotWriter();

    void write_reference(
        const reference::ContigMetadata &metadata,
        const MethylationCatalog &catalog);
    void write_diploid(
        const reference::ContigMetadata &metadata,
        const MethylationCatalog &baseline,
        const DiploidMethylationCatalog &pre_asm_catalog,
        const std::vector<variant::Variant> &variants,
        const std::vector<AsmRecord> &asm_records = {});
    void finish();

private:
    std::ostream &output_;
    std::unique_ptr<SnapshotWriterImpl> impl_;
    std::uint32_t contig_count_ = 0;
    std::uint32_t written_ = 0;
    bool finished_ = false;
};

class Snapshot {
public:
    Snapshot(
        const std::string &path,
        const crypto::Sha256Digest &expected_binding,
        const std::vector<reference::ContigMetadata> &reference_catalog);

    // Decode, validate, and compile only the requested contig. The returned
    // value owns its data and no other contig is materialized.
    SnapshotContig contig(std::uint32_t contig_index) const;
    // Planning can load only the independently framed event sub-section,
    // without inflating baseline or methylation overlays.
    std::vector<variant::Variant> variants(
        std::uint32_t contig_index) const;
    bool contig_is_diploid(std::uint32_t contig_index) const;
    bool has_diploid_contigs() const noexcept;
    const crypto::Sha256Digest &file_sha256() const noexcept {
        return file_sha256_;
    }
    const crypto::Sha256Digest &content_sha256() const noexcept {
        return content_sha256_;
    }

private:
    struct Impl;
    std::shared_ptr<const Impl> impl_;
    crypto::Sha256Digest file_sha256_ = {};
    crypto::Sha256Digest content_sha256_ = {};
};

// Legacy internal reader for the rich MethDB extended-BED inspection export.
// It is not selected by the public --methbed option, which reads a
// methylation-only profile through MethylationProfileFormat::methbed.
class MethbedSnapshot {
public:
    MethbedSnapshot(
        const std::string &path,
        const crypto::Sha256Digest &expected_binding,
        const std::vector<reference::ContigMetadata> &reference_catalog);
    ~MethbedSnapshot();

    MethbedSnapshot(const MethbedSnapshot &) = delete;
    MethbedSnapshot &operator=(const MethbedSnapshot &) = delete;
    MethbedSnapshot(MethbedSnapshot &&) = delete;
    MethbedSnapshot &operator=(MethbedSnapshot &&) = delete;

    SnapshotContig contig(const reference::Contig &contig) const;
    std::vector<variant::Variant> variants(
        std::uint32_t contig_index) const;
    bool contig_is_diploid(std::uint32_t contig_index) const;
    bool has_diploid_contigs() const noexcept;
    const crypto::Sha256Digest &file_sha256() const noexcept;

private:
    class Impl;
    std::shared_ptr<Impl> impl_;
};

// Validate legacy MethDB extended-BED rows against the reference and embedded
// variants, then restore derived runtime state.
void normalize_methbed_contig(
    const reference::Contig &reference_contig,
    SnapshotContig &methbed_contig);

// Decode every stored row without requiring the original run configuration.
// Reference-backed origins receive half-open coordinates. Insertion origins
// retain their event ordinal and insertion offset because an inserted base has
// no one-base reference interval; #variant rows expose the event anchor.
void export_snapshot_bed(const std::string &path, std::ostream &sink);

} // namespace htsim::methdb

#endif // HTSIM_METHDB_H
