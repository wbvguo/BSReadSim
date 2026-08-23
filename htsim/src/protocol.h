#ifndef HTSIM_PROTOCOL_H
#define HTSIM_PROTOCOL_H

#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include "utilities.h"
#include "types.h"

// ---- wire --------------------------------------------------------

namespace htsim::protocol {

inline constexpr std::uint16_t protocol_major = 2;
inline constexpr std::uint16_t protocol_minor = 1;
inline constexpr std::uint32_t preamble_flags = 0;
inline constexpr std::size_t maximum_string_bytes = 1024U * 1024U;
inline constexpr std::size_t maximum_frame_payload = 64U * 1024U * 1024U;
inline constexpr std::uint32_t no_reference_position = UINT32_C(0xffffffff);
inline constexpr std::uint8_t details_present = UINT8_C(0x01);
inline constexpr std::string_view config_schema_version = "1.1";
inline constexpr std::string_view rng_contract = rng::contract_id;

using Digest = crypto::Sha256Digest;
using Bases = std::vector<std::uint8_t>;

// CRC32C protects each frame independently; SHA-256 below protects ordered
// stream identity. Keeping CRC here makes the wire integrity primitive
// explicit without exposing a general checksum utility module.
std::uint32_t crc32c(const std::uint8_t *data, std::size_t size) noexcept;
std::uint32_t crc32c(const std::vector<std::uint8_t> &data) noexcept;

class ProtocolError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class FrameType : std::uint8_t {
    header = 1,
    fragment_batch = 2,
    trailer = 3,
    error = 255,
};

enum class Technology : std::uint8_t {wgbs = 1, rrbs = 2, tbs = 3};
enum class BaseEncoding : std::uint8_t {acgtn_u8 = 1};
enum class AmbiguityPolicy : std::uint8_t {preserve_n = 0, resolve_once = 1};
enum class CaptureStrand : std::uint8_t {unknown = 0, forward = 1, reverse = 2};
enum class VariantKind : std::uint8_t {snv = 1, insertion = 2, deletion = 3};
enum class VariantSource : std::uint8_t {vcf = 1, de_novo = 2};
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
    asm_profile = 2,
    beta = 3,
    pooled_cgmap = 4,
};
enum class MethylationAllele : std::uint8_t {
    shared = 0,
    reference_haplotype = 1,
    alternate_haplotype = 2,
};

struct Contig {
    std::string name;
    std::uint32_t length = 0;
    Digest reference_sha256 = {};
};

struct Header {
    std::string run_id;
    std::string core_version;
    std::string config_schema_version;
    std::string rng_contract;
    std::uint64_t master_seed = 0;
    Digest normalized_config_sha256 = {};
    Technology technology = Technology::wgbs;
    bool has_details = false;
    std::uint8_t mates_per_fragment = 1;
    BaseEncoding base_encoding = BaseEncoding::acgtn_u8;
    AmbiguityPolicy ambiguity_policy = AmbiguityPolicy::preserve_n;
    std::uint32_t read_length_r1 = 0;
    std::uint32_t read_length_r2 = 0;
    std::vector<Contig> contigs;
};

struct FragmentDetails {
    std::vector<std::uint32_t> projection_offsets;
    std::vector<std::uint32_t> variant_offsets;
    std::vector<std::uint32_t> original_n_offsets;
    std::vector<std::uint32_t> projection_template_starts;
    std::vector<std::uint32_t> projection_template_ends;
    std::vector<std::uint32_t> projection_reference_starts;
    std::vector<std::uint32_t> variant_indices;
    std::vector<std::uint32_t> variant_id_offsets;
    std::vector<std::uint32_t> variant_reference_starts;
    std::vector<std::uint32_t> variant_reference_ends;
    std::vector<std::uint32_t> variant_template_starts;
    std::vector<std::uint32_t> variant_template_ends;
    std::vector<std::uint32_t> variant_ref_offsets;
    std::vector<std::uint32_t> variant_alt_offsets;
    std::vector<std::uint32_t> site_reference_positions;
    std::vector<std::uint32_t> original_n_template_offsets;
    std::vector<std::uint8_t> variant_sources;
    std::vector<std::uint8_t> variant_kinds;
    std::vector<std::uint8_t> variant_phased_haplotypes;
    std::vector<std::uint8_t> variant_ids;
    Bases variant_ref_bases;
    Bases variant_alt_bases;
};

struct FragmentBatch {
    std::uint32_t first_fragment_ordinal = 0;
    std::vector<std::uint32_t> contig_indices;
    std::vector<std::uint32_t> reference_starts;
    std::vector<std::uint32_t> reference_ends;
    std::vector<std::uint32_t> template_offsets;
    std::vector<std::uint32_t> mate_offsets;
    std::vector<std::uint32_t> site_offsets;
    std::vector<std::uint32_t> mate_template_starts;
    std::vector<std::uint32_t> mate_template_ends;
    std::vector<std::uint32_t> site_template_offsets;
    std::vector<float> site_probabilities;
    std::vector<std::uint8_t> haplotypes;
    std::vector<std::uint8_t> capture_strands;
    std::vector<std::uint8_t> mate_indices;
    std::vector<std::uint8_t> mate_reverse_complements;
    std::vector<std::uint8_t> site_contexts;
    std::vector<std::uint8_t> methylation_sources;
    std::vector<std::uint8_t> site_alleles;
    Bases template_bases;
    std::optional<FragmentDetails> details;

    std::uint32_t fragment_count() const;
    std::uint32_t template_base_count() const;
    std::uint32_t mate_count() const;
    std::uint32_t methylation_site_count() const;
};

struct Trailer {
    std::uint64_t fragment_count = 0;
    std::uint64_t fragment_batch_count = 0;
    std::uint64_t mate_count = 0;
    std::uint64_t template_base_count = 0;
    std::uint64_t methylation_site_count = 0;
    std::uint64_t skipped_fragment_count = 0;
    std::vector<std::uint64_t> per_contig_fragment_counts;
    Digest stream_sha256 = {};
};

struct ErrorFrame {
    std::uint32_t error_code = 0;
    std::string message;
};

class PreparedFragmentBatch {
public:
    PreparedFragmentBatch(PreparedFragmentBatch &&) noexcept = default;
    PreparedFragmentBatch &operator=(PreparedFragmentBatch &&) noexcept = default;

    PreparedFragmentBatch(const PreparedFragmentBatch &) = delete;
    PreparedFragmentBatch &operator=(const PreparedFragmentBatch &) = delete;

    std::uint32_t first_ordinal() const noexcept {return first_ordinal_;}
    std::uint32_t fragment_count() const noexcept {return fragment_count_;}
    std::size_t encoded_payload_size() const noexcept {return payload_.size();}

private:
    PreparedFragmentBatch() = default;

    friend PreparedFragmentBatch prepare_fragment_batch(
        const Header &, FragmentBatch);
    friend class Writer;

    std::uint32_t first_ordinal_ = 0;
    std::uint32_t fragment_count_ = 0;
    std::uint32_t mate_count_ = 0;
    std::uint32_t template_base_count_ = 0;
    std::uint32_t methylation_site_count_ = 0;
    std::uint8_t frame_flags_ = 0;
    Digest header_payload_sha256_ = {};
    std::vector<std::uint64_t> per_contig_fragment_counts_;
    std::vector<std::uint8_t> payload_;
};

PreparedFragmentBatch prepare_fragment_batch(
    const Header &header,
    FragmentBatch batch);

class Writer {
public:
    explicit Writer(std::ostream &sink);

    Writer(const Writer &) = delete;
    Writer &operator=(const Writer &) = delete;
    Writer(Writer &&) = delete;
    Writer &operator=(Writer &&) = delete;

    void write_header(const Header &header);
    void write_batch(FragmentBatch batch);
    void write_prepared_batch(PreparedFragmentBatch &&batch);
    Trailer finish(std::uint64_t skipped_fragment_count = 0);
    void write_error(const ErrorFrame &error);

    bool complete() const noexcept {return complete_;}
    bool failed() const noexcept {return failed_;}

private:
    void ensure_open() const;
    void write_raw(const std::uint8_t *bytes, std::size_t size);
    void write_raw(const std::vector<std::uint8_t> &bytes);
    void poison() noexcept {failed_ = true;}

    std::ostream *sink_;
    std::optional<Header> header_;
    Digest header_payload_sha256_ = {};
    crypto::Sha256 digest_state_;
    std::uint64_t next_sequence_ = 0;
    std::uint64_t next_ordinal_ = 0;
    std::uint64_t fragment_count_ = 0;
    std::uint64_t fragment_batch_count_ = 0;
    std::uint64_t mate_count_ = 0;
    std::uint64_t template_base_count_ = 0;
    std::uint64_t methylation_site_count_ = 0;
    std::vector<std::uint64_t> per_contig_fragment_counts_;
    bool complete_ = false;
    bool failed_ = false;
};

} // namespace htsim::protocol

// ---- adapter --------------------------------------------------------

namespace htsim::protocol {

// Project the generator's typed fragment materialization directly onto the
// columnar wire schema. model::Fragment is the shared in-memory scientific
// boundary; this file owns the only projection from that model to transport.
FragmentBatch make_fragment_batch(
    const Header &header,
    const std::vector<model::Fragment> &fragments);

} // namespace htsim::protocol

// ---- emitter --------------------------------------------------------

namespace htsim::model {
struct Fragment;
}

namespace htsim::protocol {

class Writer;
struct Header;

// Buffer typed fragments, encode batches concurrently, and commit them in
// fragment-ordinal order. Wire validation remains the Writer's responsibility.
class BatchEmitter {
public:
    BatchEmitter(
        Writer &writer,
        const Header &header,
        std::uint32_t worker_count,
        std::uint32_t fragments_per_batch);
    ~BatchEmitter();

    BatchEmitter(const BatchEmitter &) = delete;
    BatchEmitter &operator=(const BatchEmitter &) = delete;

    void write(model::Fragment fragment);
    void finish();

private:
    class Implementation;
    std::unique_ptr<Implementation> implementation_;
};

} // namespace htsim::protocol

#endif // HTSIM_PROTOCOL_H
