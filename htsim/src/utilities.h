#ifndef HTSIM_UTILITIES_H
#define HTSIM_UTILITIES_H

#include <cstdint>
#include <stdexcept>
#include <string_view>
#include <array>
#include <cstddef>
#include <vector>
#include <functional>
#include <memory>
#include <string>

struct htsFile;

// ---- normal_sampler --------------------------------------------------------

namespace htsim::normal_sampler {

inline constexpr std::string_view algorithm_id = "box-muller-normal";

class SamplingError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Stateless standard-normal draw from one explicit Philox block.
//
// pair 0 -> U1=(top53(u64)+1)/2^53 in (0,1]
// pair 1 -> U2=top53(u64)/2^53 in [0,1)
// result -> sqrt(-2*log(U1))*cos(2*pi*U2)
//
// U1==1 returns +0 without evaluating trigonometric functions. The algorithm
// requires IEEE-754 binary64 and round-to-nearest. Counter mapping is bit-exact
// across implementations; numerical bit identity is guaranteed only for a
// frozen compiler/libm, with tolerance required across other platforms.
double standard_normal(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index);

} // namespace htsim::normal_sampler

// ---- rng --------------------------------------------------------

namespace htsim::rng {

inline constexpr std::string_view contract_id =
    "philox4x32-10+philox-domain";

enum class Stage : std::uint32_t {
    mutation = 0,
    methylation_level = 1,
    fragment = 2,
    haplotype = 3,
    site_state = 4,
    library_orientation = 5,
    conversion = 6,
    quality = 7,
    sequencing_error = 8,
    count = 9,
};

inline constexpr std::array<std::string_view, 9> stage_names = {{
    "mutation",
    "methylation-level",
    "fragment",
    "haplotype",
    "site-state",
    "library-orientation",
    "conversion",
    "quality",
    "sequencing-error",
}};

// Numeric ASCII "BSR/KEY2". This occupies the entity half of the one Philox
// block reserved exclusively for domain-key derivation.
inline constexpr std::uint64_t domain_key_entity =
    UINT64_C(0x4253522f4b455932);

class ContractError : public std::invalid_argument {
public:
    using std::invalid_argument::invalid_argument;
};

using Block = std::array<std::uint32_t, 4>;

// A complete, immutable address for one RNG block. No draw advances any state.
struct Address {
    std::uint64_t master_seed;
    Stage stage;
    std::uint32_t contig_index;
    std::uint64_t entity_ordinal;
    std::uint64_t local_index;
};

// Parse the canonical decimal-string representation used by the JSON contract.
// Signs, whitespace, leading zeroes, an empty string, and values outside
// [0, 2**64) are rejected (the single string "0" remains valid).
std::uint64_t parse_u64_decimal(std::string_view text);

bool is_valid_stage(Stage stage) noexcept;

// Frozen cross-language BSReadSim RNG identity; this is a wire/product
// contract, not the name of an implementation namespace.
//
// domain_local = (uint64(stage) << 32) | contig_index
// block = Philox4x32-10(master_seed, domain_key_entity, domain_local)
// key = uint64(block[0]) | (uint64(block[1]) << 32)
std::uint64_t derive_key(
    std::uint64_t master_seed,
    Stage stage,
    std::uint32_t contig_index);

// Low-level primitives operating on an already-derived 64-bit Philox key.
Block philox4x32_10(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index) noexcept;

std::uint32_t u32(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::size_t lane = 0);

std::uint64_t u64(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::size_t pair = 0);

double uniform01(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::size_t pair = 0);

bool bernoulli(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    double probability,
    std::size_t pair = 0);

// Fixed-cost multiply-high range reduction. The unsigned overload returns a
// value in [0, upper_exclusive); the signed overload returns [lower, upper).
std::uint64_t bounded_integer(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::uint64_t upper_exclusive);

std::int64_t bounded_integer(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::int64_t lower,
    std::int64_t upper_exclusive);

// Full-address adapters. These are pure convenience functions and derive the
// same stage/contig key on every invocation.
Block block(const Address &address);
std::uint32_t u32(const Address &address, std::size_t lane = 0);
std::uint64_t u64(const Address &address, std::size_t pair = 0);
double uniform01(const Address &address, std::size_t pair = 0);
bool bernoulli(const Address &address, double probability, std::size_t pair = 0);
std::uint64_t bounded_integer(
    const Address &address,
    std::uint64_t upper_exclusive);
std::int64_t bounded_integer(
    const Address &address,
    std::int64_t lower,
    std::int64_t upper_exclusive);

} // namespace htsim::rng

// ---- sha256 --------------------------------------------------------

namespace htsim::crypto {

using Sha256Digest = std::array<std::uint8_t, 32>;

inline constexpr std::uint64_t maximum_sha256_input_bytes =
    UINT64_MAX / UINT64_C(8);

class Sha256Error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Incremental SHA-256 with an explicit input-length boundary. digest() is a
// non-mutating snapshot, so callers may inspect an intermediate digest and
// continue updating the same instance.
class Sha256 {
public:
    Sha256() noexcept;

    void update(const std::uint8_t *data, std::size_t size);
    void update(const std::vector<std::uint8_t> &data);

    Sha256Digest digest() const noexcept;
    std::uint64_t size() const noexcept {return total_bytes_;}

private:
    std::array<std::uint32_t, 8> hash_ = {};
    std::array<std::uint8_t, 64> buffer_ = {};
    std::uint64_t total_bytes_ = 0;
    std::size_t buffered_bytes_ = 0;
};

Sha256Digest sha256(const std::uint8_t *data, std::size_t size);
Sha256Digest sha256(const std::vector<std::uint8_t> &data);

} // namespace htsim::crypto

// ---- text_snapshot --------------------------------------------------------

namespace htsim::text {

inline constexpr std::size_t maximum_line_bytes = 1024U * 1024U;

class TextSnapshotError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// The line view is valid only for the duration of the callback. Line numbers
// are one-based physical text lines after transparent gzip decompression.
using LineVisitor =
    std::function<void(std::string_view line, std::uint64_t line_number)>;
// The chunk view contains decoded bytes and is valid only during the callback.
// Chunk boundaries are deliberately unspecified.
using ChunkVisitor = std::function<void(std::string_view chunk)>;
using HtsFileVisitor = std::function<void(::htsFile *file)>;

// Stable fail-closed access to a plain or gzip-compressed regular text file.
// Construction opens the path once and computes its raw-byte SHA-256 before
// decompression. Every visit uses that descriptor, validates raw bytes while
// decoding, then performs an independent final raw pass. LF and CRLF are
// accepted; bare CR, corrupt/truncated/trailed gzip, and lines over 1 MiB are
// rejected. A failed visit permanently poisons the snapshot.
class TextSnapshot {
public:
    TextSnapshot(
        const std::string &path,
        std::size_t maximum_decoded_line_bytes = maximum_line_bytes);
    ~TextSnapshot();

    TextSnapshot(const TextSnapshot &) = delete;
    TextSnapshot &operator=(const TextSnapshot &) = delete;
    TextSnapshot(TextSnapshot &&) = delete;
    TextSnapshot &operator=(TextSnapshot &&) = delete;

    const crypto::Sha256Digest &file_sha256() const noexcept;
    void visit_lines(const LineVisitor &visitor);
    void visit_chunks(const ChunkVisitor &visitor);
    // Advanced format-specific access through the same verified descriptor.
    // The handle is owned by the snapshot and valid only during the callback.
    void visit_hts(const HtsFileVisitor &visitor);
    void verify_unchanged();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace htsim::text

#endif // HTSIM_UTILITIES_H
