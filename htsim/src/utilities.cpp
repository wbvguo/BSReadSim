#include "utilities.h"

#include <cfenv>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <algorithm>
#include <array>
#include <cstring>
#include <string>
#include <vector>
#include <atomic>
#include <cerrno>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <utility>

#include <htslib/hfile.h>
#include <htslib/bgzf.h>
#include <htslib/hts.h>

// ---- normal_sampler --------------------------------------------------------

namespace htsim::normal_sampler {
namespace {

static_assert(std::numeric_limits<double>::is_iec559
                  && std::numeric_limits<double>::radix == 2
                  && std::numeric_limits<double>::digits == 53,
              "normal sampler v1 requires IEEE-754 binary64");

constexpr double two_pi = 0x1.921fb54442d18p+2;
constexpr double inverse_uniform_scale = 0x1p-53;

double uniform_open_closed(
    std::uint64_t key,
    std::uint64_t entity,
    std::uint64_t local_index)
{
    const std::uint64_t significand =
        rng::u64(key, entity, local_index, 0) >> 11U;
    return static_cast<double>(significand + UINT64_C(1))
        * inverse_uniform_scale;
}

} // namespace

double standard_normal(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index)
{
    if (std::fegetround() != FE_TONEAREST) {
        throw SamplingError("normal sampler v1 requires round-to-nearest mode");
    }
    const double radius_uniform =
        uniform_open_closed(key, entity_ordinal, local_index);
    if (radius_uniform == 1.0) {return 0.0;}
    const double angle_uniform =
        rng::uniform01(key, entity_ordinal, local_index, 1);
    const double radius = std::sqrt(-2.0 * std::log(radius_uniform));
    const double result = radius * std::cos(two_pi * angle_uniform);
    if (!std::isfinite(result)) {
        throw SamplingError("normal sampler produced a non-finite value");
    }
    return result;
}

} // namespace htsim::normal_sampler

// ---- rng --------------------------------------------------------

namespace htsim::rng {
namespace {

constexpr std::uint32_t philox_m0 = UINT32_C(0xd2511f53);
constexpr std::uint32_t philox_m1 = UINT32_C(0xcd9e8d57);
constexpr std::uint32_t philox_w0 = UINT32_C(0x9e3779b9);
constexpr std::uint32_t philox_w1 = UINT32_C(0xbb67ae85);

struct Product64 {
    std::uint64_t low;
    std::uint64_t high;
};

Product64 multiply_u64(std::uint64_t left, std::uint64_t right) noexcept
{
    constexpr std::uint64_t mask32 = UINT64_C(0xffffffff);
    const std::uint64_t left_low = left & mask32;
    const std::uint64_t left_high = left >> 32;
    const std::uint64_t right_low = right & mask32;
    const std::uint64_t right_high = right >> 32;

    const std::uint64_t product00 = left_low * right_low;
    const std::uint64_t product01 = left_low * right_high;
    const std::uint64_t product10 = left_high * right_low;
    const std::uint64_t product11 = left_high * right_high;

    const std::uint64_t limb1 =
        (product00 >> 32) + (product01 & mask32) + (product10 & mask32);
    const std::uint64_t limb2 =
        (product01 >> 32) + (product10 >> 32) + (product11 & mask32)
        + (limb1 >> 32);
    const std::uint64_t limb3 = (product11 >> 32) + (limb2 >> 32);

    return {
        (product00 & mask32) | ((limb1 & mask32) << 32),
        (limb2 & mask32) | (limb3 << 32),
    };
}

std::uint64_t bounded_offset(const Block &random_block, std::uint64_t width)
{
    if (width == 0) {
        throw ContractError("bounded integer upper bound must be positive");
    }

    const std::uint64_t random_low =
        static_cast<std::uint64_t>(random_block[0])
        | (static_cast<std::uint64_t>(random_block[1]) << 32);
    const std::uint64_t random_high =
        static_cast<std::uint64_t>(random_block[2])
        | (static_cast<std::uint64_t>(random_block[3]) << 32);
    const Product64 low_product = multiply_u64(random_low, width);
    const Product64 high_product = multiply_u64(random_high, width);
    const std::uint64_t middle = low_product.high + high_product.low;
    const std::uint64_t carry = middle < low_product.high ? 1 : 0;
    return high_product.high + carry;
}

std::int64_t add_signed_offset(std::int64_t lower, std::uint64_t offset) noexcept
{
    if (lower >= 0) {
        return lower + static_cast<std::int64_t>(offset);
    }

    const std::uint64_t magnitude =
        static_cast<std::uint64_t>(-(lower + 1)) + 1;
    if (offset < magnitude) {
        const std::uint64_t remaining = magnitude - offset;
        if (remaining == (UINT64_C(1) << 63)) {
            return std::numeric_limits<std::int64_t>::min();
        }
        return -static_cast<std::int64_t>(remaining);
    }
    return static_cast<std::int64_t>(offset - magnitude);
}

} // namespace

std::uint64_t parse_u64_decimal(std::string_view text)
{
    if (text.empty()) {
        throw ContractError("unsigned 64-bit integer must not be empty");
    }
    if (text.size() > 1 && text.front() == '0') {
        throw ContractError("unsigned 64-bit integer must use canonical decimal form");
    }

    std::uint64_t value = 0;
    constexpr std::uint64_t maximum = std::numeric_limits<std::uint64_t>::max();
    for (char character : text) {
        if (character < '0' || character > '9') {
            throw ContractError("unsigned 64-bit integer must contain decimal digits only");
        }
        const auto digit = static_cast<std::uint64_t>(character - '0');
        if (value > (maximum - digit) / 10) {
            throw ContractError("unsigned 64-bit integer is out of range");
        }
        value = value * 10 + digit;
    }
    return value;
}

bool is_valid_stage(Stage stage) noexcept
{
    return static_cast<std::uint32_t>(stage)
        < static_cast<std::uint32_t>(Stage::count);
}

std::uint64_t derive_key(
    std::uint64_t master_seed,
    Stage stage,
    std::uint32_t contig_index)
{
    if (!is_valid_stage(stage)) {
        throw ContractError("stage is outside the RNG contract enum");
    }
    const std::uint64_t domain_local =
        (static_cast<std::uint64_t>(stage) << 32U) | contig_index;
    const Block domain_block = philox4x32_10(
        master_seed, domain_key_entity, domain_local);
    return static_cast<std::uint64_t>(domain_block[0])
        | (static_cast<std::uint64_t>(domain_block[1]) << 32U);
}

Block philox4x32_10(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index) noexcept
{
    Block counter = {{
        static_cast<std::uint32_t>(entity_ordinal),
        static_cast<std::uint32_t>(entity_ordinal >> 32),
        static_cast<std::uint32_t>(local_index),
        static_cast<std::uint32_t>(local_index >> 32),
    }};
    std::uint32_t key0 = static_cast<std::uint32_t>(key);
    std::uint32_t key1 = static_cast<std::uint32_t>(key >> 32);

    for (unsigned int round = 0; round < 10; ++round) {
        const std::uint64_t product0 =
            static_cast<std::uint64_t>(philox_m0) * counter[0];
        const std::uint64_t product1 =
            static_cast<std::uint64_t>(philox_m1) * counter[2];
        counter = {{
            static_cast<std::uint32_t>(product1 >> 32) ^ counter[1] ^ key0,
            static_cast<std::uint32_t>(product1),
            static_cast<std::uint32_t>(product0 >> 32) ^ counter[3] ^ key1,
            static_cast<std::uint32_t>(product0),
        }};
        key0 += philox_w0;
        key1 += philox_w1;
    }
    return counter;
}

std::uint32_t u32(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::size_t lane)
{
    if (lane >= 4) { throw ContractError("u32 lane must be in [0, 4)"); }
    return philox4x32_10(key, entity_ordinal, local_index)[lane];
}

std::uint64_t u64(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::size_t pair)
{
    if (pair >= 2) { throw ContractError("u64 pair must be in [0, 2)"); }
    const Block random_block = philox4x32_10(key, entity_ordinal, local_index);
    const std::size_t offset = pair * 2;
    return static_cast<std::uint64_t>(random_block[offset])
        | (static_cast<std::uint64_t>(random_block[offset + 1]) << 32);
}

double uniform01(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::size_t pair)
{
    static_assert(std::numeric_limits<double>::digits == 53,
                  "the RNG contract requires IEEE-754 binary64 precision");
    const std::uint64_t significand =
        u64(key, entity_ordinal, local_index, pair) >> 11;
    return std::ldexp(static_cast<double>(significand), -53);
}

bool bernoulli(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    double probability,
    std::size_t pair)
{
    if (!std::isfinite(probability) || probability < 0.0 || probability > 1.0) {
        throw ContractError("Bernoulli probability must be finite and in [0, 1]");
    }
    if (probability == 0.0) { return false; }
    if (probability == 1.0) { return true; }
    return uniform01(key, entity_ordinal, local_index, pair) < probability;
}

std::uint64_t bounded_integer(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::uint64_t upper_exclusive)
{
    return bounded_offset(
        philox4x32_10(key, entity_ordinal, local_index), upper_exclusive);
}

std::int64_t bounded_integer(
    std::uint64_t key,
    std::uint64_t entity_ordinal,
    std::uint64_t local_index,
    std::int64_t lower,
    std::int64_t upper_exclusive)
{
    if (upper_exclusive <= lower) {
        throw ContractError("bounded integer upper bound must exceed lower bound");
    }
    const std::uint64_t width =
        static_cast<std::uint64_t>(upper_exclusive)
        - static_cast<std::uint64_t>(lower);
    const std::uint64_t offset = bounded_offset(
        philox4x32_10(key, entity_ordinal, local_index), width);
    return add_signed_offset(lower, offset);
}

Block block(const Address &address)
{
    return philox4x32_10(
        derive_key(address.master_seed, address.stage, address.contig_index),
        address.entity_ordinal,
        address.local_index);
}

std::uint32_t u32(const Address &address, std::size_t lane)
{
    return u32(
        derive_key(address.master_seed, address.stage, address.contig_index),
        address.entity_ordinal,
        address.local_index,
        lane);
}

std::uint64_t u64(const Address &address, std::size_t pair)
{
    return u64(
        derive_key(address.master_seed, address.stage, address.contig_index),
        address.entity_ordinal,
        address.local_index,
        pair);
}

double uniform01(const Address &address, std::size_t pair)
{
    return uniform01(
        derive_key(address.master_seed, address.stage, address.contig_index),
        address.entity_ordinal,
        address.local_index,
        pair);
}

bool bernoulli(const Address &address, double probability, std::size_t pair)
{
    return bernoulli(
        derive_key(address.master_seed, address.stage, address.contig_index),
        address.entity_ordinal,
        address.local_index,
        probability,
        pair);
}

std::uint64_t bounded_integer(
    const Address &address,
    std::uint64_t upper_exclusive)
{
    return bounded_integer(
        derive_key(address.master_seed, address.stage, address.contig_index),
        address.entity_ordinal,
        address.local_index,
        upper_exclusive);
}

std::int64_t bounded_integer(
    const Address &address,
    std::int64_t lower,
    std::int64_t upper_exclusive)
{
    return bounded_integer(
        derive_key(address.master_seed, address.stage, address.contig_index),
        address.entity_ordinal,
        address.local_index,
        lower,
        upper_exclusive);
}

} // namespace htsim::rng

// ---- sha256 --------------------------------------------------------

namespace htsim::crypto {
namespace {

constexpr std::array<std::uint32_t, 64> constants = {{
    UINT32_C(0x428a2f98), UINT32_C(0x71374491), UINT32_C(0xb5c0fbcf),
    UINT32_C(0xe9b5dba5), UINT32_C(0x3956c25b), UINT32_C(0x59f111f1),
    UINT32_C(0x923f82a4), UINT32_C(0xab1c5ed5), UINT32_C(0xd807aa98),
    UINT32_C(0x12835b01), UINT32_C(0x243185be), UINT32_C(0x550c7dc3),
    UINT32_C(0x72be5d74), UINT32_C(0x80deb1fe), UINT32_C(0x9bdc06a7),
    UINT32_C(0xc19bf174), UINT32_C(0xe49b69c1), UINT32_C(0xefbe4786),
    UINT32_C(0x0fc19dc6), UINT32_C(0x240ca1cc), UINT32_C(0x2de92c6f),
    UINT32_C(0x4a7484aa), UINT32_C(0x5cb0a9dc), UINT32_C(0x76f988da),
    UINT32_C(0x983e5152), UINT32_C(0xa831c66d), UINT32_C(0xb00327c8),
    UINT32_C(0xbf597fc7), UINT32_C(0xc6e00bf3), UINT32_C(0xd5a79147),
    UINT32_C(0x06ca6351), UINT32_C(0x14292967), UINT32_C(0x27b70a85),
    UINT32_C(0x2e1b2138), UINT32_C(0x4d2c6dfc), UINT32_C(0x53380d13),
    UINT32_C(0x650a7354), UINT32_C(0x766a0abb), UINT32_C(0x81c2c92e),
    UINT32_C(0x92722c85), UINT32_C(0xa2bfe8a1), UINT32_C(0xa81a664b),
    UINT32_C(0xc24b8b70), UINT32_C(0xc76c51a3), UINT32_C(0xd192e819),
    UINT32_C(0xd6990624), UINT32_C(0xf40e3585), UINT32_C(0x106aa070),
    UINT32_C(0x19a4c116), UINT32_C(0x1e376c08), UINT32_C(0x2748774c),
    UINT32_C(0x34b0bcb5), UINT32_C(0x391c0cb3), UINT32_C(0x4ed8aa4a),
    UINT32_C(0x5b9cca4f), UINT32_C(0x682e6ff3), UINT32_C(0x748f82ee),
    UINT32_C(0x78a5636f), UINT32_C(0x84c87814), UINT32_C(0x8cc70208),
    UINT32_C(0x90befffa), UINT32_C(0xa4506ceb), UINT32_C(0xbef9a3f7),
    UINT32_C(0xc67178f2),
}};

constexpr std::array<std::uint32_t, 8> initial_hash = {{
    UINT32_C(0x6a09e667), UINT32_C(0xbb67ae85), UINT32_C(0x3c6ef372),
    UINT32_C(0xa54ff53a), UINT32_C(0x510e527f), UINT32_C(0x9b05688c),
    UINT32_C(0x1f83d9ab), UINT32_C(0x5be0cd19),
}};

std::uint32_t rotate_right(std::uint32_t value, unsigned int count) noexcept
{
    return (value >> count) | (value << (32U - count));
}

void compress(
    std::array<std::uint32_t, 8> &hash,
    const std::uint8_t *block) noexcept
{
    std::array<std::uint32_t, 64> words = {};
    for (std::size_t index = 0; index < 16; ++index) {
        const std::size_t offset = index * 4U;
        words[index] =
            (static_cast<std::uint32_t>(block[offset]) << 24U)
            | (static_cast<std::uint32_t>(block[offset + 1U]) << 16U)
            | (static_cast<std::uint32_t>(block[offset + 2U]) << 8U)
            | static_cast<std::uint32_t>(block[offset + 3U]);
    }
    for (std::size_t index = 16; index < words.size(); ++index) {
        const std::uint32_t previous_15 = words[index - 15U];
        const std::uint32_t previous_2 = words[index - 2U];
        const std::uint32_t sigma0 =
            rotate_right(previous_15, 7U) ^ rotate_right(previous_15, 18U)
            ^ (previous_15 >> 3U);
        const std::uint32_t sigma1 =
            rotate_right(previous_2, 17U) ^ rotate_right(previous_2, 19U)
            ^ (previous_2 >> 10U);
        words[index] = words[index - 16U] + sigma0
            + words[index - 7U] + sigma1;
    }

    std::uint32_t a = hash[0];
    std::uint32_t b = hash[1];
    std::uint32_t c = hash[2];
    std::uint32_t d = hash[3];
    std::uint32_t e = hash[4];
    std::uint32_t f = hash[5];
    std::uint32_t g = hash[6];
    std::uint32_t h = hash[7];
    for (std::size_t index = 0; index < words.size(); ++index) {
        const std::uint32_t sum1 =
            rotate_right(e, 6U) ^ rotate_right(e, 11U)
            ^ rotate_right(e, 25U);
        const std::uint32_t choice = (e & f) ^ ((~e) & g);
        const std::uint32_t temporary1 =
            h + sum1 + choice + constants[index] + words[index];
        const std::uint32_t sum0 =
            rotate_right(a, 2U) ^ rotate_right(a, 13U)
            ^ rotate_right(a, 22U);
        const std::uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
        const std::uint32_t temporary2 = sum0 + majority;
        h = g;
        g = f;
        f = e;
        e = d + temporary1;
        d = c;
        c = b;
        b = a;
        a = temporary1 + temporary2;
    }

    hash[0] += a;
    hash[1] += b;
    hash[2] += c;
    hash[3] += d;
    hash[4] += e;
    hash[5] += f;
    hash[6] += g;
    hash[7] += h;
}

} // namespace

Sha256::Sha256() noexcept : hash_(initial_hash) {}

void Sha256::update(const std::uint8_t *data, std::size_t size)
{
    if (size != 0U && data == nullptr) {
        throw Sha256Error("non-empty SHA-256 input has a null pointer");
    }
    if constexpr (sizeof(std::size_t) > sizeof(std::uint64_t)) {
        if (size > std::numeric_limits<std::uint64_t>::max()) {
            throw Sha256Error("SHA-256 byte count exceeds uint64");
        }
    }
    const auto increment = static_cast<std::uint64_t>(size);
    if (increment > maximum_sha256_input_bytes - total_bytes_) {
        throw Sha256Error("SHA-256 input exceeds the 2^64-bit length limit");
    }
    total_bytes_ += increment;

    std::size_t offset = 0;
    if (buffered_bytes_ != 0U) {
        const std::size_t copied =
            std::min(size, buffer_.size() - buffered_bytes_);
        if (copied != 0U) {
            std::memcpy(buffer_.data() + buffered_bytes_, data, copied);
        }
        buffered_bytes_ += copied;
        offset += copied;
        if (buffered_bytes_ == buffer_.size()) {
            compress(hash_, buffer_.data());
            buffered_bytes_ = 0;
        }
    }
    while (size - offset >= buffer_.size()) {
        compress(hash_, data + offset);
        offset += buffer_.size();
    }
    const std::size_t remaining = size - offset;
    if (remaining != 0U) {
        std::memcpy(buffer_.data(), data + offset, remaining);
        buffered_bytes_ = remaining;
    }
}

void Sha256::update(const std::vector<std::uint8_t> &data)
{
    update(data.data(), data.size());
}

Sha256Digest Sha256::digest() const noexcept
{
    auto hash = hash_;
    auto buffer = buffer_;
    std::size_t buffered_bytes = buffered_bytes_;
    const std::uint64_t total_bits = total_bytes_ * UINT64_C(8);

    buffer[buffered_bytes++] = 0x80U;
    if (buffered_bytes > 56U) {
        std::fill(
            buffer.begin() + static_cast<std::ptrdiff_t>(buffered_bytes),
            buffer.end(), 0U);
        compress(hash, buffer.data());
        buffer.fill(0U);
        buffered_bytes = 0;
    }
    std::fill(
        buffer.begin() + static_cast<std::ptrdiff_t>(buffered_bytes),
        buffer.begin() + 56, 0U);
    for (unsigned int index = 0; index < 8U; ++index) {
        buffer[63U - index] = static_cast<std::uint8_t>(
            total_bits >> (8U * index));
    }
    compress(hash, buffer.data());

    Sha256Digest result = {};
    for (std::size_t word = 0; word < hash.size(); ++word) {
        for (unsigned int byte = 0; byte < 4U; ++byte) {
            result[word * 4U + byte] = static_cast<std::uint8_t>(
                hash[word] >> (24U - 8U * byte));
        }
    }
    return result;
}

Sha256Digest sha256(const std::uint8_t *data, std::size_t size)
{
    Sha256 hasher;
    hasher.update(data, size);
    return hasher.digest();
}

Sha256Digest sha256(const std::vector<std::uint8_t> &data)
{
    return sha256(data.data(), data.size());
}

} // namespace htsim::crypto

// ---- text_snapshot --------------------------------------------------------

namespace htsim::text {
namespace {

constexpr std::size_t io_buffer_size = 64U * 1024U;

class HtsReader {
public:
    HtsReader(int descriptor, const std::string &path) : path_(path)
    {
        const int duplicate = dup(descriptor);
        if (duplicate < 0) {
            throw TextSnapshotError(
                "cannot duplicate text input " + path_ + ": "
                + std::strerror(errno));
        }
        if (lseek(duplicate, 0, SEEK_SET) != 0) {
            const int saved_errno = errno;
            (void)::close(duplicate);
            throw TextSnapshotError(
                "cannot rewind text input " + path_ + ": "
                + std::strerror(saved_errno));
        }
        hFILE *raw = hdopen(duplicate, "r");
        if (raw == nullptr) {
            const int saved_errno = errno;
            (void)::close(duplicate);
            throw TextSnapshotError(
                "HTSlib cannot adopt text input " + path_ + ": "
                + std::strerror(saved_errno));
        }
        file_ = hts_hopen(raw, path_.c_str(), "r");
        if (file_ == nullptr) {
            const int saved_errno = errno;
            const int close_status = hclose(raw);
            (void)close_status;
            throw TextSnapshotError(
                "HTSlib cannot open text input " + path_ + ": "
                + std::strerror(saved_errno));
        }
    }

    ~HtsReader()
    {
        if (file_ != nullptr) {(void)hts_close(file_);}
        std::free(line_.s);
    }

    HtsReader(const HtsReader &) = delete;
    HtsReader &operator=(const HtsReader &) = delete;

    htsFile *file() const noexcept {return file_;}

    void visit(
        const LineVisitor &visitor,
        std::size_t maximum_decoded_line_bytes)
    {
        while (true) {
            const int status = hts_getline(file_, '\n', &line_);
            if (status == -1) {break;}
            if (status < -1) {
                throw TextSnapshotError(
                    "HTSlib failed while decoding text input: " + path_);
            }
            if (line_.l > maximum_decoded_line_bytes) {
                throw TextSnapshotError("text line exceeds its configured limit");
            }
            std::string_view line(line_.s, line_.l);
            if (!line.empty() && line.back() == '\r') {line.remove_suffix(1U);}
            if (line.find('\r') != std::string_view::npos) {
                throw TextSnapshotError(
                    "text input contains a bare carriage return");
            }
            if (line_number_ == std::numeric_limits<std::uint64_t>::max()) {
                throw TextSnapshotError("text line count exceeds uint64");
            }
            ++line_number_;
            visitor(line, line_number_);
        }
    }

    void visit_chunks(const ChunkVisitor &visitor)
    {
        std::array<char, io_buffer_size> buffer = {};
        while (true) {
            ssize_t observed = -1;
            if (file_->is_cram) {
                throw TextSnapshotError(
                    "HTSlib selected a non-text CRAM backend: " + path_);
            }
            if (file_->is_bgzf) {
                observed = bgzf_read(
                    file_->fp.bgzf, buffer.data(), buffer.size());
            } else {
                observed = hread(
                    file_->fp.hfile, buffer.data(), buffer.size());
            }
            if (observed == 0) {break;}
            if (observed < 0) {
                throw TextSnapshotError(
                    "HTSlib failed while decoding text input: " + path_);
            }
            visitor(std::string_view(
                buffer.data(), static_cast<std::size_t>(observed)));
        }
    }

    void close()
    {
        if (file_ == nullptr) {return;}
        htsFile *closing = file_;
        file_ = nullptr;
        if (hts_close(closing) != 0) {
            throw TextSnapshotError(
                "HTSlib failed while closing text input: " + path_);
        }
    }

private:
    std::string path_;
    htsFile *file_ = nullptr;
    kstring_t line_ = KS_INITIALIZE;
    std::uint64_t line_number_ = 0;
};

struct FileIdentity {
    dev_t device = 0;
    ino_t inode = 0;
    mode_t mode = 0;
    off_t size = 0;
    timespec modification_time = {};
};

bool same_identity(const FileIdentity &left, const FileIdentity &right) noexcept
{
    return left.device == right.device && left.inode == right.inode
        && left.mode == right.mode && left.size == right.size
        && left.modification_time.tv_sec == right.modification_time.tv_sec
        && left.modification_time.tv_nsec == right.modification_time.tv_nsec;
}

std::string system_message(const std::string &operation, const std::string &path)
{
    return operation + " text input " + path + ": " + std::strerror(errno);
}

FileIdentity read_identity(int descriptor, const std::string &path)
{
    struct stat status = {};
    if (fstat(descriptor, &status) != 0) {
        throw TextSnapshotError(system_message("cannot inspect", path));
    }
    if (!S_ISREG(status.st_mode)) {
        throw TextSnapshotError("text input is not a regular file: " + path);
    }
    if (status.st_size < 0) {
        throw TextSnapshotError("text input has a negative byte size: " + path);
    }
    if (static_cast<std::uintmax_t>(status.st_size)
        > crypto::maximum_sha256_input_bytes) {
        throw TextSnapshotError("text input exceeds the SHA-256 limit: " + path);
    }
    return {
        status.st_dev,
        status.st_ino,
        status.st_mode,
        status.st_size,
        status.st_mtim,
    };
}

void require_identity(
    int descriptor,
    const std::string &path,
    const FileIdentity &expected)
{
    if (!same_identity(read_identity(descriptor, path), expected)) {
        throw TextSnapshotError("text input identity changed: " + path);
    }
}

void rewind_descriptor(int descriptor, const std::string &path)
{
    if (lseek(descriptor, 0, SEEK_SET) != 0) {
        throw TextSnapshotError(system_message("cannot rewind", path));
    }
}

ssize_t read_descriptor(
    int descriptor,
    std::uint8_t *buffer,
    std::size_t size,
    const std::string &path)
{
    while (true) {
        const ssize_t count = read(descriptor, buffer, size);
        if (count >= 0) {return count;}
        if (errno != EINTR) {
            throw TextSnapshotError(system_message("cannot read", path));
        }
    }
}

crypto::Sha256Digest raw_hash_pass(
    int descriptor,
    const std::string &path,
    const FileIdentity &identity)
{
    require_identity(descriptor, path, identity);
    rewind_descriptor(descriptor, path);
    crypto::Sha256 hash;
    std::uint64_t byte_count = 0;
    std::array<std::uint8_t, io_buffer_size> buffer = {};
    while (true) {
        const ssize_t observed =
            read_descriptor(descriptor, buffer.data(), buffer.size(), path);
        if (observed == 0) {break;}
        const auto size = static_cast<std::size_t>(observed);
        hash.update(buffer.data(), size);
        byte_count += static_cast<std::uint64_t>(size);
    }
    require_identity(descriptor, path, identity);
    if (byte_count != static_cast<std::uint64_t>(identity.size)) {
        throw TextSnapshotError("text input byte size changed: " + path);
    }
    return hash.digest();
}

crypto::Sha256Digest decoded_pass(
    int descriptor,
    const std::string &path,
    const FileIdentity &identity,
    const LineVisitor &visitor,
    std::size_t maximum_decoded_line_bytes)
{
    require_identity(descriptor, path, identity);
    const crypto::Sha256Digest before = raw_hash_pass(
        descriptor, path, identity);
    HtsReader reader(descriptor, path);
    reader.visit(visitor, maximum_decoded_line_bytes);
    reader.close();
    require_identity(descriptor, path, identity);
    return before;
}

crypto::Sha256Digest hts_pass(
    int descriptor,
    const std::string &path,
    const FileIdentity &identity,
    const HtsFileVisitor &visitor)
{
    require_identity(descriptor, path, identity);
    const crypto::Sha256Digest before = raw_hash_pass(
        descriptor, path, identity);
    HtsReader reader(descriptor, path);
    visitor(reader.file());
    reader.close();
    require_identity(descriptor, path, identity);
    return before;
}

crypto::Sha256Digest decoded_chunk_pass(
    int descriptor,
    const std::string &path,
    const FileIdentity &identity,
    const ChunkVisitor &visitor)
{
    require_identity(descriptor, path, identity);
    const crypto::Sha256Digest before = raw_hash_pass(
        descriptor, path, identity);
    HtsReader reader(descriptor, path);
    reader.visit_chunks(visitor);
    reader.close();
    require_identity(descriptor, path, identity);
    return before;
}

class VisitGuard {
public:
    explicit VisitGuard(std::atomic_flag &active) : active_(&active) {}
    ~VisitGuard() {active_->clear(std::memory_order_release);}

    VisitGuard(const VisitGuard &) = delete;
    VisitGuard &operator=(const VisitGuard &) = delete;

private:
    std::atomic_flag *active_;
};

} // namespace

class TextSnapshot::Impl {
public:
    Impl(
        const std::string &source_path,
        const crypto::Sha256Digest &expected,
        std::size_t line_limit)
        : path(source_path), file_sha256(expected), maximum_line_bytes(line_limit)
    {
        if (maximum_line_bytes == 0U) {
            throw TextSnapshotError("text snapshot line limit must be positive");
        }
        descriptor = open(path.c_str(), O_RDONLY | O_CLOEXEC | O_NONBLOCK);
        if (descriptor < 0) {
            throw TextSnapshotError(system_message("cannot open", path));
        }
        try {
            identity = read_identity(descriptor, path);
            if (raw_hash_pass(descriptor, path, identity) != file_sha256) {
                throw TextSnapshotError(
                    "text input SHA-256 does not match the expected digest");
            }
        } catch (...) {
            (void)close(descriptor);
            descriptor = -1;
            throw;
        }
    }

    ~Impl()
    {
        if (descriptor >= 0) {(void)close(descriptor);}
    }

    int descriptor = -1;
    std::string path;
    FileIdentity identity;
    crypto::Sha256Digest file_sha256 = {};
    std::size_t maximum_line_bytes = 0;
    std::atomic_flag visit_active = ATOMIC_FLAG_INIT;
    std::atomic<bool> poisoned{false};
};

TextSnapshot::TextSnapshot(
    const std::string &path,
    const crypto::Sha256Digest &expected_file_sha256,
    std::size_t maximum_decoded_line_bytes)
    : impl_(std::make_unique<Impl>(
          path, expected_file_sha256, maximum_decoded_line_bytes))
{
}

TextSnapshot::~TextSnapshot() = default;

const crypto::Sha256Digest &TextSnapshot::file_sha256() const noexcept
{
    return impl_->file_sha256;
}

void TextSnapshot::visit_lines(const LineVisitor &visitor)
{
    if (impl_->poisoned.load(std::memory_order_acquire)) {
        throw TextSnapshotError("text snapshot is poisoned after a prior failure");
    }
    if (!visitor) {throw TextSnapshotError("text snapshot requires a line visitor");}
    if (impl_->visit_active.test_and_set(std::memory_order_acquire)) {
        throw TextSnapshotError("text snapshot visits may not overlap or re-enter");
    }
    VisitGuard guard(impl_->visit_active);
    try {
        if (decoded_pass(
                impl_->descriptor,
                impl_->path,
                impl_->identity,
                visitor,
                impl_->maximum_line_bytes)
            != impl_->file_sha256) {
            throw TextSnapshotError("text input changed during line decoding");
        }
        if (raw_hash_pass(impl_->descriptor, impl_->path, impl_->identity)
            != impl_->file_sha256) {
            throw TextSnapshotError("text input changed after line decoding");
        }
    } catch (...) {
        impl_->poisoned.store(true, std::memory_order_release);
        throw;
    }
}

void TextSnapshot::visit_hts(const HtsFileVisitor &visitor)
{
    if (impl_->poisoned.load(std::memory_order_acquire)) {
        throw TextSnapshotError("text snapshot is poisoned after a prior failure");
    }
    if (!visitor) {
        throw TextSnapshotError("text snapshot requires an HTSlib visitor");
    }
    if (impl_->visit_active.test_and_set(std::memory_order_acquire)) {
        throw TextSnapshotError("text snapshot visits may not overlap or re-enter");
    }
    VisitGuard guard(impl_->visit_active);
    try {
        if (hts_pass(
                impl_->descriptor,
                impl_->path,
                impl_->identity,
                visitor)
            != impl_->file_sha256) {
            throw TextSnapshotError("text input changed during HTSlib decoding");
        }
        if (raw_hash_pass(impl_->descriptor, impl_->path, impl_->identity)
            != impl_->file_sha256) {
            throw TextSnapshotError("text input changed after HTSlib decoding");
        }
    } catch (...) {
        impl_->poisoned.store(true, std::memory_order_release);
        throw;
    }
}

void TextSnapshot::visit_chunks(const ChunkVisitor &visitor)
{
    if (impl_->poisoned.load(std::memory_order_acquire)) {
        throw TextSnapshotError("text snapshot is poisoned after a prior failure");
    }
    if (!visitor) {
        throw TextSnapshotError("text snapshot requires a chunk visitor");
    }
    if (impl_->visit_active.test_and_set(std::memory_order_acquire)) {
        throw TextSnapshotError("text snapshot visits may not overlap or re-enter");
    }
    VisitGuard guard(impl_->visit_active);
    try {
        if (decoded_chunk_pass(
                impl_->descriptor,
                impl_->path,
                impl_->identity,
                visitor)
            != impl_->file_sha256) {
            throw TextSnapshotError("text input changed during chunk decoding");
        }
        if (raw_hash_pass(impl_->descriptor, impl_->path, impl_->identity)
            != impl_->file_sha256) {
            throw TextSnapshotError("text input changed after chunk decoding");
        }
    } catch (...) {
        impl_->poisoned.store(true, std::memory_order_release);
        throw;
    }
}

void TextSnapshot::verify_unchanged()
{
    if (impl_->poisoned.load(std::memory_order_acquire)) {
        throw TextSnapshotError("text snapshot is poisoned after a prior failure");
    }
    if (impl_->visit_active.test_and_set(std::memory_order_acquire)) {
        throw TextSnapshotError("text snapshot visits may not overlap or re-enter");
    }
    VisitGuard guard(impl_->visit_active);
    try {
        if (raw_hash_pass(impl_->descriptor, impl_->path, impl_->identity)
            != impl_->file_sha256) {
            throw TextSnapshotError("text input changed after verified access");
        }
    } catch (...) {
        impl_->poisoned.store(true, std::memory_order_release);
        throw;
    }
}

} // namespace htsim::text
