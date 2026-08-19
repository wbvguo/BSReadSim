#include "utilities.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using htsim::crypto::Sha256;
using htsim::crypto::Sha256Digest;

[[noreturn]] void fail(const std::string &message)
{
    throw std::runtime_error(message);
}

void require(bool condition, const std::string &message)
{
    if (!condition) {fail(message);}
}

std::uint8_t hex_nibble(char value)
{
    if (value >= '0' && value <= '9') {
        return static_cast<std::uint8_t>(value - '0');
    }
    if (value >= 'a' && value <= 'f') {
        return static_cast<std::uint8_t>(value - 'a' + 10);
    }
    fail("invalid digest fixture");
}

Sha256Digest from_hex(const std::string &text)
{
    require(text.size() == 64U, "digest fixture must have 64 hex characters");
    Sha256Digest result = {};
    for (std::size_t index = 0; index < result.size(); ++index) {
        result[index] = static_cast<std::uint8_t>(
            (hex_nibble(text[index * 2U]) << 4U)
            | hex_nibble(text[index * 2U + 1U]));
    }
    return result;
}

void test_standard_vectors()
{
    require(
        htsim::crypto::sha256(nullptr, 0U)
            == from_hex("e3b0c44298fc1c149afbf4c8996fb924"
                        "27ae41e4649b934ca495991b7852b855"),
        "empty SHA-256 vector changed");
    const std::string abc = "abc";
    require(
        htsim::crypto::sha256(
            reinterpret_cast<const std::uint8_t *>(abc.data()), abc.size())
            == from_hex("ba7816bf8f01cfea414140de5dae2223"
                        "b00361a396177a9cb410ff61f20015ad"),
        "abc SHA-256 vector changed");
}

void test_incremental_boundaries()
{
    std::vector<std::uint8_t> bytes(131075U);
    for (std::size_t index = 0; index < bytes.size(); ++index) {
        bytes[index] = static_cast<std::uint8_t>(index * 37U + 11U);
    }
    const Sha256Digest expected = htsim::crypto::sha256(bytes);
    for (const std::size_t chunk : std::array<std::size_t, 9>{
             1U, 2U, 7U, 55U, 56U, 63U, 64U, 65U, 65536U}) {
        Sha256 hasher;
        std::size_t offset = 0;
        while (offset < bytes.size()) {
            const std::size_t size = std::min(chunk, bytes.size() - offset);
            hasher.update(bytes.data() + offset, size);
            offset += size;
        }
        require(hasher.size() == bytes.size(), "incremental byte count changed");
        require(hasher.digest() == expected, "incremental digest changed");
    }
}

void test_snapshot_and_validation()
{
    Sha256 hasher;
    const std::vector<std::uint8_t> prefix = {'a', 'b'};
    const std::vector<std::uint8_t> suffix = {'c'};
    hasher.update(prefix);
    const Sha256Digest intermediate = hasher.digest();
    require(intermediate == htsim::crypto::sha256(prefix),
            "digest() mutated or misreported intermediate state");
    hasher.update(suffix);
    require(hasher.digest() == from_hex(
        "ba7816bf8f01cfea414140de5dae2223"
        "b00361a396177a9cb410ff61f20015ad"),
        "updates after digest() changed");

    bool threw = false;
    try {
        hasher.update(nullptr, 1U);
    } catch (const htsim::crypto::Sha256Error &) {
        threw = true;
    }
    require(threw, "non-empty null input was accepted");
    hasher.update(nullptr, 0U);
}

} // namespace

int main()
{
    try {
        test_standard_vectors();
        test_incremental_boundaries();
        test_snapshot_and_validation();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
