#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <string>

#include "utilities.h"

namespace {

using htsim::rng::Address;
using htsim::rng::Block;
using htsim::rng::ContractError;

void check(bool condition, const char *message)
{
    if (!condition) {
        std::fprintf(stderr, "rng_test: %s\n", message);
        std::exit(EXIT_FAILURE);
    }
}

template <typename Function>
void check_contract_error(Function function, const char *message)
{
    try {
        function();
    } catch (const ContractError &) {
        return;
    }
    check(false, message);
}

} // namespace

int main()
{
    using namespace htsim::rng;
    constexpr std::uint64_t maximum_u64 =
        std::numeric_limits<std::uint64_t>::max();

    check(contract_id == "philox4x32-10+philox-domain",
          "RNG contract identifier changed");
    const std::array<std::string_view, 9> expected_stages = {{
        "mutation", "methylation-level", "fragment", "haplotype", "site-state",
        "library-orientation", "conversion", "quality", "sequencing-error",
    }};
    check(stage_names == expected_stages, "RNG stage whitelist changed");
    check(is_valid_stage(Stage::mutation), "valid stage was rejected");
    check(!is_valid_stage(Stage::count), "stage count sentinel was accepted");

    check(parse_u64_decimal("0") == 0, "zero seed parsing failed");
    check(parse_u64_decimal("18446744073709551615") == maximum_u64,
          "maximum u64 seed parsing failed");
    check_contract_error([] { (void)parse_u64_decimal(""); },
                         "empty u64 was accepted");
    check_contract_error([] { (void)parse_u64_decimal("-1"); },
                         "negative u64 was accepted");
    check_contract_error([] { (void)parse_u64_decimal("00"); },
                         "non-canonical u64 was accepted");
    check_contract_error([] { (void)parse_u64_decimal("18446744073709551616"); },
                         "overflowing u64 was accepted");

    check(derive_key(0, Stage::mutation, 0U)
              == UINT64_C(0x1cd75ac86fb9d4fd),
          "Philox domain seed-zero vector changed");
    check(derive_key(maximum_u64, Stage::sequencing_error, UINT32_MAX)
              == UINT64_C(0x8a25ab3f5a158aa5),
          "Philox domain maximum vector changed");
    check(derive_key(42, Stage::quality, 3U)
              == UINT64_C(0xba6bdd075c69a508),
          "Philox stage/contig-index vector changed");
    check(derive_key(99, Stage::fragment, 93U)
              == UINT64_C(0xc599c66331fff182),
          "Philox numeric contig vector changed");
    check(derive_key(99, Stage::fragment, 94U)
              == UINT64_C(0x1cbf5d41b9a046fd),
          "Philox adjacent-contig vector changed");
    check(derive_key(42, Stage::quality, 3U)
              != derive_key(42, Stage::quality, 0U),
          "contig domains are not isolated");
    check(derive_key(42, Stage::quality, 3U)
              != derive_key(42, Stage::conversion, 3U),
          "stage domains are not isolated");
    check_contract_error([] {
        (void)derive_key(7, static_cast<Stage>(UINT32_MAX), 0U);
    },
                         "unknown stage was accepted");

    const Block domain_block = philox4x32_10(
        42U,
        domain_key_entity,
        (static_cast<std::uint64_t>(Stage::quality) << 32U) | UINT64_C(3));
    check((static_cast<std::uint64_t>(domain_block[0])
               | (static_cast<std::uint64_t>(domain_block[1]) << 32U))
              == derive_key(42U, Stage::quality, 3U),
          "domain-key mapping is not the reserved Philox block");

    const Block zero_vector = {{
        UINT32_C(0x6627e8d5), UINT32_C(0xe169c58d),
        UINT32_C(0xbc57ac4c), UINT32_C(0x9b00dbd8),
    }};
    check(philox4x32_10(0, 0, 0) == zero_vector,
          "Random123 Philox zero vector changed");

    const std::uint64_t site_key = derive_key(0, Stage::site_state, 0U);
    check(site_key == UINT64_C(0x98ce57ddcaf9036f),
          "site-state derived key changed");
    const Block site_vector = {{
        UINT32_C(0xad602429), UINT32_C(0x806eae81),
        UINT32_C(0x3886458d), UINT32_C(0x5c257f25),
    }};
    check(philox4x32_10(site_key, 0, 0) == site_vector,
          "derived-key Philox vector changed");

    const std::uint64_t max_key =
        derive_key(maximum_u64, Stage::sequencing_error, UINT32_MAX);
    const Block max_vector = {{
        UINT32_C(0x5a0278ff), UINT32_C(0x29be839c),
        UINT32_C(0xd8f61380), UINT32_C(0x75af06ab),
    }};
    check(philox4x32_10(max_key, maximum_u64, maximum_u64) == max_vector,
          "maximum Philox counter vector changed");

    const Address address = {
        42, Stage::quality, 3U, UINT64_C(123456789), UINT64_C(987654321),
    };
    const std::uint64_t key = derive_key(
        address.master_seed, address.stage, address.contig_index);
    const Block address_vector = {{
        UINT32_C(0x9edd88a7), UINT32_C(0x95098494),
        UINT32_C(0xeddc2146), UINT32_C(0x6fc075ec),
    }};
    check(block(address) == address_vector, "full-address block vector changed");
    check(u32(address, 0) == UINT32_C(0x9edd88a7), "u32 lane zero changed");
    check(u32(address, 3) == UINT32_C(0x6fc075ec), "u32 lane three changed");
    check(u64(address, 0) == UINT64_C(0x950984949edd88a7),
          "u64 pair zero changed");
    check(u64(address, 1) == UINT64_C(0x6fc075eceddc2146),
          "u64 pair one changed");
    check(uniform01(address) == 0x1.2a1309293dbb1p-1,
          "binary64 uniform conversion changed");
    check(!bernoulli(address, 0.0), "Bernoulli p=0 must be false");
    check(bernoulli(address, 1.0), "Bernoulli p=1 must be true");
    check(!bernoulli(address, 0.5), "Bernoulli p=0.5 vector changed");
    check(!bernoulli(address, 0.25), "Bernoulli p=0.25 vector changed");

    check(bounded_integer(address, UINT64_C(10)) == 4,
          "unsigned bounded-integer vector changed");
    check(bounded_integer(address, INT64_C(-5), INT64_C(6)) == -1,
          "signed bounded-integer vector changed");
    check(bounded_integer(address, INT64_C(10), INT64_C(20)) == 14,
          "shifted bounded-integer vector changed");
    check(bounded_integer(address, maximum_u64)
              == UINT64_C(0x6fc075eceddc2146),
          "maximum-width bounded-integer vector changed");
    check(bounded_integer(
              address,
              std::numeric_limits<std::int64_t>::min(),
              std::numeric_limits<std::int64_t>::max())
              == INT64_C(-1170806242652970682),
          "full signed-range bounded-integer vector changed");

    check(philox4x32_10(key, address.entity_ordinal, address.local_index)
              == block(address),
          "full-address and pre-derived-key APIs disagree");
    check_contract_error([&] { (void)u32(address, 4); },
                         "invalid u32 lane was accepted");
    check_contract_error([&] { (void)u64(address, 2); },
                         "invalid u64 pair was accepted");
    check_contract_error([&] {
        (void)bernoulli(address, std::numeric_limits<double>::quiet_NaN());
    }, "NaN Bernoulli probability was accepted");
    check_contract_error([&] { (void)bernoulli(address, -0.1); },
                         "negative Bernoulli probability was accepted");
    check_contract_error([&] { (void)bounded_integer(address, UINT64_C(0)); },
                         "zero bounded-integer width was accepted");
    check_contract_error([&] {
        (void)bounded_integer(address, INT64_C(5), INT64_C(5));
    }, "empty signed bounded-integer range was accepted");

    return EXIT_SUCCESS;
}
