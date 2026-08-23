#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

#include "types.h"
#include "methdb.h"
#include "protocol.h"

namespace {

using htsim::model::HaplotypeMask;
using htsim::methdb::EntityError;
using htsim::methdb::SiteEntityKind;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const EntityError &) {
        return;
    }
    throw std::runtime_error(message);
}

void test_exact_layout_and_decode()
{
    using namespace htsim::methdb;
    const auto reference = reference_site_entity(UINT32_C(0xffffffff));
    const auto shared_variant = variant_reference_site_entity(
        17U, HaplotypeMask::both, 0U);
    const auto haplotype_0 = variant_reference_site_entity(
        17U, HaplotypeMask::haplotype_1, 0U);
    const auto haplotype_1 = variant_reference_site_entity(
        17U, HaplotypeMask::haplotype_2, 1U);
    const auto shared_insertion = insertion_site_entity(
        7U, 2U, HaplotypeMask::both, 1U);
    const auto insertion_0 = insertion_site_entity(
        7U, 2U, HaplotypeMask::haplotype_1, 0U);
    const auto insertion_1 = insertion_site_entity(
        7U, 2U, HaplotypeMask::haplotype_2, 1U);

    require(reference.value() == UINT64_C(0x00000000ffffffff)
                && shared_variant.value() == UINT64_C(0x1000000000000011)
                && haplotype_0.value() == UINT64_C(0x2000000000000011)
                && haplotype_1.value() == UINT64_C(0x3000000000000011)
                && shared_insertion.value() == UINT64_C(0x400000000000001e)
                && insertion_0.value() == UINT64_C(0x500000000000001e)
                && insertion_1.value() == UINT64_C(0x600000000000001e),
            "64-bit methylation entity layout changed");

    const auto decoded_reference = decode_site_entity(haplotype_1);
    const auto decoded_insertion = decode_site_entity(shared_insertion);
    require(decoded_reference.kind
                == SiteEntityKind::variant_reference_haplotype_1
                && decoded_reference.reference_position == 17U,
            "variant reference entity did not decode");
    require(decoded_insertion.kind == SiteEntityKind::insertion_shared
                && decoded_insertion.event_ordinal == 7U
                && decoded_insertion.insertion_offset == 2U,
            "insertion entity did not decode");
    require(!entity_is_insertion(haplotype_1)
                && entity_is_insertion(insertion_0),
            "insertion entity classification changed");
}

void test_invalid_values_fail_closed()
{
    using namespace htsim::methdb;
    require_error(
        [] {
            (void)variant_reference_site_entity(
                1U, HaplotypeMask::haplotype_1, 1U);
        },
        "reference site absent from the selected haplotype was accepted");
    require_error(
        [] {
            (void)insertion_site_entity(
                1U, 0U, static_cast<HaplotypeMask>(0), 0U);
        },
        "invalid haplotype mask was accepted");
    require_error(
        [] {
            (void)insertion_site_entity(
                htsim::model::no_variant_index,
                0U,
                HaplotypeMask::both,
                0U);
        },
        "no-event insertion ordinal was accepted");
    require_error(
        [] {
            (void)insertion_site_entity(
                1U, 4U, HaplotypeMask::both, 0U);
        },
        "fifth insertion base was accepted");
    require_error(
        [] {decode_site_entity(UINT64_C(0x7000000000000000));},
        "reserved entity kind was accepted");
    require_error(
        [] {decode_site_entity(UINT64_C(0x1000000100000000));},
        "non-zero reference reserved bits were accepted");
    require_error(
        [] {decode_site_entity(UINT64_C(0x4000000400000000));},
        "non-zero insertion reserved bits were accepted");
    require_error(
        [] {decode_site_entity(UINT64_C(0x40000003fffffffc));},
        "no-event insertion payload was accepted");
}

} // namespace

int main()
{
    try {
        test_exact_layout_and_decode();
        test_invalid_values_fail_closed();
    } catch (const std::exception &error) {
        std::cerr << "site_entity_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
