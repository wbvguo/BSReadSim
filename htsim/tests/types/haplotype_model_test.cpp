#include "types.h"

#include <cstdint>
#include <iostream>
#include <stdexcept>

namespace {

void require(bool condition, const char *message)
{
    if (!condition) {throw std::runtime_error(message);}
}

} // namespace

int main()
{
    using htsim::model::HaplotypeMask;
    using htsim::model::is_haplotype_mask;
    using htsim::model::mask_contains;

    require(!is_haplotype_mask(0), "zero is not a haplotype mask");
    require(is_haplotype_mask(1), "haplotype 1 mask is valid");
    require(is_haplotype_mask(2), "haplotype 2 mask is valid");
    require(is_haplotype_mask(3), "shared mask is valid");
    require(!is_haplotype_mask(4), "reserved mask is invalid");
    require(
        mask_contains(HaplotypeMask::haplotype_1, 0)
            && !mask_contains(HaplotypeMask::haplotype_1, 1),
        "haplotype 1 membership changed");
    require(
        mask_contains(HaplotypeMask::haplotype_2, 1)
            && !mask_contains(HaplotypeMask::haplotype_2, 0),
        "haplotype 2 membership changed");
    require(
        mask_contains(HaplotypeMask::both, 0)
            && mask_contains(HaplotypeMask::both, 1)
            && !mask_contains(HaplotypeMask::both, 2),
        "shared haplotype membership changed");
    require(
        htsim::model::maximum_insertion_bases == 4,
        "insertion representation limit changed");
    std::cout << "haplotype model tests passed\n";
    return 0;
}
