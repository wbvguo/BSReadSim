#ifndef HTSIM_FRAGMENT_H
#define HTSIM_FRAGMENT_H

#include <cstdint>
#include <stdexcept>
#include <vector>
#include <string_view>

#include "variant.h"
#include "methdb.h"
#include "types.h"
#include "reference.h"

// ---- allocation --------------------------------------------------------

namespace htsim::allocation {

class AllocationError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Allocate an exact positive number of fragments in input order. Each weight
// is a non-negative uint32 measure owned by the caller (for example eligible
// starts or an exact target-score sum) for the corresponding FASTA contig.
// Hamilton's largest-remainder method is used: floor(total * weight / sum)
// first, followed by descending remainder and ascending contig-index ties.
// Zero-weight entries are never eligible for a remainder seat.
std::vector<std::uint32_t> largest_remainder(
    const std::vector<std::uint32_t> &weights,
    std::uint32_t total_fragments);

// Floating weights are used only for reference-calibrated target profiles.
// The same Hamilton allocation is applied after validating a finite positive
// total. Ties remain stable by input index.
std::vector<std::uint32_t> largest_remainder_real(
    const std::vector<double> &weights,
    std::uint32_t total_fragments);

} // namespace htsim::allocation

namespace htsim::depth_count {

class DepthCountError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Convert mean sequencing depth to an exact uint32 fragment count by first
// resolving the required number of reads, then rounding up to a complete
// fragment bundle. The frozen evaluation order is:
//   raw_reads = (double(effective_reference_bases) * depth) / read_length
//   ceil(raw_reads / (paired_end ? 2 : 1))
// under round-to-nearest floating point.
std::uint32_t fragments(
    double depth,
    std::uint64_t effective_reference_bases,
    std::uint32_t read_length,
    bool paired_end);

} // namespace htsim::depth_count

// ---- insert_length --------------------------------------------------------

namespace htsim::insert_length {

class InsertLengthError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct Parameters {
    std::uint32_t minimum = 0;
    std::uint32_t mean = 0;
    std::uint32_t maximum = 0;
    double standard_deviation = 0.0;
};

class Sampler {
public:
    Sampler(
        std::uint64_t master_seed,
        std::uint32_t contig_index,
        const Parameters &parameters);

    std::uint32_t sample(std::uint64_t candidate_ordinal) const;

private:
    std::uint64_t key_ = 0;
    Parameters parameters_;
};

// Addressed clamped-normal insert length for one per-contig fragment candidate.
// RNG stage=Stage::fragment, entity=candidate_ordinal, local_index=0. The returned
// length and every configured boundary are uint32; the ordinal/key are uint64.
std::uint32_t sample(
    std::uint64_t master_seed,
    std::uint32_t contig_index,
    std::uint64_t candidate_ordinal,
    const Parameters &parameters);

} // namespace htsim::insert_length

// ---- fragment_builder --------------------------------------------------------

namespace htsim::fragment_builder {

class FragmentBuilderError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

struct ReadLayout {
    enum class InsertCoordinate : std::uint8_t {
        reference = 0,
        haplotype = 1,
    };

    std::uint32_t insert_length = 0;
    std::uint32_t read_length = 0;
    bool paired_end = false;
    InsertCoordinate insert_coordinate = InsertCoordinate::reference;
};

// Select the typed state retained after fragment construction. Full Details
// requires the complete projection. Common columns need only the template,
// methylation rows, and structural mate slices; omitting the remaining state
// avoids constructing data that the production wire path does not carry.
enum class FragmentDetail : std::uint8_t {
    full = 0,
    common_columns = 1,
};

// Exact upper bound for a variant-free protocol fragment payload when every
// template base is methylatable and every mate base references a site.
std::uint64_t maximum_payload_bytes(const ReadLayout &layout);
void require_payload_fits_protocol(const ReadLayout &layout);

// Exact upper bound for an already projected variant-aware interval. The
// bound includes every typed variant event and assumes every projected base is
// methylatable. Coordinates/counts remain uint32 at this component boundary;
// only the byte-size arithmetic is widened to uint64.
std::uint64_t maximum_payload_bytes(
    const haplotype::ProjectedInterval &projection,
    const ReadLayout &layout);
void require_payload_fits_protocol(
    const haplotype::ProjectedInterval &projection,
    const ReadLayout &layout);

// Pure construction boundary. The caller owns sampling, haplotype RNG,
// fragment ordinal assignment, protocol writing, and all I/O.
model::Fragment build_fragment(
    const reference::Contig &contig,
    const methdb::MethylationCatalog &catalog,
    std::uint64_t fragment_ordinal,
    std::uint32_t reference_start,
    std::uint8_t haplotype,
    model::CaptureStrand capture_strand,
    const ReadLayout &layout,
    FragmentDetail detail = FragmentDetail::full,
    bool reverse_molecule = false);

// Consume one typed haplotype projection and attach its diploid methylation
// sites. This boundary performs no sampling or I/O. Passing by value allows a
// temporary projection to transfer its potentially large arrays without a
// copy.
model::Fragment build_fragment(
    haplotype::ProjectedInterval projection,
    const methdb::DiploidMethylationCatalog &catalog,
    std::uint64_t fragment_ordinal,
    model::CaptureStrand capture_strand,
    const ReadLayout &layout,
    FragmentDetail detail = FragmentDetail::full,
    bool reverse_molecule = false);

} // namespace htsim::fragment_builder

#endif // HTSIM_FRAGMENT_H
