#include "protocol.h"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using htsim::protocol::make_fragment_batch;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Callable>
void require_error(Callable &&callable, const std::string &message)
{
    try {
        callable();
    } catch (const htsim::protocol::ProtocolError &) {
        return;
    }
    throw std::runtime_error(message);
}

htsim::protocol::Header make_header(bool has_details)
{
    using namespace htsim::protocol;
    Header header;
    header.run_id = "12345678-1234-4234-8234-123456789abc";
    header.core_version = "0.2.0";
    header.config_schema_version = std::string(config_schema_version);
    header.rng_contract = std::string(rng_contract);
    header.master_seed = 9U;
    header.normalized_config_sha256.fill(0x11U);
    header.technology = Technology::wgbs;
    header.has_details = has_details;
    header.mates_per_fragment = 2U;
    header.base_encoding = BaseEncoding::acgtn_u8;
    header.ambiguity_policy = AmbiguityPolicy::preserve_n;
    header.read_length_r1 = 3U;
    header.read_length_r2 = 3U;
    Contig contig{"chr1", 200U, {}};
    contig.reference_sha256.fill(0x22U);
    header.contigs.push_back(contig);
    return header;
}

htsim::model::Fragment make_fragment(std::uint64_t ordinal = 0U)
{
    using namespace htsim::model;
    Fragment fragment;
    fragment.fragment_ordinal = ordinal;
    fragment.contig_index = 0U;
    fragment.haplotype = 1U;
    fragment.capture_strand = CaptureStrand::reverse;
    fragment.reference_start = 100U;
    fragment.reference_end = 108U;
    fragment.template_bases = {0U, 3U, 4U, 3U, 2U, 0U, 3U};
    fragment.reference_positions = {100, 101, 102, 103, -1, 104, 107};
    fragment.base_variant_indices = {
        no_variant_index, 2U, no_variant_index, no_variant_index,
        5U, no_variant_index, no_variant_index,
    };
    fragment.variants = {
        {2U, "v2", VariantSource::vcf,
         VariantKind::snv, 1U, 101U, 102U, {1U}, {3U}},
        {5U, "v5", VariantSource::vcf,
         VariantKind::insertion, 255U, 104U, 104U, {}, {2U}},
        {7U, "v7", VariantSource::vcf,
         VariantKind::deletion, 1U, 105U, 107U, {1U, 2U}, {}},
    };
    fragment.methylation_sites = {
        {0U, 4U, -1, MethylationContext::cg_g,
         MethylationSource::beta, MethylationAllele::shared, 0.25F},
    };
    fragment.mates = {
        {0U, false, 0U, 3U, 100U, 103U, {}},
        {1U, true, 4U, 7U, 104U, 108U, {{2U, 0U}}},
    };
    return fragment;
}

htsim::model::Fragment make_compact_fragment(
    std::uint64_t ordinal = 0U)
{
    auto fragment = make_fragment(ordinal);
    fragment.reference_positions.clear();
    fragment.base_variant_indices.clear();
    fragment.variants.clear();
    for (auto &mate : fragment.mates) {
        mate.reference_start = 0U;
        mate.reference_end = 0U;
        mate.site_refs.clear();
    }
    return fragment;
}

std::string encode_no_annotations(const htsim::model::Fragment &fragment)
{
    using namespace htsim::protocol;
    const Header header = make_header(false);
    std::ostringstream sink(std::ios::binary);
    Writer writer(sink);
    writer.write_header(header);
    writer.write_batch(make_fragment_batch(header, {fragment}));
    (void)writer.finish();
    return sink.str();
}

void test_full_annotation_projection_is_canonical()
{
    using namespace htsim::protocol;
    const Header header = make_header(true);
    FragmentBatch batch = make_fragment_batch(header, {make_fragment()});
    require(batch.first_fragment_ordinal == 0U,
            "first fragment ordinal changed");
    require(batch.template_offsets == std::vector<std::uint32_t>({0U, 7U}),
            "template prefix changed");
    require(batch.mate_offsets == std::vector<std::uint32_t>({0U, 2U}),
            "mate prefix changed");
    require(batch.site_offsets == std::vector<std::uint32_t>({0U, 1U}),
            "site prefix changed");
    require(batch.details.has_value(), "Full Details columns were omitted");
    const FragmentDetails &details = *batch.details;
    require(
        details.projection_offsets == std::vector<std::uint32_t>({0U, 3U})
        && details.projection_template_starts
            == std::vector<std::uint32_t>({0U, 5U, 6U})
        && details.projection_template_ends
            == std::vector<std::uint32_t>({4U, 6U, 7U})
        && details.projection_reference_starts
            == std::vector<std::uint32_t>({100U, 104U, 107U}),
        "dense projection was not compressed into maximal runs");
    require(
        details.variant_indices == std::vector<std::uint32_t>({2U, 5U, 7U})
        && details.variant_template_starts
            == std::vector<std::uint32_t>({1U, 4U, 6U})
        && details.variant_template_ends
            == std::vector<std::uint32_t>({2U, 5U, 6U}),
        "typed events did not retain physical template spans");
    require(
        details.site_reference_positions
            == std::vector<std::uint32_t>({no_reference_position}),
        "inserted methylation site lost its reference sentinel");
    require(
        details.original_n_offsets == std::vector<std::uint32_t>({0U, 1U})
        && details.original_n_template_offsets
            == std::vector<std::uint32_t>({2U}),
        "preserved-N provenance changed");

    std::ostringstream sink(std::ios::binary);
    Writer writer(sink);
    writer.write_header(header);
    writer.write_batch(std::move(batch));
    const Trailer trailer = writer.finish();
    require(trailer.fragment_count == 1U
                && trailer.fragment_batch_count == 1U
                && trailer.mate_count == 2U,
            "adapter output failed writer aggregate accounting");
}

void test_no_annotation_omits_only_sparse_columns()
{
    using namespace htsim::protocol;
    const Header header = make_header(false);
    const auto full_fragment = make_fragment();
    const auto compact_fragment = make_compact_fragment();
    FragmentBatch batch = make_fragment_batch(header, {compact_fragment});
    require(!batch.details.has_value(), "no-Details batch retained Details columns");
    require(batch.template_bases == full_fragment.template_bases,
            "no-Details mode changed common template bases");
    require(encode_no_annotations(full_fragment) == encode_no_annotations(compact_fragment),
            "compact fragment changed columnar protocol common wire bytes");
    std::ostringstream sink(std::ios::binary);
    Writer writer(sink);
    writer.write_header(header);
    writer.write_batch(std::move(batch));
    require(writer.finish().fragment_count == 1U,
            "no-Details adapter output was rejected");
}

void test_adapter_rejects_noncanonical_typed_input()
{
    const auto header = make_header(true);
    auto fragment = make_fragment();
    fragment.reference_positions[1] = 100;
    require_error(
        [&] {make_fragment_batch(header, {fragment});},
        "non-increasing dense positions were accepted");

    auto second = make_fragment(2U);
    require_error(
        [&] {make_fragment_batch(header, {make_fragment(), second});},
        "non-consecutive fragment ordinals were accepted");

    fragment = make_fragment();
    fragment.base_variant_indices[4] = htsim::model::no_variant_index;
    require_error(
        [&] {make_fragment_batch(header, {fragment});},
        "incomplete insertion span was accepted");

    const auto compact = make_compact_fragment();
    require_error(
        [&] {make_fragment_batch(header, {compact});},
        "compact fragment was accepted for Full Details");

    const auto without_annotations_header = make_header(false);
    fragment = make_compact_fragment();
    fragment.reference_positions = make_fragment().reference_positions;
    require_error(
        [&] {make_fragment_batch(without_annotations_header, {fragment});},
        "partially compact projection arrays were accepted");

    fragment = make_compact_fragment();
    fragment.variants = make_fragment().variants;
    require_error(
        [&] {make_fragment_batch(without_annotations_header, {fragment});},
        "compact fragment retaining variant events was accepted");
}

std::string encode_through_emitter(std::uint32_t workers)
{
    using namespace htsim::protocol;
    const Header header = make_header(false);
    std::ostringstream sink(std::ios::binary);
    Writer writer(sink);
    writer.write_header(header);
    {
        BatchEmitter emitter(writer, header, workers, 2U);
        for (std::uint64_t ordinal = 0U; ordinal < 5U; ++ordinal) {
            emitter.write(make_compact_fragment(ordinal));
        }
        emitter.finish();
    }
    const Trailer trailer = writer.finish();
    require(
        trailer.fragment_count == 5U
            && trailer.fragment_batch_count == 3U
            && trailer.mate_count == 10U,
        "batch emitter aggregate accounting changed");
    return sink.str();
}

void test_batch_emitter_preserves_order()
{
    using namespace htsim::protocol;
    require(
        encode_through_emitter(1U) == encode_through_emitter(3U),
        "parallel batch encoding changed protocol bytes");

    const Header header = make_header(false);
    std::ostringstream sink(std::ios::binary);
    Writer writer(sink);
    writer.write_header(header);
    require_error(
        [&] {BatchEmitter invalid(writer, header, 0U, 1U);},
        "batch emitter accepted zero workers");
    BatchEmitter emitter(writer, header, 1U, 1U);
    emitter.write(make_compact_fragment());
    emitter.finish();
    require_error(
        [&] {emitter.write(make_compact_fragment(1U));},
        "batch emitter accepted a fragment after finish");
    require_error(
        [&] {emitter.finish();},
        "batch emitter accepted a second finish");
}

} // namespace

int main()
{
    try {
        test_full_annotation_projection_is_canonical();
        test_no_annotation_omits_only_sparse_columns();
        test_adapter_rejects_noncanonical_typed_input();
        test_batch_emitter_preserves_order();
    } catch (const std::exception &error) {
        std::cerr << "protocol_adapter_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
