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

htsim::protocol::Header make_header(
    htsim::protocol::TruthMode truth_mode)
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
    header.truth_columns = truth_mode;
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
    fragment.base_event_ids = {
        no_variant_event, 2U, no_variant_event, no_variant_event,
        5U, no_variant_event, no_variant_event,
    };
    fragment.variant_events = {
        {2U, VariantKind::snv, 1U, 101U, 102U, {1U}, {3U}},
        {5U, VariantKind::insertion, 255U, 104U, 104U, {}, {2U}},
        {7U, VariantKind::deletion, 1U, 105U, 107U, {1U, 2U}, {}},
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
    fragment.base_event_ids.clear();
    fragment.variant_events.clear();
    for (auto &mate : fragment.mates) {
        mate.reference_start = 0U;
        mate.reference_end = 0U;
        mate.site_refs.clear();
    }
    return fragment;
}

std::string encode_no_truth(const htsim::model::Fragment &fragment)
{
    using namespace htsim::protocol;
    const Header header = make_header(TruthMode::none);
    std::ostringstream sink(std::ios::binary);
    Writer writer(sink);
    writer.write_header(header);
    writer.write_batch(make_fragment_batch(header, {fragment}));
    (void)writer.finish();
    return sink.str();
}

void test_full_truth_projection_is_canonical()
{
    using namespace htsim::protocol;
    const Header header = make_header(TruthMode::full);
    FragmentBatch batch = make_fragment_batch(header, {make_fragment()});
    require(batch.first_fragment_ordinal == 0U,
            "first fragment ordinal changed");
    require(batch.template_offsets == std::vector<std::uint32_t>({0U, 7U}),
            "template prefix changed");
    require(batch.mate_offsets == std::vector<std::uint32_t>({0U, 2U}),
            "mate prefix changed");
    require(batch.site_offsets == std::vector<std::uint32_t>({0U, 1U}),
            "site prefix changed");
    require(batch.truth.has_value(), "Full Truth columns were omitted");
    const TruthColumns &truth = *batch.truth;
    require(
        truth.projection_offsets == std::vector<std::uint32_t>({0U, 3U})
        && truth.projection_template_begins
            == std::vector<std::uint32_t>({0U, 5U, 6U})
        && truth.projection_template_ends
            == std::vector<std::uint32_t>({4U, 6U, 7U})
        && truth.projection_reference_begins
            == std::vector<std::uint32_t>({100U, 104U, 107U}),
        "dense projection was not compressed into maximal runs");
    require(
        truth.event_ids == std::vector<std::uint32_t>({2U, 5U, 7U})
        && truth.event_template_begins
            == std::vector<std::uint32_t>({1U, 4U, 6U})
        && truth.event_template_ends
            == std::vector<std::uint32_t>({2U, 5U, 6U}),
        "typed events did not retain physical template spans");
    require(
        truth.site_reference_positions
            == std::vector<std::uint32_t>({no_reference_position}),
        "inserted methylation site lost its reference sentinel");
    require(
        truth.original_n_offsets == std::vector<std::uint32_t>({0U, 1U})
        && truth.original_n_template_offsets
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

void test_no_truth_omits_only_sparse_columns()
{
    using namespace htsim::protocol;
    const Header header = make_header(TruthMode::none);
    const auto full_fragment = make_fragment();
    const auto compact_fragment = make_compact_fragment();
    FragmentBatch batch = make_fragment_batch(header, {compact_fragment});
    require(!batch.truth.has_value(), "no-Truth batch retained Truth columns");
    require(batch.template_bases == full_fragment.template_bases,
            "no-Truth mode changed common template bases");
    require(encode_no_truth(full_fragment) == encode_no_truth(compact_fragment),
            "compact fragment changed columnar protocol common wire bytes");
    std::ostringstream sink(std::ios::binary);
    Writer writer(sink);
    writer.write_header(header);
    writer.write_batch(std::move(batch));
    require(writer.finish().fragment_count == 1U,
            "no-Truth adapter output was rejected");
}

void test_adapter_rejects_noncanonical_typed_input()
{
    const auto header = make_header(htsim::protocol::TruthMode::full);
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
    fragment.base_event_ids[4] = htsim::model::no_variant_event;
    require_error(
        [&] {make_fragment_batch(header, {fragment});},
        "incomplete insertion span was accepted");

    const auto compact = make_compact_fragment();
    require_error(
        [&] {make_fragment_batch(header, {compact});},
        "compact fragment was accepted for Full Truth");

    const auto no_truth_header = make_header(
        htsim::protocol::TruthMode::none);
    fragment = make_compact_fragment();
    fragment.reference_positions = make_fragment().reference_positions;
    require_error(
        [&] {make_fragment_batch(no_truth_header, {fragment});},
        "partially compact projection arrays were accepted");

    fragment = make_compact_fragment();
    fragment.variant_events = make_fragment().variant_events;
    require_error(
        [&] {make_fragment_batch(no_truth_header, {fragment});},
        "compact fragment retaining variant events was accepted");
}

std::string encode_through_emitter(std::uint32_t workers)
{
    using namespace htsim::protocol;
    const Header header = make_header(TruthMode::none);
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

    const Header header = make_header(TruthMode::none);
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
        test_full_truth_projection_is_canonical();
        test_no_truth_omits_only_sparse_columns();
        test_adapter_rejects_noncanonical_typed_input();
        test_batch_emitter_preserves_order();
    } catch (const std::exception &error) {
        std::cerr << "protocol_adapter_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
