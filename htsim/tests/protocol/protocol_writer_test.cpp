#include "protocol.h"

#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>

namespace {

using namespace htsim::protocol;

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

template <typename Callable>
void require_protocol_error(Callable &&callable, const std::string &message)
{
    try {
        callable();
    } catch (const ProtocolError &) {
        return;
    }
    throw std::runtime_error(message);
}

std::uint32_t load_le32(const std::string &bytes, std::size_t offset)
{
    require(offset + 4U <= bytes.size(), "load_le32 offset is out of range");
    std::uint32_t result = 0;
    for (unsigned int byte = 0; byte < 4; ++byte) {
        result |= static_cast<std::uint32_t>(
                      static_cast<unsigned char>(bytes[offset + byte]))
            << (8U * byte);
    }
    return result;
}

std::uint64_t load_le64(const std::string &bytes, std::size_t offset)
{
    require(offset + 8U <= bytes.size(), "load_le64 offset is out of range");
    std::uint64_t result = 0;
    for (unsigned int byte = 0; byte < 8; ++byte) {
        result |= static_cast<std::uint64_t>(
                      static_cast<unsigned char>(bytes[offset + byte]))
            << (8U * byte);
    }
    return result;
}

Digest sha256_text(const std::string &value)
{
    return htsim::crypto::sha256(
        reinterpret_cast<const std::uint8_t *>(value.data()), value.size());
}

Header make_header(bool details)
{
    Header header;
    header.run_id = "00000000-0000-4000-8000-000000000002";
    header.core_version = "2.0.0-alpha.1";
    header.rng_contract = std::string(rng_contract);
    header.master_seed = UINT64_C(0x0123456789abcdef);
    header.normalized_config_sha256 = sha256_text("protocol-config");
    header.technology = Technology::wgbs;
    header.has_details = details;
    header.mates_per_fragment = 1;
    header.base_encoding = BaseEncoding::acgtn_u8;
    header.ambiguity_policy = AmbiguityPolicy::preserve_n;
    header.read_length_r1 = 3;
    header.read_length_r2 = 0;
    std::string reference;
    for (unsigned int index = 0; index < 25; ++index) {reference += "ACGT";}
    header.contigs = {{"chrMini", 100, sha256_text(reference)}};
    return header;
}

FragmentBatch make_no_annotation_batch()
{
    FragmentBatch batch;
    batch.first_fragment_ordinal = 0;
    batch.contig_indices = {0, 0};
    batch.reference_starts = {10, 20};
    batch.reference_ends = {14, 24};
    batch.template_offsets = {0, 4, 8};
    batch.mate_offsets = {0, 1, 2};
    batch.site_offsets = {0, 1, 2};
    batch.mate_template_starts = {0, 1};
    batch.mate_template_ends = {3, 4};
    batch.site_template_offsets = {1, 1};
    batch.site_probabilities = {0.25F, 0.75F};
    batch.haplotypes = {0, 1};
    batch.capture_strands = {
        static_cast<std::uint8_t>(CaptureStrand::unknown),
        static_cast<std::uint8_t>(CaptureStrand::forward),
    };
    batch.mate_indices = {0, 0};
    batch.mate_reverse_complements = {0, 0};
    batch.site_contexts = {
        static_cast<std::uint8_t>(MethylationContext::cg_c),
        static_cast<std::uint8_t>(MethylationContext::cg_g),
    };
    batch.methylation_sources = {
        static_cast<std::uint8_t>(MethylationSource::beta),
        static_cast<std::uint8_t>(MethylationSource::cgmap),
    };
    batch.site_alleles = {
        static_cast<std::uint8_t>(MethylationAllele::shared),
        static_cast<std::uint8_t>(MethylationAllele::alternate_haplotype),
    };
    batch.template_bases = {0, 1, 2, 3, 3, 2, 1, 0};
    return batch;
}

FragmentBatch make_full_annotation_batch()
{
    FragmentBatch batch;
    batch.first_fragment_ordinal = 0;
    batch.contig_indices = {0, 0, 0};
    batch.reference_starts = {10, 20, 30};
    batch.reference_ends = {14, 24, 34};
    batch.template_offsets = {0, 4, 9, 12};
    batch.mate_offsets = {0, 1, 2, 3};
    batch.site_offsets = {0, 1, 3, 4};
    batch.mate_template_starts = {0, 0, 0};
    batch.mate_template_ends = {3, 3, 3};
    batch.site_template_offsets = {2, 1, 2, 1};
    batch.site_probabilities = {0.5F, 0.25F, 0.75F, 1.0F};
    batch.haplotypes = {1, 1, 1};
    batch.capture_strands = {0, 0, 0};
    batch.mate_indices = {0, 0, 0};
    batch.mate_reverse_complements = {0, 0, 0};
    batch.site_contexts = {
        static_cast<std::uint8_t>(MethylationContext::cg_g),
        static_cast<std::uint8_t>(MethylationContext::cg_c),
        static_cast<std::uint8_t>(MethylationContext::cg_g),
        static_cast<std::uint8_t>(MethylationContext::cg_c),
    };
    batch.methylation_sources.assign(
        4, static_cast<std::uint8_t>(MethylationSource::beta));
    batch.site_alleles = {
        static_cast<std::uint8_t>(MethylationAllele::alternate_haplotype),
        static_cast<std::uint8_t>(MethylationAllele::shared),
        static_cast<std::uint8_t>(MethylationAllele::alternate_haplotype),
        static_cast<std::uint8_t>(MethylationAllele::shared),
    };
    batch.template_bases = {0, 3, 2, 3, 0, 1, 2, 3, 0, 4, 1, 3};

    FragmentDetails details;
    details.projection_offsets = {0, 1, 3, 5};
    details.variant_offsets = {0, 1, 2, 3};
    details.original_n_offsets = {0, 0, 0, 1};
    details.projection_template_starts = {0, 0, 3, 0, 2};
    details.projection_template_ends = {4, 2, 5, 2, 3};
    details.projection_reference_starts = {10, 20, 22, 30, 33};
    details.variant_indices = {1, 2, 3};
    details.variant_id_offsets = {0, 2, 4, 6};
    details.variant_reference_starts = {11, 22, 32};
    details.variant_reference_ends = {12, 22, 33};
    details.variant_template_starts = {1, 2, 2};
    details.variant_template_ends = {2, 3, 2};
    details.variant_ref_offsets = {0, 1, 1, 2};
    details.variant_alt_offsets = {0, 1, 2, 2};
    details.site_reference_positions = {12, 21, no_reference_position, 31};
    details.original_n_template_offsets = {0};
    details.variant_kinds = {
        static_cast<std::uint8_t>(VariantKind::snv),
        static_cast<std::uint8_t>(VariantKind::insertion),
        static_cast<std::uint8_t>(VariantKind::deletion),
    };
    details.variant_sources = {
        static_cast<std::uint8_t>(VariantSource::vcf),
        static_cast<std::uint8_t>(VariantSource::vcf),
        static_cast<std::uint8_t>(VariantSource::de_novo),
    };
    details.variant_phased_haplotypes = {1, 1, 1};
    details.variant_ids = {'v', '1', 'v', '2', 'v', '3'};
    details.variant_ref_bases = {1, 2};
    details.variant_alt_bases = {3, 2};
    batch.details = std::move(details);
    return batch;
}

std::pair<std::string, Trailer> make_no_annotation_stream(bool prepared)
{
    std::ostringstream output(std::ios::out | std::ios::binary);
    Writer writer(output);
    const Header header = make_header(false);
    writer.write_header(header);
    if (prepared) {
        writer.write_prepared_batch(
            prepare_fragment_batch(header, make_no_annotation_batch()));
    } else {
        writer.write_batch(make_no_annotation_batch());
    }
    const Trailer trailer = writer.finish(7);
    require(writer.complete() && !writer.failed(), "valid writer state is invalid");
    return {output.str(), trailer};
}

std::string make_full_annotation_stream()
{
    std::ostringstream output(std::ios::out | std::ios::binary);
    Writer writer(output);
    writer.write_header(make_header(true));
    writer.write_batch(make_full_annotation_batch());
    writer.finish();
    return output.str();
}

std::string make_error_stream()
{
    std::ostringstream output(std::ios::out | std::ios::binary);
    Writer writer(output);
    writer.write_header(make_header(false));
    writer.write_error({1204, "batch exceeds limit"});
    return output.str();
}

void test_no_annotation_envelope_and_counts()
{
    const auto result = make_no_annotation_stream(false);
    const std::string &stream = result.first;
    const Trailer &trailer = result.second;
    require(stream.size() == 500U, "no-Details stream size changed");
    require(stream.substr(0, 8) == std::string("BSRSTRM\0", 8), "magic changed");
    require(static_cast<unsigned char>(stream[8]) == 2U,
            "protocol major version changed");
    const std::size_t batch_start = 232U;
    require(load_le32(stream, batch_start) == 136U, "batch payload size changed");
    require(static_cast<unsigned char>(stream[batch_start + 4U]) == 2U,
            "batch frame type changed");
    require(load_le64(stream, batch_start + 8U) == 1U,
            "batch frame sequence changed");
    const std::size_t batch_payload = batch_start + 16U;
    require(load_le32(stream, batch_payload) == 0U
                && load_le32(stream, batch_payload + 4U) == 2U
                && load_le32(stream, batch_payload + 8U) == 8U
                && load_le32(stream, batch_payload + 12U) == 2U
                && load_le32(stream, batch_payload + 16U) == 2U,
            "batch authenticated counts changed");
    require(trailer.fragment_count == 2U
                && trailer.fragment_batch_count == 1U
                && trailer.mate_count == 2U
                && trailer.template_base_count == 8U
                && trailer.methylation_site_count == 2U
                && trailer.skipped_fragment_count == 7U
                && trailer.per_contig_fragment_counts
                    == std::vector<std::uint64_t>{2U},
            "trailer aggregates changed");
}

void test_prepared_batch_is_byte_identical()
{
    const auto scalar = make_no_annotation_stream(false);
    const auto prepared = make_no_annotation_stream(true);
    require(scalar.first == prepared.first, "prepared batch changed stream bytes");
    require(scalar.second.stream_sha256 == prepared.second.stream_sha256,
            "prepared batch changed stream digest");
}

void test_full_annotation_shape()
{
    const std::string stream = make_full_annotation_stream();
    require(stream.size() == 840U,
            "Full-Details stream size changed: " + std::to_string(stream.size()));
    require(load_le32(stream, 232U) == 476U,
            "Full-Details payload size changed: "
                + std::to_string(load_le32(stream, 232U)));
    require(static_cast<unsigned char>(stream[237U]) == details_present,
            "Full-Details frame flag changed");
}

void test_invalid_header_writes_zero_bytes_and_poisons()
{
    std::ostringstream output(std::ios::out | std::ios::binary);
    Writer writer(output);
    Header invalid = make_header(false);
    invalid.read_length_r1 = 0;
    require_protocol_error(
        [&] {writer.write_header(invalid);},
        "invalid header was accepted");
    require(output.str().empty(), "invalid header emitted a preamble");
    require(writer.failed(), "invalid header did not poison writer");
    require_protocol_error(
        [&] {writer.write_header(make_header(false));},
        "poisoned writer accepted a header");
}

void test_invalid_batch_and_order_poison_writer()
{
    {
        std::ostringstream output(std::ios::out | std::ios::binary);
        Writer writer(output);
        writer.write_header(make_header(false));
        FragmentBatch invalid = make_no_annotation_batch();
        invalid.template_offsets.front() = 1;
        require_protocol_error(
            [&] {writer.write_batch(std::move(invalid));},
            "invalid prefix was accepted");
        require(writer.failed(), "invalid batch did not poison writer");
    }
    {
        std::ostringstream output(std::ios::out | std::ios::binary);
        Writer writer(output);
        const Header header = make_header(false);
        writer.write_header(header);
        FragmentBatch out_of_order = make_no_annotation_batch();
        out_of_order.first_fragment_ordinal = 1;
        PreparedFragmentBatch prepared =
            prepare_fragment_batch(header, std::move(out_of_order));
        require_protocol_error(
            [&] {writer.write_prepared_batch(std::move(prepared));},
            "out-of-order prepared batch was accepted");
        require(writer.failed(), "out-of-order batch did not poison writer");
    }
}

void test_prepared_header_mismatch_fails_closed()
{
    Header prepared_header = make_header(false);
    prepared_header.master_seed += 1U;
    PreparedFragmentBatch prepared =
        prepare_fragment_batch(prepared_header, make_no_annotation_batch());
    std::ostringstream output(std::ios::out | std::ios::binary);
    Writer writer(output);
    writer.write_header(make_header(false));
    require_protocol_error(
        [&] {writer.write_prepared_batch(std::move(prepared));},
        "prepared batch for another header was accepted");
}

void test_sparse_annotation_failures_are_rejected()
{
    FragmentBatch invalid = make_full_annotation_batch();
    invalid.details->projection_template_ends[1] = 3;
    require_protocol_error(
        [&] {(void)prepare_fragment_batch(make_header(true), std::move(invalid));},
        "overlapping details projection was accepted");

    invalid = make_full_annotation_batch();
    invalid.details->variant_reference_starts[1] = 21;
    invalid.details->variant_reference_ends[1] = 21;
    require_protocol_error(
        [&] {(void)prepare_fragment_batch(make_header(true), std::move(invalid));},
        "incorrect insertion anchor was accepted");

    invalid = make_full_annotation_batch();
    invalid.details->original_n_offsets = {0, 0, 0, 0};
    require_protocol_error(
        [&] {(void)prepare_fragment_batch(make_header(true), std::move(invalid));},
        "incomplete N provenance was accepted");
}

void test_terminal_writer_rejects_reuse()
{
    std::ostringstream output(std::ios::out | std::ios::binary);
    Writer writer(output);
    writer.write_header(make_header(false));
    writer.finish();
    require_protocol_error(
        [&] {writer.write_batch(make_no_annotation_batch());},
        "complete writer accepted another batch");
}

void write_binary(const std::string &path, const std::string &bytes)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {throw std::runtime_error("failed to open fixture output");}
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    if (!output) {throw std::runtime_error("failed to write fixture output");}
}

void write_fixtures(const std::string &prefix)
{
    write_binary(prefix + "-none.bin", make_no_annotation_stream(false).first);
    write_binary(prefix + "-full.bin", make_full_annotation_stream());
    write_binary(prefix + "-error.bin", make_error_stream());
}

} // namespace

int main(int argc, char *argv[])
{
    try {
        test_no_annotation_envelope_and_counts();
        test_prepared_batch_is_byte_identical();
        test_full_annotation_shape();
        test_invalid_header_writes_zero_bytes_and_poisons();
        test_invalid_batch_and_order_poison_writer();
        test_prepared_header_mismatch_fails_closed();
        test_sparse_annotation_failures_are_rejected();
        test_terminal_writer_rejects_reuse();
        if (argc == 3 && std::string(argv[1]) == "--write-fixtures") {
            write_fixtures(argv[2]);
        } else if (argc != 1) {
            throw std::runtime_error(
                "usage: protocol_writer_test [--write-fixtures PREFIX]");
        }
    } catch (const std::exception &error) {
        std::cerr << "protocol_writer_test: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
