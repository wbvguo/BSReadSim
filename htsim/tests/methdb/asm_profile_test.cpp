#include "methdb.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <unistd.h>
#include <zlib.h>

#include "protocol.h"
#include "reference.h"
#include "utilities.h"

namespace {

using htsim::methdb::AsmProfile;
using htsim::methdb::AsmProfileError;
using htsim::methdb::AsmProfileFormat;
using htsim::methdb::AsmRecord;
using htsim::model::Bases;
using htsim::model::MethylationContext;
using htsim::reference::Contig;
using htsim::reference::ContigMetadata;

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-asm-profile-XXXXXX";
        const int descriptor = mkstemp(pattern);
        if (descriptor < 0) {throw std::runtime_error("mkstemp failed");}
        path_ = pattern;
        if (close(descriptor) != 0) {
            throw std::runtime_error("temporary descriptor close failed");
        }
    }

    ~TempFile() {if (!path_.empty()) {(void)unlink(path_.c_str());}}
    const std::string &path() const noexcept {return path_;}

private:
    std::string path_;
};

void require(bool condition, const std::string &message)
{
    if (!condition) {throw std::runtime_error(message);}
}

std::uint16_t q(float probability)
{
    return htsim::methdb::probability_to_u16(probability);
}

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const AsmProfileError &) {
        return;
    }
    throw std::runtime_error(message);
}

std::vector<std::uint8_t> bytes_of(const std::string &text)
{
    return {text.begin(), text.end()};
}

void write_bytes(const std::string &path, const std::vector<std::uint8_t> &bytes)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output.write(
        reinterpret_cast<const char *>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
    if (!output) {throw std::runtime_error("temporary file write failed");}
}

std::vector<std::uint8_t> gzip_bytes(const std::string &text)
{
    z_stream stream = {};
    if (deflateInit2(
            &stream,
            Z_BEST_SPEED,
            Z_DEFLATED,
            MAX_WBITS + 16,
            8,
            Z_DEFAULT_STRATEGY) != Z_OK) {
        throw std::runtime_error("deflateInit2 failed");
    }
    std::vector<std::uint8_t> output(128U + text.size());
    stream.next_in = reinterpret_cast<Bytef *>(
        const_cast<char *>(text.data()));
    stream.avail_in = static_cast<uInt>(text.size());
    stream.next_out = output.data();
    stream.avail_out = static_cast<uInt>(output.size());
    const int status = deflate(&stream, Z_FINISH);
    if (status != Z_STREAM_END) {
        (void)deflateEnd(&stream);
        throw std::runtime_error("deflate failed");
    }
    output.resize(output.size() - stream.avail_out);
    if (deflateEnd(&stream) != Z_OK) {
        throw std::runtime_error("deflateEnd failed");
    }
    return output;
}

Bases encode(const std::string &sequence)
{
    Bases result;
    result.reserve(sequence.size());
    for (const char base : sequence) {
        switch (base) {
        case 'A': result.push_back(0U); break;
        case 'C': result.push_back(1U); break;
        case 'G': result.push_back(2U); break;
        case 'T': result.push_back(3U); break;
        case 'N': result.push_back(4U); break;
        default: throw std::runtime_error("invalid fixture base");
        }
    }
    return result;
}

htsim::crypto::Sha256Digest sequence_digest(const std::string &sequence)
{
    return htsim::crypto::sha256(bytes_of(sequence));
}

std::vector<ContigMetadata> reference_catalog()
{
    return {
        {"chr1", 13U, sequence_digest("ACGACAGCAATCG")},
        {"chr2", 3U, sequence_digest("TCG")},
    };
}

Contig make_contig(
    std::uint32_t index,
    const std::string &name,
    const std::string &sequence)
{
    Contig contig;
    contig.index = index;
    contig.name = name;
    contig.bases = encode(sequence);
    contig.length = contig.bases.size();
    contig.reference_sha256 = sequence_digest(sequence);
    return contig;
}

std::string row(
    const std::string &chromosome,
    const std::string &nucleotide,
    const std::string &position,
    const std::string &context,
    const std::string &dinucleotide,
    const std::string &total,
    const std::string &variant_position,
    const std::string &reference,
    const std::string &alternate,
    const std::string &reference_methylation,
    const std::string &alternate_methylation,
    const std::string &fold_change = "2",
    const std::string &p_value = "0.01",
    const std::string &comment = "1/2;1/2;1/2")
{
    return chromosome + '\t' + nucleotide + '\t' + position + '\t'
        + context + '\t' + dinucleotide + '\t' + total + '\t'
        + variant_position + '\t' + reference + '\t' + alternate + '\t'
        + reference_methylation + '\t' + alternate_methylation + '\t'
        + fold_change + '\t' + p_value + '\t' + comment + '\n';
}

std::string valid_text()
{
    return
        "# htsim ASM fixture\n\n"
        + row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "0.25", "0.75")
        + row("chr1", "G", "3", "CG", "CG", "na", "4", "A", "T", "0", "1", "na", "na")
        + row("chr1", "C", "5", "CHG", "CA", "0.4", "10", "A", "G", "0.2", "0.8")
        + row("chr1", "C", "8", "CHH", "CA", "0.3", "10", "A", "C", "0.1", "0.9")
        + row("chr2", "G", "3", "CG", "CG", "0.6", "1", "T", "A", "0.4", "0.6");
}

void verify_valid_profile(const std::vector<std::uint8_t> &input)
{
    TempFile file;
    write_bytes(file.path(), input);
    const auto digest = htsim::crypto::sha256(input);
    AsmProfile profile(file.path(), reference_catalog());
    require(profile.file_sha256() == digest, "ASM digest changed");
    require(profile.row_count() == 5U, "ASM row count changed");

    const Contig chr1_contig = make_contig(0U, "chr1", "ACGACAGCAATCG");
    const Contig chr2_contig = make_contig(1U, "chr2", "TCG");
    const auto chr1 = profile.records(chr1_contig);
    require(chr1.size() == 4U, "chr1 ASM spool count changed");
    require(
        chr1[0].target_reference_position == 1U
            && chr1[0].linked_variant_position == 3U
            && chr1[0].context == MethylationContext::cg_c
            && chr1[0].reference_probability_u16 == q(0.25F)
            && chr1[0].alternate_probability_u16 == q(0.75F)
            && chr1[0].linked_reference_base == 0U
            && chr1[0].linked_alternate_base == 3U,
        "first ASM record changed");
    require(
        chr1[1].context == MethylationContext::cg_g
            && chr1[1].reference_probability_u16 == q(0.0F)
            && chr1[1].alternate_probability_u16 == q(1.0F),
        "reverse ASM record changed");
    require(
        chr1[2].context == MethylationContext::chg_c
            && chr1[3].context == MethylationContext::chh_c,
        "ASM non-CpG contexts changed");

    const auto chr2 = profile.records(chr2_contig);
    require(
        chr2.size() == 1U
            && chr2[0].target_reference_position == 2U
            && chr2[0].linked_variant_position == 0U
            && chr2[0].context == MethylationContext::cg_g,
        "chr2 ASM record changed");
    require(
        profile.records(chr1_contig).size() == chr1.size(),
        "repeated per-contig ASM spool read changed");

    profile.validate_contig(chr1_contig);
    profile.validate_contig(chr2_contig);
    require_error(
        [&] {(void)profile.records(make_contig(2U, "chr3", "CGA"));},
        "out-of-range ASM contig index was accepted");
}

void test_plain_and_gzip_profile()
{
    verify_valid_profile(bytes_of(valid_text()));
    verify_valid_profile(gzip_bytes(valid_text()));
}

std::string valid_asm_bed_text()
{
    return
        "track type=bed name=asm-fixture\n"
        "chr1\t1\t2\tasm-1\t0\t+\t3\t4\tA\tT\t0.25\t0.75\n"
        "chr1\t2\t3\tasm-2\t0\t-\t3\t4\tA\tT\t0\t1\n"
        "chr1\t4\t5\tasm-3\t0\t+\t9\t10\tA\tG\t0.2\t0.8\n"
        "chr1\t7\t8\tasm-4\t0\t+\t9\t10\tA\tC\t0.1\t0.9\n"
        "chr2\t2\t3\tasm-5\t0\t.\t0\t1\tT\tA\t0.4\t0.6\n";
}

void verify_valid_asm_bed(const std::vector<std::uint8_t> &input)
{
    TempFile file;
    write_bytes(file.path(), input);
    const auto digest = htsim::crypto::sha256(input);
    AsmProfile profile(
        file.path(), reference_catalog(), AsmProfileFormat::bed);
    require(profile.file_sha256() == digest, "ASM BED digest changed");
    require(profile.row_count() == 5U, "ASM BED row count changed");

    const auto chr1 = profile.records(
        make_contig(0U, "chr1", "ACGACAGCAATCG"));
    require(chr1.size() == 4U, "ASM BED chr1 row count changed");
    require(
        chr1[0].target_reference_position == 1U
            && chr1[0].linked_variant_position == 3U
            && chr1[0].context == MethylationContext::cg_c
            && chr1[0].dinucleotide_second == 2U
            && chr1[0].reference_probability_u16 == q(0.25F)
            && chr1[0].alternate_probability_u16 == q(0.75F)
            && chr1[1].context == MethylationContext::cg_g
            && chr1[2].context == MethylationContext::chg_c
            && chr1[3].context == MethylationContext::chh_c,
        "ASM BED rows did not normalize to typed ASM records");
    const auto chr2 = profile.records(make_contig(1U, "chr2", "TCG"));
    require(
        chr2.size() == 1U
            && chr2[0].target_reference_position == 2U
            && chr2[0].linked_variant_position == 0U
            && chr2[0].context == MethylationContext::cg_g,
        "unstranded ASM BED row did not resolve from the reference");
}

void test_asm_bed_plain_and_gzip()
{
    verify_valid_asm_bed(bytes_of(valid_asm_bed_text()));
    verify_valid_asm_bed(gzip_bytes(valid_asm_bed_text()));
}

void require_asm_bed_error(
    const std::string &text,
    bool validate_reference,
    const std::string &message)
{
    TempFile file;
    const auto bytes = bytes_of(text);
    write_bytes(file.path(), bytes);
    require_error(
        [&] {
            AsmProfile profile(
                file.path(),
                reference_catalog(),
                AsmProfileFormat::bed);
            if (validate_reference) {
                (void)profile.records(
                    make_contig(0U, "chr1", "ACGACAGCAATCG"));
            }
        },
        message);
}

void test_asm_bed_fail_closed_boundaries()
{
    const std::vector<std::string> parse_invalid = {
        "# no data\n",
        "chr1\t1\t2\tasm\t0\t+\t3\t4\tA\tT\t0.2\n",
        "chr1\t1\t3\tasm\t0\t+\t3\t4\tA\tT\t0.2\t0.8\n",
        "chr1\t1\t2\tasm\t1001\t+\t3\t4\tA\tT\t0.2\t0.8\n",
        "chr1\t1\t2\tasm\t0\tx\t3\t4\tA\tT\t0.2\t0.8\n",
        "chr1\t1\t2\tasm\t0\t+\t3\t4\tA\tA\t0.2\t0.8\n",
        "chr1\t1\t2\tasm\t0\t+\t3\t4\tA\tT\tna\t0.8\n",
        "chr1\t1\t2\tasm\t0\t+\t3\t5\tA\tT\t0.2\t0.8\n",
    };
    for (const std::string &text : parse_invalid) {
        require_asm_bed_error(
            text, false, "invalid ASM BED syntax was accepted");
    }
    require_asm_bed_error(
        "chr1\t2\t3\tasm\t0\t+\t3\t4\tA\tT\t0.2\t0.8\n",
        true,
        "ASM BED strand/reference mismatch was accepted");
    require_asm_bed_error(
        "chr1\t1\t2\tasm\t0\t+\t3\t4\tC\tT\t0.2\t0.8\n",
        true,
        "ASM BED linked REF/reference mismatch was accepted");
}

void require_parse_error(const std::string &text, const std::string &message)
{
    TempFile file;
    const auto bytes = bytes_of(text);
    write_bytes(file.path(), bytes);
    require_error(
        [&] {
            (void)AsmProfile(
                file.path(), reference_catalog());
        },
        message);
}

void test_strict_text_boundaries()
{
    const std::vector<std::pair<std::string, std::string>> invalid = {
        {"# no data\n\n", "empty ASM was accepted"},
        {"chr1\tC\t2\tCG\tCG\t0.5\t4\tA\tT\t0.2\t0.8\t2\t0.01\n",
         "thirteen fields were accepted"},
        {row("missing", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8"),
         "unknown contig was accepted"},
        {row("chr1", "C", "0", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8"),
         "target position zero was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "0", "A", "T", "0.2", "0.8"),
         "linked position zero was accepted"},
        {row("chr1", "C", "14", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8"),
         "outside target was accepted"},
        {row("chr1", "A", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8"),
         "non-C/G target was accepted"},
        {row("chr1", "C", "2", "CX", "CG", "0.5", "4", "A", "T", "0.2", "0.8"),
         "unknown context was accepted"},
        {row("chr1", "C", "2", "CG", "CA", "0.5", "4", "A", "T", "0.2", "0.8"),
         "wrong context dinucleotide was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "NaN", "4", "A", "T", "0.2", "0.8"),
         "NaN total was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "4", "N", "T", "0.2", "0.8"),
         "non-ACGT REF was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "A", "0.2", "0.8"),
         "equal REF and ALT were accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "na", "0.8"),
         "missing REF_METH was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "1.1"),
         "ALT_METH over one was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8", "inf"),
         "infinite fold change was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8", "2", "1.1"),
         "p-value over one was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8", "2", "0.1", ""),
         "empty comment was accepted"},
        {row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8")
             + row("chr1", "G", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8"),
         "duplicate target was accepted"},
        {row("chr2", "G", "3", "CG", "CG", "0.5", "1", "T", "A", "0.2", "0.8")
             + row("chr1", "C", "2", "CG", "CG", "0.5", "4", "A", "T", "0.2", "0.8"),
         "returning contig was accepted"},
    };
    for (const auto &fixture : invalid) {
        require_parse_error(fixture.first, fixture.second);
    }

    TempFile file;
    const auto bytes = bytes_of(valid_text());
    write_bytes(file.path(), bytes);
    require_error(
        [&] {(void)AsmProfile(file.path(), {});},
        "empty reference catalog was accepted");
}

void test_reference_validation_boundaries()
{
    TempFile file;
    const auto bytes = bytes_of(valid_text());
    write_bytes(file.path(), bytes);
    AsmProfile profile(file.path(), reference_catalog());
    require_error(
        [&] {profile.validate_contig(make_contig(0U, "chr1", "AGGACAGCAATCG"));},
        "ASM target/reference mismatch was accepted");
    require_error(
        [&] {profile.validate_contig(make_contig(0U, "chr1", "ACGTCAGCAATCG"));},
        "ASM linked REF mismatch was accepted");
    require_error(
        [&] {profile.validate_contig(make_contig(0U, "wrong", "ACGACAGCAATCG"));},
        "ASM catalog identity mismatch was accepted");

    const Bases bases = encode("ACGACAGCAATCG");
    require_error(
        [&] {
            htsim::methdb::validate_asm_records(
                bases,
                {
                    {2U, 3U, q(0.2F), q(0.8F),
                     MethylationContext::cg_g, 2U, 0U, 3U},
                    {1U, 3U, q(0.2F), q(0.8F),
                     MethylationContext::cg_c, 2U, 0U, 3U},
                });
        },
        "unsorted normalized ASM records were accepted");
    require_error(
        [&] {
            htsim::methdb::validate_asm_records(
                bases,
                {{1U, 3U, q(0.2F), q(0.8F),
                  MethylationContext::chh_c, 0U, 0U, 3U}});
        },
        "wrong normalized ASM context was accepted");
    require_error(
        [&] {
            htsim::methdb::validate_asm_records(
                bases,
                {{1U, 3U, q(0.2F), q(0.8F),
                  MethylationContext::cg_c, 2U, 1U, 3U}});
        },
        "wrong normalized linked REF was accepted");
}

} // namespace

int main()
{
    try {
        test_plain_and_gzip_profile();
        test_asm_bed_plain_and_gzip();
        test_asm_bed_fail_closed_boundaries();
        test_strict_text_boundaries();
        test_reference_validation_boundaries();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "asm_profile_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
