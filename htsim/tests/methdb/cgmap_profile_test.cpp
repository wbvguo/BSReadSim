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

using htsim::methdb::CgmapProfile;
using htsim::methdb::CgmapProfileError;
using htsim::methdb::CgmapRecord;
using htsim::methdb::MethylationProfileFormat;
using htsim::model::Bases;
using htsim::model::MethylationContext;
using htsim::reference::Contig;
using htsim::reference::ContigMetadata;

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-cgmap-profile-XXXXXX";
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
    } catch (const CgmapProfileError &) {
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

std::string valid_text()
{
    return
        "# CGmap fixture\n"
        "\n"
        "chr1\tC\t2\tCG\tCG\t0.25\t1\t4\n"
        "chr1\tG\t3\tCG\tCG\tna\t0\t0\n"
        "chr1\tC\t5\tCHG\tCA\t0.5\t2\t4\n"
        "chr1\tC\t8\tCHH\tCA\t1\t4\t4\n"
        "chr2\tG\t3\tCG\tCG\t0\t0\t10\n";
}

void verify_valid_profile(const std::vector<std::uint8_t> &input)
{
    TempFile file;
    write_bytes(file.path(), input);
    const auto digest = htsim::crypto::sha256(input);
    CgmapProfile profile(file.path(), reference_catalog());
    require(profile.file_sha256() == digest, "CGmap digest changed");
    require(profile.row_count() == 5U, "CGmap row count changed");
    require(
        profile.defined_probability_count() == 4U,
        "CGmap defined-probability count changed");

    const Contig chr1_contig = make_contig(0U, "chr1", "ACGACAGCAATCG");
    const Contig chr2_contig = make_contig(1U, "chr2", "TCG");
    const auto chr1 = profile.records(chr1_contig);
    require(chr1.size() == 4U, "chr1 CGmap spool count changed");
    require(
        chr1[0].reference_position == 1U
            && chr1[0].context == MethylationContext::cg_c
            && chr1[0].has_probability
            && chr1[0].probability_u16 == q(0.25F),
        "first CGmap record changed");
    require(
        chr1[1].reference_position == 2U
            && chr1[1].context == MethylationContext::cg_g
            && !chr1[1].has_probability
            && chr1[1].probability_u16 == 0U,
        "CGmap na record changed");
    require(
        chr1[2].context == MethylationContext::chg_c
            && chr1[3].context == MethylationContext::chh_c,
        "CGmap non-CpG contexts changed");

    const auto chr2 = profile.records(chr2_contig);
    require(
        chr2.size() == 1U && chr2[0].reference_position == 2U
            && chr2[0].context == MethylationContext::cg_g
            && chr2[0].has_probability
            && chr2[0].probability_u16 == 0U,
        "chr2 CGmap record changed");
    require(
        profile.records(chr1_contig).size() == chr1.size(),
        "repeated per-contig spool read changed");

    profile.validate_contig(chr1_contig);
    profile.validate_contig(chr2_contig);
    require_error(
        [&] {(void)profile.records(make_contig(2U, "chr3", "CGA"));},
        "out-of-range CGmap contig index was accepted");
}

void test_plain_and_gzip_profile()
{
    verify_valid_profile(bytes_of(valid_text()));
    verify_valid_profile(gzip_bytes(valid_text()));
}

std::string valid_bed_methyl_text(bool extended)
{
    if (extended) {
        return
            "track type=bedMethyl name=fixture\n"
            "chr1\t1\t2\tm\t4\t+\t1\t2\t255,0,0\t4\t25\t1\t3\t0\t0\t0\t0\t0\n"
            "chr1\t2\t3\tm\t0\t-\t2\t3\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
            "chr1\t4\t5\th\t4\t+\t4\t5\t255,0,0\t4\t50\t2\t2\t0\t0\t0\t0\t0\n"
            "chr1\t7\t8\th\t4\t+\t7\t8\t255,0,0\t4\t100\t4\t0\t0\t0\t0\t0\t0\n"
            "chr2\t2\t3\tm\t10\t.\t2\t3\t255,0,0\t10\t0\t0\t10\t0\t0\t0\t0\t0\n";
    }
    return
        "browser position chr1:1-13\n"
        "chr1\t1\t2\tm\t4\t+\t1\t2\t255,0,0\t4\t25\n"
        "chr1\t2\t3\tm\t0\t-\t2\t3\t0\t0\t0\n"
        "chr1\t4\t5\th\t4\t+\t4\t5\t255,0,0\t4\t50\n"
        "chr1\t7\t8\th\t4\t+\t7\t8\t255,0,0\t4\t100\n"
        "chr2\t2\t3\tm\t10\t.\t2\t3\t255,0,0\t10\t0\n";
}

void verify_valid_bed_methyl(const std::vector<std::uint8_t> &input)
{
    TempFile file;
    write_bytes(file.path(), input);
    const auto digest = htsim::crypto::sha256(input);
    CgmapProfile profile(
        file.path(),
        reference_catalog(),
        MethylationProfileFormat::bed_methyl);
    require(profile.file_sha256() == digest, "bedMethyl digest changed");
    require(profile.row_count() == 5U, "bedMethyl row count changed");
    require(
        profile.defined_probability_count() == 5U,
        "bedMethyl defined-probability count changed");

    const auto chr1 = profile.records(
        make_contig(0U, "chr1", "ACGACAGCAATCG"));
    require(chr1.size() == 4U, "bedMethyl chr1 row count changed");
    require(
        chr1[0].reference_position == 1U
            && chr1[0].context == MethylationContext::cg_c
            && chr1[0].dinucleotide_second == 2U
            && chr1[0].probability_u16 == q(0.25F)
            && chr1[1].context == MethylationContext::cg_g
            && chr1[2].context == MethylationContext::chg_c
            && chr1[3].context == MethylationContext::chh_c,
        "bedMethyl rows did not normalize to typed CGmap records");
    const auto chr2 = profile.records(make_contig(1U, "chr2", "TCG"));
    require(
        chr2.size() == 1U
            && chr2[0].context == MethylationContext::cg_g
            && chr2[0].probability_u16 == 0U,
        "unstranded bedMethyl row did not resolve from the reference");
}

void test_bed_methyl_plain_gzip_and_widths()
{
    for (const bool extended : {false, true}) {
        const std::string text = valid_bed_methyl_text(extended);
        verify_valid_bed_methyl(bytes_of(text));
        verify_valid_bed_methyl(gzip_bytes(text));
    }
}

void require_bed_methyl_error(
    const std::string &text,
    bool validate_reference,
    const std::string &message)
{
    TempFile file;
    const auto bytes = bytes_of(text);
    write_bytes(file.path(), bytes);
    require_error(
        [&] {
            CgmapProfile profile(
                file.path(),
                reference_catalog(),
                MethylationProfileFormat::bed_methyl);
            if (validate_reference) {
                (void)profile.records(
                    make_contig(0U, "chr1", "ACGACAGCAATCG"));
            }
        },
        message);
}

void test_bed_methyl_fail_closed_boundaries()
{
    const std::vector<std::string> parse_invalid = {
        "# no data\n",
        "chr1\t1\t2\tm\t4\t+\t1\t2\t255,0,0\t4\n",
        "chr1\t1\t3\tm\t4\t+\t1\t3\t255,0,0\t4\t25\n",
        "chr1\t1\t2\tm\t4\t+\t1\t2\t256,0,0\t4\t25\n",
        "chr1\t1\t2\tm\t4\t+\t1\t2\t255,0,0\t4\t101\n",
        "chr1\t1\t2\tm\t4\t+\t1\t2\t255,0,0\t4\t25\n"
        "chr1\t2\t3\tm\t4\t-\t2\t3\t255,0,0\t4\t25\t1\t3\t0\t0\t0\t0\t0\n",
        "chr1\t1\t2\tm\t4\t+\t1\t2\t255,0,0\t4\t25\t1\t2\t0\t0\t0\t0\t0\n",
    };
    for (const std::string &text : parse_invalid) {
        require_bed_methyl_error(
            text, false, "invalid bedMethyl syntax was accepted");
    }
    require_bed_methyl_error(
        "chr1\t2\t3\tm\t4\t+\t2\t3\t255,0,0\t4\t25\n",
        true,
        "bedMethyl strand/reference mismatch was accepted");
    require_bed_methyl_error(
        "chr1\t3\t4\tm\t4\t.\t3\t4\t255,0,0\t4\t25\n",
        true,
        "bedMethyl non-C/G target was accepted");
}

void verify_methbg(const std::vector<std::uint8_t> &input)
{
    TempFile file;
    write_bytes(file.path(), input);
    CgmapProfile profile(
        file.path(), reference_catalog(), MethylationProfileFormat::methbg);
    require(profile.row_count() == 4U, "MethBG row count changed");
    require(
        profile.defined_probability_count() == 4U,
        "MethBG defined-probability count changed");
    const auto records = profile.records(
        make_contig(0U, "chr1", "ACGACAGCAATCG"));
    require(
        records.size() == 3U
            && records[0].context == MethylationContext::cg_c
            && records[0].probability_u16 == q(0.25F)
            && records[1].context == MethylationContext::chg_c
            && records[2].context == MethylationContext::chh_c,
        "MethBG rows did not normalize from the reference");
}

void test_methbg_plain_gzip_and_boundaries()
{
    const std::string text =
        "track type=bedGraph name=fixture\n"
        "chr1\t1\t2\t0.25\n"
        "chr1\t4\t5\t0.5\n"
        "chr1\t7\t8\t1\n"
        "chr2\t2\t3\t0\n";
    verify_methbg(bytes_of(text));
    verify_methbg(gzip_bytes(text));

    TempFile invalid;
    write_bytes(invalid.path(), bytes_of("chr1\t1\t3\t0.5\n"));
    require_error(
        [&] {
            (void)CgmapProfile(
                invalid.path(),
                reference_catalog(),
                MethylationProfileFormat::methbg);
        },
        "invalid MethBG interval was accepted");
}

void verify_methbed(const std::string &text)
{
    TempFile file;
    write_bytes(file.path(), bytes_of(text));
    CgmapProfile profile(
        file.path(), reference_catalog(), MethylationProfileFormat::methbed);
    require(profile.row_count() == 4U, "MethBED row count changed");
    const auto records = profile.records(
        make_contig(0U, "chr1", "ACGACAGCAATCG"));
    require(
        records.size() == 3U
            && records[0].context == MethylationContext::cg_c
            && records[0].probability_u16 == q(0.25F)
            && records[1].context == MethylationContext::chg_c
            && records[1].probability_u16 == q(0.5F)
            && records[2].context == MethylationContext::chh_c
            && records[2].probability_u16 == q(1.0F),
        "MethBED rows did not normalize from the reference");
}

void test_methbed_widths_and_boundaries()
{
    verify_methbed(
        "chr1\t1\t2\t.\t250\t+\n"
        "chr1\t4\t5\t.\t500\t+\n"
        "chr1\t7\t8\t.\t1000\t+\n"
        "chr2\t2\t3\t.\t0\t-\n");
    verify_methbed(
        "browser position chr1:1-13\n"
        "chr1\t1\t2\t.\t250\t+\t1\t4\tC\tCG\n"
        "chr1\t4\t5\t.\t500\t+\t2\t4\tC\tCHG\n"
        "chr1\t7\t8\t.\t1000\t+\t4\t4\tC\tCHH\n"
        "chr2\t2\t3\t.\t0\t-\t0\t10\tG\tCG\n");

    TempFile invalid;
    write_bytes(
        invalid.path(),
        bytes_of("chr1\t1\t2\t.\t250\t+\t1\t4\tC\tCHH\n"));
    CgmapProfile profile(
        invalid.path(), reference_catalog(), MethylationProfileFormat::methbed);
    require_error(
        [&] {
            (void)profile.records(
                make_contig(0U, "chr1", "ACGACAGCAATCG"));
        },
        "MethBED context/reference mismatch was accepted");
}

void require_parse_error(const std::string &text, const std::string &message)
{
    TempFile file;
    const auto bytes = bytes_of(text);
    write_bytes(file.path(), bytes);
    require_error(
        [&] {
            (void)CgmapProfile(
                file.path(), reference_catalog());
        },
        message);
}

void test_strict_text_boundaries()
{
    const std::vector<std::pair<std::string, std::string>> invalid = {
        {"# no data\n\n", "empty CGmap was accepted"},
        {"chr1\tC\t2\tCG\tCG\t0.5\t1\n", "seven fields were accepted"},
        {"missing\tC\t2\tCG\tCG\t0.5\t1\t2\n", "unknown contig was accepted"},
        {"chr1\tC\t0\tCG\tCG\t0.5\t1\t2\n", "position zero was accepted"},
        {"chr1\tC\t14\tCG\tCG\t0.5\t1\t2\n", "outside position was accepted"},
        {"chr1\tA\t2\tCG\tCG\t0.5\t1\t2\n", "non-C/G nucleotide was accepted"},
        {"chr1\tC\t2\tCX\tCG\t0.5\t1\t2\n", "unknown context was accepted"},
        {"chr1\tC\t2\tCG\tGA\t0.5\t1\t2\n", "unknown dinucleotide was accepted"},
        {"chr1\tC\t2\tCG\tCA\t0.5\t1\t2\n",
         "non-CG dinucleotide for CG context was accepted"},
        {"chr1\tC\t5\tCHG\tCG\t0.5\t1\t2\n",
         "CG dinucleotide for CHG context was accepted"},
        {"chr1\tC\t2\tCG\tCG\tNaN\t1\t2\n", "NaN was accepted"},
        {"chr1\tC\t2\tCG\tCG\tNA\t1\t2\n", "uppercase NA was accepted"},
        {"chr1\tC\t2\tCG\tCG\t-0.1\t1\t2\n", "negative level was accepted"},
        {"chr1\tC\t2\tCG\tCG\t1.1\t1\t2\n", "level over one was accepted"},
        {"chr1\tC\t2\tCG\tCG\t0.5\t-1\t2\n", "negative count was accepted"},
        {"chr1\tC\t2\tCG\tCG\t0.5\t4294967296\t4294967296\n",
         "count over uint32 was accepted"},
        {"chr1\tC\t2\tCG\tCG\t0.5\t3\t2\n", "MC greater than NC was accepted"},
        {"chr1\tC\t2\tCG\tCG\t0.5\t1\t2\n"
         "chr1\tG\t2\tCG\tCG\t0.5\t1\t2\n",
         "duplicate position was accepted"},
        {"chr2\tG\t3\tCG\tCG\t0.5\t1\t2\n"
         "chr1\tC\t2\tCG\tCG\t0.5\t1\t2\n",
         "returning contig was accepted"},
    };
    for (const auto &fixture : invalid) {
        require_parse_error(fixture.first, fixture.second);
    }

    TempFile file;
    const auto bytes = bytes_of(valid_text());
    write_bytes(file.path(), bytes);
    require_error(
        [&] {(void)CgmapProfile(file.path(), {});},
        "empty reference catalog was accepted");
}

void test_reference_validation_boundaries()
{
    TempFile file;
    const auto bytes = bytes_of(valid_text());
    write_bytes(file.path(), bytes);
    CgmapProfile profile(file.path(), reference_catalog());
    require_error(
        [&] {profile.validate_contig(make_contig(0U, "chr1", "AGGACAGCAATCG"));},
        "CGmap nucleotide/reference mismatch was accepted");
    require_error(
        [&] {profile.validate_contig(make_contig(0U, "chr1", "ACGACAACAATCG"));},
        "CGmap context/reference mismatch was accepted");
    require_error(
        [&] {profile.validate_contig(make_contig(0U, "wrong", "ACGACAGCAATCG"));},
        "CGmap catalog identity mismatch was accepted");

    const Bases bases = encode("ACGACAGCAATCG");
    require_error(
        [&] {
            htsim::methdb::validate_cgmap_records(
                bases,
                {
                    {2U, q(0.2F), MethylationContext::cg_g, true, 2U},
                    {1U, q(0.2F), MethylationContext::cg_c, true, 2U},
                });
        },
        "unsorted normalized CGmap records were accepted");
    require_error(
        [&] {
            htsim::methdb::validate_cgmap_records(
                bases,
                {{1U, q(0.2F), MethylationContext::chh_c, true, 0U}});
        },
        "wrong normalized CGmap context was accepted");
    require_error(
        [&] {
            htsim::methdb::validate_cgmap_records(
                bases,
                {{1U, q(0.2F), MethylationContext::cg_c, false, 2U}});
        },
        "nonzero undefined CGmap probability was accepted");

    const std::string wrong_dinucleotide =
        "chr1\tC\t8\tCHH\tCT\t0.5\t1\t2\n";
    TempFile dinucleotide_file;
    const auto dinucleotide_bytes = bytes_of(wrong_dinucleotide);
    write_bytes(dinucleotide_file.path(), dinucleotide_bytes);
    CgmapProfile dinucleotide_profile(
        dinucleotide_file.path(),
        reference_catalog());
    require_error(
        [&] {
            dinucleotide_profile.validate_contig(
                make_contig(0U, "chr1", "ACGACAGCAATCG"));
        },
        "CGmap dinucleotide/reference mismatch was accepted");
}

} // namespace

int main()
{
    try {
        test_plain_and_gzip_profile();
        test_bed_methyl_plain_gzip_and_widths();
        test_bed_methyl_fail_closed_boundaries();
        test_methbg_plain_gzip_and_boundaries();
        test_methbed_widths_and_boundaries();
        test_strict_text_boundaries();
        test_reference_validation_boundaries();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "cgmap_profile_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
