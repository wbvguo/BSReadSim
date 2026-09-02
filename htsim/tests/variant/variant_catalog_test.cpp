#include "variant.h"

#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>
#include <zlib.h>

namespace {

using htsim::model::HaplotypeMask;
using htsim::model::VariantKind;
using htsim::variant::ContigVariants;
using htsim::variant::Variant;
using htsim::variant::VariantCatalogError;
using htsim::variant::VariantFile;

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-variant-catalog-XXXXXX";
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

template <typename Operation>
void require_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const VariantCatalogError &) {
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
    if (!output) {throw std::runtime_error("temporary VCF write failed");}
}

std::vector<std::uint8_t> gzip_bytes(const std::string &text)
{
    z_stream stream = {};
    require(deflateInit2(
                &stream,
                Z_BEST_SPEED,
                Z_DEFLATED,
                MAX_WBITS + 16,
                8,
                Z_DEFAULT_STRATEGY) == Z_OK,
            "deflateInit2 failed");
    std::vector<std::uint8_t> result(text.size() + 128U);
    stream.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(text.data()));
    stream.avail_in = static_cast<uInt>(text.size());
    stream.next_out = result.data();
    stream.avail_out = static_cast<uInt>(result.size());
    const int status = deflate(&stream, Z_FINISH);
    require(status == Z_STREAM_END, "deflate failed");
    result.resize(result.size() - stream.avail_out);
    require(deflateEnd(&stream) == Z_OK, "deflateEnd failed");
    return result;
}

std::vector<htsim::reference::ContigMetadata> reference_catalog()
{
    return {
        {"chr1", 20, {}},
        {"chr2", 8, {}},
    };
}

htsim::model::Bases encode(const std::string &text)
{
    htsim::model::Bases result;
    for (const char base : text) {
        switch (base) {
        case 'A': result.push_back(0U); break;
        case 'C': result.push_back(1U); break;
        case 'G': result.push_back(2U); break;
        case 'T': result.push_back(3U); break;
        default: throw std::runtime_error("invalid test base");
        }
    }
    return result;
}

std::string valid_vcf()
{
    return
        "##fileformat=VCFv4.3\r\n"
        "##source=fixture\r\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\r\n"
        "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\r\n"
        "chr1\t5\t.\tA\tAGT\t.\tPASS\t.\tDP:GT\t9:1/0\r\n"
        "chr1\t9\t.\tACG\tA\t.\tPASS\t.\tGT\t1|1\r\n"
        "chr2\t1\t.\tT\tC\t.\tPASS\t.\tGT\t0/0";
}

VariantFile load(
    const std::vector<std::uint8_t> &bytes,
    TempFile &file,
    std::uint64_t seed = 123)
{
    write_bytes(file.path(), bytes);
    return VariantFile(
        file.path(),
        reference_catalog(),
        seed);
}

void test_normalization_phasing_and_reference_validation()
{
    TempFile file;
    const auto bytes = bytes_of(valid_vcf());
    const VariantFile variants = load(bytes, file);
    require(variants.variant_count() == 3U
                && variants.variants(0).size() == 3U
                && variants.variants(1).empty()
                && variants.file_sha256() == htsim::crypto::sha256(bytes),
            "VCF event counts or raw identity changed");

    const Variant &snv = variants.variants(0)[0];
    const Variant &insertion = variants.variants(0)[1];
    const Variant &deletion = variants.variants(0)[2];
    require(snv.kind == VariantKind::snv
                && snv.reference_start == 1U && snv.reference_end == 2U
                && snv.ref_bases == encode("C") && snv.alt_bases == encode("T")
                && snv.alt_haplotypes == HaplotypeMask::haplotype_2,
            "phased SNV projection changed");
    require(insertion.kind == VariantKind::insertion
                && insertion.reference_start == 5U
                && insertion.reference_end == 5U
                && insertion.ref_bases.empty()
                && insertion.alt_bases == encode("GT")
                && insertion.alt_haplotypes == HaplotypeMask::haplotype_1,
            "insertion normalization or frozen unphased vector changed");
    require(deletion.kind == VariantKind::deletion
                && deletion.reference_start == 9U
                && deletion.reference_end == 11U
                && deletion.ref_bases == encode("CG")
                && deletion.alt_bases.empty()
                && deletion.alt_haplotypes == HaplotypeMask::both,
            "deletion normalization or homozygous mask changed");

    const ContigVariants checked(
        encode("ACGTACGTACGTACGTACGT"), variants.variants(0), 0);
    require(checked.variants().size() == 3U && checked.contig_index() == 0U
                && checked.reference_length() == 20U,
            "reference-validated VCF catalog lost events");
    require_error(
        [&] {
            (void)ContigVariants(
                encode("AAGTACGTACGTACGTACGT"), variants.variants(0), 0);
        },
        "VCF REF mismatch was accepted");

    TempFile other_seed_file;
    const VariantFile other_seed = load(bytes, other_seed_file, 126);
    require(other_seed.variants(0)[1].alt_haplotypes
                == HaplotypeMask::haplotype_2,
            "master seed did not control unphased heterozygote assignment");
}

void test_gzip_snapshot_and_format_versions()
{
    TempFile compressed_file;
    const auto compressed = gzip_bytes(valid_vcf());
    const VariantFile compressed_variants = load(compressed, compressed_file);
    require(compressed_variants.variant_count() == 3U
                && compressed_variants.file_sha256()
                    == htsim::crypto::sha256(compressed),
            "gzip VCF snapshot changed event projection or raw identity");

    std::string version_42 = valid_vcf();
    version_42.replace(
        version_42.find("VCFv4.3"), std::string("VCFv4.3").size(), "VCFv4.2");
    TempFile version_file;
    require(load(bytes_of(version_42), version_file).variant_count() == 3U,
            "VCF 4.2 was not accepted");
}

void test_equal_vcf_positions_are_allowed_for_distinct_normalized_events()
{
    const std::string text =
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t2\tsnv\tC\tT\t.\tPASS\t.\tGT\t1|0\n"
        "chr1\t2\tdel\tCG\tC\t.\tPASS\t.\tGT\t0|1\n";
    TempFile file;
    const VariantFile variants = load(bytes_of(text), file);
    require(variants.variant_count() == 2U,
            "distinct events sharing an anchored VCF POS were rejected");
    require(variants.variants(0)[0].reference_start == 1U
                && variants.variants(0)[1].reference_start == 2U,
            "equal-POS VCF rows did not normalize into canonical order");
}

void test_unsupported_variant_shapes_are_skipped()
{
    const std::string text =
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t2\tkept-1\tC\tT\t.\tPASS\t.\tGT\t0|1\n"
        "chr1\t5\tcomplex-insertion\tC\tTCGA\t.\tPASS\t.\tGT\t0|1\n"
        "chr1\t9\tmnp\tAC\tGT\t.\tPASS\t.\tGT\t1|1\n"
        "chr1\t13\tlong-insertion\tA\tACGTAC\t.\tPASS\t.\tGT\t0|1\n"
        "chr1\t14\tlong-deletion\tCAAAAA\tC\t.\tPASS\t.\tGT\t1|1\n"
        "chr1\t20\tkept-2\tT\tA\t.\tPASS\t.\tGT\t1|0\n";
    TempFile file;
    const VariantFile variants = load(bytes_of(text), file);
    require(variants.variant_count() == 2U
                && variants.row_count() == 6U
                && variants.input_contig_count() == 1U
                && variants.reference_genotype_count() == 0U
                && variants.skipped_mnp_count() == 1U
                && variants.skipped_complex_replacement_count() == 1U
                && variants.skipped_long_indel_count() == 2U
                && variants.skipped_unsupported_count() == 4U
                && variants.variants(0).size() == 2U
                && variants.variants(0)[0].id == "kept-1"
                && variants.variants(0)[1].id == "kept-2",
            "unsupported VCF accounting or filtering changed");
}

void test_four_base_indels_are_retained()
{
    const std::string text =
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
        "chr1\t1\tins-4\tA\tAACGT\t.\tPASS\t.\tGT\t0|1\n"
        "chr1\t5\tdel-4\tACGTA\tA\t.\tPASS\t.\tGT\t1|0\n";
    TempFile file;
    const VariantFile variants = load(bytes_of(text), file);
    require(variants.variant_count() == 2U
                && variants.variants(0)[0].kind == VariantKind::insertion
                && variants.variants(0)[0].alt_bases.size() == 4U
                && variants.variants(0)[1].kind == VariantKind::deletion
                && variants.variants(0)[1].ref_bases.size() == 4U,
            "four-base VCF indel was skipped");
    (void)ContigVariants(
        encode("ACGTACGTACGTACGTACGT"), variants.variants(0), 0);
}

void test_file_rejections()
{
    const std::string header =
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n";
    const std::vector<std::string> invalid = {
        "",
        "##fileformat=VCFv4.1\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n",
        "##fileformat=VCFv4.3\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n",
        header + "missing\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\n",
        header + "chr1\t0\t.\tA\tC\t.\tPASS\t.\tGT\t0|1\n",
        header + "chr1\t1\t.\tA\tC,G\t.\tPASS\t.\tGT\t0|1\n",
        header + "chr1\t1\t.\ta\tC\t.\tPASS\t.\tGT\t0|1\n",
        header + "chr1\t1\t.\tA\tC\t.\tPASS\t.\tDP\t4\n",
        header + "chr1\t1\t.\tA\tC\t.\tPASS\t.\tGT\t0|2\n",
        header
            + "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\n"
              "chr1\t2\t.\tC\tA\t.\tPASS\t.\tGT\t0|1\n",
        header
            + "chr2\t1\t.\tT\tC\t.\tPASS\t.\tGT\t0|1\n"
              "chr1\t2\t.\tC\tT\t.\tPASS\t.\tGT\t0|1\n",
        header
            + "chr1\t2\t.\tCG\tC\t.\tPASS\t.\tGT\t0|1\n"
              "chr1\t3\t.\tG\tT\t.\tPASS\t.\tGT\t0|1\n",
        header + "\n",
    };
    for (const std::string &text : invalid) {
        TempFile file;
        const auto bytes = bytes_of(text);
        write_bytes(file.path(), bytes);
        require_error(
            [&] {
                (void)VariantFile(
                    file.path(),
                    reference_catalog(),
                    123);
            },
            "invalid VCF was accepted: " + text);
    }
}

void test_typed_catalog_rejections()
{
    const auto bases = encode("ACGTACGT");
    std::vector<Variant> wrong_contig = {
        {1, 1, 2, VariantKind::snv, encode("C"), encode("T"),
         HaplotypeMask::haplotype_1},
    };
    require_error(
        [&] {(void)ContigVariants(bases, wrong_contig, 0);},
        "event from another contig was accepted");

    std::vector<Variant> invalid_mask = wrong_contig;
    invalid_mask[0].contig_index = 0;
    invalid_mask[0].alt_haplotypes = static_cast<HaplotypeMask>(0);
    require_error(
        [&] {(void)ContigVariants(bases, invalid_mask, 0);},
        "zero haplotype mask was accepted");

    std::vector<Variant> unsorted = {
        {0, 5, 6, VariantKind::snv, encode("C"), encode("A"),
         HaplotypeMask::both},
        {0, 1, 2, VariantKind::snv, encode("C"), encode("T"),
         HaplotypeMask::both},
    };
    require_error(
        [&] {(void)ContigVariants(bases, unsorted, 0);},
        "non-canonical typed events were accepted");

    const std::vector<std::vector<Variant>> malformed = {
        {{0, 1, 2, VariantKind::snv, encode("C"), encode("TA"),
          HaplotypeMask::both}},
        {{0, 2, 2, VariantKind::insertion, {}, {},
          HaplotypeMask::both}},
        {{0, 1, 2, VariantKind::deletion, encode("C"), encode("T"),
          HaplotypeMask::both}},
        {{0, 1, 2, static_cast<VariantKind>(99), encode("C"), encode("T"),
          HaplotypeMask::both}},
        {{0, 1, 2, VariantKind::snv, {9}, encode("T"),
          HaplotypeMask::both}},
    };
    for (const auto &events : malformed) {
        require_error(
            [&] {(void)ContigVariants(bases, events, 0);},
            "malformed typed VCF event was accepted");
    }
}

} // namespace

int main()
{
    try {
        test_normalization_phasing_and_reference_validation();
        test_gzip_snapshot_and_format_versions();
        test_equal_vcf_positions_are_allowed_for_distinct_normalized_events();
        test_unsupported_variant_shapes_are_skipped();
        test_four_base_indels_are_retained();
        test_file_rejections();
        test_typed_catalog_rejections();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "variant_catalog_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
