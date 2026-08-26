#include "methdb.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

namespace {

using htsim::methdb::CatalogSite;
using htsim::methdb::DiploidMethylationCatalog;
using htsim::methdb::DiploidSite;
using htsim::methdb::MethylationCatalog;
using htsim::methdb::Snapshot;
using htsim::methdb::SnapshotWriter;
using htsim::model::MethylationAllele;
using htsim::model::MethylationContext;
using htsim::model::MethylationSource;

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-methdb-snapshot-XXXXXX";
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
void require_snapshot_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const htsim::methdb::SnapshotError &) {
        return;
    }
    throw std::runtime_error(message);
}

std::vector<std::uint8_t> bytes_of(const std::string &text)
{
    return {text.begin(), text.end()};
}

void write_bytes(const std::string &path, const std::string &bytes)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    output.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
    if (!output) {throw std::runtime_error("temporary MethDB write failed");}
}

std::vector<CatalogSite> reference_sites()
{
    std::vector<CatalogSite> sites;
    sites.reserve(50000U);
    for (std::uint32_t index = 0U; index < 50000U; ++index) {
        const std::uint32_t bits = index * UINT32_C(2654435761);
        const float probability = static_cast<float>(
            static_cast<double>(bits & UINT32_C(0x00ffffff))
            / static_cast<double>(UINT32_C(0x01000000)));
        sites.push_back({
            index * 2U,
            probability,
            index % 2U == 0U
                ? MethylationContext::cg_c
                : MethylationContext::chh_g,
            index % 3U == 0U
                ? MethylationSource::beta
                : MethylationSource::cgmap,
        });
    }
    sites.back().methylation_probability = 1.0F;
    return sites;
}

float unorm16_probability(float value)
{
    constexpr double scale = 65535.0;
    const auto encoded = static_cast<std::uint32_t>(
        static_cast<double>(value) * scale + 0.5);
    return static_cast<float>(encoded) / static_cast<float>(scale);
}

bool same_unorm16_probability(float left, float right)
{
    constexpr float maximum_quantization_error = 0.5F / 65535.0F;
    return std::fabs(left - right)
        <= maximum_quantization_error + std::numeric_limits<float>::epsilon();
}

bool same_diploid_sites(
    const std::vector<DiploidSite> &left,
    const std::vector<DiploidSite> &right)
{
    if (left.size() != right.size()) {return false;}
    for (std::size_t index = 0U; index < left.size(); ++index) {
        if (left[index].origin_id != right[index].origin_id
            || left[index].context != right[index].context
            || left[index].methylation_source
                != right[index].methylation_source
            || left[index].allele != right[index].allele
            || !same_unorm16_probability(
                left[index].methylation_probability,
                right[index].methylation_probability)) {
            return false;
        }
    }
    return true;
}

void test_methdb_is_deterministic_compact_and_uses_unorm16_probabilities()
{
    const auto binding = htsim::crypto::sha256(bytes_of("binding"));
    const std::vector<CatalogSite> sites = reference_sites();
    const std::vector<htsim::reference::ContigMetadata> metadata = {
        {"chrReference", 100001U,
         htsim::crypto::sha256(bytes_of("reference-1"))},
        {"chrDiploid", 100U,
         htsim::crypto::sha256(bytes_of("reference-2"))},
    };
    const MethylationCatalog reference_catalog(100001U, sites);
    const std::vector<DiploidSite> shared = {
        {0U, MethylationContext::cg_c, MethylationSource::beta,
         MethylationAllele::shared, 0.25F},
        {7U, MethylationContext::chg_c, MethylationSource::cgmap,
         MethylationAllele::shared, 0.75F},
        {UINT64_C(0x8000000000000004), MethylationContext::chh_g,
         MethylationSource::beta, MethylationAllele::shared, 0.5F},
    };
    const std::array<std::vector<DiploidSite>, 2> haplotypes = {{
        {{2U, MethylationContext::cg_c, MethylationSource::asm_source,
          MethylationAllele::alternate_haplotype, 0.125F}},
        {{2U, MethylationContext::cg_c, MethylationSource::asm_source,
          MethylationAllele::reference_haplotype, 0.875F}},
    }};
    const DiploidMethylationCatalog diploid_catalog(
        1U, 100U, shared, haplotypes);

    const auto serialize = [&]() {
        std::ostringstream output(std::ios::binary);
        SnapshotWriter writer(output, binding, 2U);
        writer.write_reference(metadata[0], reference_catalog);
        writer.write_diploid(metadata[1], diploid_catalog);
        writer.finish();
        return output.str();
    };
    const std::string first = serialize();
    const std::string repeated = serialize();
    require(first == repeated, "MethDB bytes changed across repeated writes");
    constexpr std::size_t magic_size = sizeof(htsim::methdb::methdb_magic) - 1U;
    require(first.substr(0U, magic_size) == htsim::methdb::methdb_magic,
            "MethDB writer did not emit the methdb magic");
    require(first.size() > magic_size
                && static_cast<std::uint8_t>(first[magic_size])
                    == htsim::methdb::methdb_version,
            "MethDB writer did not emit the current version");
    const std::size_t float32_fixed_row_bytes = sites.size() * 12U
        + (shared.size() + haplotypes[0].size() + haplotypes[1].size()) * 16U;
    require(first.size() < float32_fixed_row_bytes * 2U / 3U,
            "MethDB did not materially reduce float32 fixed-row storage");

    TempFile file;
    write_bytes(file.path(), first);
    const Snapshot loaded(file.path(), binding, metadata);
    require(loaded.file_sha256() == htsim::crypto::sha256(bytes_of(first)),
            "MethDB reader did not retain its computed digest");
    const auto &loaded_reference = loaded.contig(0U);
    const auto &loaded_diploid = loaded.contig(1U);
    require(!loaded_reference.diploid
                && loaded_reference.reference_sites.size() == sites.size(),
            "MethDB reference row count changed");
    bool observed_unorm16_quantization = false;
    for (std::size_t index = 0U; index < sites.size(); ++index) {
        const CatalogSite &expected = sites[index];
        const CatalogSite &observed = loaded_reference.reference_sites[index];
        require(observed.reference_position == expected.reference_position
                    && observed.context == expected.context
                    && observed.methylation_source
                        == expected.methylation_source
                    && same_unorm16_probability(
                        observed.methylation_probability,
                        expected.methylation_probability),
                "MethDB reference row exceeded uint16 UNORM precision");
        require(observed.methylation_probability
                    == unorm16_probability(expected.methylation_probability),
                "MethDB reference row did not use uint16 UNORM encoding");
        observed_unorm16_quantization = observed_unorm16_quantization
            || observed.methylation_probability
                != expected.methylation_probability;
    }
    require(loaded_reference.reference_sites.front().methylation_probability == 0.0F
                && loaded_reference.reference_sites.back().methylation_probability
                    == 1.0F,
            "MethDB uint16 UNORM endpoints changed");
    require(observed_unorm16_quantization,
            "MethDB probabilities were not stored as uint16 UNORM");
    require(loaded_diploid.diploid
                && same_diploid_sites(loaded_diploid.shared_sites, shared)
                && same_diploid_sites(
                    loaded_diploid.haplotype_sites[0], haplotypes[0])
                && same_diploid_sites(
                    loaded_diploid.haplotype_sites[1], haplotypes[1]),
            "MethDB diploid rows exceeded uint16 UNORM precision");

    std::ostringstream table;
    htsim::methdb::export_snapshot_bed(file.path(), table);
    const std::string exported = table.str();
    require(exported.find("#format\tmethdb-bed\n") == 0U,
            "MethDB BED format marker changed");
    require(exported.find(
                "#contig\t0\tchrReference\t100001\t") != std::string::npos,
            "MethDB BED lost reference contig metadata");
    require(exported.find(
                "chrReference\t0\t1\tmethdb:reference:0\t0\t+\treference\t"
                "0\treference\t.\t.\tCG-C\tbeta\tshared\t0\t0\n")
                != std::string::npos,
            "MethDB BED lost a reference-backed site");
    require(exported.find(
                "#insertion\tchrDiploid\tshared\t9223372036854775812\t1\t0\t"
                "CHH-G\tbeta\tshared\t32768\t") != std::string::npos,
            "MethDB BED did not preserve an insertion origin");
    std::ostringstream repeated_table;
    htsim::methdb::export_snapshot_bed(file.path(), repeated_table);
    require(repeated_table.str() == exported,
            "MethDB BED changed across repeated exports");

    std::string wrong_version = first;
    wrong_version[magic_size] = static_cast<char>(
        htsim::methdb::methdb_version + 1U);
    write_bytes(file.path(), wrong_version);
    require_snapshot_error(
        [&] {(void)Snapshot(file.path(), binding, metadata);},
        "unsupported MethDB version was accepted");

    for (const std::string &invalid : {
             first.substr(0U, first.size() - 1U), first + "x"}) {
        write_bytes(file.path(), invalid);
        require_snapshot_error(
            [&] {
                (void)Snapshot(
                    file.path(),
                    binding,
                metadata);
            },
            "truncated or trailing MethDB bytes were accepted");
        require_snapshot_error(
            [&] {
                std::ostringstream invalid_table;
                htsim::methdb::export_snapshot_bed(
                    file.path(), invalid_table);
            },
            "truncated or trailing MethDB bytes were exported");
    }
}

} // namespace

int main()
{
    try {
        test_methdb_is_deterministic_compact_and_uses_unorm16_probabilities();
    } catch (const std::exception &error) {
        std::cerr << "methdb_snapshot_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
