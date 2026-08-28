#include "methdb.h"

#include <array>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>

namespace {

using htsim::methdb::AsmRecord;
using htsim::methdb::CatalogSite;
using htsim::methdb::DiploidMethylationCatalog;
using htsim::methdb::DiploidRuntimeArrays;
using htsim::methdb::MethylationCatalog;
using htsim::methdb::RuntimeSite;
using htsim::methdb::Snapshot;
using htsim::methdb::SnapshotWriter;
using htsim::model::HaplotypeMask;
using htsim::model::MethylationAllele;
using htsim::model::MethylationContext;
using htsim::model::MethylationSource;
using htsim::model::VariantKind;
using htsim::model::VariantSource;

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

std::uint32_t read_u32(const std::string &bytes, std::size_t &offset)
{
    require(offset + 4U <= bytes.size(), "test MethDB directory is truncated");
    std::uint32_t value = 0U;
    for (unsigned shift = 0U; shift < 32U; shift += 8U) {
        value |= static_cast<std::uint32_t>(
            static_cast<std::uint8_t>(bytes[offset++])) << shift;
    }
    return value;
}

std::uint64_t read_u64(const std::string &bytes, std::size_t &offset)
{
    require(offset + 8U <= bytes.size(), "test MethDB directory is truncated");
    std::uint64_t value = 0U;
    for (unsigned shift = 0U; shift < 64U; shift += 8U) {
        value |= static_cast<std::uint64_t>(
            static_cast<std::uint8_t>(bytes[offset++])) << shift;
    }
    return value;
}

std::vector<std::size_t> first_section_payloads(const std::string &bytes)
{
    constexpr std::size_t magic_size = sizeof(htsim::methdb::methdb_magic) - 1U;
    std::size_t offset = magic_size + 1U + 1U + 1U + 2U + 32U;
    const std::uint32_t contig_count = read_u32(bytes, offset);
    std::vector<std::size_t> result;
    result.reserve(contig_count);
    for (std::uint32_t contig = 0U; contig < contig_count; ++contig) {
        const std::uint32_t name_size = read_u32(bytes, offset);
        offset += name_size + 4U + 32U + 1U;
        require(offset < bytes.size(), "test MethDB contig header is truncated");
        const std::uint8_t section_count =
            static_cast<std::uint8_t>(bytes[offset++]);
        offset += 2U;
        for (std::uint8_t section = 0U; section < section_count; ++section) {
            offset += 1U + 1U + 2U + 4U + 4U;
            (void)read_u64(bytes, offset);
            const std::uint64_t compressed_size = read_u64(bytes, offset);
            offset += 32U;
            require(offset < bytes.size()
                        && compressed_size <= bytes.size() - offset,
                    "test MethDB section payload is truncated");
            if (section == 0U) {result.push_back(offset);}
            offset += static_cast<std::size_t>(compressed_size);
        }
    }
    return result;
}

std::vector<CatalogSite> large_reference_sites()
{
    std::vector<CatalogSite> sites;
    sites.reserve(50000U);
    for (std::uint32_t index = 0U; index < 50000U; ++index) {
        const std::uint32_t bits = index * UINT32_C(2654435761);
        sites.push_back({
            index * 2U,
            static_cast<std::uint16_t>(bits & UINT32_C(0xffff)),
            index % 2U == 0U
                ? MethylationContext::cg_c
                : MethylationContext::chh_g,
            index % 3U == 0U
                ? MethylationSource::beta
                : MethylationSource::cgmap,
        });
    }
    sites.front().probability_u16 = 0U;
    sites.back().probability_u16 = UINT16_MAX;
    return sites;
}

RuntimeSite runtime_site(
    std::uint32_t position,
    std::uint16_t probability,
    MethylationContext context,
    MethylationSource source,
    MethylationAllele allele = MethylationAllele::shared,
    bool reference_equivalent = true)
{
    return htsim::methdb::pack_runtime_site(
        position,
        probability,
        context,
        source,
        allele,
        reference_equivalent);
}

bool same_arrays(
    const DiploidRuntimeArrays &left,
    const DiploidRuntimeArrays &right)
{
    return left.reference_shared == right.reference_shared
        && left.reference_haplotypes == right.reference_haplotypes
        && left.insertion_shared == right.insertion_shared
        && left.insertion_haplotypes == right.insertion_haplotypes;
}

void test_v2_is_canonical_compact_and_lazy()
{
    const auto binding = htsim::crypto::sha256(bytes_of("binding-v2"));
    const std::vector<htsim::reference::ContigMetadata> metadata = {
        {"chrReference", 100001U,
         htsim::crypto::sha256(bytes_of("reference-1"))},
        {"chrDiploid", 100U,
         htsim::crypto::sha256(bytes_of("reference-2"))},
    };

    const std::vector<CatalogSite> reference_sites = large_reference_sites();
    const MethylationCatalog reference_catalog(100001U, reference_sites);
    const std::vector<CatalogSite> baseline_sites = {
        {0U, 0U, MethylationContext::cg_c, MethylationSource::beta},
        {2U, 32768U, MethylationContext::cg_c, MethylationSource::beta},
        {7U, 49151U, MethylationContext::chg_c,
         MethylationSource::cgmap},
    };
    const MethylationCatalog diploid_baseline(100U, baseline_sites);

    std::vector<htsim::variant::Variant> variants;
    htsim::variant::Variant linked;
    linked.contig_index = 1U;
    linked.reference_start = 50U;
    linked.reference_end = 51U;
    linked.kind = VariantKind::snv;
    linked.ref_bases = {3U};
    linked.alt_bases = {0U};
    linked.alt_haplotypes = HaplotypeMask::haplotype_1;
    linked.id = "linked";
    linked.source = VariantSource::vcf;
    variants.push_back(linked);
    htsim::variant::Variant insertion;
    insertion.contig_index = 1U;
    insertion.reference_start = 60U;
    insertion.reference_end = 60U;
    insertion.kind = VariantKind::insertion;
    insertion.alt_bases = {1U};
    insertion.alt_haplotypes = HaplotypeMask::both;
    insertion.id = "inserted";
    insertion.source = VariantSource::de_novo;
    variants.push_back(insertion);

    DiploidRuntimeArrays pre_asm_sites;
    pre_asm_sites.reference_shared = {
        runtime_site(0U, 0U, MethylationContext::cg_c,
                     MethylationSource::beta),
        runtime_site(2U, 32768U, MethylationContext::cg_c,
                     MethylationSource::beta),
        runtime_site(7U, 49151U, MethylationContext::chg_c,
                     MethylationSource::cgmap),
    };
    pre_asm_sites.insertion_shared = {
        runtime_site(
            (1U << 2U) | 0U,
            24576U,
            MethylationContext::chh_g,
            MethylationSource::beta,
            MethylationAllele::shared,
            false),
    };
    const DiploidMethylationCatalog pre_asm_catalog(
        1U, 100U, pre_asm_sites);

    AsmRecord asm_record;
    asm_record.target_reference_position = 2U;
    asm_record.linked_variant_position = 50U;
    asm_record.reference_probability_u16 = 8192U;
    asm_record.alternate_probability_u16 = 57343U;
    asm_record.context = MethylationContext::cg_c;
    asm_record.dinucleotide_second = 2U;
    asm_record.linked_reference_base = 3U;
    asm_record.linked_alternate_base = 0U;
    const std::vector<AsmRecord> asm_records = {asm_record};

    const auto serialize = [&]() {
        std::ostringstream output(std::ios::binary);
        SnapshotWriter writer(output, binding, 2U);
        writer.write_reference(metadata[0], reference_catalog);
        writer.write_diploid(
            metadata[1],
            diploid_baseline,
            pre_asm_catalog,
            variants,
            asm_records);
        writer.finish();
        return output.str();
    };
    const std::string first = serialize();
    const std::string repeated = serialize();
    require(first == repeated, "MethDB v2 bytes changed across writes");
    constexpr std::size_t magic_size = sizeof(htsim::methdb::methdb_magic) - 1U;
    require(first.substr(0U, magic_size) == htsim::methdb::methdb_magic,
            "MethDB v2 magic changed");
    require(first.size() > magic_size
                && static_cast<std::uint8_t>(first[magic_size])
                    == htsim::methdb::methdb_version,
            "MethDB v2 version byte changed");
    require(first.size() < reference_sites.size() * 8U,
            "MethDB v2 sections were not materially compressed");

    TempFile file;
    write_bytes(file.path(), first);
    const Snapshot snapshot(file.path(), binding, metadata);
    require(snapshot.file_sha256() == htsim::crypto::sha256(bytes_of(first)),
            "MethDB v2 file digest changed");
    require(snapshot.content_sha256()
                != htsim::crypto::Sha256Digest{},
            "MethDB v2 canonical content digest is empty");

    // Construction indexes only metadata. Each value below independently owns
    // the one contig that was requested.
    const auto loaded_reference = snapshot.contig(0U);
    require(!loaded_reference.diploid
                && loaded_reference.reference_sites.size()
                    == reference_sites.size(),
            "MethDB v2 reference contig changed");
    for (std::size_t index = 0U; index < reference_sites.size(); ++index) {
        const auto &left = loaded_reference.reference_sites[index];
        const auto &right = reference_sites[index];
        require(left.reference_position == right.reference_position
                    && left.probability_u16 == right.probability_u16
                    && left.context == right.context
                    && left.methylation_source == right.methylation_source,
                "MethDB v2 baseline did not round-trip exactly");
    }

    const auto loaded_diploid = snapshot.contig(1U);
    require(loaded_diploid.diploid
                && loaded_diploid.variants.size() == variants.size(),
            "MethDB v2 diploid/event layer changed");
    for (std::size_t index = 0U; index < variants.size(); ++index) {
        const auto &left = loaded_diploid.variants[index];
        const auto &right = variants[index];
        require(left.contig_index == right.contig_index
                    && left.reference_start == right.reference_start
                    && left.reference_end == right.reference_end
                    && left.kind == right.kind
                    && left.ref_bases == right.ref_bases
                    && left.alt_bases == right.alt_bases
                    && left.alt_haplotypes == right.alt_haplotypes
                    && left.id == right.id
                    && left.source == right.source,
                "MethDB v2 embedded event changed");
    }

    DiploidRuntimeArrays expected_sites;
    expected_sites.reference_shared = {
        runtime_site(0U, 0U, MethylationContext::cg_c,
                     MethylationSource::beta),
        runtime_site(7U, 49151U, MethylationContext::chg_c,
                     MethylationSource::cgmap),
    };
    expected_sites.insertion_shared = pre_asm_sites.insertion_shared;
    expected_sites.reference_haplotypes[0].push_back(runtime_site(
        2U,
        57343U,
        MethylationContext::cg_c,
        MethylationSource::asm_source,
        MethylationAllele::alternate_haplotype,
        true));
    expected_sites.reference_haplotypes[1].push_back(runtime_site(
        2U,
        8192U,
        MethylationContext::cg_c,
        MethylationSource::asm_source,
        MethylationAllele::reference_haplotype,
        true));
    require(
        same_arrays(loaded_diploid.diploid_sites, expected_sites),
        "MethDB v2 packed diploid runtime arrays changed");

    std::ostringstream table;
    htsim::methdb::export_snapshot_bed(file.path(), table);
    const std::string exported = table.str();
    require(exported.find("#format\tmethdb-bed-v2\n") == 0U,
            "MethDB BED v2 marker changed");
    require(exported.find("#content_sha256\t") != std::string::npos,
            "MethDB BED lost canonical content identity");
    require(exported.find("#contig\t1\tchrDiploid\t100\t")
                != std::string::npos,
            "MethDB BED lost diploid contig metadata");
    require(exported.find(
                "#variant\tchrDiploid\t0\t50\t51\tSNV\tT\tA\t1\t"
                "linked\tvcf\n") != std::string::npos,
            "MethDB BED lost embedded variant authority");
    require(exported.find("#insertion\tchrDiploid\tshared\t")
                != std::string::npos,
            "MethDB BED lost insertion identity");

    const std::vector<std::size_t> payloads = first_section_payloads(first);
    require(payloads.size() == 2U,
            "test could not locate independently framed contigs");
    std::string corrupt_second_contig = first;
    corrupt_second_contig[payloads[1]] ^= 1;
    write_bytes(file.path(), corrupt_second_contig);
    const Snapshot lazy(file.path(), binding, metadata);
    require(lazy.contig(0U).reference_sites.size() == reference_sites.size(),
            "corrupt inactive contig prevented lazy first-contig loading");
    require_snapshot_error(
        [&] {(void)lazy.contig(1U);},
        "corrupt requested contig section was accepted");

    std::string corrupt_content_root = first;
    corrupt_content_root.back() ^= 1;
    write_bytes(file.path(), corrupt_content_root);
    require_snapshot_error(
        [&] {(void)Snapshot(file.path(), binding, metadata);},
        "corrupt MethDB canonical content root was accepted");

    write_bytes(file.path(), first);

    const auto wrong_binding = htsim::crypto::sha256(bytes_of("wrong"));
    require_snapshot_error(
        [&] {(void)Snapshot(file.path(), wrong_binding, metadata);},
        "MethDB v2 accepted a wrong binding");
    std::vector<htsim::reference::ContigMetadata> wrong_metadata = metadata;
    wrong_metadata[1].length += 1U;
    require_snapshot_error(
        [&] {(void)Snapshot(file.path(), binding, wrong_metadata);},
        "MethDB v2 accepted wrong reference metadata");

    for (const std::uint8_t version : {
             std::uint8_t{1U},
             static_cast<std::uint8_t>(
                 htsim::methdb::methdb_version + 1U)}) {
        std::string wrong_version = first;
        wrong_version[magic_size] = static_cast<char>(version);
        write_bytes(file.path(), wrong_version);
        require_snapshot_error(
            [&] {(void)Snapshot(file.path(), binding, metadata);},
            "unsupported MethDB version was accepted");
    }

    for (const std::string &invalid : {
             first.substr(0U, first.size() - 1U), first + "x"}) {
        write_bytes(file.path(), invalid);
        require_snapshot_error(
            [&] {(void)Snapshot(file.path(), binding, metadata);},
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
        test_v2_is_canonical_compact_and_lazy();
    } catch (const std::exception &error) {
        std::cerr << "methdb_snapshot_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
