#include "tbs.h"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <unistd.h>
#include <zlib.h>

#include "utilities.h"

namespace {

using htsim::model::Bases;
using htsim::tbs::CandidateCatalog;
using htsim::tbs::DiploidCandidateCatalog;
using htsim::tbs::SamplingMode;
using htsim::tbs::Target;
using htsim::tbs::TargetFile;
using htsim::tbs::TbsCatalogError;

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-tbs-catalog-XXXXXX";
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
    } catch (const TbsCatalogError &) {
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
    if (!output) {throw std::runtime_error("temporary BED write failed");}
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
    return {{"chr1", 30, {}}, {"chr2", 20, {}}};
}

Bases encode(const std::string &text)
{
    Bases result;
    for (const char base : text) {
        switch (base) {
        case 'A': result.push_back(0); break;
        case 'C': result.push_back(1); break;
        case 'G': result.push_back(2); break;
        case 'T': result.push_back(3); break;
        case 'N': result.push_back(4); break;
        default: throw std::runtime_error("invalid test base");
        }
    }
    return result;
}

TargetFile load(const std::vector<std::uint8_t> &bytes, TempFile &file)
{
    write_bytes(file.path(), bytes);
    return TargetFile(
        file.path(), htsim::crypto::sha256(bytes), reference_catalog());
}

void test_bed6_projection_and_canonical_order()
{
    const std::string bed =
        "track name=x\r\n"
        "chr2\t7\t8\tminus\t2.5\t-\r\n"
        "# comment\r\n"
        "\r\n"
        "chr1\t10\t12\tunknown\t0\t.\r\n"
        "chr1\t2\t3\tplus\t1000\t+";
    TempFile plain;
    TargetFile targets = load(bytes_of(bed), plain);
    require(targets.target_count() == 3U, "BED metadata filtering changed");
    require(targets.targets(0).size() == 2U
                && targets.targets(0)[0].name == "plus"
                && targets.targets(0)[0].capture_strand
                    == htsim::model::CaptureStrand::forward
                && targets.targets(0)[1].name == "unknown"
                && targets.targets(1).size() == 1U
                && targets.targets(1)[0].score == 2.5
                && targets.targets(1)[0].capture_strand
                    == htsim::model::CaptureStrand::reverse,
            "BED rows were not projected in canonical reference order");

    TempFile compressed_file;
    const auto compressed = gzip_bytes(bed);
    TargetFile compressed_targets = load(compressed, compressed_file);
    require(compressed_targets.target_count() == targets.target_count()
                && compressed_targets.targets(0)[0].name == "plus"
                && compressed_targets.targets(1)[0].name == "minus",
            "gzip BED projection changed");
    require(compressed_targets.file_sha256()
                == htsim::crypto::sha256(compressed),
            "raw gzip digest was not retained");
}

void test_bed_rejections()
{
    const auto catalog = reference_catalog();
    const std::vector<std::string> invalid = {
        "chr1\t1\t2\tname\t0\n",
        "missing\t1\t2\tname\t0\t+\n",
        "chr1\t2\t2\tname\t0\t+\n",
        "chr1\t1\t31\tname\t0\t+\n",
        "chr1\t-1\t2\tname\t0\t+\n",
        "chr1\t1\t2\t\t0\t+\n",
        "chr1\t1\t2\tname\t-1\t+\n",
        "chr1\t1\t2\tname\tnan\t+\n",
        "chr1\t1\t2\tname\t0\t?\n",
        "chr1\t1\t2\tone\t0\t+\nchr1\t1\t2\ttwo\t1\t+\n",
        "# only metadata\n",
    };
    for (const std::string &bed : invalid) {
        TempFile file;
        const auto bytes = bytes_of(bed);
        write_bytes(file.path(), bytes);
        require_error(
            [&] {
                (void)TargetFile(
                    file.path(), htsim::crypto::sha256(bytes), catalog);
            },
            "invalid BED input was accepted: " + bed);
    }

    TempFile file;
    const auto valid = bytes_of("chr1\t1\t2\tname\t0\t+\n");
    write_bytes(file.path(), valid);
    require_error(
        [&] {(void)TargetFile(file.path(), {}, catalog);},
        "BED digest mismatch was accepted");
}

void test_fixed_center_candidates_and_ambiguity()
{
    const std::vector<Target> targets = {
        {0, 0, 1, "left-edge", 1.0,
         htsim::model::CaptureStrand::forward},
        {0, 5, 6, "forward", 2.0,
         htsim::model::CaptureStrand::forward},
        {0, 8, 10, "reverse", 3.0,
         htsim::model::CaptureStrand::reverse},
        {0, 11, 12, "right-edge", 4.0,
         htsim::model::CaptureStrand::unknown},
    };
    const CandidateCatalog catalog(
        encode("AACGTACGTAAA"), targets, 0, 0.0, 5, 2, true, 0.0);
    require(catalog.choice_count() == 2U,
            "boundary targets were not filtered from fixed-center candidates");
    const auto sampled = catalog.sample(0U, 7, 0, 64);
    const auto first_found = std::find_if(
        sampled.candidates.begin(), sampled.candidates.end(),
        [](const htsim::tbs::Candidate &candidate) {
            return candidate.target_ordinal == 1U;
        });
    const auto second_found = std::find_if(
        sampled.candidates.begin(), sampled.candidates.end(),
        [](const htsim::tbs::Candidate &candidate) {
            return candidate.target_ordinal == 2U;
        });
    require(first_found != sampled.candidates.end()
                && second_found != sampled.candidates.end(),
            "fixed candidates were not both reachable through uniform sampling");
    const auto &first = *first_found;
    const auto &second = *second_found;
    require(first.reference_start == 3U && first.reference_end == 8U
                && first.target_start == 5U && first.target_ordinal == 1U
                && first.capture_strand
                    == htsim::model::CaptureStrand::forward
                && second.reference_start == 7U && second.reference_end == 12U
                && second.target_start == 8U && second.target_ordinal == 2U
                && second.capture_strand
                    == htsim::model::CaptureStrand::reverse,
            "fixed-center placement or strand projection changed");
    require(first.haplotype_mask
                == htsim::model::HaplotypeMask::both,
            "TBS candidate lost the frozen availability bit mask");

    const std::vector<Target> n_target = {
        {0, 3, 4, "n", 0.0, htsim::model::CaptureStrand::unknown}};
    const CandidateCatalog rejected(
        encode("ANAAAA"), n_target, 0, 0.0, 4, 2, false, 0.0);
    require(rejected.choice_count() == 0U,
            "ambiguous TBS mate was accepted at zero threshold");
    const CandidateCatalog accepted(
        encode("ANAAAA"), n_target, 0, 0.0, 4, 2, false, 0.5);
    require(accepted.choice_count() == 1U,
            "ambiguous TBS fraction boundary changed");

    auto wrong_contig = n_target;
    wrong_contig[0].contig_index = 1;
    require_error(
        [&] {
            (void)CandidateCatalog(
                encode("ANAAAA"), wrong_contig, 0, 0.0, 4, 2, false, 0.5);
        },
        "target from a different contig was accepted");
    auto unsorted = targets;
    std::swap(unsorted[1], unsorted[2]);
    require_error(
        [&] {
            (void)CandidateCatalog(
                encode("AACGTACGTAAA"), unsorted, 0, 0.0, 5, 2, true, 0.0);
        },
        "non-canonical target order was accepted");
}

void test_sampling_is_chunk_independent()
{
    const std::vector<Target> targets = {
        {0, 3, 4, "a", 0.0, htsim::model::CaptureStrand::unknown},
        {0, 5, 6, "b", 0.0, htsim::model::CaptureStrand::unknown},
        {0, 7, 8, "c", 0.0, htsim::model::CaptureStrand::unknown},
    };
    const CandidateCatalog catalog(
        encode("AAAAAAAAAAAA"), targets, 0, 0.0, 3, 2, false, 0.0);
    const auto all = catalog.sample(0U, 91, 0, 11);
    auto joined = catalog.sample(0U, 91, 0, 4);
    const auto tail = catalog.sample(0U, 91, 4, 7);
    joined.candidates.insert(
        joined.candidates.end(), tail.candidates.begin(), tail.candidates.end());
    joined.skipped_count += tail.skipped_count;
    require(all.candidates.size() == joined.candidates.size()
                && all.skipped_count == joined.skipped_count,
            "TBS sampling count changed across chunk boundaries");
    for (std::size_t index = 0; index < all.candidates.size(); ++index) {
        const auto &left = all.candidates[index];
        const auto &right = joined.candidates[index];
        require(left.reference_start == right.reference_start
                    && left.reference_end == right.reference_end
                    && left.target_ordinal == right.target_ordinal,
                "TBS sample changed across chunk boundaries");
    }
    require_error(
        [&] {
            (void)catalog.sample(
                0U, 91, std::numeric_limits<std::uint64_t>::max(), 2);
        },
        "overflowing TBS ordinal range was accepted");
}

void test_target_score_sampling_is_exact_and_chunk_independent()
{
    const std::vector<Target> targets = {
        {0, 3, 4, "zero", 0.0,
         htsim::model::CaptureStrand::unknown},
        {0, 5, 6, "one", 1.0,
         htsim::model::CaptureStrand::forward},
        {0, 7, 8, "three", 3.0,
         htsim::model::CaptureStrand::reverse},
    };
    const CandidateCatalog catalog(
        encode("AAAAAAAAAAAA"),
        targets,
        0,
        0.0,
        3,
        2,
        false,
        0.0,
        SamplingMode::output_weight);
    require(catalog.choice_count() == 3U
                && catalog.allocation_weight() == 4U,
            "TBS target-score catalog lost exact eligible weights");

    const auto all = catalog.sample(0U, 19U, 0U, 16U);
    const std::vector<std::uint32_t> frozen = {
        1U, 2U, 2U, 1U, 2U, 2U, 2U, 2U,
        2U, 1U, 2U, 2U, 2U, 2U, 1U, 2U,
    };
    std::vector<std::uint32_t> observed;
    observed.reserve(all.candidates.size());
    const std::uint64_t key = htsim::rng::derive_key(
        19U, htsim::rng::Stage::fragment, 0U);
    for (std::size_t index = 0; index < all.candidates.size(); ++index) {
        const auto &candidate = all.candidates[index];
        observed.push_back(candidate.target_ordinal);
        const std::uint64_t draw = htsim::rng::bounded_integer(
            key, index, UINT64_C(1), UINT64_C(4));
        const std::uint32_t expected = draw == 0U ? 1U : 2U;
        require(candidate.target_ordinal == expected,
                "TBS target-score selection left its frozen RNG address");
    }
    require(observed == frozen,
            "TBS target-score exact v1 selection vector changed");
    require(std::find(observed.begin(), observed.end(), 0U) == observed.end(),
            "zero-weight TBS target was sampled");

    auto joined = catalog.sample(0U, 19U, 0U, 5U);
    const auto tail = catalog.sample(0U, 19U, 5U, 11U);
    joined.candidates.insert(
        joined.candidates.end(), tail.candidates.begin(), tail.candidates.end());
    require(joined.candidates.size() == all.candidates.size(),
            "TBS target-score chunking changed the sample count");
    for (std::size_t index = 0; index < all.candidates.size(); ++index) {
        require(joined.candidates[index].target_ordinal
                    == all.candidates[index].target_ordinal,
                "TBS target-score chunking changed target selection");
    }

    const CandidateCatalog uniform(
        encode("AAAAAAAAAAAA"), targets, 0, 0.0, 3, 2, false, 0.0);
    require(uniform.allocation_weight() == uniform.choice_count(),
            "uniform TBS allocation stopped counting choices");
}

void test_target_score_boundaries_fail_closed()
{
    const auto build = [](std::vector<Target> targets) {
        return CandidateCatalog(
            encode("AAAAAAAAAAAA"),
            targets,
            0,
            0.0,
            3,
            2,
            false,
            0.0,
            SamplingMode::output_weight);
    };
    require_error(
        [&] {
            (void)build({
                // This target is outside the fixed-center fragment domain.
                // Its score must still be validated before eligibility.
                {0, 0, 1, "fractional", 0.5,
                 htsim::model::CaptureStrand::unknown},
            });
        },
        "fractional TBS target-score weight was accepted");
    require_error(
        [&] {
            (void)build({
                {0, 3, 4, "wide", 4294967296.0,
                 htsim::model::CaptureStrand::unknown},
            });
        },
        "TBS target-score weight above uint32 was accepted");
    require_error(
        [&] {
            (void)build({
                {0, 3, 4, "max", 4294967295.0,
                 htsim::model::CaptureStrand::unknown},
                {0, 5, 6, "overflow", 1.0,
                 htsim::model::CaptureStrand::unknown},
            });
        },
        "TBS per-contig target-score sum overflow was accepted");

    const CandidateCatalog zero = build({
        {0, 3, 4, "zero-a", 0.0,
         htsim::model::CaptureStrand::unknown},
        {0, 5, 6, "zero-b", 0.0,
         htsim::model::CaptureStrand::unknown},
    });
    require(zero.choice_count() == 2U && zero.allocation_weight() == 0U,
            "all-zero target scores changed eligibility metadata");
    require_error(
        [&] {(void)zero.sample(0U, 0U, 0U, 1U);},
        "all-zero TBS target scores were sampled");
    require_error(
        [&] {
            (void)CandidateCatalog(
                encode("AAAAAAAAAAAA"),
                {},
                0,
                0.0,
                3,
                2,
                false,
                0.0,
                static_cast<SamplingMode>(2U));
        },
        "unknown TBS sampling mode was accepted");
}

void test_normal_center_sampling_is_chunk_independent()
{
    const std::vector<Target> targets = {
        {0, 1, 2, "left", 0.0, htsim::model::CaptureStrand::forward},
        {0, 6, 7, "middle", 0.0, htsim::model::CaptureStrand::reverse},
        {0, 11, 12, "right", 0.0, htsim::model::CaptureStrand::unknown},
    };
    const CandidateCatalog catalog(
        encode("AAAAAAAAAAAA"), targets, 0, 4.0, 5, 2, true, 0.0);
    require(catalog.choice_count() == 3U,
            "normal-center catalog lost a target that can displace inward");
    const auto all = catalog.sample(0U, 123, 0, 19);
    auto joined = catalog.sample(0U, 123, 0, 7);
    const auto tail = catalog.sample(0U, 123, 7, 12);
    joined.candidates.insert(
        joined.candidates.end(), tail.candidates.begin(), tail.candidates.end());
    joined.skipped_count += tail.skipped_count;
    require(all.candidates.size() == 19U
                && all.skipped_count > 0U
                && all.skipped_count == joined.skipped_count,
            "normal-center sampling count/skips changed across chunks");
    for (std::size_t index = 0; index < all.candidates.size(); ++index) {
        const auto &left = all.candidates[index];
        const auto &right = joined.candidates[index];
        require(left.reference_start == right.reference_start
                    && left.reference_end == right.reference_end
                    && left.target_ordinal == right.target_ordinal
                    && left.capture_strand == right.capture_strand,
                "normal-center candidate changed across chunks");
        require(left.reference_end - left.reference_start == 5U,
                "normal-center candidate changed fixed insert length");
    }

    require_error(
        [&] {
            (void)CandidateCatalog(
                encode("AAAAAAAAAAAA"),
                targets,
                0,
                std::numeric_limits<double>::infinity(),
                5,
                2,
                true,
                0.0);
        },
        "non-finite center standard deviation was accepted");
    require_error(
        [&] {
            (void)CandidateCatalog(
                encode("AAAAAAAAAAAA"),
                targets,
                0,
                -1.0,
                5,
                2,
                true,
                0.0);
        },
        "negative center standard deviation was accepted");

    const CandidateCatalog no_sequenceable_start(
        encode("NNNNNNNNNNNN"), targets, 0, 4.0, 5, 2, true, 0.0);
    require(no_sequenceable_start.choice_count() == 0U,
            "positive dispersion retained targets on an unusable contig");
    require_error(
        [&] {(void)no_sequenceable_start.sample(0U, 123, 0, 1);},
        "positive dispersion sampled a contig without a usable start");
}

htsim::reference::Contig make_contig(const std::string &sequence)
{
    htsim::reference::Contig contig;
    contig.index = 0U;
    contig.name = "chrDiploidTbs";
    contig.bases = encode(sequence);
    contig.length = contig.bases.size();
    return contig;
}

void test_diploid_tbs_centers_on_constructed_haplotypes()
{
    const auto contig = make_contig("AAAAAAAAAAAA");
    const std::vector<htsim::variant::Event> events = {
        {0U, 6U, 6U, htsim::model::VariantKind::insertion,
         {}, encode("TT"),
         htsim::model::HaplotypeMask::haplotype_1},
    };
    const htsim::variant::ContigVariants variants(
        contig.bases, events, contig.index);
    const std::vector<Target> targets = {
        {0U, 6U, 7U, "center", 1.0,
         htsim::model::CaptureStrand::forward},
    };
    const DiploidCandidateCatalog catalog(
        contig, variants, targets, 0.0, 5U, 2U, true, 0.0);
    require(catalog.choice_count() == 1U
                && catalog.allocation_weight() == 1U,
            "diploid TBS stopped treating one BED row as one target choice");
    const auto batch = catalog.sample(0U, 41U, 0U, 0U, 64U);
    const auto haplotype_0 = std::find_if(
        batch.candidates.begin(), batch.candidates.end(),
        [](const htsim::tbs::Candidate &candidate) {
            return candidate.haplotype_mask
                == htsim::model::HaplotypeMask::haplotype_1;
        });
    const auto haplotype_1 = std::find_if(
        batch.candidates.begin(), batch.candidates.end(),
        [](const htsim::tbs::Candidate &candidate) {
            return candidate.haplotype_mask
                == htsim::model::HaplotypeMask::haplotype_2;
        });
    require(haplotype_0 != batch.candidates.end()
                && haplotype_0->reference_start == 4U
                && haplotype_0->reference_end == 7U
                && haplotype_0->template_length == 5U,
            "TBS insertion-bearing haplotype was not fragmented physically");
    require(haplotype_1 != batch.candidates.end()
                && haplotype_1->reference_start == 4U
                && haplotype_1->reference_end == 9U
                && haplotype_1->template_length == 5U,
            "TBS reference haplotype center moved beside an insertion");
}

void test_deleted_tbs_center_removes_only_its_haplotype_bit()
{
    const auto contig = make_contig("AAAAAAAAAAAA");
    const std::vector<htsim::variant::Event> events = {
        {0U, 5U, 7U, htsim::model::VariantKind::deletion,
         encode("AA"), {},
         htsim::model::HaplotypeMask::haplotype_1},
    };
    const htsim::variant::ContigVariants variants(
        contig.bases, events, contig.index);
    const std::vector<Target> targets = {
        {0U, 6U, 7U, "deleted-center", 1.0,
         htsim::model::CaptureStrand::reverse},
    };
    const DiploidCandidateCatalog catalog(
        contig, variants, targets, 0.0, 5U, 2U, true, 0.0);
    require(catalog.choice_count() == 1U,
            "a BED center deleted from one haplotype remained eligible");
    const auto batch = catalog.sample(0U, 9U, 0U, 0U, 1U);
    require(batch.candidates.size() == 1U
                && batch.candidates[0].haplotype_mask
                    == htsim::model::HaplotypeMask::haplotype_2
                && batch.candidates[0].capture_strand
                    == htsim::model::CaptureStrand::reverse,
            "deleted TBS center cleared the wrong availability bit");
}

void test_diploid_tbs_selects_target_before_haplotype()
{
    const auto contig = make_contig("AAAAAAAAAAAAAAAAAAAA");
    const std::vector<htsim::variant::Event> events = {
        {0U, 5U, 7U, htsim::model::VariantKind::deletion,
         encode("AA"), {},
         htsim::model::HaplotypeMask::haplotype_1},
    };
    const htsim::variant::ContigVariants variants(
        contig.bases, events, contig.index);
    const std::vector<Target> targets = {
        {0U, 6U, 7U, "one-haplotype", 1.0,
         htsim::model::CaptureStrand::forward},
        {0U, 14U, 15U, "two-haplotypes", 3.0,
         htsim::model::CaptureStrand::reverse},
    };
    const DiploidCandidateCatalog weighted(
        contig,
        variants,
        targets,
        0.0,
        5U,
        2U,
        true,
        0.0,
        SamplingMode::output_weight);
    require(weighted.choice_count() == 2U,
            "diploid TBS expanded one target into haplotype choices");
    require(weighted.allocation_weight() == 4U,
            "diploid TBS multiplied target output weight by haplotype count");

    constexpr std::uint64_t seed = 271U;
    constexpr std::uint64_t first_fragment_ordinal = 17U;
    constexpr std::uint32_t output_count = 128U;
    const auto all = weighted.sample(
        contig.index,
        seed,
        0U,
        first_fragment_ordinal,
        output_count);
    const std::uint64_t fragment_key = htsim::rng::derive_key(
        seed, htsim::rng::Stage::fragment, contig.index);
    const std::uint64_t haplotype_key = htsim::rng::derive_key(
        seed, htsim::rng::Stage::haplotype, contig.index);
    bool found_haplotype_1 = false;
    bool found_haplotype_2 = false;
    for (std::uint32_t index = 0U; index < output_count; ++index) {
        const auto &candidate = all.candidates[index];
        const std::uint64_t target_draw = htsim::rng::bounded_integer(
            fragment_key, index, UINT64_C(1), UINT64_C(4));
        const std::uint32_t expected_target = target_draw == 0U ? 0U : 1U;
        require(candidate.target_ordinal == expected_target,
                "diploid TBS target draw depends on haplotype availability");
        if (expected_target == 0U) {
            require(candidate.haplotype_mask
                        == htsim::model::HaplotypeMask::haplotype_2,
                    "single-copy TBS target did not select its only haplotype");
            continue;
        }
        const bool expected_second = htsim::rng::bernoulli(
            haplotype_key,
            first_fragment_ordinal + index,
            UINT64_C(0),
            0.5);
        const auto expected_mask = expected_second
            ? htsim::model::HaplotypeMask::haplotype_2
            : htsim::model::HaplotypeMask::haplotype_1;
        require(candidate.haplotype_mask == expected_mask,
                "diploid TBS haplotype draw left its independent RNG stage");
        found_haplotype_1 = found_haplotype_1
            || expected_mask == htsim::model::HaplotypeMask::haplotype_1;
        found_haplotype_2 = found_haplotype_2
            || expected_mask == htsim::model::HaplotypeMask::haplotype_2;
    }
    require(found_haplotype_1 && found_haplotype_2,
            "two-copy TBS target did not exercise both conditional haplotypes");

    auto joined = weighted.sample(
        contig.index, seed, 0U, first_fragment_ordinal, 37U);
    const auto tail = weighted.sample(
        contig.index,
        seed,
        37U,
        first_fragment_ordinal + 37U,
        output_count - 37U);
    joined.candidates.insert(
        joined.candidates.end(),
        tail.candidates.begin(),
        tail.candidates.end());
    require(joined.candidates.size() == all.candidates.size(),
            "diploid TBS chunking changed the target-first sample count");
    for (std::size_t index = 0U; index < all.candidates.size(); ++index) {
        const auto &left = all.candidates[index];
        const auto &right = joined.candidates[index];
        require(left.target_ordinal == right.target_ordinal
                    && left.haplotype_mask == right.haplotype_mask
                    && left.reference_start == right.reference_start
                    && left.reference_end == right.reference_end,
                "diploid TBS target/haplotype selection changed across chunks");
    }

    const DiploidCandidateCatalog uniform(
        contig, variants, targets, 0.0, 5U, 2U, true, 0.0);
    require(uniform.choice_count() == 2U
                && uniform.allocation_weight() == 2U,
            "uniform diploid TBS stopped allocating one unit per BED target");

    const DiploidCandidateCatalog dispersed(
        contig,
        variants,
        targets,
        2.0,
        5U,
        2U,
        true,
        0.0,
        SamplingMode::output_weight);
    const auto displaced = dispersed.sample(
        contig.index, seed, 0U, first_fragment_ordinal, 32U);
    require(displaced.candidates.size() == 32U,
            "displaced diploid TBS did not emit its requested targets");
    for (std::size_t index = 0U; index < displaced.candidates.size(); ++index) {
        require(displaced.candidates[index].target_ordinal
                    == all.candidates[index].target_ordinal
                    && displaced.candidates[index].haplotype_mask
                        == all.candidates[index].haplotype_mask,
                "center displacement changed the prior target/haplotype draw");
    }
    auto displaced_joined = dispersed.sample(
        contig.index, seed, 0U, first_fragment_ordinal, 11U);
    const auto displaced_tail = dispersed.sample(
        contig.index, seed, 11U, first_fragment_ordinal + 11U, 21U);
    displaced_joined.candidates.insert(
        displaced_joined.candidates.end(),
        displaced_tail.candidates.begin(),
        displaced_tail.candidates.end());
    displaced_joined.skipped_count += displaced_tail.skipped_count;
    require(displaced_joined.skipped_count == displaced.skipped_count,
            "displaced diploid TBS skip accounting changed across chunks");
    for (std::size_t index = 0U; index < displaced.candidates.size(); ++index) {
        const auto &left = displaced.candidates[index];
        const auto &right = displaced_joined.candidates[index];
        require(left.target_ordinal == right.target_ordinal
                    && left.haplotype_mask == right.haplotype_mask
                    && left.reference_start == right.reference_start
                    && left.reference_end == right.reference_end,
                "displaced diploid TBS changed across chunks");
    }

    require_error(
        [&] {
            (void)weighted.sample(
                contig.index,
                seed,
                0U,
                std::numeric_limits<std::uint64_t>::max(),
                2U);
        },
        "overflowing diploid TBS fragment ordinal range was accepted");
}

} // namespace

int main()
{
    try {
        test_bed6_projection_and_canonical_order();
        test_bed_rejections();
        test_fixed_center_candidates_and_ambiguity();
        test_sampling_is_chunk_independent();
        test_target_score_sampling_is_exact_and_chunk_independent();
        test_target_score_boundaries_fail_closed();
        test_normal_center_sampling_is_chunk_independent();
        test_diploid_tbs_centers_on_constructed_haplotypes();
        test_deleted_tbs_center_removes_only_its_haplotype_bit();
        test_diploid_tbs_selects_target_before_haplotype();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "tbs_catalog_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
