#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <future>
#include <iterator>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include "protocol.h"
#include "reference.h"

namespace {

using htsim::crypto::Sha256Digest;
using htsim::model::Bases;
using htsim::reference::Contig;
using htsim::reference::ReferenceError;
using htsim::reference::ReferenceSnapshot;

static_assert(
    std::is_same_v<
        decltype(htsim::reference::ContigMetadata{}.length),
        std::uint64_t>,
    "reference catalog lengths must match the reference input range");
static_assert(
    std::is_same_v<decltype(Contig{}.index), std::uint32_t>,
    "materialized contig indices must match the protocol u32 field");

[[noreturn]] void fail(const std::string &message)
{
    std::fprintf(stderr, "reference_test: %s\n", message.c_str());
    std::exit(EXIT_FAILURE);
}

void require(bool condition, const std::string &message)
{
    if (!condition) {fail(message);}
}

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-reference-XXXXXX";
        const int descriptor = mkstemp(pattern);
        if (descriptor < 0) {fail("failed to create a temporary file");}
        path_ = pattern;
        if (close(descriptor) != 0) {fail("failed to close a temporary file");}
    }

    ~TempFile()
    {
        if (!path_.empty()) {(void)unlink(path_.c_str());}
    }

    TempFile(const TempFile &) = delete;
    TempFile &operator=(const TempFile &) = delete;

    const std::string &path() const noexcept {return path_;}

private:
    std::string path_;
};

std::vector<std::uint8_t> bytes_of(std::string_view text)
{
    return std::vector<std::uint8_t>(text.begin(), text.end());
}

std::vector<std::uint8_t> read_bytes(const std::string &path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {fail("failed to open a temporary file for reading");}
    const std::vector<char> raw(
        (std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    if (input.bad()) {fail("failed while reading a temporary file");}
    return std::vector<std::uint8_t>(raw.begin(), raw.end());
}

void write_bytes(const std::string &path, const std::vector<std::uint8_t> &bytes)
{
    std::ofstream output(path, std::ios::binary | std::ios::trunc);
    if (!output) {fail("failed to open a temporary file for writing");}
    output.write(
        reinterpret_cast<const char *>(bytes.data()),
        static_cast<std::streamsize>(bytes.size()));
    if (!output) {fail("failed while writing a temporary file");}
}

void write_text(const std::string &path, std::string_view text)
{
    write_bytes(path, bytes_of(text));
}

std::vector<std::uint8_t> gzip_bytes(std::string_view text)
{
    TempFile file;
    gzFile output = gzopen(file.path().c_str(), "wb");
    if (output == nullptr) {fail("failed to open a gzip fixture");}
    const int written = gzwrite(
        output, text.data(), static_cast<unsigned int>(text.size()));
    if (written != static_cast<int>(text.size())) {
        (void)gzclose(output);
        fail("failed to write a gzip fixture");
    }
    if (gzclose(output) != Z_OK) {fail("failed to close a gzip fixture");}
    return read_bytes(file.path());
}

std::vector<std::uint8_t> concatenate(
    std::vector<std::uint8_t> left,
    const std::vector<std::uint8_t> &right)
{
    left.insert(left.end(), right.begin(), right.end());
    return left;
}

// Add a gzip FEXTRA field.  Header bytes are outside the uncompressed-data
// CRC, so this creates a valid member whose exact raw size is controllable.
std::vector<std::uint8_t> pad_gzip_member(
    std::vector<std::uint8_t> member,
    std::size_t target_size)
{
    require(member.size() >= 18 && member[0] == 0x1fU && member[1] == 0x8bU,
            "gzip padding fixture has no canonical header");
    require(member[3] == 0 && target_size >= member.size() + 2,
            "gzip padding target is invalid");
    const std::size_t extra_size = target_size - member.size() - 2;
    require(extra_size <= 0xffffU, "gzip padding exceeds FEXTRA capacity");
    member[3] = 0x04U;
    std::vector<std::uint8_t> extra(2 + extra_size, 0);
    extra[0] = static_cast<std::uint8_t>(extra_size);
    extra[1] = static_cast<std::uint8_t>(extra_size >> 8);
    member.insert(member.begin() + 10, extra.begin(), extra.end());
    require(member.size() == target_size, "gzip member padding size is wrong");
    return member;
}

Sha256Digest digest(const std::vector<std::uint8_t> &bytes)
{
    return htsim::crypto::sha256(bytes);
}

Sha256Digest digest(std::string_view text)
{
    return digest(bytes_of(text));
}

template <typename Operation>
std::string require_reference_error(Operation operation, const std::string &message)
{
    try {
        operation();
    } catch (const ReferenceError &error) {
        return error.what();
    } catch (const std::exception &error) {
        fail(message + "; wrong exception: " + error.what());
    }
    fail(message);
}

ReferenceSnapshot snapshot_for(const TempFile &file)
{
    return ReferenceSnapshot(file.path());
}

std::vector<Contig> collect_contigs(ReferenceSnapshot &snapshot)
{
    std::vector<Contig> result;
    snapshot.visit_contigs([&](const Contig &contig) {result.push_back(contig);});
    return result;
}

void verify_two_contigs(ReferenceSnapshot &snapshot, const Sha256Digest &file_digest)
{
    require(snapshot.file_sha256() == file_digest, "file digest was not retained");
    require(snapshot.catalog().size() == 2, "wrong number of catalog contigs");
    require(snapshot.catalog()[0].name == "chr1"
                && snapshot.catalog()[0].length == UINT64_C(5)
                && snapshot.catalog()[0].reference_sha256 == digest("ACGTN"),
            "first catalog entry is wrong");
    require(snapshot.catalog()[1].name == "chr2"
                && snapshot.catalog()[1].length == UINT64_C(4)
                && snapshot.catalog()[1].reference_sha256 == digest("NNTT"),
            "second catalog entry is wrong");

    const std::vector<Contig> first = collect_contigs(snapshot);
    require(first.size() == 2, "wrong number of visited contigs");
    require(first[0].index == 0 && first[0].name == "chr1"
                && first[0].length == UINT64_C(5)
                && first[0].bases == Bases({0, 1, 2, 3, 4}),
            "first materialized contig is wrong");
    require(first[1].index == 1 && first[1].name == "chr2"
                && first[1].length == UINT64_C(4)
                && first[1].bases == Bases({4, 4, 3, 3}),
            "second materialized contig is wrong");

    const std::vector<Contig> second = collect_contigs(snapshot);
    require(second.size() == first.size() && second[0].bases == first[0].bases
                && second[1].bases == first[1].bases,
            "a second validated visit changed the reference");
}

void test_plain_gzip_catalog_and_repeatable_visits()
{
    const std::string fasta = ">chr1 description\r\nacgtn\r\n>chr2\nNNtt\n";

    TempFile plain;
    write_text(plain.path(), fasta);
    const auto plain_bytes = read_bytes(plain.path());
    ReferenceSnapshot plain_snapshot(plain.path());
    verify_two_contigs(plain_snapshot, digest(plain_bytes));

    TempFile gzip;
    const auto compressed = gzip_bytes(fasta);
    write_bytes(gzip.path(), compressed);
    ReferenceSnapshot gzip_snapshot(gzip.path());
    verify_two_contigs(gzip_snapshot, digest(compressed));
}

void test_path_replacement_does_not_redirect_snapshot()
{
    TempFile original;
    write_text(original.path(), ">old\nACGT\n");
    ReferenceSnapshot snapshot = snapshot_for(original);

    TempFile replacement;
    write_text(replacement.path(), ">new\nNNNN\n");
    require(rename(replacement.path().c_str(), original.path().c_str()) == 0,
            "failed to atomically replace a snapshot path");

    const std::vector<Contig> contigs = collect_contigs(snapshot);
    require(contigs.size() == 1 && contigs[0].name == "old"
                && contigs[0].bases == Bases({0, 1, 2, 3}),
            "a path replacement redirected the open snapshot descriptor");
}

void test_non_regular_input_is_rejected_without_blocking()
{
    TempFile fifo;
    require(unlink(fifo.path().c_str()) == 0, "failed to prepare FIFO fixture");
    require(mkfifo(fifo.path().c_str(), 0600) == 0, "failed to create FIFO fixture");
    require_reference_error(
        [&] {(void)ReferenceSnapshot(fifo.path());},
        "a FIFO was accepted as a reference");
}

void test_gzip_member_and_trailing_data_rules()
{
    const auto continuation = concatenate(
        gzip_bytes(">chr1\nAC"), gzip_bytes("GT\n>chr2\nNN\n"));
    TempFile continued;
    write_bytes(continued.path(), continuation);
    ReferenceSnapshot continued_snapshot(continued.path());
    const auto continued_contigs = collect_contigs(continued_snapshot);
    require(continued_contigs.size() == 2
                && continued_contigs[0].bases == Bases({0, 1, 2, 3}),
            "concatenated gzip members did not form one decompressed FASTA stream");

    const auto empty_then_fasta = concatenate(
        gzip_bytes(""), gzip_bytes(">chr1\nA\n"));
    TempFile empty_member;
    write_bytes(empty_member.path(), empty_then_fasta);
    ReferenceSnapshot empty_snapshot(empty_member.path());
    require(empty_snapshot.catalog().size() == 1,
            "an empty leading gzip member broke the FASTA stream");

    std::vector<std::uint8_t> many_members;
    const auto empty_gzip = gzip_bytes("");
    for (std::size_t index = 0; index < 4096; ++index) {
        many_members.insert(
            many_members.end(), empty_gzip.begin(), empty_gzip.end());
    }
    many_members = concatenate(
        std::move(many_members), gzip_bytes(">chr1\nA\n"));
    TempFile many_member_file;
    write_bytes(many_member_file.path(), many_members);
    ReferenceSnapshot many_member_snapshot(many_member_file.path());
    require(many_member_snapshot.catalog().size() == 1,
            "a long iterative gzip-member chain was rejected");

    const auto boundary_member =
        pad_gzip_member(gzip_bytes(">chr1\nAC\n"), 65535);
    const auto split_magic = concatenate(
        boundary_member, gzip_bytes(">chr2\nGT\n"));
    TempFile split;
    write_bytes(split.path(), split_magic);
    ReferenceSnapshot split_snapshot(split.path());
    require(split_snapshot.catalog().size() == 2,
            "gzip member magic split across input buffers was rejected");

    const auto valid_member = gzip_bytes(">chr1\nAC\n");
    std::vector<std::vector<std::uint8_t>> trailing_cases;
    trailing_cases.push_back(concatenate(valid_member, {'x'}));
    trailing_cases.push_back(concatenate(valid_member, {0x1fU}));
    trailing_cases.push_back(concatenate(valid_member, {0x1fU, 0x8bU}));
    trailing_cases.push_back(concatenate(valid_member, bytes_of(">chr2\nA\n")));
    for (std::size_t index = 0; index < trailing_cases.size(); ++index) {
        TempFile trailing;
        write_bytes(trailing.path(), trailing_cases[index]);
        require_reference_error(
            [&] {
                (void)ReferenceSnapshot(trailing.path());
            },
            "gzip trailing data case was accepted");
    }

    auto corrupt_second = gzip_bytes(">chr2\nGT\n");
    require(corrupt_second.size() > 12, "second gzip member is too short");
    corrupt_second[corrupt_second.size() - 6] ^= 0x80U;
    const auto corrupt = concatenate(valid_member, corrupt_second);
    TempFile corrupt_file;
    write_bytes(corrupt_file.path(), corrupt);
    require_reference_error(
        [&] {(void)ReferenceSnapshot(corrupt_file.path());},
        "a corrupt second gzip member was accepted");

    auto truncated = valid_member;
    truncated.resize(truncated.size() - 4);
    TempFile truncated_file;
    write_bytes(truncated_file.path(), truncated);
    require_reference_error(
        [&] {(void)ReferenceSnapshot(truncated_file.path());},
        "a truncated gzip stream was accepted");
}

void test_fasta_line_and_header_semantics()
{
    TempFile no_final_newline;
    write_text(no_final_newline.path(), ">chr1\nACGT");
    ReferenceSnapshot valid = snapshot_for(no_final_newline);
    require(valid.catalog().size() == 1 && valid.catalog()[0].length == 4,
            "FASTA without a final newline was rejected");

    const std::vector<std::string> invalid_documents = {
        "ACGT\n",
        "@chr1\nACGT\n+\n!!!!\n",
        ">chr1\nAC\n\nGT\n",
        ">chr1\rAC\n",
        ">chr1\nAC>GT\n",
        ">chr1\nAC GT\n",
        ">chr1\nACGU\n",
        ";comment\n>chr1\nA\n",
        ">empty\n>next\nA\n",
        ">\nA\n",
        "> chr1\nA\n",
        ">chr1\x01" "description\nA\n",
    };
    for (const std::string &document : invalid_documents) {
        TempFile invalid;
        write_text(invalid.path(), document);
        const auto bytes = read_bytes(invalid.path());
        require_reference_error(
            [&] {(void)ReferenceSnapshot(invalid.path());},
            "an invalid FASTA line/header document was accepted");
    }

    TempFile invalid_utf8;
    write_bytes(invalid_utf8.path(), {'>', 0xffU, '\n', 'A', '\n'});
    const auto invalid_name = read_bytes(invalid_utf8.path());
    require_reference_error(
        [&] {(void)ReferenceSnapshot(invalid_utf8.path());},
        "an invalid UTF-8 contig name was accepted");

    // Descriptions are discarded at the input boundary, so only the retained
    // contig-name token is required to be UTF-8.
    TempFile ignored_description;
    write_bytes(
        ignored_description.path(), {'>', 'x', ' ', 0xffU, '\n', 'A', '\n'});
    const auto ignored_bytes = read_bytes(ignored_description.path());
    ReferenceSnapshot ignored(ignored_description.path());
    require(ignored.catalog()[0].name == "x",
            "discarded header description changed the retained contig name");

    std::string maximum_name(htsim::protocol::maximum_string_bytes, 'x');
    TempFile maximum;
    write_text(maximum.path(), ">" + maximum_name + "\nA\n");
    ReferenceSnapshot maximum_snapshot = snapshot_for(maximum);
    require(maximum_snapshot.catalog()[0].name.size() == maximum_name.size(),
            "a protocol-limit contig name was rejected");

    maximum_name.push_back('x');
    TempFile oversized;
    write_text(oversized.path(), ">" + maximum_name + "\nA\n");
    const auto oversized_bytes = read_bytes(oversized.path());
    require_reference_error(
        [&] {(void)ReferenceSnapshot(oversized.path());},
        "an oversized contig name was accepted");
}

void test_sha_and_io_boundaries()
{
    for (const std::size_t total_size : {
             std::size_t{55}, std::size_t{56}, std::size_t{63},
             std::size_t{64}, std::size_t{65}, std::size_t{65535},
             std::size_t{65536}, std::size_t{65537}, std::size_t{131075}}) {
        require(total_size >= 5, "invalid SHA boundary fixture size");
        const std::string fasta =
            ">x\n" + std::string(total_size - 4, 'a') + "\n";
        require(fasta.size() == total_size, "SHA boundary fixture length is wrong");
        TempFile file;
        write_text(file.path(), fasta);
        const auto bytes = read_bytes(file.path());
        ReferenceSnapshot snapshot(file.path());
        require(snapshot.file_sha256() == digest(bytes)
                    && snapshot.catalog()[0].reference_sha256
                        == digest(std::string(total_size - 4, 'A')),
                "streaming SHA-256 failed at a block or I/O boundary");
    }
}

void test_duplicate_empty_and_invalid_base_inputs()
{
    const std::vector<std::string> documents = {
        ">same\nA\n>same description\nC\n",
        "",
        ">chr1\nACGX\n",
    };
    for (const std::string &document : documents) {
        TempFile file;
        write_text(file.path(), document);
        const auto bytes = read_bytes(file.path());
        require_reference_error(
            [&] {(void)ReferenceSnapshot(file.path());},
            "duplicate, empty, or invalid-base FASTA was accepted");
    }
}

void test_visit_failure_poisoning_and_reentry()
{
    TempFile callback_file;
    write_text(callback_file.path(), ">chr1\nAC\n>chr2\nGT\n");
    ReferenceSnapshot callback_snapshot = snapshot_for(callback_file);
    try {
        callback_snapshot.visit_contigs([](const Contig &) {
            throw std::runtime_error("visitor sentinel");
        });
        fail("a visitor exception was swallowed");
    } catch (const std::runtime_error &error) {
        require(std::string(error.what()) == "visitor sentinel",
                "visitor exception was not propagated unchanged");
    }
    require_reference_error(
        [&] {callback_snapshot.visit_contigs([](const Contig &) {});},
        "a snapshot remained usable after a visitor failure");

    TempFile reentry_file;
    write_text(reentry_file.path(), ">chr1\nAC\n");
    ReferenceSnapshot reentry_snapshot = snapshot_for(reentry_file);
    reentry_snapshot.visit_contigs([&](const Contig &) {
        require_reference_error(
            [&] {reentry_snapshot.visit_contigs([](const Contig &) {});},
            "a reentrant snapshot visit was accepted");
    });
    require(collect_contigs(reentry_snapshot).size() == 1,
            "a rejected reentrant call poisoned a successful outer visit");

    require_reference_error(
        [&] {reentry_snapshot.visit_contigs({});},
        "an empty visitor was accepted");
    require(collect_contigs(reentry_snapshot).size() == 1,
            "an empty visitor poisoned an otherwise valid snapshot");
}

void test_concurrent_visit_is_rejected()
{
    TempFile file;
    write_text(file.path(), ">chr1\nAC\n");
    ReferenceSnapshot snapshot = snapshot_for(file);
    std::promise<void> entered;
    std::promise<void> release_promise;
    std::shared_future<void> release = release_promise.get_future().share();
    std::exception_ptr worker_error;
    std::thread worker([&] {
        try {
            snapshot.visit_contigs([&](const Contig &) {
                entered.set_value();
                release.wait();
            });
        } catch (...) {
            worker_error = std::current_exception();
        }
    });
    entered.get_future().wait();
    require_reference_error(
        [&] {snapshot.visit_contigs([](const Contig &) {});},
        "a concurrent snapshot visit was accepted");
    release_promise.set_value();
    worker.join();
    if (worker_error) {std::rethrow_exception(worker_error);}
    require(collect_contigs(snapshot).size() == 1,
            "a rejected concurrent visit poisoned the completed visit");
}

void test_in_place_change_is_detected_and_poisons()
{
    TempFile file;
    write_text(file.path(), ">chr1\nAAAA\n>chr2\nCCCC\n");
    ReferenceSnapshot snapshot = snapshot_for(file);
    bool changed = false;
    require_reference_error(
        [&] {
            snapshot.visit_contigs([&](const Contig &) {
                if (!changed) {
                    changed = true;
                    write_text(file.path(), ">chr1\nTTTT\n>chr2\nGGGG\n");
                }
            });
        },
        "an in-place reference change was accepted");
    require_reference_error(
        [&] {snapshot.visit_contigs([](const Contig &) {});},
        "an in-place reference change did not poison the snapshot");
}

} // namespace

int main()
{
    test_plain_gzip_catalog_and_repeatable_visits();
    test_path_replacement_does_not_redirect_snapshot();
    test_non_regular_input_is_rejected_without_blocking();
    test_gzip_member_and_trailing_data_rules();
    test_fasta_line_and_header_semantics();
    test_sha_and_io_boundaries();
    test_duplicate_empty_and_invalid_base_inputs();
    test_visit_failure_poisoning_and_reentry();
    test_concurrent_visit_is_rejected();
    test_in_place_change_is_detected_and_poisons();
    return EXIT_SUCCESS;
}
