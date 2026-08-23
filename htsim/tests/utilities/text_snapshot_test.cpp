#include "utilities.h"

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

namespace {

using htsim::text::TextSnapshot;
using htsim::text::TextSnapshotError;

class TempFile {
public:
    TempFile()
    {
        char pattern[] = "/tmp/htsim-text-snapshot-XXXXXX";
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
    } catch (const TextSnapshotError &) {
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

std::vector<std::pair<std::uint64_t, std::string>> collect(TextSnapshot &snapshot)
{
    std::vector<std::pair<std::uint64_t, std::string>> result;
    snapshot.visit_lines([&](std::string_view line, std::uint64_t number) {
        result.emplace_back(number, line);
    });
    return result;
}

void test_plain_and_gzip_lines()
{
    const std::string text = "first\r\n\nthird";
    TempFile plain;
    const auto plain_bytes = bytes_of(text);
    write_bytes(plain.path(), plain_bytes);
    TextSnapshot plain_snapshot(
        plain.path(), htsim::crypto::sha256(plain_bytes));
    const auto expected = std::vector<std::pair<std::uint64_t, std::string>>{
        {1, "first"}, {2, ""}, {3, "third"}};
    require(collect(plain_snapshot) == expected, "plain line decoding changed");
    require(collect(plain_snapshot) == expected, "repeat visit changed lines");

    TempFile gzip;
    const auto compressed = gzip_bytes(text);
    write_bytes(gzip.path(), compressed);
    TextSnapshot gzip_snapshot(
        gzip.path(), htsim::crypto::sha256(compressed));
    require(collect(gzip_snapshot) == expected, "gzip line decoding changed");
}

void test_fail_closed_boundaries()
{
    TempFile input;
    const auto text = bytes_of("one\ntwo\n");
    write_bytes(input.path(), text);
    require_error(
        [&] {TextSnapshot snapshot(input.path(), {});},
        "digest mismatch was accepted");

    const auto bare_cr = bytes_of("one\rtwo");
    write_bytes(input.path(), bare_cr);
    TextSnapshot bare_snapshot(
        input.path(), htsim::crypto::sha256(bare_cr));
    require_error(
        [&] {(void)collect(bare_snapshot);},
        "bare carriage return was accepted");
    require_error(
        [&] {(void)collect(bare_snapshot);},
        "failed visit did not poison the snapshot");

    auto trailed = gzip_bytes("one\n");
    trailed.push_back(0U);
    write_bytes(input.path(), trailed);
    TextSnapshot trailing_snapshot(
        input.path(), htsim::crypto::sha256(trailed));
    require_error(
        [&] {(void)collect(trailing_snapshot);},
        "gzip trailing data was accepted");

    std::vector<std::uint8_t> long_line(
        htsim::text::maximum_line_bytes + 1U,
        static_cast<std::uint8_t>('A'));
    write_bytes(input.path(), long_line);
    TextSnapshot long_snapshot(
        input.path(), htsim::crypto::sha256(long_line));
    require_error(
        [&] {(void)collect(long_snapshot);},
        "oversized line was accepted");
}

void test_reentry_and_mutation_poison()
{
    TempFile input;
    const auto text = bytes_of("one\ntwo\n");
    write_bytes(input.path(), text);
    TextSnapshot snapshot(input.path(), htsim::crypto::sha256(text));
    bool reentry_rejected = false;
    snapshot.visit_lines([&](std::string_view, std::uint64_t line) {
        if (line != 1U) {return;}
        try {
            snapshot.visit_lines([](std::string_view, std::uint64_t) {});
        } catch (const TextSnapshotError &) {
            reentry_rejected = true;
        }
    });
    require(reentry_rejected, "reentrant visit was accepted");

    TextSnapshot mutating(input.path(), htsim::crypto::sha256(text));
    require_error(
        [&] {
            mutating.visit_lines([&](std::string_view, std::uint64_t line) {
                if (line == 1U) {write_bytes(input.path(), bytes_of("changed\n"));}
            });
        },
        "same-path mutation during visit was accepted");
}

} // namespace

int main()
{
    try {
        test_plain_and_gzip_lines();
        test_fail_closed_boundaries();
        test_reentry_and_mutation_poison();
        return EXIT_SUCCESS;
    } catch (const std::exception &error) {
        std::cerr << "text_snapshot_test: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
