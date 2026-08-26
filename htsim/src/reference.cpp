#include "reference.h"

#include <array>
#include <atomic>
#include <cstddef>
#include <limits>
#include <string_view>
#include <unordered_set>
#include <utility>

#include "utilities.h"

namespace htsim::reference {
namespace {

constexpr std::size_t maximum_contig_name_bytes = 1024U * 1024U;
constexpr std::size_t sequence_hash_buffer_bytes = 64U * 1024U;

bool valid_utf8(std::string_view text) noexcept
{
    std::size_t index = 0;
    while (index < text.size()) {
        const auto first = static_cast<unsigned char>(text[index]);
        if (first <= 0x7fU) {
            ++index;
            continue;
        }

        std::size_t length = 0;
        unsigned char second_min = 0x80U;
        unsigned char second_max = 0xbfU;
        if (first >= 0xc2U && first <= 0xdfU) {
            length = 2;
        } else if (first >= 0xe0U && first <= 0xefU) {
            length = 3;
            if (first == 0xe0U) {second_min = 0xa0U;}
            if (first == 0xedU) {second_max = 0x9fU;}
        } else if (first >= 0xf0U && first <= 0xf4U) {
            length = 4;
            if (first == 0xf0U) {second_min = 0x90U;}
            if (first == 0xf4U) {second_max = 0x8fU;}
        } else {
            return false;
        }
        if (index + length > text.size()) {return false;}
        const auto second = static_cast<unsigned char>(text[index + 1U]);
        if (second < second_min || second > second_max) {return false;}
        for (std::size_t continuation = 2; continuation < length;
             ++continuation) {
            const auto value = static_cast<unsigned char>(
                text[index + continuation]);
            if (value < 0x80U || value > 0xbfU) {return false;}
        }
        index += length;
    }
    return true;
}

struct ParsedContig {
    std::string name;
    model::Bases bases;
    std::uint64_t length = 0;
    crypto::Sha256Digest reference_sha256 = {};
};

using ParsedContigVisitor = std::function<void(ParsedContig &&)>;

enum class LineKind {
    start,
    header,
    sequence,
};

// FASTA semantics remain deliberately narrower than HTSlib's generic text
// support. HTSlib owns decompression and descriptor I/O; this parser preserves
// BSReadSim's ACGTN-only, UTF-8-name, ordering, hashing, and memory contracts.
class FastaParser {
public:
    FastaParser(
        bool materialize_bases,
        ParsedContigVisitor visitor,
        const std::vector<ContigMetadata> *expected_catalog = nullptr)
        : materialize_bases_(materialize_bases), visitor_(std::move(visitor)),
          expected_catalog_(expected_catalog)
    {
        if (!visitor_) {
            throw ReferenceError("FASTA parser requires a contig visitor");
        }
    }

    void feed(std::string_view chunk)
    {
        for (const char byte : chunk) {
            consume(static_cast<std::uint8_t>(
                static_cast<unsigned char>(byte)));
        }
    }

    void finish()
    {
        if (pending_carriage_return_) {
            throw ReferenceError("FASTA contains a bare carriage return");
        }
        if (line_kind_ == LineKind::header) {finish_header();}
        if (active_contig_) {finish_contig();}
        if (contig_count_ == 0) {
            throw ReferenceError("FASTA contains no contigs");
        }
    }

    std::uint64_t contig_count() const noexcept {return contig_count_;}

private:
    void consume(std::uint8_t value)
    {
        if (pending_carriage_return_) {
            if (value != static_cast<std::uint8_t>('\n')) {
                throw ReferenceError("FASTA contains a bare carriage return");
            }
            pending_carriage_return_ = false;
            finish_line();
            return;
        }
        if (value == static_cast<std::uint8_t>('\r')) {
            pending_carriage_return_ = true;
            return;
        }
        if (value == static_cast<std::uint8_t>('\n')) {
            finish_line();
            return;
        }

        if (line_kind_ == LineKind::start) {
            if (value == static_cast<std::uint8_t>('>')) {
                begin_header();
                return;
            }
            if (!active_contig_) {
                throw ReferenceError(
                    "FASTA sequence appears before its header");
            }
            line_kind_ = LineKind::sequence;
        }

        if (line_kind_ == LineKind::header) {
            consume_header_byte(value);
        } else {
            consume_base(value);
        }
    }

    void begin_header()
    {
        if (active_contig_) {finish_contig();}
        line_kind_ = LineKind::header;
        header_name_.clear();
        collecting_name_ = true;
    }

    void consume_header_byte(std::uint8_t value)
    {
        if (value == 0
            || (value < 0x20U && value != static_cast<std::uint8_t>('\t'))
            || value == 0x7fU) {
            throw ReferenceError(
                "FASTA header contains a control character");
        }
        if (collecting_name_
            && (value == static_cast<std::uint8_t>(' ')
                || value == static_cast<std::uint8_t>('\t'))) {
            if (header_name_.empty()) {
                throw ReferenceError("FASTA contig name is empty");
            }
            collecting_name_ = false;
            return;
        }
        if (!collecting_name_) {return;}
        if (header_name_.size() >= maximum_contig_name_bytes) {
            throw ReferenceError(
                "FASTA contig name exceeds the protocol limit");
        }
        header_name_.push_back(static_cast<char>(value));
    }

    void finish_header()
    {
        if (header_name_.empty()) {
            throw ReferenceError("FASTA contig name is empty");
        }
        if (!valid_utf8(header_name_)) {
            throw ReferenceError("FASTA contig name is not valid UTF-8");
        }
        if (!names_.insert(header_name_).second) {
            throw ReferenceError("FASTA contig names must be unique");
        }

        current_ = {};
        current_.name = header_name_;
        sequence_hash_ = crypto::Sha256{};
        normalized_buffer_size_ = 0;
        if (expected_catalog_ != nullptr) {
            if (contig_count_ >= expected_catalog_->size()) {
                throw ReferenceError(
                    "reference FASTA gained a contig during a visit");
            }
            const ContigMetadata &expected = expected_catalog_->at(
                static_cast<std::size_t>(contig_count_));
            if (current_.name != expected.name) {
                throw ReferenceError(
                    "reference FASTA contig name disagrees with the verified "
                    "catalog");
            }
            if (expected.length > current_.bases.max_size()) {
                throw ReferenceError(
                    "reference FASTA contig cannot be materialized on this "
                    "platform");
            }
            current_.bases.reserve(static_cast<std::size_t>(expected.length));
        }
        active_contig_ = true;
    }

    void consume_base(std::uint8_t value)
    {
        std::uint8_t encoded = 0;
        std::uint8_t normalized = 0;
        switch (value) {
        case 'A':
        case 'a': encoded = 0; normalized = 'A'; break;
        case 'C':
        case 'c': encoded = 1; normalized = 'C'; break;
        case 'G':
        case 'g': encoded = 2; normalized = 'G'; break;
        case 'T':
        case 't': encoded = 3; normalized = 'T'; break;
        case 'N':
        case 'n': encoded = 4; normalized = 'N'; break;
        default:
            throw ReferenceError(
                "FASTA sequence contains a non-ACGTN byte");
        }
        if (current_.length == crypto::maximum_sha256_input_bytes) {
            throw ReferenceError(
                "FASTA contig exceeds the SHA-256 length limit");
        }
        ++current_.length;
        normalized_buffer_[normalized_buffer_size_++] = normalized;
        if (normalized_buffer_size_ == normalized_buffer_.size()) {
            flush_sequence_hash();
        }
        if (materialize_bases_) {current_.bases.push_back(encoded);}
    }

    void finish_line()
    {
        if (line_kind_ == LineKind::start) {
            throw ReferenceError("FASTA contains an empty line");
        }
        if (line_kind_ == LineKind::header) {finish_header();}
        line_kind_ = LineKind::start;
    }

    void finish_contig()
    {
        if (current_.length == 0) {
            throw ReferenceError("FASTA contig sequence is empty");
        }
        if (contig_count_ == std::numeric_limits<std::uint32_t>::max()) {
            throw ReferenceError("FASTA contig count exceeds UINT32_MAX");
        }
        if (materialize_bases_
            && current_.bases.size() != current_.length) {
            throw ReferenceError(
                "materialized FASTA contig length overflowed size_t");
        }
        flush_sequence_hash();
        current_.reference_sha256 = sequence_hash_.digest();
        visitor_(std::move(current_));
        current_ = {};
        active_contig_ = false;
        ++contig_count_;
    }

    void flush_sequence_hash()
    {
        if (normalized_buffer_size_ == 0) {return;}
        sequence_hash_.update(
            normalized_buffer_.data(), normalized_buffer_size_);
        normalized_buffer_size_ = 0;
    }

    bool materialize_bases_ = false;
    ParsedContigVisitor visitor_;
    const std::vector<ContigMetadata> *expected_catalog_ = nullptr;
    LineKind line_kind_ = LineKind::start;
    bool pending_carriage_return_ = false;
    bool collecting_name_ = false;
    bool active_contig_ = false;
    std::string header_name_;
    ParsedContig current_;
    crypto::Sha256 sequence_hash_;
    std::array<std::uint8_t, sequence_hash_buffer_bytes> normalized_buffer_ = {};
    std::size_t normalized_buffer_size_ = 0;
    std::unordered_set<std::string> names_;
    std::uint64_t contig_count_ = 0;
};

class VisitStateGuard {
public:
    explicit VisitStateGuard(std::atomic_flag &active) : active_(&active) {}
    ~VisitStateGuard() {active_->clear(std::memory_order_release);}

    VisitStateGuard(const VisitStateGuard &) = delete;
    VisitStateGuard &operator=(const VisitStateGuard &) = delete;

private:
    std::atomic_flag *active_;
};

[[noreturn]] void translate_text_error(const text::TextSnapshotError &error)
{
    throw ReferenceError(error.what());
}

} // namespace

class ReferenceSnapshot::Impl {
public:
    explicit Impl(const std::string &path)
    {
        try {
            snapshot = std::make_unique<text::TextSnapshot>(path);
            file_sha256 = snapshot->file_sha256();
            FastaParser parser(false, [this](ParsedContig &&contig) {
                catalog.push_back(ContigMetadata{
                    std::move(contig.name),
                    contig.length,
                    contig.reference_sha256,
                });
            });
            snapshot->visit_chunks(
                [&](std::string_view chunk) {parser.feed(chunk);});
            parser.finish();
            snapshot->verify_unchanged();
            if (catalog.size() != parser.contig_count()) {
                throw ReferenceError(
                    "reference FASTA catalog count is inconsistent");
            }
        } catch (const text::TextSnapshotError &error) {
            translate_text_error(error);
        }
    }

    crypto::Sha256Digest file_sha256 = {};
    std::unique_ptr<text::TextSnapshot> snapshot;
    std::vector<ContigMetadata> catalog;
    std::atomic_flag visit_active = ATOMIC_FLAG_INIT;
    std::atomic<bool> poisoned{false};
};

ReferenceSnapshot::ReferenceSnapshot(const std::string &path)
    : impl_(std::make_unique<Impl>(path))
{
}

ReferenceSnapshot::~ReferenceSnapshot() = default;

const std::vector<ContigMetadata> &ReferenceSnapshot::catalog() const noexcept
{
    return impl_->catalog;
}

const crypto::Sha256Digest &ReferenceSnapshot::file_sha256() const noexcept
{
    return impl_->file_sha256;
}

void ReferenceSnapshot::visit_contigs(const ContigVisitor &visitor)
{
    if (!visitor) {
        throw ReferenceError("reference contig visitor must not be empty");
    }
    if (impl_->visit_active.test_and_set(std::memory_order_acquire)) {
        throw ReferenceError("reference snapshot visit is already active");
    }
    VisitStateGuard active_guard(impl_->visit_active);
    if (impl_->poisoned.load(std::memory_order_acquire)) {
        throw ReferenceError(
            "reference snapshot is poisoned after a prior failed visit");
    }

    try {
        std::uint64_t index = 0;
        FastaParser parser(
            true,
            [&](ParsedContig &&parsed) {
                if (index >= impl_->catalog.size()) {
                    throw ReferenceError(
                        "reference FASTA gained a contig during a visit");
                }
                const ContigMetadata &expected =
                    impl_->catalog[static_cast<std::size_t>(index)];
                if (parsed.name != expected.name
                    || parsed.length != expected.length
                    || parsed.reference_sha256 != expected.reference_sha256) {
                    throw ReferenceError(
                        "reference FASTA contig disagrees with the verified "
                        "catalog");
                }
                Contig contig{
                    static_cast<std::uint32_t>(index),
                    std::move(parsed.name),
                    std::move(parsed.bases),
                    parsed.length,
                    parsed.reference_sha256,
                };
                visitor(contig);
                ++index;
            },
            &impl_->catalog);

        try {
            impl_->snapshot->visit_chunks(
                [&](std::string_view chunk) {parser.feed(chunk);});
            parser.finish();
            // The final contig visitor runs at parser.finish(), after HTSlib
            // reaches EOF, so close the mutation window with one more raw pass.
            impl_->snapshot->verify_unchanged();
        } catch (const text::TextSnapshotError &error) {
            translate_text_error(error);
        }
        if (index != impl_->catalog.size()
            || parser.contig_count() != impl_->catalog.size()) {
            throw ReferenceError(
                "reference FASTA lost a contig during a visit");
        }
    } catch (...) {
        impl_->poisoned.store(true, std::memory_order_release);
        throw;
    }
}

} // namespace htsim::reference
