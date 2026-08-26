#ifndef HTSIM_REFERENCE_H
#define HTSIM_REFERENCE_H

#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "types.h"
#include "utilities.h"

namespace htsim::reference {

class ReferenceError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Immutable metadata retained for the lifetime of a verified reference
// snapshot.  Entries are in FASTA order and contain no sequence allocation.
struct ContigMetadata {
    std::string name;
    std::uint64_t length = 0;
    // SHA-256 of the normalized uppercase ASCII sequence.
    crypto::Sha256Digest reference_sha256 = {};
};

// One materialized contig supplied temporarily by visit_contigs().  The
// reference passed to a visitor is valid only until that visitor returns.
struct Contig {
    std::uint32_t index = 0;
    std::string name;
    // Shared base encoding: 0=A, 1=C, 2=G, 3=T, 4=N.
    model::Bases bases;
    std::uint64_t length = 0;
    crypto::Sha256Digest reference_sha256 = {};
};

using ContigVisitor = std::function<void(const Contig &)>;

// A stable, fail-closed view of one regular reference file.
//
// Construction opens the path exactly once, hashes the raw file before any
// decompression, and builds a metadata-only catalog from the same descriptor.
// Each visit is another complete, independently hash- and identity-validated
// pass over that descriptor.  Successful visits may be repeated, but visits
// may not overlap or re-enter.  A parsing, validation, I/O, or visitor failure
// permanently poisons the snapshot; subsequent visits are rejected.
class ReferenceSnapshot {
public:
    explicit ReferenceSnapshot(const std::string &path);
    ~ReferenceSnapshot();

    ReferenceSnapshot(const ReferenceSnapshot &) = delete;
    ReferenceSnapshot &operator=(const ReferenceSnapshot &) = delete;
    ReferenceSnapshot(ReferenceSnapshot &&) = delete;
    ReferenceSnapshot &operator=(ReferenceSnapshot &&) = delete;

    const std::vector<ContigMetadata> &catalog() const noexcept;
    const crypto::Sha256Digest &file_sha256() const noexcept;

    void visit_contigs(const ContigVisitor &visitor);

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace htsim::reference

#endif // HTSIM_REFERENCE_H
