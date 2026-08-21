#ifndef HTSIM_CORE_H
#define HTSIM_CORE_H

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
#include <iosfwd>

#include "utilities.h"
#include "protocol.h"

// ---- config --------------------------------------------------------

namespace htsim::core {

class CoreConfigError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class Technology {
    wgbs,
    rrbs,
    tbs,
};

enum class CoverageMode {
    uniform,
    profile,
    target_score,
};

struct BetaShape {
    double alpha = 0.0;
    double beta = 0.0;
};

struct CoreConfig {
    // Transport details are execution metadata, not scientific configuration.
    bool emit_details = false;
    std::uint32_t protocol_batch_fragments = 64;
    std::string run_id;
    crypto::Sha256Digest normalized_config_sha256 = {};
    std::uint64_t master_seed = 0;
    std::uint64_t catalog_seed = 0;
    std::string reference_path;
    crypto::Sha256Digest reference_sha256 = {};
    std::optional<std::string> vcf_path;
    std::optional<crypto::Sha256Digest> vcf_sha256;
    std::optional<std::string> cgmap_path;
    std::optional<crypto::Sha256Digest> cgmap_sha256;
    std::optional<std::string> bed_methyl_path;
    std::optional<crypto::Sha256Digest> bed_methyl_sha256;
    std::optional<std::string> methdb_path;
    std::optional<crypto::Sha256Digest> methdb_sha256;
    std::optional<std::string> asm_path;
    std::optional<crypto::Sha256Digest> asm_sha256;
    std::optional<std::string> asm_bed_path;
    std::optional<crypto::Sha256Digest> asm_bed_sha256;

    Technology technology = Technology::wgbs;
    bool paired_end = false;
    std::uint32_t read_length_1 = 0;
    std::optional<std::uint32_t> read_length_2;
    std::uint32_t insert_min = 0;
    std::uint32_t insert_mean = 0;
    std::uint32_t insert_max = 0;
    double insert_stddev = 0.0;
    std::optional<double> depth;
    std::optional<std::uint32_t> read_pairs;
    double max_ambiguous_fraction = 0.0;
    std::uint32_t chunk_size = 0;
    std::uint32_t core_workers = 1;

    double mutation_rate = 0.0;
    double indel_fraction = 0.0;
    double indel_extension_probability = 0.0;
    bool homozygous_only = false;

    bool collect_non_cpg = true;
    bool cgmap_pool = false;
    bool update_variant_boundaries = true;
    BetaShape beta_cg;
    BetaShape beta_chg;
    BetaShape beta_chh;

    CoverageMode coverage = CoverageMode::uniform;
    std::optional<std::string> coverage_profile_path;
    std::optional<std::string> coverage_profile_format;
    std::optional<std::string> coverage_profile_version;
    std::optional<crypto::Sha256Digest> coverage_profile_sha256;
    std::vector<std::string> rrbs_cut_sites;
    // External RRBS candidate scores are verified by exact regenerated-row
    // matching.  Unlike immutable scientific inputs, this exchange path has
    // no paired digest option and no hash in the BED format.
    std::optional<std::string> rrbs_candidate_bed_path;
    std::optional<std::string> tbs_bed_path;
    std::optional<crypto::Sha256Digest> tbs_bed_sha256;
    std::optional<double> tbs_center_stddev;
};

// One semantic boundary is shared by argv parsing and direct generator calls.
void validate_core_config(const CoreConfig &config);
CoreConfig parse_core_config(const std::vector<std::string> &arguments);
CoreConfig parse_core_config(int argc, char *argv[]);

} // namespace htsim::core

// ---- generator --------------------------------------------------------

namespace htsim::core {

inline constexpr const char core_version[] = "0.4.0";

class CoreGeneratorError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Run the generator and write the columnar protocol bytes to sink.
// All capability, reference, planning, allocation, and header-size checks run
// before the protocol Writer is constructed, so those failures leave sink at
// zero bytes. A failure during record emission leaves a truncated stream that
// the strict Python supervisor rejects together with the nonzero process exit.
protocol::Trailer generate_core_stream(
    const CoreConfig &config,
    std::ostream &sink);

// Serialize the exact RRBS catalog owned by C++ as a plain candidate BED.
// This is the first half of the external-model exchange workflow and emits no
// binary protocol frames.
void generate_rrbs_candidate_bed(
    const CoreConfig &config,
    std::ostream &sink);

// Serialize the exact normalized methylation probability catalog. The
// snapshot is independent of fragment/state/conversion/sequencing draws.
void generate_methdb_catalog(
    const CoreConfig &config,
    std::ostream &sink);

} // namespace htsim::core

#endif // HTSIM_CORE_H
