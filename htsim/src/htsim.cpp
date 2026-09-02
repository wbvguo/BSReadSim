#include <exception>
#include <iostream>
#include <string_view>

#include "bam.h"
#include "core.h"
#include "methdb.h"

namespace {

void print_help(std::ostream &output)
{
    output
        << "Usage: htsim-core [core contract options]\n"
        << "       htsim-core rrbs-catalog [core contract options]\n"
        << "       htsim-core methdb-build [core contract options]\n"
        << "       htsim-core variant-catalog [core contract options]\n"
        << "       htsim-core validate-inputs [core contract options]\n"
        << "       htsim-core methdb-export INPUT.methdb\n"
        << "       htsim-core --sam-to-bam LEVEL THREADS\n"
        << "Run bsreadsim --help for the supported public interface.\n"
        << "Output controls: --emit-details none|full.\n";
}

} // namespace

int main(int argc, char *argv[])
{
    std::ios::sync_with_stdio(false);
    if (argc == 2 && argv != nullptr && argv[1] != nullptr) {
        const std::string_view option(argv[1]);
        if (option == "--help") {
            print_help(std::cout);
            return 0;
        }
        if (option == "--version") {
            std::cout << "htsim-core " << htsim::core::core_version << '\n';
            return 0;
        }
    }

    try {
        if (argc == 4 && argv != nullptr && argv[1] != nullptr
            && argv[2] != nullptr && argv[3] != nullptr
            && std::string_view(argv[1]) == "--sam-to-bam") {
            const std::string_view level_text(argv[2]);
            if (level_text.size() != 1U || level_text.front() < '0'
                || level_text.front() > '9') {
                throw htsim::bam::BamStreamError(
                    "--sam-to-bam LEVEL must be an integer in [0, 9]");
            }
            const std::string_view threads_text(argv[3]);
            unsigned int compression_threads = 0U;
            if (threads_text.empty()) {
                throw htsim::bam::BamStreamError(
                    "--sam-to-bam THREADS must be an integer in [0, 64]");
            }
            for (const char value : threads_text) {
                if (value < '0' || value > '9') {
                    throw htsim::bam::BamStreamError(
                        "--sam-to-bam THREADS must be an integer in [0, 64]");
                }
                compression_threads = compression_threads * 10U
                    + static_cast<unsigned int>(value - '0');
                if (compression_threads > 64U) {
                    throw htsim::bam::BamStreamError(
                        "--sam-to-bam THREADS must be an integer in [0, 64]");
                }
            }
            htsim::bam::sam_to_bam(
                "-",
                "-",
                static_cast<int>(level_text.front() - '0'),
                static_cast<int>(compression_threads));
            return 0;
        }
        if (argc == 3 && argv != nullptr && argv[1] != nullptr
            && argv[2] != nullptr
            && std::string_view(argv[1]) == "methdb-export") {
            htsim::methdb::export_snapshot_bed(argv[2], std::cout);
            std::cout.flush();
            if (!std::cout) {
                throw htsim::methdb::SnapshotError(
                    "failed while flushing MethBED");
            }
            return 0;
        }
        if (argc >= 2 && argv != nullptr && argv[1] != nullptr
            && std::string_view(argv[1]) == "rrbs-catalog") {
            const htsim::core::CoreConfig config =
                htsim::core::parse_core_config(argc - 1, argv + 1);
            htsim::core::generate_rrbs_candidate_bed(config, std::cout);
            std::cout.flush();
            if (!std::cout) {
                throw htsim::core::CoreGeneratorError(
                    "failed while flushing the RRBS candidate BED");
            }
            return 0;
        }
        if (argc >= 2 && argv != nullptr && argv[1] != nullptr
            && std::string_view(argv[1]) == "methdb-build") {
            const htsim::core::CoreConfig config =
                htsim::core::parse_core_config(argc - 1, argv + 1);
            htsim::core::build_methdb_snapshot(config, std::cout);
            std::cout.flush();
            if (!std::cout) {
                throw htsim::core::CoreGeneratorError(
                    "failed while flushing the MethDB snapshot");
            }
            return 0;
        }
        if (argc >= 2 && argv != nullptr && argv[1] != nullptr
            && std::string_view(argv[1]) == "variant-catalog") {
            const htsim::core::CoreConfig config =
                htsim::core::parse_core_config(argc - 1, argv + 1);
            htsim::core::generate_variant_catalog_vcf(config, std::cout);
            std::cout.flush();
            if (!std::cout) {
                throw htsim::core::CoreGeneratorError(
                    "failed while flushing the variant VCF");
            }
            return 0;
        }
        if (argc >= 2 && argv != nullptr && argv[1] != nullptr
            && std::string_view(argv[1]) == "validate-inputs") {
            const htsim::core::CoreConfig config =
                htsim::core::parse_core_config(argc - 1, argv + 1);
            htsim::core::validate_inputs(config, std::cout);
            std::cout.flush();
            if (!std::cout) {
                throw htsim::core::CoreGeneratorError(
                    "failed while flushing the validation summary");
            }
            return 0;
        }
        const htsim::core::CoreConfig config =
            htsim::core::parse_core_config(argc, argv);
        (void)htsim::core::generate_core_stream(config, std::cout);
        std::cout.flush();
        if (!std::cout) {
            throw htsim::core::CoreGeneratorError(
                "failed while flushing the protocol stream");
        }
        return 0;
    } catch (const std::exception &error) {
        std::cerr << "htsim-core: " << error.what() << '\n';
        return 2;
    }
}
