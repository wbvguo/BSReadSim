#include <exception>
#include <iostream>
#include <string_view>

#include "bam.h"
#include "core.h"

namespace {

void print_help(std::ostream &output)
{
    output
        << "Usage: htsim-core [core contract options]\n"
        << "       htsim-core rrbs-catalog [core contract options]\n"
        << "       htsim-core methdb-catalog [core contract options]\n"
        << "       htsim-core --sam-to-bam LEVEL\n"
        << "The complete argv contract is documented in docs/core-cli.md.\n"
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
        if (argc == 3 && argv != nullptr && argv[1] != nullptr
            && argv[2] != nullptr
            && std::string_view(argv[1]) == "--sam-to-bam") {
            const std::string_view level_text(argv[2]);
            if (level_text.size() != 1U || level_text.front() < '0'
                || level_text.front() > '9') {
                throw htsim::bam::BamStreamError(
                    "--sam-to-bam LEVEL must be an integer in [0, 9]");
            }
            htsim::bam::sam_to_bam(
                "-", "-", static_cast<int>(level_text.front() - '0'));
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
            && std::string_view(argv[1]) == "methdb-catalog") {
            const htsim::core::CoreConfig config =
                htsim::core::parse_core_config(argc - 1, argv + 1);
            htsim::core::generate_methdb_catalog(config, std::cout);
            std::cout.flush();
            if (!std::cout) {
                throw htsim::core::CoreGeneratorError(
                    "failed while flushing the MethDB snapshot");
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
