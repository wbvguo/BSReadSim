#include "bam.h"

#include <string>

#include <htslib/hts.h>
#include <htslib/sam.h>

namespace htsim::bam {
namespace {

void close_ignoring_errors(samFile *stream) noexcept
{
    if (stream != nullptr) {
        (void)sam_close(stream);
    }
}

} // namespace

void sam_to_bam(
    const std::string &input_path,
    const std::string &output_path,
    int compression_level)
{
    if (input_path.empty() || output_path.empty()) {
        throw BamStreamError("SAM/BAM stream paths must be non-empty");
    }
    if (compression_level < 0 || compression_level > 9) {
        throw BamStreamError("BAM compression level must be in [0, 9]");
    }

    samFile *input = sam_open(input_path.c_str(), "r");
    if (input == nullptr) {
        throw BamStreamError("failed to open the SAM input stream");
    }
    samFile *output = nullptr;
    sam_hdr_t *header = nullptr;
    bam1_t *record = nullptr;
    try {
        const std::string output_mode =
            std::string("wb") + static_cast<char>('0' + compression_level);
        output = sam_open(output_path.c_str(), output_mode.c_str());
        if (output == nullptr) {
            throw BamStreamError("failed to open the BAM output stream");
        }
        header = sam_hdr_read(input);
        if (header == nullptr) {
            throw BamStreamError("failed to parse the SAM header");
        }
        if (sam_hdr_write(output, header) < 0) {
            throw BamStreamError("failed to write the BAM header");
        }
        record = bam_init1();
        if (record == nullptr) {
            throw BamStreamError("failed to allocate a BAM record");
        }

        int status = 0;
        while ((status = sam_read1(input, header, record)) >= 0) {
            if (sam_write1(output, header, record) < 0) {
                throw BamStreamError("failed while writing a BAM record");
            }
        }
        if (status < -1) {
            throw BamStreamError("failed while parsing a SAM alignment record");
        }

        bam_destroy1(record);
        record = nullptr;
        sam_hdr_destroy(header);
        header = nullptr;
        if (sam_close(input) != 0) {
            input = nullptr;
            throw BamStreamError("failed while closing the SAM input stream");
        }
        input = nullptr;
        if (sam_close(output) != 0) {
            output = nullptr;
            throw BamStreamError("failed while finalizing the BAM output stream");
        }
        output = nullptr;
    } catch (...) {
        if (record != nullptr) {
            bam_destroy1(record);
        }
        if (header != nullptr) {
            sam_hdr_destroy(header);
        }
        close_ignoring_errors(input);
        close_ignoring_errors(output);
        throw;
    }
}

} // namespace htsim::bam
