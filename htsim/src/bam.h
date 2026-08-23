#ifndef HTSIM_BAM_H
#define HTSIM_BAM_H

#include <stdexcept>
#include <string>

namespace htsim::bam {

class BamStreamError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

// Convert an authenticated SAM stream to BAM through HTSlib. A path of "-"
// denotes standard input or standard output, enabling Python to stream records
// without materializing an intermediate SAM file.
void sam_to_bam(
    const std::string &input_path,
    const std::string &output_path,
    int compression_level);

} // namespace htsim::bam

#endif
