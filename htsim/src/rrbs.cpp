#include "rrbs.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <map>
#include <ostream>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "variant.h"
#include "utilities.h"

// ---- catalog --------------------------------------------------------

namespace htsim::rrbs {
namespace {

std::uint8_t encode_base(char base)
{
    switch (base) {
    case 'A': return 0;
    case 'C': return 1;
    case 'G': return 2;
    case 'T': return 3;
    case 'N': return 4;
    default: throw RrbsCatalogError("cut motif contains an invalid base");
    }
}

std::uint32_t population_count(std::uint64_t value) noexcept
{
    std::uint32_t count = 0;
    while (value != 0) {
        value &= value - 1U;
        ++count;
    }
    return count;
}

std::uint32_t physical_copy_count(model::HaplotypeMask mask)
{
    switch (mask) {
    case model::HaplotypeMask::haplotype_1:
    case model::HaplotypeMask::haplotype_2:
        return 1U;
    case model::HaplotypeMask::both:
        return 2U;
    }
    throw RrbsCatalogError("RRBS candidate has an invalid haplotype mask");
}

std::uint32_t checked_allocation_weight(
    const std::vector<Candidate> &candidates)
{
    std::uint64_t total = 0U;
    for (const Candidate &candidate : candidates) {
        total += physical_copy_count(candidate.haplotype_mask);
        if (total > std::numeric_limits<std::uint32_t>::max()) {
            throw RrbsCatalogError(
                "RRBS physical-copy allocation weight exceeds uint32");
        }
    }
    return static_cast<std::uint32_t>(total);
}

std::string base36(std::uint32_t value)
{
    constexpr char digits[] = "0123456789abcdefghijklmnopqrstuvwxyz";
    std::string result;
    do {
        result.push_back(digits[value % 36U]);
        value /= 36U;
    } while (value != 0U);
    std::reverse(result.begin(), result.end());
    return result;
}

std::uint64_t envelope_key(const Candidate &candidate) noexcept
{
    return (static_cast<std::uint64_t>(candidate.reference_start) << 32U)
        | candidate.reference_end;
}

class AmbiguousBaseIndex {
public:
    explicit AmbiguousBaseIndex(const model::Bases &bases)
    {
        if (bases.size() > std::numeric_limits<std::uint32_t>::max()) {
            throw RrbsCatalogError("contig length exceeds uint32");
        }
        length_ = static_cast<std::uint32_t>(bases.size());
        words_.assign((bases.size() + 63U) / 64U, 0);
        for (std::size_t index = 0; index < bases.size(); ++index) {
            const std::uint8_t base = bases[index];
            if (base > 4) {
                throw RrbsCatalogError(
                    "contig contains a base outside protocol encoding");
            }
            if (base == 4) {
                words_[index / 64U] |= UINT64_C(1) << (index % 64U);
            }
        }
        prefix_.reserve(words_.size() + 1U);
        prefix_.push_back(0);
        for (const std::uint64_t word : words_) {
            const std::uint32_t count = population_count(word);
            if (count > std::numeric_limits<std::uint32_t>::max()
                    - prefix_.back()) {
                throw RrbsCatalogError("ambiguous-base count exceeds uint32");
            }
            prefix_.push_back(prefix_.back() + count);
        }
    }

    std::uint32_t count(std::uint32_t begin, std::uint32_t end) const
    {
        if (begin > end || end > length_) {
            throw RrbsCatalogError("ambiguous-base query is outside the contig");
        }
        return rank(end) - rank(begin);
    }

private:
    std::uint32_t rank(std::uint32_t end) const
    {
        const std::size_t word = static_cast<std::size_t>(end / 64U);
        const unsigned bit = static_cast<unsigned>(end % 64U);
        std::uint32_t result = prefix_.at(word);
        if (bit != 0U) {
            const std::uint64_t mask = (UINT64_C(1) << bit) - 1U;
            result += population_count(words_.at(word) & mask);
        }
        return result;
    }

    std::vector<std::uint64_t> words_;
    std::vector<std::uint32_t> prefix_;
    std::uint32_t length_ = 0;
};

bool motif_matches(
    const model::Bases &bases,
    std::size_t start,
    const CutSite &site)
{
    for (std::size_t offset = 0; offset < site.motif.size(); ++offset) {
        const std::uint8_t expected = site.motif[offset];
        const std::uint8_t observed = bases[start + offset];
        if (expected == 4) {
            if (observed >= 4) {return false;}
        } else if (observed != expected) {
            return false;
        }
    }
    return true;
}

std::uint32_t maximum_ambiguous_count(
    std::uint32_t read_length,
    double fraction)
{
    if (!std::isfinite(fraction) || fraction < 0.0 || fraction > 1.0) {
        throw RrbsCatalogError(
            "maximum ambiguous fraction must be finite and in [0, 1]");
    }
    return static_cast<std::uint32_t>(std::floor(
        fraction * static_cast<double>(read_length)));
}

bool sequenceable(
    const AmbiguousBaseIndex &ambiguity,
    std::uint32_t start,
    std::uint32_t end,
    std::uint32_t read_length,
    bool paired_end,
    std::uint32_t maximum_n)
{
    if (ambiguity.count(start, start + read_length) > maximum_n) {
        return false;
    }
    return !paired_end
        || ambiguity.count(end - read_length, end) <= maximum_n;
}

std::vector<Candidate> build_candidates(
    const model::Bases &bases,
    const std::vector<CutSite> &cut_sites,
    std::uint32_t minimum_insert_length,
    std::uint32_t maximum_insert_length,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction,
    model::HaplotypeMask haplotype_mask,
    const haplotype::HaplotypeLayout *layout)
{
    if (minimum_insert_length == 0U
        || maximum_insert_length < minimum_insert_length) {
        throw RrbsCatalogError("RRBS insert-length range is invalid");
    }
    if (read_length == 0U || read_length > minimum_insert_length) {
        throw RrbsCatalogError("read length must fit every RRBS candidate");
    }
    if (!model::is_haplotype_mask(
            static_cast<std::uint8_t>(haplotype_mask))) {
        throw RrbsCatalogError("RRBS candidate haplotype mask is invalid");
    }
    if (layout != nullptr
        && (layout->length() != bases.size()
            || !model::mask_contains(
                haplotype_mask, layout->haplotype()))) {
        throw RrbsCatalogError("RRBS haplotype layout identity is invalid");
    }
    const std::uint32_t maximum_n = maximum_ambiguous_count(
        read_length, max_ambiguous_fraction);
    const AmbiguousBaseIndex ambiguity(bases);
    const std::vector<CutPosition> cuts = find_cut_positions(bases, cut_sites);
    if (cuts.size() < 2U) {return {};}

    std::vector<std::uint64_t> recognition_prefix(cuts.size() + 1U, 0U);
    for (std::size_t index = 0U; index < cuts.size(); ++index) {
        recognition_prefix[index + 1U] = recognition_prefix[index]
            + cuts[index].recognition_count;
    }
    std::vector<std::uint32_t> gc_prefix(bases.size() + 1U, 0U);
    for (std::size_t index = 0U; index < bases.size(); ++index) {
        const std::uint32_t increment =
            bases[index] == 1U || bases[index] == 2U ? 1U : 0U;
        gc_prefix[index + 1U] = gc_prefix[index] + increment;
    }

    std::vector<Candidate> candidates;
    for (std::size_t left = 0U; left + 1U < cuts.size(); ++left) {
        for (std::size_t right = left + 1U; right < cuts.size(); ++right) {
            const std::uint32_t length = cuts[right].position
                - cuts[left].position;
            if (length < minimum_insert_length) {continue;}
            if (length > maximum_insert_length) {break;}
            if (!sequenceable(
                    ambiguity,
                    cuts[left].position,
                    cuts[right].position,
                    read_length,
                    paired_end,
                    maximum_n)) {
                continue;
            }

            std::uint32_t reference_start = cuts[left].position;
            std::uint32_t reference_end = cuts[right].position;
            bool include_start_insertion = true;
            bool include_end_insertion = false;
            if (layout != nullptr) {
                const auto start = layout->boundary(cuts[left].position);
                const auto end = layout->boundary(cuts[right].position);
                // The wire contract carries each insertion as one whole typed
                // event. A cut inside its ALT bases would lose event identity,
                // so that physical candidate is excluded rather than truncated.
                if (!start || !end) {continue;}
                reference_start = start->right_reference_start;
                reference_end = end->left_reference_end;
                include_start_insertion =
                    start->include_insertion_in_right_fragment;
                include_end_insertion =
                    end->include_insertion_in_left_fragment;
                if (reference_start > reference_end) {continue;}
            }

            const std::uint64_t count = recognition_prefix[right + 1U]
                - recognition_prefix[left];
            if (count > std::numeric_limits<std::uint32_t>::max()) {
                throw RrbsCatalogError(
                    "candidate restriction count exceeds uint32");
            }
            if (candidates.size()
                == std::numeric_limits<std::uint32_t>::max()) {
                throw RrbsCatalogError("RRBS candidate count exceeds uint32");
            }
            candidates.push_back(Candidate{
                reference_start,
                reference_end,
                length,
                gc_prefix[cuts[right].position]
                    - gc_prefix[cuts[left].position],
                static_cast<std::uint32_t>(count),
                haplotype_mask,
                include_start_insertion,
                include_end_insertion,
            });
        }
    }
    return candidates;
}

std::vector<std::uint32_t> sample_indices(
    const std::vector<Candidate> &candidates,
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count)
{
    if (output_count == 0U) {return {};}
    if (candidates.empty()) {
        throw RrbsCatalogError("cannot sample from zero RRBS candidates");
    }
    if (static_cast<std::uint64_t>(output_count) - 1U
        > std::numeric_limits<std::uint64_t>::max()
            - first_candidate_ordinal) {
        throw RrbsCatalogError("candidate ordinal range exceeds uint64");
    }
    const std::uint64_t key = rng::derive_key(
        master_seed, rng::Stage::fragment, contig_index);
    const std::uint32_t first_mass =
        physical_copy_count(candidates.front().haplotype_mask);
    const bool equal_mass = std::all_of(
        candidates.begin(), candidates.end(),
        [first_mass](const Candidate &candidate) {
            return physical_copy_count(candidate.haplotype_mask) == first_mass;
        });
    std::vector<std::uint64_t> cumulative;
    std::uint64_t total_mass = 0U;
    if (!equal_mass) {
        cumulative.reserve(candidates.size());
        for (const Candidate &candidate : candidates) {
            total_mass += physical_copy_count(candidate.haplotype_mask);
            cumulative.push_back(total_mass);
        }
    }
    std::vector<std::uint32_t> indices;
    indices.reserve(output_count);
    for (std::uint32_t index = 0U; index < output_count; ++index) {
        const std::uint64_t ordinal = first_candidate_ordinal + index;
        if (equal_mass) {
            const std::uint64_t rank = rng::bounded_integer(
                key, ordinal, UINT64_C(1), candidates.size());
            indices.push_back(static_cast<std::uint32_t>(rank));
        } else {
            const std::uint64_t draw = rng::bounded_integer(
                key, ordinal, UINT64_C(1), total_mass);
            const auto found = std::upper_bound(
                cumulative.begin(), cumulative.end(), draw);
            if (found == cumulative.end()) {
                throw RrbsCatalogError(
                    "RRBS physical-copy selection exceeded its cumulative mass");
            }
            indices.push_back(static_cast<std::uint32_t>(
                found - cumulative.begin()));
        }
    }
    return indices;
}

} // namespace

std::vector<CutSite> parse_cut_sites(
    const std::vector<std::string> &declarations)
{
    if (declarations.empty()) {
        throw RrbsCatalogError("RRBS requires at least one cut motif");
    }
    std::set<std::string> unique;
    std::vector<CutSite> result;
    result.reserve(declarations.size());
    for (const std::string &declaration : declarations) {
        if (!unique.insert(declaration).second) {
            throw RrbsCatalogError("duplicate RRBS cut motif: " + declaration);
        }
        const std::size_t separator = declaration.find('|');
        if (separator == std::string::npos
            || separator + 1U == declaration.size()
            || declaration.find('|', separator + 1U) != std::string::npos) {
            throw RrbsCatalogError(
                "cut motif must contain one | and a non-empty right side");
        }
        const std::size_t motif_length = declaration.size() - 1U;
        if (motif_length == 0 || motif_length > maximum_motif_length) {
            throw RrbsCatalogError("cut motif length is outside [1, 1024]");
        }
        CutSite site;
        site.cut_offset = static_cast<std::uint32_t>(separator);
        site.motif.reserve(motif_length);
        for (const char base : declaration) {
            if (base != '|') {site.motif.push_back(encode_base(base));}
        }
        result.push_back(std::move(site));
    }
    return result;
}

std::vector<CutPosition> find_cut_positions(
    const model::Bases &contig_bases,
    const std::vector<CutSite> &cut_sites)
{
    if (contig_bases.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw RrbsCatalogError("contig length exceeds uint32");
    }
    if (cut_sites.empty()) {
        throw RrbsCatalogError("RRBS requires at least one cut motif");
    }
    for (const std::uint8_t base : contig_bases) {
        if (base > 4) {
            throw RrbsCatalogError(
                "contig contains a base outside protocol encoding");
        }
    }

    std::vector<std::uint32_t> raw_positions;
    for (const CutSite &site : cut_sites) {
        if (site.motif.empty() || site.motif.size() > maximum_motif_length
            || site.cut_offset > site.motif.size()) {
            throw RrbsCatalogError("cut motif representation is invalid");
        }
        if (site.motif.size() > contig_bases.size()) {continue;}
        const std::size_t last = contig_bases.size() - site.motif.size();
        for (std::size_t start = 0; start <= last; ++start) {
            if (!motif_matches(contig_bases, start, site)) {continue;}
            const std::size_t position = start + site.cut_offset;
            if (position > std::numeric_limits<std::uint32_t>::max()) {
                throw RrbsCatalogError("cut position exceeds uint32");
            }
            raw_positions.push_back(static_cast<std::uint32_t>(position));
        }
    }
    std::sort(raw_positions.begin(), raw_positions.end());
    std::vector<CutPosition> result;
    for (const std::uint32_t position : raw_positions) {
        if (result.empty() || result.back().position != position) {
            result.push_back({position, 1});
        } else {
            if (result.back().recognition_count
                == std::numeric_limits<std::uint32_t>::max()) {
                throw RrbsCatalogError("recognition count exceeds uint32");
            }
            ++result.back().recognition_count;
        }
    }
    return result;
}

std::vector<std::string> candidate_ids(
    std::string_view contig_name,
    const std::vector<Candidate> &candidates)
{
    if (contig_name.empty()) {
        throw RrbsCatalogError("RRBS candidate contig name is empty");
    }
    for (const char value : contig_name) {
        const auto byte = static_cast<unsigned char>(value);
        if (byte == 0U || byte <= 0x20U || byte == 0x7fU) {
            throw RrbsCatalogError(
                "RRBS candidate contig name contains whitespace or a control byte");
        }
    }
    if (candidates.size() > std::numeric_limits<std::uint32_t>::max()) {
        throw RrbsCatalogError("RRBS candidate count exceeds uint32");
    }

    std::unordered_map<std::uint64_t, std::uint32_t> totals;
    totals.reserve(candidates.size());
    for (const Candidate &candidate : candidates) {
        std::uint32_t &count = totals[envelope_key(candidate)];
        if (count == std::numeric_limits<std::uint32_t>::max()) {
            throw RrbsCatalogError(
                "RRBS reference-envelope multiplicity exceeds uint32");
        }
        ++count;
    }

    std::unordered_map<std::uint64_t, std::uint32_t> observed;
    observed.reserve(totals.size());
    std::vector<std::string> result;
    result.reserve(candidates.size());
    std::unordered_set<std::string> unique;
    unique.reserve(candidates.size());
    for (const Candidate &candidate : candidates) {
        const std::uint64_t key = envelope_key(candidate);
        std::string id = std::string(contig_name) + ":"
            + std::to_string(candidate.reference_start) + "-"
            + std::to_string(candidate.reference_end);
        if (totals.at(key) > 1U) {
            id += "~" + base36(observed[key]++);
        }
        if (!unique.insert(id).second) {
            throw RrbsCatalogError("RRBS candidate ID is not unique");
        }
        result.push_back(std::move(id));
    }
    return result;
}

void write_candidate_bed_header(std::ostream &output)
{
    output
        << "#chrom\tstart\tend\tcandidate_id\tscore\tstrand"
        << "\thaplotype_mask\ttemplate_length\tgc_count"
        << "\trestriction_site_count\n";
    if (!output) {
        throw RrbsCatalogError("failed while writing the RRBS candidate header");
    }
}

void write_candidate_bed_contig(
    std::ostream &output,
    std::string_view contig_name,
    const std::vector<Candidate> &candidates)
{
    const std::vector<std::string> ids = candidate_ids(
        contig_name, candidates);
    for (std::size_t index = 0U; index < candidates.size(); ++index) {
        const Candidate &candidate = candidates[index];
        (void)physical_copy_count(candidate.haplotype_mask);
        output
            << contig_name << '\t'
            << candidate.reference_start << '\t'
            << candidate.reference_end << '\t'
            << ids[index] << "\t1\t.\t"
            << static_cast<unsigned int>(
                   static_cast<std::uint8_t>(candidate.haplotype_mask))
            << '\t' << candidate.template_length
            << '\t' << candidate.gc_count
            << '\t' << candidate.restriction_site_count << '\n';
    }
    if (!output) {
        throw RrbsCatalogError("failed while writing RRBS candidate rows");
    }
}

namespace {

bool candidate_metadata_line(std::string_view line) noexcept
{
    return line.empty() || line.front() == '#';
}

std::vector<std::string_view> split_candidate_row(std::string_view line)
{
    std::vector<std::string_view> fields;
    fields.reserve(10U);
    std::size_t begin = 0U;
    while (true) {
        const std::size_t end = line.find('\t', begin);
        fields.push_back(line.substr(
            begin,
            end == std::string_view::npos
                ? line.size() - begin
                : end - begin));
        if (end == std::string_view::npos) {break;}
        begin = end + 1U;
    }
    if (fields.size() != 10U) {
        throw RrbsCatalogError(
            "RRBS candidate row must contain exactly ten tab-separated fields");
    }
    return fields;
}

std::uint32_t parse_candidate_u32(
    std::string_view text,
    const char *field)
{
    if (text.empty()) {
        throw RrbsCatalogError(std::string("RRBS candidate ") + field
                               + " is empty");
    }
    std::uint32_t value = 0U;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value, 10);
    if (parsed.ec != std::errc{}
        || parsed.ptr != text.data() + text.size()) {
        throw RrbsCatalogError(std::string("RRBS candidate ") + field
                               + " is not a uint32 decimal");
    }
    return value;
}

std::optional<double> parse_candidate_score(std::string_view text)
{
    if (text == ".") {return std::nullopt;}
    if (text.empty()) {
        throw RrbsCatalogError("RRBS candidate score is empty");
    }
    double value = 0.0;
    const auto parsed = std::from_chars(
        text.data(), text.data() + text.size(), value,
        std::chars_format::general);
    if (parsed.ec != std::errc{}
        || parsed.ptr != text.data() + text.size()
        || !std::isfinite(value) || value < 0.0) {
        throw RrbsCatalogError(
            "RRBS candidate score must be . or a finite non-negative number");
    }
    return value;
}

void validate_candidate_id(std::string_view id)
{
    if (id.empty()) {
        throw RrbsCatalogError("RRBS candidate ID is empty");
    }
    for (const char value : id) {
        const auto byte = static_cast<unsigned char>(value);
        if (byte == 0U || byte < 0x21U || byte == 0x7fU) {
            throw RrbsCatalogError(
                "RRBS candidate ID contains whitespace or a control byte");
        }
    }
}

} // namespace

CandidateBed::CandidateBed(
    const std::string &path,
    const std::vector<reference::ContigMetadata> &reference_catalog)
    : rows_by_contig_(reference_catalog.size())
{
    try {
        if (path.empty()) {
            throw RrbsCatalogError("RRBS candidate BED path is empty");
        }
        std::unordered_map<std::string, std::uint32_t> contig_indices;
        contig_indices.reserve(reference_catalog.size());
        for (std::size_t index = 0U; index < reference_catalog.size(); ++index) {
            if (index > std::numeric_limits<std::uint32_t>::max()) {
                throw RrbsCatalogError(
                    "reference contig index exceeds uint32");
            }
            if (!contig_indices.emplace(
                    reference_catalog[index].name,
                    static_cast<std::uint32_t>(index)).second) {
                throw RrbsCatalogError(
                    "reference catalog contains duplicate contig names");
            }
        }
        std::vector<std::unordered_set<std::string>> ids(
            reference_catalog.size());
        std::ifstream input(path, std::ios::binary);
        if (!input) {
            throw RrbsCatalogError(
                "cannot open RRBS candidate BED: " + path);
        }
        std::string line;
        std::uint64_t line_number = 0U;
        while (std::getline(input, line)) {
            if (line_number == std::numeric_limits<std::uint64_t>::max()) {
                throw RrbsCatalogError(
                    "RRBS candidate BED line count exceeds uint64");
            }
            ++line_number;
            if (!line.empty() && line.back() == '\r') {line.pop_back();}
            if (line.find('\r') != std::string::npos) {
                throw RrbsCatalogError(
                    "RRBS candidate BED contains a bare carriage return");
            }
            if (line.size() > text::maximum_line_bytes) {
                throw RrbsCatalogError(
                    "RRBS candidate BED line exceeds the 1 MiB limit");
            }
            if (candidate_metadata_line(line)) {continue;}
            try {
                const auto fields = split_candidate_row(line);
                const auto found = contig_indices.find(std::string(fields[0]));
                if (found == contig_indices.end()) {
                    throw RrbsCatalogError(
                        "RRBS candidate names an unknown reference contig");
                }
                const std::uint32_t contig_index = found->second;
                const std::uint32_t start = parse_candidate_u32(
                    fields[1], "start");
                const std::uint32_t end = parse_candidate_u32(
                    fields[2], "end");
                if (start > end
                    || end > reference_catalog[contig_index].length) {
                    throw RrbsCatalogError(
                        "RRBS candidate reference envelope is outside its contig");
                }
                validate_candidate_id(fields[3]);
                const std::string id(fields[3]);
                if (!ids[contig_index].insert(id).second) {
                    throw RrbsCatalogError(
                        "RRBS candidate BED contains a duplicate candidate ID");
                }
                const std::optional<double> score =
                    parse_candidate_score(fields[4]);
                if (fields[5] != ".") {
                    throw RrbsCatalogError(
                        "RRBS candidate strand must be .");
                }
                const std::uint32_t mask_value = parse_candidate_u32(
                    fields[6], "haplotype_mask");
                if (mask_value > std::numeric_limits<std::uint8_t>::max()
                    || !model::is_haplotype_mask(
                        static_cast<std::uint8_t>(mask_value))) {
                    throw RrbsCatalogError(
                        "RRBS candidate haplotype_mask must be 1, 2, or 3");
                }
                const std::uint32_t template_length = parse_candidate_u32(
                    fields[7], "template_length");
                const std::uint32_t gc_count = parse_candidate_u32(
                    fields[8], "gc_count");
                const std::uint32_t restriction_count = parse_candidate_u32(
                    fields[9], "restriction_site_count");
                if (template_length == 0U || gc_count > template_length) {
                    throw RrbsCatalogError(
                        "RRBS candidate length or GC count is invalid");
                }
                if (rows_by_contig_[contig_index].size()
                    == std::numeric_limits<std::uint32_t>::max()) {
                    throw RrbsCatalogError(
                        "RRBS candidate count for one contig exceeds uint32");
                }
                if (row_count_ == std::numeric_limits<std::uint64_t>::max()) {
                    throw RrbsCatalogError(
                        "RRBS candidate BED row count exceeds uint64");
                }
                rows_by_contig_[contig_index].push_back(CandidateBedRow{
                    start,
                    end,
                    std::move(id),
                    score,
                    static_cast<model::HaplotypeMask>(mask_value),
                    template_length,
                    gc_count,
                    restriction_count,
                });
                ++row_count_;
            } catch (const RrbsCatalogError &error) {
                throw RrbsCatalogError(
                    "RRBS candidate BED line "
                    + std::to_string(line_number) + ": " + error.what());
            }
        }
        if (input.bad()) {
            throw RrbsCatalogError(
                "failed while reading RRBS candidate BED: " + path);
        }
        if (row_count_ == 0U) {
            throw RrbsCatalogError(
                "RRBS candidate BED contains no candidate rows");
        }
    } catch (const RrbsCatalogError &) {
        throw;
    } catch (const std::exception &error) {
        throw RrbsCatalogError(error.what());
    }
}

std::vector<double> CandidateBed::match_scores(
    std::uint32_t contig_index,
    std::string_view contig_name,
    const std::vector<Candidate> &candidates,
    bool require_scores) const
{
    if (contig_index >= rows_by_contig_.size()) {
        throw RrbsCatalogError(
            "RRBS candidate BED contig index is out of range");
    }
    const auto &rows = rows_by_contig_[contig_index];
    if (rows.size() != candidates.size()) {
        throw RrbsCatalogError(
            "RRBS candidate BED row count disagrees with the regenerated "
            "candidate set for "
            + std::string(contig_name));
    }
    std::unordered_map<std::string_view, const CandidateBedRow *> by_id;
    by_id.reserve(rows.size());
    for (const CandidateBedRow &row : rows) {
        if (!by_id.emplace(row.candidate_id, &row).second) {
            throw RrbsCatalogError(
                "RRBS candidate BED contains a duplicate candidate ID");
        }
    }
    const std::vector<std::string> ids = candidate_ids(
        contig_name, candidates);
    std::vector<double> scores;
    scores.reserve(candidates.size());
    for (std::size_t index = 0U; index < candidates.size(); ++index) {
        const auto found = by_id.find(ids[index]);
        if (found == by_id.end()) {
            throw RrbsCatalogError(
                "RRBS candidate BED is missing regenerated candidate "
                + ids[index]);
        }
        const CandidateBedRow &row = *found->second;
        const Candidate &candidate = candidates[index];
        if (row.reference_start != candidate.reference_start
            || row.reference_end != candidate.reference_end
            || row.haplotype_mask != candidate.haplotype_mask
            || row.template_length != candidate.template_length
            || row.gc_count != candidate.gc_count
            || row.restriction_site_count
                != candidate.restriction_site_count) {
            throw RrbsCatalogError(
                "RRBS candidate BED changed a fixed field for " + ids[index]);
        }
        if (require_scores && !row.score) {
            throw RrbsCatalogError(
                "RRBS profile candidate is missing a score: " + ids[index]);
        }
        scores.push_back(row.score.value_or(0.0));
    }
    return scores;
}

ProfileSampler::ProfileSampler(
    const std::vector<Candidate> &candidates,
    const std::vector<double> &scores)
{
    if (scores.size() != candidates.size()) {
        throw RrbsCatalogError(
            "RRBS profile score count disagrees with its candidate set");
    }
    cumulative_weights_.reserve(candidates.size());
    for (std::size_t index = 0U; index < candidates.size(); ++index) {
        const double score = scores[index];
        if (!std::isfinite(score) || score < 0.0) {
            throw RrbsCatalogError(
                "RRBS profile score must be finite and non-negative");
        }
        const long double weight = static_cast<long double>(score)
            * physical_copy_count(candidates[index].haplotype_mask);
        if (!std::isfinite(weight)) {
            throw RrbsCatalogError(
                "RRBS profile candidate weight is outside the floating domain");
        }
        const long double next = total_weight_ + weight;
        if (!std::isfinite(next)
            || (weight > 0.0L && next == total_weight_)) {
            throw RrbsCatalogError(
                "RRBS profile cumulative weight lost numeric precision");
        }
        total_weight_ = next;
        cumulative_weights_.push_back(total_weight_);
    }
    allocation_weight_ = static_cast<double>(total_weight_);
    if (!std::isfinite(allocation_weight_)) {
        throw RrbsCatalogError(
            "RRBS profile allocation weight exceeds binary64");
    }
}

std::vector<std::uint32_t> ProfileSampler::sample_indices(
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count) const
{
    if (output_count == 0U) {return {};}
    if (total_weight_ <= 0.0L || cumulative_weights_.empty()) {
        throw RrbsCatalogError(
            "cannot sample an RRBS contig with zero profile weight");
    }
    if (static_cast<std::uint64_t>(output_count) - 1U
        > std::numeric_limits<std::uint64_t>::max()
            - first_candidate_ordinal) {
        throw RrbsCatalogError("candidate ordinal range exceeds uint64");
    }
    const std::uint64_t key = rng::derive_key(
        master_seed, rng::Stage::fragment, contig_index);
    std::vector<std::uint32_t> result;
    result.reserve(output_count);
    for (std::uint32_t index = 0U; index < output_count; ++index) {
        const std::uint64_t ordinal = first_candidate_ordinal + index;
        long double draw = static_cast<long double>(rng::uniform01(
            key, ordinal, UINT64_C(1))) * total_weight_;
        if (draw >= total_weight_) {
            draw = std::nextafter(total_weight_, 0.0L);
        }
        const auto found = std::upper_bound(
            cumulative_weights_.begin(), cumulative_weights_.end(), draw);
        if (found == cumulative_weights_.end()) {
            throw RrbsCatalogError(
                "RRBS profile selection exceeded its cumulative weight");
        }
        result.push_back(static_cast<std::uint32_t>(
            found - cumulative_weights_.begin()));
    }
    return result;
}

CandidateCatalog::CandidateCatalog(
    const model::Bases &contig_bases,
    const std::vector<CutSite> &cut_sites,
    std::uint32_t minimum_insert_length,
    std::uint32_t maximum_insert_length,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction)
{
    candidates_ = build_candidates(
        contig_bases,
        cut_sites,
        minimum_insert_length,
        maximum_insert_length,
        read_length,
        paired_end,
        max_ambiguous_fraction,
        model::HaplotypeMask::both,
        nullptr);
}

std::uint32_t CandidateCatalog::candidate_count() const noexcept
{
    return static_cast<std::uint32_t>(candidates_.size());
}

std::uint32_t CandidateCatalog::allocation_weight() const
{
    return checked_allocation_weight(candidates_);
}

const Candidate &CandidateCatalog::candidate(std::uint32_t index) const
{
    if (index >= candidates_.size()) {
        throw RrbsCatalogError("RRBS candidate index is out of range");
    }
    return candidates_[index];
}

std::vector<std::uint32_t> CandidateCatalog::sample_indices(
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count) const
{
    return rrbs::sample_indices(
        candidates_,
        contig_index,
        master_seed,
        first_candidate_ordinal,
        output_count);
}

DiploidCandidateCatalog::DiploidCandidateCatalog(
    const reference::Contig &contig,
    const variant::ContigVariants &variants,
    const std::vector<CutSite> &cut_sites,
    std::uint32_t minimum_insert_length,
    std::uint32_t maximum_insert_length,
    std::uint32_t read_length,
    bool paired_end,
    double max_ambiguous_fraction)
{
    for (std::uint8_t haplotype_index = 0U;
         haplotype_index < 2U;
         ++haplotype_index) {
        const haplotype::HaplotypeLayout layout(
            contig, variants, haplotype_index, true);
        const auto mask = haplotype_index == 0U
            ? model::HaplotypeMask::haplotype_1
            : model::HaplotypeMask::haplotype_2;
        std::vector<Candidate> haplotype_candidates = build_candidates(
            layout.bases(),
            cut_sites,
            minimum_insert_length,
            maximum_insert_length,
            read_length,
            paired_end,
            max_ambiguous_fraction,
            mask,
            &layout);
        if (haplotype_candidates.size()
            > std::numeric_limits<std::uint32_t>::max() - candidates_.size()) {
            throw RrbsCatalogError(
                "diploid RRBS candidate count exceeds uint32");
        }
        candidates_.insert(
            candidates_.end(),
            std::make_move_iterator(haplotype_candidates.begin()),
            std::make_move_iterator(haplotype_candidates.end()));
    }
}

std::uint32_t DiploidCandidateCatalog::candidate_count() const noexcept
{
    return static_cast<std::uint32_t>(candidates_.size());
}

std::uint32_t DiploidCandidateCatalog::allocation_weight() const
{
    return checked_allocation_weight(candidates_);
}

const Candidate &DiploidCandidateCatalog::candidate(std::uint32_t index) const
{
    if (index >= candidates_.size()) {
        throw RrbsCatalogError("RRBS candidate index is out of range");
    }
    return candidates_[index];
}

std::vector<std::uint32_t> DiploidCandidateCatalog::sample_indices(
    std::uint32_t contig_index,
    std::uint64_t master_seed,
    std::uint64_t first_candidate_ordinal,
    std::uint32_t output_count) const
{
    return rrbs::sample_indices(
        candidates_,
        contig_index,
        master_seed,
        first_candidate_ordinal,
        output_count);
}

} // namespace htsim::rrbs
