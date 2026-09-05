#include <stdint.h>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <htslib/hts.h>
#include "struct.h"
#include "utility.h"
#include "mode.h"
#include "haplo.h"


static const int MAX_ATTEMPTS = 100000;    /* tries to place one fragment */


// for BED
void parse_bed(const char *fname, const std::vector<chr_rec>& chr_vec, std::vector<std::vector<probe_rec>>& probe_vec, uint64_t *rows)
{
    std::map<std::string, int> chr_idx = chr_index(chr_vec);
    std::vector<std::vector<std::pair<probe_rec, std::string>>> named(chr_vec.size());
    text_in *in = in_open(fname, "target BED");
    uint64_t count = 0;
    while (in_next(in)) {
        if (bed_header(in->line.s)) continue;
        if (in_split(in) != 6) IN_DIE(in, "target rows must have exactly six BED fields");
        std::vector<char*> &f = in->f;
        int chr = in_chr(in, chr_idx, f[0]);
        probe_rec tmp_probe;
        uint32_t start = in_u32(in, f[1], "start"), end = in_u32(in, f[2], "end");
        if (start >= end || end > chr_vec[chr].chr_len) IN_DIE(in, "target must satisfy 0 <= start < end <= contig length");
        if (!*f[3]) IN_DIE(in, "target name is empty");
        tmp_probe.pos_l = (int)start;
        tmp_probe.pos_r = (int)end;
        tmp_probe.score = in_double(in, f[4], "score");
        if (tmp_probe.score < 0) IN_DIE(in, "score must be non-negative");
        tmp_probe.strand = in_strand(in, f[5]);
        named[chr].push_back(std::make_pair(tmp_probe, std::string(f[3])));
        ++count;
    }
    in_close(in);
    if (count == 0) die("target BED contains no target rows");

    // sort by start, end, strand, name; one target per interval and strand
    probe_vec.assign(chr_vec.size(), std::vector<probe_rec>());
    for (size_t c = 0; c < named.size(); ++c) {
        std::sort(named[c].begin(), named[c].end(), [](const std::pair<probe_rec, std::string> &a, const std::pair<probe_rec, std::string> &b) {
            const probe_rec &x = a.first, &y = b.first;
            if (x.pos_l != y.pos_l) return x.pos_l < y.pos_l;
            if (x.pos_r != y.pos_r) return x.pos_r < y.pos_r;
            if (x.strand != y.strand) return x.strand < y.strand;
            return a.second < b.second;
        });
        for (size_t i = 0; i < named[c].size(); ++i) {
            const probe_rec &probe = named[c][i].first;
            const std::vector<probe_rec> &probes = probe_vec[c];
            if (!probes.empty() && probes.back().pos_l == probe.pos_l && probes.back().pos_r == probe.pos_r && probes.back().strand == probe.strand) {
                die("target BED repeats the target %s", named[c][i].second.c_str());
            }
            probe_vec[c].push_back(probe);
        }
    }
    if (rows) *rows = count;
}


// for GC-bias
void parse_bias_file(const char *fname, std::vector<double>& eff_vec)
{
    text_in *in = in_open(fname, "GC profile");
    double sum = 0;
    eff_vec.clear();
    while (in_next(in)) {
        eff_vec.push_back(in_prob(in, in->line.s, "GC profile probability"));
        sum += eff_vec.back();
    }
    in_close(in);
    if (eff_vec.size() < 2) die("GC profile needs at least two bins");
    if (fabs(sum - 1) > 1e-6) die("GC profile probabilities must sum to 1");
}


// for length and score calculation

// round(gc / len * (bins - 1))
static int gc_bin(int gc, int len, int bins) {return (int)((2 * (uint64_t)gc * (bins - 1) + len) / (2 * (uint64_t)len));}

// the merged length of sorted, disjoint intervals
static uint64_t union_length(const std::vector<std::pair<int, int>>& merged)
{
    uint64_t total = 0;
    for (size_t i = 0; i < merged.size(); ++i) total += merged[i].second - merged[i].first;
    return total;
}

// where the fragments of a target fall: within the longest insert plus 6 centre SDs of
// its centre. A fragment beyond (a few in a billion) is drawn again, so the methylome
// of a targeted contig is needed only there.
static std::pair<int, int> target_reach(const probe_rec &probe, int l, const expt_param *expt_set)
{
    int centre = probe.pos_l + (probe.pos_r - probe.pos_l) / 2;
    double margin = expt_set->max_insert + 6 * expt_set->sd_center;
    return std::make_pair((int)std::max(0.0, floor(centre - margin)), (int)std::min((double)l, ceil(centre + margin)));
}

// WGBS/WGS: reference starts of a mean-length fragment whose reads (read 1 at the left
// end, and read 2 at the right end if paired) have few enough N bases
void wg_domain(const ref_seq *ref, const expt_param *expt_set, int bins, chr_rec *tmp_chr)
{
    int l = (int)ref->len, len = expt_set->mean_insert;
    int read1 = expt_set->read_length[0], read2 = expt_set->read_length[1];
    int limit1 = max_n(expt_set, 0), limit2 = max_n(expt_set, 1);
    tmp_chr->gc_starts.assign(bins, 0);
    uint64_t starts = 0;
    if (len <= l && read1 <= len && read2 <= len) {
        const char *s = ref->seq;
        #define IS_N(i) (nst_nt4_table[(int)s[i]] > 3)
        #define IS_GC(i) (cg_table[nst_nt4_table[(int)s[i]]])
        int left = 0, right = 0, gc = 0;
        for (int i = 0; i < read1; ++i) left += IS_N(i);
        for (int i = 0; i < read2; ++i) right += IS_N(len - read2 + i);
        for (int i = 0; i < len; ++i) gc += IS_GC(i);
        for (int start = 0;; ++start) {
            if (left <= limit1 && (!expt_set->paired_end || right <= limit2)) {
                ++starts;
                if (bins) ++tmp_chr->gc_starts[gc_bin(gc, len, bins)];
            }
            int end = start + len;
            if (end == l) break;
            left += IS_N(start + read1) - IS_N(start);
            right += IS_N(end) - IS_N(end - read2);
            gc += IS_GC(end) - IS_GC(start);
        }
        #undef IS_N
        #undef IS_GC
    }
    tmp_chr->score = (double)starts;
    tmp_chr->eff_len = starts ? l : 0;
}

// RRBS: copies (times score) of each candidate; a candidate of both haplotypes counts twice
void rrbs_domain(const std::vector<frag_rrbs_rec>& cands, frag_domain *dom, chr_rec *tmp_chr)
{
    double total = 0;
    std::vector<std::pair<int, int>> envelopes;
    dom->cands.clear();
    dom->cum.clear();
    for (size_t i = 0; i < cands.size(); ++i) {
        const frag_rrbs_rec &cand = cands[i];
        frag_rec frag;
        frag.pos_l = cand.pos_l; frag.skip = cand.skip; frag.len = cand.len;
        frag.haplo = cand.haplo == 3 ? -1 : cand.haplo - 1;
        dom->cands.push_back(frag);
        total += cand.score * (cand.haplo == 3 ? 2 : 1);
        dom->cum.push_back(total);
        envelopes.push_back(std::make_pair(cand.env_l, cand.env_r));
    }
    tmp_chr->score = total;
    dom->reach = merge_intervals(envelopes);
    tmp_chr->eff_len = union_length(dom->reach);
}

// targeted: a target is usable if the reads of a mean-length reference fragment
// centred on it (moved inside the contig if needed) can have few enough N bases (read 1
// at the left end and read 2 at the right; a single read at either end).
void probe_domain(const ref_seq *ref, const expt_param *expt_set, std::vector<probe_rec>& probe_vec, frag_domain *dom, chr_rec *tmp_chr)
{
    int l = (int)ref->len, len = expt_set->mean_insert, limit1 = max_n(expt_set, 0), limit2 = max_n(expt_set, 1);
    int read1 = expt_set->read_length[0], read2 = expt_set->paired_end ? expt_set->read_length[1] : read1;
    double total = 0;
    std::vector<std::pair<int, int>> intervals, reach;
    dom->usable.clear();
    dom->cum.clear();
    for (size_t i = 0; i < probe_vec.size(); ++i) {
        probe_rec &probe = probe_vec[i];
        intervals.push_back(std::make_pair(probe.pos_l, probe.pos_r));
        double weight = expt_set->sampling == SCORE ? probe.score : 1.0;
        if (weight <= 0 || len > l) continue;
        int centre = probe.pos_l + (probe.pos_r - probe.pos_l) / 2;
        int begin = std::min(centre - std::min(centre, len / 2), l - len);
        int n_left = count_n(ref, begin, begin + read1), n_right = count_n(ref, begin + len - read2, begin + len);
        bool ok = expt_set->paired_end ? n_left <= limit1 && n_right <= limit2 : n_left <= limit1 || n_right <= limit1;
        if (!ok) continue;
        total += weight;
        dom->usable.push_back((int)i);
        dom->cum.push_back(total);
        reach.push_back(target_reach(probe, l, expt_set));
    }
    dom->reach = merge_intervals(reach);
    tmp_chr->score = total;
    tmp_chr->eff_len = union_length(merge_intervals(intervals));
}

void calibrate_gc(std::vector<double>& eff_vec, std::vector<chr_rec>& chr_vec, bool drop_unreachable)
{
    // accepting a start of bin i with probability proportional to p_i / starts_i makes accepted fragments follow p
    size_t bins = eff_vec.size();
    std::vector<double> starts(bins, 0), acceptance(bins, 0);
    for (size_t c = 0; c < chr_vec.size(); ++c) {
        for (size_t i = 0; i < bins; ++i) starts[i] += (double)chr_vec[c].gc_starts[i];
    }
    double highest = 0;
    for (size_t i = 0; i < bins; ++i) {
        if (eff_vec[i] == 0) continue;
        if (starts[i] == 0) {
            if (drop_unreachable) continue;
            die("GC profile bin %zu has probability %g but no fragment has that GC content", i, eff_vec[i]);
        }
        acceptance[i] = eff_vec[i] / starts[i];
        highest = std::max(highest, acceptance[i]);
    }
    if (highest == 0) die("no fragment matches the GC profile");
    for (size_t i = 0; i < bins; ++i) acceptance[i] /= highest;
    for (size_t c = 0; c < chr_vec.size(); ++c) {
        chr_vec[c].score = 0;
        for (size_t i = 0; i < bins; ++i) chr_vec[c].score += (double)chr_vec[c].gc_starts[i] * acceptance[i];
    }
    eff_vec.swap(acceptance);
}


// for count calculation
void cal_chr_count(const expt_param *expt_set, std::vector<chr_rec>& chr_vec)
{
    double sum = 0;
    uint64_t bases = 0;
    for (size_t i = 0; i < chr_vec.size(); ++i) {
        sum += chr_vec[i].score;
        if (chr_vec[i].score > 0) bases += chr_vec[i].eff_len;
    }
    uint64_t N = expt_set->fragments;
    if (N == 0) {   // fragments for a mean depth over the effective bases (read bases per fragment), rounded up
        if (bases == 0) die("depth needs a reference region that can produce fragments");
        int read_bases = expt_set->read_length[0] + (expt_set->paired_end ? expt_set->read_length[1] : 0);
        double n = ceil((double)bases * expt_set->depth / read_bases);
        if (n < 1) die("depth yields no fragments");
        N = (uint64_t)n;
    }
    if (!(sum > 0)) die("no contig can produce a sequenceable fragment");

    std::vector<std::pair<double, size_t>> remainders;
    uint64_t assigned = 0;
    for (size_t i = 0; i < chr_vec.size(); ++i) {
        long double share = (long double)N * chr_vec[i].score / sum;
        chr_vec[i].count = (uint64_t)share;
        assigned += chr_vec[i].count;
        if (chr_vec[i].score > 0) remainders.push_back(std::make_pair(-(double)(share - chr_vec[i].count), i));
    }
    std::sort(remainders.begin(), remainders.end());
    for (size_t k = 0; assigned < N; ++k, ++assigned) ++chr_vec[remainders[k % remainders.size()].second].count;
}


// for fragment generation
static int draw_insert(const expt_param *expt_set, Rng &rng)
{
    if (expt_set->sd_insert == 0 || expt_set->min_insert == expt_set->max_insert) return expt_set->mean_insert;
    double len = round(rng.norm(expt_set->mean_insert, expt_set->sd_insert));
    return (int)std::min<double>(std::max<double>(len, expt_set->min_insert), expt_set->max_insert);
}

// the strand read 1 comes from: bisulfite libraries follow the capture strand ('+', '-',
// or '.' for either), and non-directional ones also yield the complementary strands
static int8_t draw_strand(const expt_param *expt_set, Rng &rng, int capture)
{
    if (!expt_set->bisulfite) return rng.unif() < 0.5 ? OB : OT;
    bool bottom = capture == '-' || (capture == '.' && rng.unif() < 0.5);
    bool complementary = !expt_set->directional && rng.unif() < 0.5;
    if (bottom) return complementary ? CTOB : OB;
    return complementary ? CTOT : OT;
}

// every read of the fragment has few enough N bases: read 1 at the left end (the right
// one when it is reverse), read 2 at the other
static bool sequenceable(const expt_param *expt_set, frag_seq *fs, int strand)
{
    int len = (int)fs->seq.size(), left_mate = strand == OB || strand == CTOT ? 1 : 0;
    for (int mate = 0; mate < (expt_set->paired_end ? 2 : 1); ++mate) {
        int read = expt_set->read_length[mate];
        int n = mate == left_mate ? count_n(fs->seq, 0, read) : count_n(fs->seq, len - read, len);
        if (n > max_n(expt_set, mate)) return false;
    }
    return true;
}

static size_t pick(std::vector<double>& cum, Rng &rng)
{
    double x = rng.unif() * cum.back();
    size_t found = std::upper_bound(cum.begin(), cum.end(), x) - cum.begin();
    return std::min(found, cum.size() - 1);
}

bool gen_frag(int chr, uint64_t ordinal, const ref_seq *ref, mutseq_t *rseq, const expt_param *expt_set, std::vector<double>& eff_vec,
              std::vector<probe_rec>& probe_vec, frag_domain *dom, frag_rec *frag, frag_seq *fs)
{
    Rng rng(seed_of(expt_set->seed, chr, ordinal));
    int l = (int)ref->len;
    frag->skip = 0;

    if (expt_set->tech_mode == RRBS) {      // a candidate as it is
        *frag = dom->cands[pick(dom->cum, rng)];
        if (frag->haplo < 0) frag->haplo = rng.unif() < 0.5 ? 1 : 0;
        frag->strand= draw_strand(expt_set, rng, '.');
        gen_frag_seq(&rseq[frag->haplo], frag->pos_l, l, frag->skip, frag->len, fs);
        return true;
    }

    const probe_rec *probe = 0;
    if (expt_set->tech_mode != WGBS && expt_set->tech_mode != WGS) {
        probe = &probe_vec[dom->usable[pick(dom->cum, rng)]];
        frag->strand = draw_strand(expt_set, rng, probe->strand);
    }
    for (int attempt = 0; attempt < MAX_ATTEMPTS; ++attempt) {
        frag->haplo = rng.unif() < 0.5 ? 1 : 0;
        mutseq_t *hap = &rseq[frag->haplo];
        frag->len = draw_insert(expt_set, rng);
        if (probe) {        // centred on the target, off by a normal shift
            int shift = expt_set->sd_center > 0 ? (int)llround(rng.norm(0, expt_set->sd_center)) : 0;
            int centre = probe->pos_l + (probe->pos_r - probe->pos_l) / 2 + shift;
            if (centre < 0 || centre >= l) continue;
            frag->pos_l = centre - frag->len / 2;
            if (frag->pos_l < 0) continue;
            while (frag->pos_l < l && (hap->s[frag->pos_l] & mutmsk) == DELETE) ++frag->pos_l;
        } else {            // anywhere
            frag->pos_l = (int)rng.below(l);
            if ((hap->s[frag->pos_l] & mutmsk) == DELETE) continue;
            frag->strand = draw_strand(expt_set, rng, '.');
        }
        if (gen_frag_seq(hap, frag->pos_l, l, 0, frag->len, fs) < frag->len) continue;
        if (probe) {
            std::pair<int, int> span = target_reach(*probe, l, expt_set);
            if (fs->pos.front() < span.first || fs->pos.back() >= span.second) continue;
        }
        if (!sequenceable(expt_set, fs, frag->strand)) continue;
        if (!eff_vec.empty()) {     // GC sampling
            int gc = 0;
            for (int k = 0; k < frag->len; ++k) gc += cg_table[fs->seq[k]];
            if (!(rng.unif() < eff_vec[gc_bin(gc, frag->len, (int)eff_vec.size())])) continue;
        }
        return true;
    }
    return false;
}
