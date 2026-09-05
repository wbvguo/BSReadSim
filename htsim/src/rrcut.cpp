#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <math.h>
#include <string.h>
#include <htslib/hts.h>
#include "struct.h"
#include "utility.h"
#include "rrcut.h"


void parse_cut_site(const std::vector<std::string>& sites, std::vector<cut_rec>& cut_vec)
{
    cut_vec.clear();
    for (size_t i = 0; i < sites.size(); ++i) {
        cut_rec tmp_site;
        for (size_t k = 0; k < sites[i].size(); ++k) {
            if (sites[i][k] == '|') {tmp_site.idx = (int)tmp_site.seq.size(); continue;}
            tmp_site.seq.push_back(nst_nt4_table[(int)sites[i][k]]);
        }
        cut_vec.push_back(tmp_site);
    }
}


// a cut of a haplotype: the boundary before haplotype base off
typedef struct {
    int off;            /* haplotype offset of the base after the cut */
    int recog;          /* motif recognitions that cut here */
    int pos, skip;      /* that base: base skip of reference position pos */
    int env_l;          /* first reference base at or after the cut */
    int prev_pos;       /* reference position of the base before the cut (the anchor if inserted) */
    int gc;             /* C/G bases before the cut */
} cutpos_rec;

// the haplotype as one sequence
static void hap_flat(mutseq_t *hap, std::vector<uint8_t>& seq)
{
    seq.clear();
    seq.reserve(hap->l);
    for (int i = 0; i < hap->l; ++i) {
        int c = hap->s[i];
        if ((c & mutmsk) == DELETE) continue;
        seq.push_back(c & 0xf);
        for (int j = 0; j < mut_ins(c); ++j) seq.push_back(mut_ins_base(c, j));
    }
}

static void gen_cut_pos(std::vector<uint8_t>& seq, std::vector<cut_rec>& cut_vec, std::vector<cutpos_rec>& cutpos_vec)
{
    std::vector<int> offs;
    for (size_t s = 0; s < cut_vec.size(); ++s) {
        const cut_rec &site = cut_vec[s];
        int len = (int)site.seq.size();
        for (size_t start = 0; start + len <= seq.size(); ++start) {
            int k = 0;
            for (; k < len; ++k) {
                uint8_t base = seq[start + k];
                if (site.seq[k] == 4 ? base == 4 : base != site.seq[k]) break;
            }
            if (k == len) offs.push_back((int)start + site.idx);
        }
    }
    std::sort(offs.begin(), offs.end());
    cutpos_vec.clear();
    for (size_t i = 0; i < offs.size(); ++i) {
        if (cutpos_vec.empty() || cutpos_vec.back().off != offs[i]) {
            cutpos_rec tmp_cut = {};
            tmp_cut.off = offs[i];
            cutpos_vec.push_back(tmp_cut);
        }
        ++cutpos_vec.back().recog;
    }
    // C/G bases before each cut
    int gc = 0;
    size_t scanned = 0;
    for (size_t k = 0; k < cutpos_vec.size(); ++k) {
        for (; scanned < (size_t)cutpos_vec[k].off; ++scanned) gc += cg_table[seq[scanned]];
        cutpos_vec[k].gc = gc;
    }
}

// where on the reference each cut is
static void map_cut_pos(mutseq_t *hap, std::vector<cutpos_rec>& cutpos_vec)
{
    size_t next = 0, wait = 0;  // next cut to place; first cut waiting for its env_l
    int off = 0, last_pos = -1;
    for (int i = 0; i < hap->l && wait < cutpos_vec.size(); ++i) {
        int c = hap->s[i], mut_type = c & mutmsk;
        if (mut_type == DELETE) continue;
        int num_ins = mut_ins(c);
        for (int j = -1; j < num_ins; ++j, ++off) {
            for (; next < cutpos_vec.size() && cutpos_vec[next].off == off; ++next) {
                cutpos_vec[next].pos = i;
                cutpos_vec[next].skip = j + 1;
                cutpos_vec[next].prev_pos = last_pos;
            }
            if (j < 0) for (; wait < next; ++wait) cutpos_vec[wait].env_l = i;
            last_pos = i;
        }
    }
    for (; next < cutpos_vec.size(); ++next) {
        cutpos_vec[next].pos = hap->l;
        cutpos_vec[next].skip = 0;
        cutpos_vec[next].prev_pos = last_pos;
    }
    for (; wait < cutpos_vec.size(); ++wait) cutpos_vec[wait].env_l = hap->l;
}

// restriction fragments between any two cuts of one haplotype, size-selected, with sequenceable reads
static void add_cut_frag(mutseq_t *hap, int haplo, const expt_param *expt_set, std::vector<cut_rec>& cut_vec,
                         std::vector<frag_rrbs_rec>& frag_vec)
{
    std::vector<uint8_t> seq;
    std::vector<cutpos_rec> cutpos_vec;
    hap_flat(hap, seq);
    gen_cut_pos(seq, cut_vec, cutpos_vec);
    map_cut_pos(hap, cutpos_vec);

    int size = (int)seq.size(), n = (int)cutpos_vec.size();
    // per cut: whether the read it starts, and the read it ends, has few enough N bases (as
    // either mate if paired: the strand is drawn later); recognitions up to it
    int mates = expt_set->paired_end && expt_set->read_length[1] != expt_set->read_length[0] ? 2 : 1;
    std::vector<char> ok_after(n, 1), ok_before(n, 1);
    std::vector<int> recog(n + 1, 0);
    for (int k = 0; k < n; ++k) {
        int off = cutpos_vec[k].off;
        for (int mate = 0; mate < mates; ++mate) {
            int read = expt_set->read_length[mate], limit = max_n(expt_set, mate);
            ok_after[k] &= off + read <= size && count_n(seq, off, off + read) <= limit;
            ok_before[k] &= off >= read && count_n(seq, off - read, off) <= limit;
        }
        recog[k + 1] = recog[k] + cutpos_vec[k].recog;
    }
    for (int left = 0; left < n; ++left) {
        if (!ok_after[left]) continue;
        for (int right = left + 1; right < n; ++right) {
            int len = cutpos_vec[right].off - cutpos_vec[left].off;
            if (len < expt_set->min_insert) continue;
            if (len > expt_set->max_insert) break;
            if (expt_set->paired_end && !ok_before[right]) continue;
            frag_rrbs_rec tmp_frag;
            tmp_frag.haplo = haplo;
            tmp_frag.pos_l = cutpos_vec[left].pos;
            tmp_frag.skip  = cutpos_vec[left].skip;
            tmp_frag.env_l = cutpos_vec[left].env_l;
            tmp_frag.env_r = cutpos_vec[right].prev_pos + 1;
            if (tmp_frag.env_l >= tmp_frag.env_r) continue;    // inserted bases only
            tmp_frag.len   = len;
            tmp_frag.gc    = cutpos_vec[right].gc - cutpos_vec[left].gc;
            tmp_frag.n_cuts= recog[right + 1] - recog[left];
            frag_vec.push_back(tmp_frag);
        }
    }
}

void gen_cut_frag(mutseq_t *hap1, mutseq_t *hap2, bool has_variants, const expt_param *expt_set,
                  std::vector<cut_rec>& cut_vec, std::vector<frag_rrbs_rec>& frag_vec)
{
    frag_vec.clear();
    if (!has_variants) {
        add_cut_frag(hap1, 3, expt_set, cut_vec, frag_vec);
        return;
    }
    add_cut_frag(hap1, 1, expt_set, cut_vec, frag_vec);
    add_cut_frag(hap2, 2, expt_set, cut_vec, frag_vec);
}


void rrbs_ids(const char *chr_id, std::vector<frag_rrbs_rec>& frag_vec, std::vector<std::string>& ids)
{
    std::map<uint64_t, int> total, seen;
    for (size_t i = 0; i < frag_vec.size(); ++i) ++total[(uint64_t)frag_vec[i].env_l << 32 | frag_vec[i].env_r];
    ids.clear();
    for (size_t i = 0; i < frag_vec.size(); ++i) {
        uint64_t key = (uint64_t)frag_vec[i].env_l << 32 | frag_vec[i].env_r;
        std::string id = std::string(chr_id) + ":" + std::to_string(frag_vec[i].env_l) + "-" + std::to_string(frag_vec[i].env_r);
        if (total[key] > 1) id += "~" + std::to_string(seen[key]++);
        ids.push_back(id);
    }
}

void save_rrcut_header(out_t *out)
{
    out_printf(out, "#chrom\tstart\tend\tcandidate_id\tscore\tstrand\thaplotype_mask\ttemplate_length\tgc_count\trestriction_site_count\n");
}

void save_rrcut_bed(out_t *out, const char *chr_id, std::vector<frag_rrbs_rec>& frag_vec)
{
    std::vector<std::string> ids;
    rrbs_ids(chr_id, frag_vec, ids);
    for (size_t i = 0; i < frag_vec.size(); ++i) {
        frag_rrbs_rec &f = frag_vec[i];
        out_printf(out, "%s\t%d\t%d\t%s\t1\t.\t%d\t%d\t%d\t%d\n", chr_id, f.env_l, f.env_r, ids[i].c_str(), f.haplo, f.len, f.gc, f.n_cuts);
    }
}


void parse_rrcut_bed(const char *fname, const std::vector<chr_rec>& chr_vec,
                     std::vector<std::map<std::string, rrbs_row>>& rows_vec, uint64_t *rows)
{
    std::map<std::string, int> chr_idx = chr_index(chr_vec);
    rows_vec.assign(chr_vec.size(), std::map<std::string, rrbs_row>());
    text_in *in = in_open(fname, "RRBS candidate BED");
    uint64_t count = 0;
    while (in_next(in)) {
        if (bed_header(in->line.s)) continue;
        if (in_split(in) != 10) IN_DIE(in, "candidate rows must have exactly ten fields");
        std::vector<char*> &f = in->f;
        int chr = in_chr(in, chr_idx, f[0]);
        rrbs_row row;
        row.env_l = (int)in_u32(in, f[1], "start");
        row.env_r = (int)in_u32(in, f[2], "end");
        if (row.env_l > row.env_r || (uint32_t)row.env_r > chr_vec[chr].chr_len) IN_DIE(in, "candidate envelope is outside its contig");
        if (!*f[3]) IN_DIE(in, "candidate ID is empty");
        if (strcmp(f[4], ".")) {
            row.score = in_double(in, f[4], "score");
            if (row.score < 0) IN_DIE(in, "score must be . or non-negative");
        }
        if (strcmp(f[5], ".")) IN_DIE(in, "candidate strand must be .");
        row.haplo  = (int)in_u32(in, f[6], "haplotype_mask");
        row.len    = (int)in_u32(in, f[7], "template_length");
        row.gc     = (int)in_u32(in, f[8], "gc_count");
        row.n_cuts = (int)in_u32(in, f[9], "restriction_site_count");
        if (row.haplo < 1 || row.haplo > 3) IN_DIE(in, "haplotype_mask must be 1, 2, or 3");
        if (!rows_vec[chr].insert(std::make_pair(std::string(f[3]), row)).second) IN_DIE(in, "duplicate candidate ID %s", f[3]);
        ++count;
    }
    in_close(in);
    if (count == 0) die("RRBS candidate BED contains no candidate rows");
    if (rows) *rows = count;
}

void match_rrcut_bed(const char *chr_id, std::map<std::string, rrbs_row>& rows,
                     std::vector<frag_rrbs_rec>& frag_vec, bool use_score)
{
    if (rows.size() != frag_vec.size()) {
        die("RRBS candidate BED has %zu rows for %s but the reference yields %zu candidates", rows.size(), chr_id, frag_vec.size());
    }
    std::vector<std::string> ids;
    rrbs_ids(chr_id, frag_vec, ids);
    for (size_t i = 0; i < frag_vec.size(); ++i) {
        std::map<std::string, rrbs_row>::iterator found = rows.find(ids[i]);
        if (found == rows.end()) die("RRBS candidate BED lacks candidate %s", ids[i].c_str());
        rrbs_row &row = found->second;
        frag_rrbs_rec &f = frag_vec[i];
        if (row.env_l != f.env_l || row.env_r != f.env_r || row.haplo != f.haplo || row.len != f.len || row.gc != f.gc || row.n_cuts != f.n_cuts) {
            die("RRBS candidate BED changed a fixed field of %s", ids[i].c_str());
        }
        if (use_score) {
            if (row.score < 0) die("RRBS candidate %s has no score", ids[i].c_str());
            f.score = row.score;
        }
    }
}
