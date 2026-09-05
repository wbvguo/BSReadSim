#ifndef RRCUT_H
#define RRCUT_H

#include <vector>
#include <map>
#include <string>
#include "struct.h"
#include "utility.h"

// RRBS candidates: the restriction fragments of each haplotype (of the reference
// when a contig has no variants) whose length is in [min_insert, max_insert] and
// whose reads have few enough N bases.

// cut sites such as "C|CGG"; N in a motif matches any of A/C/G/T
void parse_cut_site(const std::vector<std::string>& sites, std::vector<cut_rec>& cut_vec);

void gen_cut_frag(mutseq_t *hap1, mutseq_t *hap2, bool has_variants, const expt_param *expt_set,
                  std::vector<cut_rec>& cut_vec, std::vector<frag_rrbs_rec>& frag_vec);

// candidate IDs: "contig:start-end" of the envelope, with "~N" when envelopes repeat
void rrbs_ids(const char *chr_id, std::vector<frag_rrbs_rec>& frag_vec, std::vector<std::string>& ids);

// candidate BED: written by rrbs-catalog, read back with scores by --rrbs-candidates
void save_rrcut_header(out_t *out);
void save_rrcut_bed(out_t *out, const char *chr_id, std::vector<frag_rrbs_rec>& frag_vec);

typedef struct {
    int env_l = 0, env_r = 0, haplo = 0, len = 0, gc = 0, n_cuts = 0;
    double score = -1;          /* -1 for "." */
} rrbs_row;
void parse_rrcut_bed(const char *fname, const std::vector<chr_rec>& chr_vec,
                     std::vector<std::map<std::string, rrbs_row>>& rows_vec, uint64_t *rows);
// check the rows of a contig against its candidates; take their scores when use_score
void match_rrcut_bed(const char *chr_id, std::map<std::string, rrbs_row>& rows,
                     std::vector<frag_rrbs_rec>& frag_vec, bool use_score);

#endif
