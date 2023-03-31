#ifndef RRCUT_H
#define RRCUT_H

#include <vector>
#include <map>
#include <algorithm>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"


void parse_cut_site(char *cut_str, std::vector<cut_rec>& cut_vec);

void gen_cut_pos(mutseq_t *hap1, mutseq_t *hap2, std::vector<cutpos_rec>& cutpos_vec, std::vector<cut_rec>& cut_vec);

void gen_cut_frag(const kseq_t *ks, expt_param *expt_set, std::vector<frag_rrbs_rec> &frag_vec, std::vector<cutpos_rec>& cutpos_vec, std::vector<cut_rec>& cut_vec);

void output_rrcut_bed(const char *fname, const char *chr_id, std::vector<frag_rrbs_rec> &frag_vec);

#endif