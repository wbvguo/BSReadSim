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

void gen_cut_pos(mutseq_t *hap1, mutseq_t *hap2, std::vector<cut_rec>& cut_vec, std::vector<cutpos_rec>& cutpos_vec);

void gen_cut_frag(int chr_len, expt_param *expt_set, std::vector<cut_rec>& cut_vec, std::vector<cutpos_rec>& cutpos_vec, std::vector<frag_rrbs_rec> &frag_vec);

void save_rrcut_bed(char *fname, char *chr_id, bool to_stdout, std::vector<frag_rrbs_rec> &frag_vec);

#endif
