#ifndef MODE_H
#define MODE_H

#include <stdint.h>
#include <vector>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"


// for BED
int parse_bed_line(char *line, char *chr_id, frag_rec *tmp_probe, int ncol_coll, int ncol_skip);

void parse_bed_chr(char *fname, char *chr_id, std::vector<frag_rec>& probe_vec, int tech_mode);


// for GC-bias
void parse_bias_file(char *fname, std::vector<float>& eff_vec);


// for length and count calculation
void collect_len_score_chr(const kseq_t *ks, chr_rec *tmp_len, char *bed_file, int tech_mode, std::vector<frag_rec>& probe_vec);

void cal_chr_count(const char *fn, char *chr_id, char *bed_file, uint64_t N, uint64_t chr_N, 
                    expt_param *expt_set, std::vector<frag_rec>& probe_vec, std::map<std::string, chr_rec> &chr_count);


// for fragment generation
void gen_frag_vec(std::uniform_int_distribution<int> *dis_ud, std::discrete_distribution<int> *dis_dd, uint32_t *posidx_arr, int chr_len, int chunk_size,
                  std::vector<frag_rec> &frag_vec, std::vector<frag_rec> &probe_vec, std::vector<float> &eff_vec, expt_param *expt_set);

//void check_frag_vec(std::vector<frag_rec> &frag_vec, mutseq_t *hap1, mutseq_t *hap2, expt_param *expt_set);


#endif

