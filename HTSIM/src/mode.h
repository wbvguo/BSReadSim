#ifndef MODE_H
#define MODE_H

#include <stdint.h>
#include <vector>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"


// for BED
int parse_bed_line(char *line, char *chr_id, frag_rec *tmp_probe);

int parse_bed_line_rrbs(char *line, char *chr_id, frag_rec *tmp_probe);

void parse_bed_chr(char *fname, char *chr_id, std::vector<frag_rec>& probe_vec);


// for GC-bias
void parse_bias_file(char *fname, std::vector<float>& eff_vec);


// for length and count calculation
void collect_len_score_chr(const kseq_t *ks, chr_rec *tmp_len, char *bed_file, std::vector<frag_rec>& probe_vec);

void cal_chr_count(const char *fn, char *chr_id, char *bed_file, uint64_t N, uint64_t chr_N, 
                    expt_param *expt_set, std::vector<frag_rec>& probe_vec, std::map<std::string, chr_rec> &chr_count);


#endif
