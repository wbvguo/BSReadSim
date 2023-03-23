#ifndef MODE_H
#define MODE_H

#include <algorithm>
#include <stdint.h>
#include <vector>
#include "htsim2.h"
#include "kseq.h"
#include "vcf.h"


// for VCF
void parse_vcf_chr(char *fname, char *chr_id, std::vector<snp_rec>& snp_vec);

// for BED
int parse_bed_line(char *line, char *chr_id, probe_rec *tmp_probe, probe_meta *tmp_probe_meta);

void parse_bed_chr(char *fname, char *chr_id, std::vector<probe_rec>& probe_vec);


// for RRBS
void parse_cut_rec(char *cut_str, std::vector<cut_rec>& cut_vec);

void gen_cut_pos(const kseq_t *ks, std::vector<cut_pos>& cutpos_vec, std::vector<cut_rec>& cut_vec);


//gen cut fragment


// for GC-bias
void parse_bias_file(char *fname, std::vector<float>& eff_vec);


// for length calculation
void cal_length_chr(const kseq_t *ks, chr_rec *tmp_len, int tech_mode, char *bed_file, int min_insert, int max_insert, 
                    std::vector<probe_rec>& probe_vec, std::vector<cut_pos>& cutpos_vec, std::vector<cut_rec>& cut_vec);

#endif
