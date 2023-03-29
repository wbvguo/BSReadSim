#ifndef MODE_H
#define MODE_H

#include <stdint.h>
#include <vector>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"

// for VCF
void parse_vcf_chr(char *fname, char *chr_id, std::vector<snp_rec>& snp_vec);

// for BED
int parse_bed_line(char *line, char *chr_id, probe_rec *tmp_probe, probe_meta *tmp_probe_meta);

void parse_bed_chr(char *fname, char *chr_id, std::vector<probe_rec>& probe_vec);


// for GC-bias
void parse_bias_file(char *fname, std::vector<float>& eff_vec);


// for length calculation
void collect_len_score_chr(const kseq_t *ks, chr_rec *tmp_len, char *bed_file, std::vector<probe_rec>& probe_vec);

#endif
