#ifndef HAPLO_H
#define HAPLO_H

#include <stdint.h>
#include <vector>
#include <string>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"

void parse_vcf_chr(char *fname, char *chr_id, std::vector<snp_rec>& snp_vec);

void sim_mut_vcf(const kseq_t *ks, char * vcf_file, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, std::vector<snp_rec> &snp_vec);

void sim_mut_diref(const kseq_t *ks, bool is_hap, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, mut_params *tmp_params);

void sim_print_mutref(const char *name, const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2, int output_fmt);

#endif