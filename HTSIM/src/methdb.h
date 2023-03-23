#ifndef METHDB_H
#define METHDB_H

#include <vector>
#include <random>
#include <gsl/gsl_randist.h>
#include "htsim2.h"
#include "kseq.h"
#include "vcf.h"


// create/stream MethDB
void create_methdb(const kseq_t *ks, uint32_t *pos_idx_arr, std::vector<meth_rec>& meth_vec);

void save_methdb(std::vector<meth_rec>& meth_vec, char *fname);

void parse_methdb_line(char *line, meth_rec *tmp_meth);

void load_methdb(uint32_t *pos_idx_arr, std::vector<meth_rec>& meth_vec, char *fname);


// for CGmap
int parse_cgmap_line(char *line, char *chr_id, meth_rec *tmp_meth);

void pool_cgmap(std::vector<meth_rec>& meth_vec, int seed);

void fill_cgmap_chr(char *fname, char *chr_id, uint32_t *pos_idx_arr, std::vector<meth_rec>& meth_vec, bool cgmap_pool, int seed);


// for ASM
int parse_asm_line(char *line, char *chr_id, meth_rec *tmp_meth);

void fill_asm_chr(char *fname, char *chr_id, uint32_t *pos_idx_arr, std::vector<meth_rec>& meth_vec);


// fill with distribution
void parse_param(char *param_str, std::vector<param_rec>& param_vec);

float gen_beta(gsl_rng *rng, uint8_t context, std::vector<param_rec>& param_vec);

void fill_beta(std::vector<meth_rec>& meth_vec, std::vector<param_rec>& param_vec, int seed);

#endif
