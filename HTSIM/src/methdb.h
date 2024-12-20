#ifndef METHDB_H
#define METHDB_H

#include <vector>
#include <map>
#include <zlib.h>
#include <random>
#include <gsl/gsl_randist.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"


// create/save/load MethDB
void create_methdb(const kseq_t *ks, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec);

void update_variant(const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec, 
                    meth_param *meth_set, std::map<int, param_rec>& params_map, std::map<int, snpmeth_rec>& snpmeth_map);

void save_methdb(char *fname, char *chr_id, std::vector<meth_rec>& meth_vec);

int parse_methdb_line(char *line, char *chr_id, meth_rec *tmp_meth);

void load_methdb(char *fname, char *chr_id, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec);


// for CGmap
int parse_cgmap_line(char *line, char *chr_id, meth_rec *tmp_meth);

void pool_cgmap(std::vector<meth_rec>& meth_vec, int seed);

void fill_cgmap_chr(char *fname, char *chr_id, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec, meth_param *meth_set);


// for ASM
int parse_asm_line(char *line, char *chr_id, meth_rec *tmp_meth);

void fill_asm_chr(char *fname, char *chr_id, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec);


// fill with distribution
void parse_param(char *param_str, std::map<int, param_rec>& params_map);

float gen_beta(gsl_rng *rng, uint8_t context, std::map<int, param_rec>& params_map);

void fill_beta(std::vector<meth_rec>& meth_vec, std::map<int, param_rec>& params_map, int seed);


// for SNP meth
void collect_snpmeth(const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, std::map<int, snpmeth_rec>& snpmeth_map);

void save_snpmeth(char *fname, char *chr_id, std::map<int, snpmeth_rec>& snpmeth_map);

int parse_snpmeth_line(char *line, char *chr_id, snpmeth_rec *tmp_snpmeth);

void load_snpmeth(char *fname, char *chr_id, uint32_t *posidx_arr, std::map<int, snpmeth_rec>& snpmeth_map);


#endif
