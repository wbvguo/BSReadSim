#ifndef HAPLO_H
#define HAPLO_H

#include <stdint.h>
#include <vector>
#include <string>
#include "struct.h"
#include "utility.h"

// counts reported by validate-inputs
typedef struct {
    uint64_t rows = 0, contigs = 0, retained = 0, reference_genotypes = 0;
    uint64_t mnp = 0, complex_replacement = 0, long_indel = 0;
} vcf_summary;

// read a one-sample VCF into one variant map per reference contig
void parse_vcf(const char *fname, const std::vector<chr_rec>& chr_vec, uint64_t seed_phase,
               std::vector<std::map<int, snpmeth_rec>>& vcf_vec, vcf_summary *summary);

// generate random variants of a contig
void sim_mut_diref(const ref_seq *ref, const mut_param *mut_set, int chr_idx, std::map<int, snpmeth_rec>& snpmeth_map);

// build hap1/hap2 from the variants of a contig (checking them against the reference)
void sim_mut_map(const ref_seq *ref, std::map<int, snpmeth_rec>& snpmeth_map, mutseq_t *hap1, mutseq_t *hap2);

// read haplotype bases from reference position pos_l on: those of reference interval
// [pos_l, pos_r), without the first skip of them, at most len; returns their number
int gen_frag_seq(mutseq_t *hap, int pos_l, int pos_r, int skip, int len, frag_seq *fs);
// context of every base of a fragment (0 for none; non-CpG dropped unless collect_non_cpg)
void frag_contexts(frag_seq *fs, std::vector<uint8_t>& ctx, bool collect_non_cpg);
// 7-mer index centred on base k (NO_KMER if it contains N or passes an end)
uint16_t frag_kmer(frag_seq *fs, int k);

// variants as VCF
void save_vcf_header(out_t *out);
void save_vcf_chr(out_t *out, const ref_seq *ref, std::map<int, snpmeth_rec>& snpmeth_map);

#endif
