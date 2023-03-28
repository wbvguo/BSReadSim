/* The MIT License
   Copyright (c) 2008 Genome Research Ltd (GRL).
                 2011 Heng Li <lh3@live.co.uk>
                 2022 Wenbin Guo <wbguo@ucla.edu>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* This program is based on WGSIM(v0.3.1-r13)[https://github.com/lh3/wgsim.git], with heavy 
 * modifications to simulate WGS/RRS/TS or WGBS/RRBS/TBS reads in BSReadSim for diploid organisms */

#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <zlib.h>
#include <random>
#include <vector>
#include <map>
#include <algorithm>
#include "kseq.h"
#include "vcf.h"

#include "mode.h"
#include "methdb.h"
#include "htsim2.h"
KSEQ_INIT(gzFile, gzread)

#define PACKAGE_VERSION "1.0.3"

extern const uint8_t nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

extern std::map<std::string, int> base_map = {{"C", 0}, {"G",1}};
extern std::map<std::string, int> context_map = {{"CG",1}, {"CHG",3}, {"CHH",7}};
extern std::map<int, int> params_map = {{1,0}, {3,1}, {7,2}, {9,0}, {11,1}, {15,2}}; //idx in param_vec

param_rec param_cg = {0.5,0.5};
param_rec param_chg= {0.5,0.5};
param_rec param_chh= {0.5,0.5};


static double ERR_RATE = 0.005;
static double MUT_RATE = 0.01;
static double INDEL_FRAC = 0.15;
static double INDEL_EXTN = 0.3;
static double MAX_N_RATIO= 0.05;


static uint8_t MATCH  = 0x00;
static uint8_t SNV    = 0x01;
static uint8_t INSR   = 0x03;
static uint8_t CONVT  = 0x05; // reserved, not used
static uint8_t SEQERR = 0x09;
const  uint8_t mut_table[16] = {
    0, 1, 0, 2, 
    0, 0, 0, 0, 
    0, 3, 3, 3, 
    3, 3, 3, 3
}; // MXIE


static uint8_t CG = 0x01; //5to3
static uint8_t CHG= 0x03;
static uint8_t CHH= 0x07;
static uint8_t GC = 0x09; //3to5
static uint8_t GDC= 0x0b;
static uint8_t GDD= 0x0f;
//0110**: 24-27; 01**10: 18, 22, 30; 01****: the rest of 16-31
//1001**: 36-39; 10**01: 33, 41, 45; 10****: the rest of 32-47
//encode not as 1,3,5; have problem with print 5 (or 13) when putcns
extern const uint8_t cg_context_table[64] = {
    0,   0,   0,   0,    0,   0,   0,   0, 
    0,   0,   0,   0,    0,   0,   0,   0, 
    CHH, CHH, CHG, CHH,  CHH, CHH, CHG, CHH, 
    CG,  CG,  CG,  CG,	 CHH, CHH, CHG, CHH, 
    GDD, GDC, GDD, GDD,  GC,  GC,  GC,  GC, 
    GDD, GDC, GDD, GDD,  GDD, GDC, GDD, GDD,
    0,   0,   0,   0,    0,   0,   0,   0, 
    0,   0,   0,   0,    0,   0,   0,   0, 
};

extern const uint8_t cg_table[5] = {0, 1, 1, 0, 0}; // for C/G check



//  global variables, only changed at program start.
int mean_insert=500, sd_insert=50, min_insert=100, max_insert=1000, size_l=70, size_r=70, sd_center=50;

// if the leftmost 4 bit is non-zero, then it must be snp or indel
enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
static mut_t mutmsk = (mut_t)0xf000;

typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
} mutseq_t;


std::vector<float> eff_vec;
std::vector<snp_rec> snp_vec;
// std::vector<cut_rec> cut_vec;
// std::vector<cut_pos>  cutpos_vec;
std::vector<meth_rec> meth_vec;
std::vector<param_rec> param_vec;
std::vector<frag_rec> frag_vec;
std::vector<probe_rec> probe_vec;
std::map<std::string, chr_rec> chr_count;

// initialize random generator for general distributions
std::random_device rd;
std::mt19937 gen(rd());
// initialize random generator for standard normal distribution 
std::random_device rn;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen_rn(rn()); //Standard mersenne_twister_engine seeded with rd()
std::normal_distribution<float> dis_rn(0.0, 1.0);
// initialize random generator for uniform distribution between [0,1]
std::random_device ru;
std::mt19937 gen_ru(ru());
std::uniform_real_distribution<float> dis_ru(0.0,1.0);


void sim_mut_vcf(const kseq_t *ks, char * vcf_file, mutseq_t *hap1, mutseq_t *hap2, uint32_t *pos_idx_arr) 
{
    // initiate
    mutseq_t *ret[2];

    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    
    // parse VCF
    parse_vcf_chr(vcf_file, ks->name.s, snp_vec);

    int vec_ptr = 0;
    int i, deleting = 0;
    int deletion_count = 0;
    int c;

    for (i = 0; i != ks->seq.l; ++i){
        c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (cg_table[c]) {pos_idx_arr[i]= 2;}
        if (snp_vec.size() == 0){ continue;} // ignore the rest if there is no SNP

        if (deleting){
            if(deletion_count > 0){
                if (deleting & 1){ ret[0]->s[i] |= DELETE;}
                if (deleting & 2){ ret[1]->s[i] |= DELETE;}
                deletion_count--;
                pos_idx_arr[i] |= 1;
                continue;
            } else {deleting = 0;}
        }

        if(vec_ptr < snp_vec.size() && i == snp_vec[vec_ptr].pos && c < 4){
            int geno_int = snp_vec[vec_ptr].geno;
            int is_phased= geno_int & 0x000f;
            int snp_hap1 = (geno_int & 0x003f) >> 4;
            int snp_hap2 = (geno_int & 0x00ff) >> 6;
            int ref_len  = (geno_int & 0x0fff) >> 8;
            int base_offset = geno_int >> 12;
            pos_idx_arr[i] |= 1;
            //fprintf(stderr, "%d,%d,%d,%d,%d,%d\n", i, is_pahsed,snp_hap1, snp_hap2, ref_len, base_offset);

            if (!is_phased){ // for unphased genotype, randomly swap the haplotype
                if(drand48() < 0.5){
                    int tmp_hap = snp_hap1;
                    snp_hap1 = snp_hap2;
                    snp_hap2 = tmp_hap;
                }
            }

            if(base_offset == 0 && ref_len == 1){ // SNP substitution
                c = snp_vec[vec_ptr].alt;

                if (snp_hap1 == 1 && snp_hap2 == 1){
                    ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                } else if (snp_hap1 == 1 && snp_hap2 == 0){
                    ret[0]->s[i] = SUBSTITUTE|c;
                } else if (snp_hap1 == 0 && snp_hap2 == 1){
                    ret[1]->s[i] = SUBSTITUTE|c;
                } else{continue;}
            } else if (base_offset < 0 ) { // deletion
                c = snp_vec[vec_ptr].ref;
                deletion_count = abs(base_offset) - 1; //minus one because here already delete one base

                if (snp_hap1 == 1 && snp_hap2 == 1){
                    ret[0]->s[i] = ret[1]->s[i] =  DELETE;
                    deleting = 3;
                } else if (snp_hap1 == 1 && snp_hap2 == 0){
                    ret[0]->s[i] =  DELETE;
                    deleting = 1;
                } else if (snp_hap1 == 0 && snp_hap2 == 1){
                    ret[1]->s[i] =  DELETE;
                    deleting = 2;
                } else{continue;}
            } else if (base_offset > 0){ // inserstion
                int num_ins = base_offset;
                int ins_msk = (1 << (num_ins*2)) - 1;
                int ins = snp_vec[vec_ptr].alt & ins_msk;
                //fprintf(stderr, "%d,%d,%d,%d\n", num_ins, ins_msk, alt_vec[vec_ptr], ins);

                if (snp_hap1 == 1 && snp_hap2 == 1){
                    ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else if (snp_hap1 == 1 && snp_hap2 == 0){
                    ret[0]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else if (snp_hap1 == 0 && snp_hap2 == 1){
                    ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else{continue;}
            }
            ++vec_ptr;
        }
    }
}

void sim_mut_diref(const kseq_t *ks, bool is_hap, mutseq_t *hap1, mutseq_t *hap2, uint32_t *pos_idx_arr)
{
    int i, deleting = 0;
    mutseq_t *ret[2];

    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    for (i = 0; i != ks->seq.l; ++i) {
        int c;
        c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (cg_table[c]) {pos_idx_arr[i]= 2;}
        if (deleting) {
            if (drand48() < INDEL_EXTN) {
                if (deleting & 1) ret[0]->s[i] |= DELETE;
                if (deleting & 2) ret[1]->s[i] |= DELETE;
                pos_idx_arr[i] |= 1;
                continue;
            } else deleting = 0;
        }
        if (c < 4 && drand48() < MUT_RATE) { // mutation
            if (drand48() >= INDEL_FRAC) { // substitution
                double r = drand48();
                c = (c + (int)(r * 3.0 + 1)) & 3;
                if (is_hap || drand48() < 0.333333) { // hom
                    ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                } else { // het
                    ret[drand48()<0.5?0:1]->s[i] = SUBSTITUTE|c;
                }
            } else { // indel
                if (drand48() < 0.5) { // deletion
                    if (is_hap || drand48() < 0.333333) { // hom-del
                        ret[0]->s[i] = ret[1]->s[i] = DELETE;
                        deleting = 3;
                    } else { // het-del
                        deleting = drand48()<0.5?1:2;
                        ret[deleting-1]->s[i] = DELETE;
                    }
                } else { // insertion
                    int num_ins = 0, ins = 0;
                    do {
                        num_ins++;
                        ins = (ins << 2) | (int)(drand48() * 4.0);
                    } while (num_ins < 4 && drand48() < INDEL_EXTN);

                    if (is_hap || drand48() < 0.333333) { // hom-ins
                        ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                    } else { // het-ins
                        ret[drand48()<0.5?0:1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                    }
                }
            }
            pos_idx_arr[i] |= 1;
        }
    }
}

void sim_print_mutref(const char *name, const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2)
{
    int i, j = 0; // j keeps the end of the last deletion
    for (i = 0; i != ks->seq.l; ++i) {
        int c[3];
        c[0] = nst_nt4_table[(int)ks->seq.s[i]];
        c[1] = hap1->s[i]; c[2] = hap2->s[i];
        if (c[0] >= 4) continue;
        if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
            if (c[1] == c[2]) { // hom
                if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t-\n", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]); // coordinate is 1-based
                } else if ((c[1]&mutmsk) == DELETE) { // del
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+1);
                        for (j = i; j < ks->seq.l && hap1->s[j] == hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (((c[1] & mutmsk) >> 12) <= 4) { // ins
                    printf("%s\t%d\t-\t", name, i+1);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t-\n");
                } // else: deleted base in a long deletion
            } else { // het
                if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t+\n", name, i+1, "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
                } else if ((c[1]&mutmsk) == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+1);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if ((c[2]&mutmsk) == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+1);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap2->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (((c[1] & mutmsk) >> 12) <= 4 && ((c[1] & mutmsk) >> 12) > 0) { // ins1
                    printf("%s\t%d\t-\t", name, i+1);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t+\n");
                } else if (((c[2] & mutmsk) >> 12) <= 4 || ((c[2] & mutmsk) >> 12) > 0) { // ins2
                    printf("%s\t%d\t-\t", name, i+1);
                    int n = (c[2]&mutmsk) >> 12, ins = c[2] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t+\n");
                } // else: deleted base in a long deletion
            }
        }
    }
}

void sim_print_mutref_0base(const char *name, const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2)
{
    int i, j = 0; // j keeps the end of the last deletion
    for (i = 0; i != ks->seq.l; ++i) {
        int c[3];
        c[0] = nst_nt4_table[(int)ks->seq.s[i]];
        c[1] = hap1->s[i]; c[2] = hap2->s[i];
        if (c[0] >= 4) continue;
        if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
            if (c[1] == c[2]) { // hom
                if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t-\n", name, i, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]); // coordinate is 1-based
                } else if ((c[1]&mutmsk) == DELETE) { // del
                    if (i >= j) {
                        printf("%s\t%d\t", name, i);
                        for (j = i; j < ks->seq.l && hap1->s[j] == hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (((c[1] & mutmsk) >> 12) <= 4) { // ins
                    printf("%s\t%d\t-\t", name, i);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t-\n");
                } // else: deleted base in a long deletion
            } else { // het
                if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t+\n", name, i, "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
                } else if ((c[1]&mutmsk) == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if ((c[2]&mutmsk) == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap2->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (((c[1] & mutmsk) >> 12) <= 4 && ((c[1] & mutmsk) >> 12) > 0) { // ins1
                    printf("%s\t%d\t-\t", name, i);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t+\n");
                } else if (((c[2] & mutmsk) >> 12) <= 4 || ((c[2] & mutmsk) >> 12) > 0) { // ins2
                    printf("%s\t%d\t-\t", name, i);
                    int n = (c[2]&mutmsk) >> 12, ins = c[2] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        n--;
                    }
                    printf("\t+\n");
                } // else: deleted base in a long deletion
            }
        }
    }
}

bool compare_frag_rec(const frag_rec &a, const frag_rec &b){
    return a.pos_l < b.pos_l;
}

bool compare_frag_rec2(const frag_rec &a, const frag_rec &b){
    if (a.haplo < b.haplo) {
        return true;
    } else if (a.haplo > b.haplo) {
        return false;
    } else {
        return a.pos_l < b.pos_r;
    }
}

void gen_frag_vec(std::uniform_int_distribution<int> *dis_ud, std::discrete_distribution<int> *dis_dd, 
                    uint32_t *pos_idx_arr, std::vector<frag_rec> &frag_vec, int chunk_size, int tech_mode, bool is_uniform)
{
    frag_rec tmp_frag;
    int pos_l, pos_r, insert_dev, insert_len, frag_idx, probe_center, frag_center, i;

    frag_vec.clear();
    if(tech_mode ==2){
        for(i = 0; i < chunk_size; ++i){
            frag_idx = is_uniform ? (*dis_ud)(gen) : (*dis_dd)(gen);
            probe_rec tmp_probe = probe_vec[frag_idx];
            probe_center = (tmp_probe.pos_l + tmp_probe.pos_r) >> 1;
            frag_center= probe_center + (int)(sd_center * dis_rn(gen_rn));
            insert_dev = (int)(sd_insert * dis_rn(gen_rn));
            insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));
            
            tmp_frag.pos_l = frag_center - insert_len>>1;
            tmp_frag.pos_r = frag_center + insert_len>>1; 
            tmp_frag.strand= tmp_probe.strand;     // denotes the strand
            tmp_frag.haplo = drand48()<0.5?0:1;
            frag_vec.push_back(tmp_frag);
            tmp_frag = {};
        }
    }else if (tech_mode == 1){
        for(i = 0; i < chunk_size; ++i){
            frag_idx = is_uniform ? (*dis_ud)(gen) : (*dis_dd)(gen);
            probe_rec tmp_probe = probe_vec[frag_idx];
            tmp_frag.pos_l = tmp_probe.pos_l;
            tmp_frag.pos_r = tmp_frag.pos_r;
            tmp_frag.haplo = tmp_probe.strand;     // denotes the haplotype
            frag_vec.push_back(tmp_frag);
            tmp_frag = {};
        }
    }else{
        if(is_uniform){
            for(i = 0; i < chunk_size; ++i){
                pos_l = (*dis_ud)(gen);
                insert_dev = (int)(sd_insert * dis_rn(gen_rn));
                insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));
                //pos_r = std::min(pos_l + insert_len, tot_size -2); //will not pass boundary
                tmp_frag.pos_l = pos_l;
                tmp_frag.pos_r = pos_l+insert_len;
                tmp_frag.haplo = drand48()<0.5?0:1;
                frag_vec.push_back(tmp_frag);
                tmp_frag = {};
            }
        }else{
            int bin_size = eff_vec.size();
            while (frag_vec.size()< chunk_size){
                pos_l = (*dis_ud)(gen);
                insert_dev = (int)(sd_insert * dis_rn(gen_rn));
                insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));
                pos_r = pos_l + insert_len;
                int gc_count =0;
                for(int kk = pos_l; kk <= pos_r; ++kk){gc_count += (pos_idx_arr[kk] & 0x2)>>1;}
                float gc_prob = eff_vec[(int)(gc_count*bin_size/insert_len+0.5)];
                if(dis_ru(gen_ru) > gc_prob){   // when initiate eff_vec, judge the value in case it's too small
                    tmp_frag.pos_l = pos_l;
                    tmp_frag.pos_r = pos_r;
                    tmp_frag.haplo = drand48()<0.5?0:1;
                    frag_vec.push_back(tmp_frag);
                    tmp_frag= {};
                }
            }
        }
    }
    if(tech_mode==1){
        std::sort(frag_vec.begin(), frag_vec.end(), compare_frag_rec2);
    }else{
        std::sort(frag_vec.begin(), frag_vec.end(), compare_frag_rec);
    }
}

void check_frag_vec(std::vector<frag_rec> &frag_vec, mutseq_t *hap1, mutseq_t *hap2, int tech_mode){
    // find out the start position for read2, check if the boundaries satisfy
    if(tech_mode == 1){
        mutseq_t *ret[2];
        ret[0] = hap1; ret[1] = hap2;
        for(size_t i =0; i < frag_vec.size(); ++i){
            int start2 = frag_vec[i].pos_r;
            int haplo  = frag_vec[i].haplo;
            for(int k=0; k < size_r; k++){
                int c = ret[haplo]->s[start2], mut_type = c & mutmsk;
                if(mutmsk == DELETE){
                    --start2;
                    --k;
                }else if(mutmsk == INSERT){
                    int num_ins = mut_type>>12;
                    if(k + num_ins > size_r){
                        start2 -= size_r - k;
                    }else{
                        start2 -= num_ins;
                        k += num_ins;
                    }
                }else{
                    --start2;
                }
            }
        }
    }else{
        for(size_t i =0; i < frag_vec.size(); ++i){
            frag_vec[i].start2 = frag_vec[i].pos_r - size_r;
        }
    }
}

void sim_core(const char *fn, bool is_hap, 
                bool is_uniform, int tech_mode, int output_fmt, int chunk_size,
                char *vcf_file, char *bed_file, char *chr_id,
                bool methdb_save, char *methdb_file, char *cgmap_file, bool cgmap_pool, char *asm_file, int meth_seed)
{
    kseq_t *ks;
    mutseq_t rseq[2];
    gzFile   fp_fa;
    uint64_t contig_eff_len, ii = 0;
    uint64_t tot_sub = 0, tot_indel = 0, tot_err = 0, tot_pairs = 0, n_pairs =0;
    uint8_t *tmp_seq[2];    	// save sequence
    int8_t  *tmp_offset[2]; 	// save offset per base
    uint8_t *tmp_context[2];	// save context (CG, CHG, CHH)
    uint8_t *tmp_mutation[2];	// save mutation status
    mut_t *target;
    char *qstr;
    int size[2], Q, max_size, l = 0;

    max_size = std::max(size_l, size_r);
    qstr = (char*)calloc(max_size+1, 1);
    tmp_seq[0] = (uint8_t*)calloc(max_size+2, 1);
    tmp_seq[1] = (uint8_t*)calloc(max_size+2, 1);
    tmp_offset[0]= (int8_t*)calloc(max_size+2, 1);
    tmp_offset[1]= (int8_t*)calloc(max_size+2, 1);
    tmp_context[0] = (uint8_t*)calloc(max_size+2, 1);
    tmp_context[1] = (uint8_t*)calloc(max_size+2, 1);
    tmp_mutation[0]= (uint8_t*)calloc(max_size+2, 1);
    tmp_mutation[1]= (uint8_t*)calloc(max_size+2, 1);
    size[0] = size_l; size[1] = size_r;

    // TODO: thouroughly need to check i++ and ++i

    Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;

    bool bool_chr_set   = strcmp(chr_id, "None") && strlen(chr_id);
    bool bool_vcf_set   = strcmp(vcf_file,"None") && strlen(vcf_file);
    bool bool_methdb_set= strcmp(methdb_file,"None") && strlen(methdb_file);
    bool bool_cgmap_set = strcmp(cgmap_file,"None") && strlen(cgmap_file);
    bool bool_asm_set   = strcmp(asm_file,"None") && strlen(asm_file);
    bool bool_cgmap_pool= cgmap_pool;
    bool bool_methdb_save=methdb_save;

    // start simulate
    fp_fa = gzopen(fn, "r");
    ks = kseq_init(fp_fa);

    while ((l = kseq_read(ks)) >= 0) {  //here l is the chromosome length
        if (bool_chr_set) {if (strcmp(chr_id, ks->name.s)!=0){continue;}}
        if (l < mean_insert + 3 * sd_insert) {continue;} //already output skip message in the main
        n_pairs = chr_count[std::string(ks->name.s)].count;
        if(!n_pairs){continue;} // skip if no reads simulated from this contig

        // print out simulation information
        tot_pairs += n_pairs;
        fprintf(stderr, "[%s] contig '%s': simulate %ld reads...\n", __func__, ks->name.s, n_pairs);

        // create the pos_idx array, last 2 bits records, 1. whether it's C/G, 2. if it's a SNP position
        uint32_t* pos_idx_arr = (uint32_t*)malloc(ks->seq.l * sizeof(uint32_t));
        if (pos_idx_arr == NULL) {fprintf(stderr, "ERROR: could not allocate memory\n");exit(EXIT_FAILURE);} else {
            // Initialize array as 0, record if a site is SNP/INDEL position (base can be either REF/ALT)
            // uint8_t pos_idx_arr[ks->seq.l] = {0}; // only work when length is small, otherwise overflow
            memset(pos_idx_arr, 0, ks->seq.l * sizeof(uint32_t));
        }

        // introduce mutations and print them to stdout
        fprintf(stdout, "Contig Variant Start\n");
        if(bool_vcf_set){
            sim_mut_vcf(ks, vcf_file, rseq, rseq+1, pos_idx_arr);
            if(snp_vec.size() == 0){fprintf(stdout, "%s\n", ks->name.s);}   //if no variants, print chromosome id
        } else {
            sim_mut_diref(ks, is_hap, rseq, rseq+1, pos_idx_arr);
            if(MUT_RATE == 0.0){fprintf(stdout, "%s\n", ks->name.s);}       //if no variants, print chromosome id            
        }
        if(output_fmt){sim_print_mutref_0base(ks->name.s, ks, rseq, rseq+1);}else{sim_print_mutref(ks->name.s, ks, rseq, rseq+1);}
        fprintf(stdout, "Contig Variant End\n");

        // load or create the methdb
        if(bool_methdb_set){
            // can check methdb_file's filename is the same as ks.name
            load_methdb(pos_idx_arr, meth_vec, methdb_file); 
        }else{
            create_methdb(ks, pos_idx_arr, meth_vec);
            if(bool_cgmap_set){fill_cgmap_chr(cgmap_file, ks->name.s, pos_idx_arr, meth_vec, bool_cgmap_pool, meth_seed);}
            if(bool_asm_set){fill_asm_chr(asm_file, ks->name.s, pos_idx_arr, meth_vec);}
            fill_beta(meth_vec, param_vec, meth_seed);
            if(bool_methdb_save){save_methdb(meth_vec, methdb_file);}
        }

        // initialize random number generator to generate read positions
        uint64_t unif_begin = 2, unif_end = ks->seq.l-max_insert-2;     // ensure read doesn't pass boundary with 2 base offset
        std::vector<float> weights;

        if(tech_mode){
            parse_bed_chr(bed_file, ks->name.s, probe_vec);
            if(probe_vec.size() == 0){continue;} // parse probe, skip if empty
            if(is_uniform){unif_begin = 0; unif_end = probe_vec.size()-1;}else{
                std::vector<float> weights(probe_vec.size());
                for(int i=0; i < probe_vec.size(); ++i){ weights.push_back(probe_vec[i].score);}
            }
        }
        std::uniform_int_distribution<int> dis_ud(unif_begin, unif_end);
        std::discrete_distribution<int> dis_dd(weights.begin(), weights.end()); //TODO: undefined behavior if weights empty
        // std::vector<float>().swap(weights); //  wights should remain valid for the lifetime of dis_dd

        frag_rec tmp_frag;
        while(ii < n_pairs){// the core loop
            // first generate #chunk_size fragments
            gen_frag_vec(&dis_ud, &dis_dd, pos_idx_arr, frag_vec, chunk_size, tech_mode, is_uniform);
            //target = rseq[drand48()<0.5?0:1].s; // haplotype from which the reads are generated
            // generate the read sequences

            int tmp_size = std::min(chunk_size, (int)(n_pairs - ii));
            ii += tmp_size;
            for(int idx; idx < tmp_size; ++idx){
                tmp_frag = frag_vec[idx];
                //cover_pos hold if the read covers a snp *position* (the read don't have to contain the ALT allele)
                //j hold read1/read2, k hold the length of read, ix hold the cursor transversing read
                int n_sub[2]={0,0}, n_indel[2]={0,0}, n_err[2]={0,0}, cover_pos[2]={0,0}; 
                int ext_coor[2], i, j, k, ix;
                int start[2] = {tmp_frag.pos_l, tmp_frag.start2};
                int end[2] = {start[0], start[1]};
                int offset[2] = {0, 0};

                // x: select read1 or read2; ext_coor: extend corrdinates;
                #define __gen_read(x, start_pos, iter) do {                     \
                    /* generate reads assign mutation flag; */                  \
                    for (i = (start_pos), k = 0, ext_coor[x] = -10; i >= 0 && i < ks->seq.l && k < size[x]; iter) { \
                        int c = target[i], mut_type = c & mutmsk;               \
                        if (ext_coor[x] < 0) {                                  \
                            /* avoid indel as the first base */                 \
                            if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue; \
                            start[x] = i;                                       \
                            end[x] = i;                                         \
                            ext_coor[x] = i;                                    \
                        }                                                       \
                        if (mut_type == DELETE){                                \
                            ++offset[x];                                        \
                            ++end[x];                                           \
                            ++n_indel[x];                                       \
                        }                                                       \
                        else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
                            /* context: 0x00 Match, 0x01 SNP, 0x03 INSERT       \
                                        0x01 CG, 0x03 CHG, 0x07 CHH (>>)        \
                                        0x09 GC, 0x0b GDC, 0x0f GDD (<<) */     \
                            tmp_seq[x][k] = c & 0xf;                            \
                            tmp_offset[x][k] = offset[x];                       \
                            if (mut_type == SUBSTITUTE) {                       \
                                ++n_sub[x];                                     \
                                tmp_mutation[x][k] = SNV;                       \
                            } else {                                            \
                                tmp_mutation[x][k] = MATCH;                     \
                            }                                                   \
                            ++end[x];                                           \
                            ++k;                                                \
                        } else {                                                \
                            tmp_seq[x][k] = c & 0xf;                            \
                            tmp_offset[x][k] = offset[x];                       \
                            tmp_mutation[x][k] = MATCH;/*The base is ref*/      \
                            ++n_indel[x];                                       \
                            ++end[x];                                           \
                            ++k;                                                \
                            int num_ins, ins;                                   \
                            for (num_ins = mut_type>>12, ins = c>>4; num_ins > 0 && k < size[x]; --num_ins, ins >>= 2){ \
                                --offset[x];                                    \
                                tmp_seq[x][k] = ins & 0x3;                      \
                                tmp_offset[x][k] = offset[x];                   \
                                tmp_mutation[x][k] = INSR;                      \
                                ++k;                                            \
                            }                                                   \
                        }                                                       \
                        cover_pos[x] |= pos_idx_arr[i];                         \
                    }                                                           \
                    /* append CG context flag;                                  \
                    currently not handling bounday context by indel*/           \
                    for (ix=0; ix < k; ++ix) {                                  \
                        int c_d1, c_d2;                                         \
                        int c = tmp_seq[x][ix];                                 \
                        uint8_t context = 0;                                    \
                        if (cg_table[(uint8_t) c]){                             \
                            if (c == 1) {                                       \
                                /*handle the last 2 base*/                      \
                                if(ix > k-3){                                   \
                                    int ix_ext = k - ix; /*think if ix=k-1*/    \
                                    c_d1 = target[end[x]+ix_ext];               \
                                    c_d2 = target[end[x]+ix_ext+1];             \
                                } else {                                        \
                                    c_d1 = tmp_seq[x][ix+1];                    \
                                    c_d2 = tmp_seq[x][ix+2];                    \
                                }                                               \
                            } else {                                            \
                                /*handle the first 2 base*/                     \
                                if(ix < 2){                                     \
                                    int ix_ext = 2 - ix; /*think if ix=1 */     \
                                    c_d1 = target[start[x]-ix_ext+1];           \
                                    c_d2 = target[start[x]-ix_ext];             \
                                } else {                                        \
                                    c_d1 = tmp_seq[x][ix-1];                    \
                                    c_d2 = tmp_seq[x][ix-2];                    \
                                }                                               \
                            }                                                   \
                            uint8_t context_idx = c << 4 | c_d1 <<2 | c_d2;     \
                            context = cg_context_table[context_idx];            \
                        }                                                       \
                    tmp_context[x][ix] = context;                               \
                    }                                                           \
                    if (k != size[x]) {ext_coor[x] = -10;}                      \
                } while (0)

                __gen_read(0, start[0], ++i);
                __gen_read(1, start[1], ++i);

                if (ext_coor[0] < 0 || ext_coor[1] < 0) { // failed to generate the read(s)
                    --ii;
                    continue;
                }
                for(j = 0; j < 2; ++j){ //check the number of Ns
                    int n_n =0;
                    for (i = 0; i < size[j]; ++i) {
                        int c = tmp_seq[j][i];
                        if (c >= 4) { // actually c should be never larger than 4 if everything is correct
                            ++n_n; 
                            tmp_seq[j][i] = 4;
                        }
                        qstr[i] = Q; // generate the quality score
                    }
                    qstr[i] = 0;
                    if ((double)n_n / size[j] > MAX_N_RATIO) break;
                }
                if (j < 2) { // too many ambiguous bases on one of the reads
                    --ii;
                    continue;
                }
                
                int flag_pos, flag_mut; //flag_pos: whether has mutation position; flag_mut: whether has mutation
                for(i=tmp_frag.pos_l; i < tmp_frag.pos_r; ++i){
                    flag_pos |= pos_idx_arr[i];
                    flag_mut |= target[i] & mutmsk;
                }
                int flag_indel= n_indel[0] | n_indel[1];
                
                //TODO: also output the methylation value

                // print reads to stdout: mode 0 print string (WGS), else print numbers (WGBS)
                if(output_fmt == 0){
                    // flip and get the reverse complementary
                    int is_flip = drand48() < 0.5? 0 : 1;
                    for (k = 0; k < size[1]; ++k) { 
                        if (k <= int(size[1]/2)) { // reverse
                            int tmp_base  = tmp_seq[1][k];
                            tmp_seq[1][k] = tmp_seq[1][size[1]-k];
                            tmp_seq[1][size[1]-k] = tmp_base;
                            int tmp_cigar = tmp_mutation[1][k];
                            tmp_mutation[1][k] = tmp_mutation[1][size[1]-k];
                            tmp_mutation[1][size[1]-k] = tmp_cigar;
                        }
                        tmp_seq[1][k] = tmp_seq[1][k] < 4? 3 - tmp_seq[1][k] : 4; // complement
                    }
                    for (j = 0; j < 2; ++j) {
                        int jj = j ^ is_flip; // x^0=x; x^1=!x (when x is binary)
                        // header: 1-based coordinates for string output
                        fprintf(stdout, "@%s:%d:%d:%llx:%d:%d:%d/%d\n", ks->name.s, start[0]+1, end[1]+1, (long long)ii, flag_pos, flag_mut, flag_indel, j+1);
                        // sequence (introduce random sequencing error)
                        for (i = 0; i < size[jj]; ++i) {
                            int c = tmp_seq[jj][i];
                            if (drand48() < ERR_RATE){
                                // c = (c + (int)(drand48() * 3.0 + 1)) & 3; // random sequencing errors
                                c = (c + 1) & 3; // recurrent sequencing errors
                                ++n_err[jj];
                                tmp_mutation[jj][i] |= SEQERR;
                                tmp_seq[jj][i] = c;
                            }
                            fputc("ACGTN"[c], stdout);
                        }
                        fprintf(stdout, "\n");
                        // comment
                        fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:%d:", start[jj]+1, end[jj]+1, cover_pos[jj], n_sub[jj], n_indel[jj], n_err[jj], end[1]-start[0], start[1]-end[0]);
                        for (i = 0; i < size[jj]; ++i) {
                            int c = (tmp_mutation[jj][i] & 0x0f);
                            fputc("MXIE"[mut_table[c]], stdout);
                        }
                        fprintf(stdout, "\n");
                        // quality
                        fprintf(stdout, "%s\n", qstr);
                    }
                } else {
                    for (j = 0; j < 2; ++j) {
                        // header: 0-based coordinates for number output
                        fprintf(stdout, "@%s:%d:%d:%llx %d %d %d %d %d ", ks->name.s, start[0]+1, end[1]+1, (long long)ii, j, flag_pos, flag_mut, flag_indel, tmp_frag.strand); 
                        for (i = 0; i < size[j]; ++i) {
                            //fprintf(stdout, "%x", tmp_mutation[j][i]); //this will output the hex number, below line will output the ascii
                            fputc(tmp_mutation[j][i], stdout);
                        }
                        fprintf(stdout, "\n");
                        // sequence (no sequencing error, represented by 0-4)
                        for (i = 0; i < size[j]; ++i) {
                            //fprintf(stdout, "%d", tmp_seq[j][i]);
                            fputc(tmp_seq[j][i], stdout);
                        }
                        fprintf(stdout, "\n");
                        // comment
                        fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:", start[j], end[j], cover_pos[j], n_sub[j], n_indel[j], end[1]-start[0], start[1]-end[0]);
                        const char *pad = "";
                        for (i = 0; i < size[j]; ++i) {
                            fprintf(stdout, "%s%d", pad, tmp_offset[j][i]);
                            pad = ",";
                        }
                        fprintf(stdout, "\n");
                        // quality
                        for (i = 0; i < size[j]; ++i) {
                            //fprintf(stdout, "%x", tmp_context[j][i]);
                            fputc(tmp_context[j][i], stdout);
                        }
                        fprintf(stdout, "\n");
                    }
                }
                tot_sub   += (int)(n_sub[0]  + n_sub[1] > 0);
                tot_indel += (int)(n_indel[0]+ n_indel[1] > 0);
                tot_err   += (int)(n_err[0]  + n_err[1] > 0);
            }
        }
        free(rseq[0].s); free(rseq[1].s);
        free(pos_idx_arr);
    }

    fprintf(stderr, "[%s] Generated %lu read pairs, with %lu contain SNP, %lu contain INDEL", __func__, tot_pairs, tot_sub, tot_indel);
    if (output_fmt == 0){fprintf(stderr, " and %lu contain sequencing errors\n", tot_err);} else{fprintf(stderr, "\n");}
    kseq_destroy(ks);
    gzclose(fp_fa);
    free(qstr);
    free(tmp_seq[0]); free(tmp_seq[1]);
    free(tmp_offset[0]); free(tmp_offset[1]);
    free(tmp_context[0]); free(tmp_context[1]);
    free(tmp_mutation[0]); free(tmp_mutation[1]);
}


static int simu_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "htsim (high throughput reads simulator) for WGS or WGBS/RRBS/TBS reads simulation\n");
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Contact: Wenbin Guo <wbguo@ucla.edu>; \n\n");
    fprintf(stderr, "Usage:   htsim [options] <ref.fa> \n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "general setting:\n");
    fprintf(stderr, "         -I INT        mean insert size (outer distance between 2 ends) [500]\n");
    fprintf(stderr, "         -J INT        standard deviation of insert size [50]\n");
    fprintf(stderr, "         -K INT        minimum insert size [100]\n");
    fprintf(stderr, "         -L INT        maximum insert size [1000]\n");
    fprintf(stderr, "         -c STRING     contig name, only output reads from this contig, default is to output all contigs [None]\n");
    fprintf(stderr, "         -C INT        chunk size for read generation [5000]\n");
    fprintf(stderr, "         -n INT        number of read pairs to generate for specified contig, disabled by default [0]\n");
    fprintf(stderr, "         -N INT        total number of read pairs to generate, will caculate from depth when N is 0 [0]\n");
    fprintf(stderr, "         -d INT        average sequencing depth, only used when n/N is not specified [10]\n");
    fprintf(stderr, "         -1 INT        length of the first read [70]\n");
    fprintf(stderr, "         -2 INT        length of the second read [70]\n");
    fprintf(stderr, "         -E FLOAT      base error rate (set to be 0 for bisulfite sequencing) [%.3f]\n", ERR_RATE);
    fprintf(stderr, "         -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [%.2f]\n", MAX_N_RATIO);
    fprintf(stderr, "         -O INT        output format: 0 for letters; 1 for ascii numbers (for python module) [0]\n");
    fprintf(stderr, "mutation setting:\n");
    fprintf(stderr, "         -v STRING     path to the genetic variant file (.vcf/vcf.gz) [None]\n");
    fprintf(stderr, "         -R FLOAT      rate of mutations [%.4f]\n", MUT_RATE);
    fprintf(stderr, "         -F FLOAT      fraction of indels [%.2f]\n", INDEL_FRAC);
    fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", INDEL_EXTN);
    fprintf(stderr, "         -H INT        haplotype mode: 0 for disable, nonzero for enable (all variants are homozygotes) [0]\n");
    fprintf(stderr, "         -s INT        seed for random generator [-1]\n");
    fprintf(stderr, "MethDB setting:\n");
    fprintf(stderr, "         -a STRING     ASM file for allelic-specific methylation simulation (.asm/.asm.gz) [None]\n");
    fprintf(stderr, "         -m STRING     methylation reference (.CGmap/.CGmap.gz) [None]\n");
    fprintf(stderr, "         -p INT        pool CGmap: 0 for disable, nonzero for enable (randomly draw value based on context) [0]\n");
    fprintf(stderr, "         -W INT        write MethDB to disk: 0 for disable, nonzero for enable [0]\n");
    fprintf(stderr, "         -M STRING     MethDB intermediate file [None]\n");
    fprintf(stderr, "         -S INT        seed for methylation value generation (for cgmap pool or beta distribution) [-1]\n");
    fprintf(stderr, "         -P STRING     Parameter string for beta distribution [-1]\n");
    fprintf(stderr, "technology setting:\n");
    fprintf(stderr, "         -T INT        technology: 0 for Whole genome; 1 for Reduced representation; 2 for Targeted [0]\n");
    fprintf(stderr, "         -u INT        uniform coverage: 0 for diable, nonzero for enable [1]\n");
    fprintf(stderr, "         -B STRING     GC Bias reference for WGS/WGBS, only used when -u set to be 0 [None]\n");
    fprintf(stderr, "         -x STRING     enzyme cutting site string for reduced representation sequencing [None]\n");
    fprintf(stderr, "         -b STRING     BED file for reduced-representation / targeted sequencing (.bed/.bed.gz) [None]\n");
    fprintf(stderr, "         -D INT        fragment center's deviaiton from the probe center [50]\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[])
{
    //default
    uint64_t N  = 0, chr_N  = 0;
    bool is_hap = false, is_uniform = true;
    int tech_mode = 0, output_fmt = 0, seed_snp = -1, seed_meth = -1, depth = 10, methdb_save =0, cgmap_pool=0, chunk_size=5000;

    char none_default[] = "None";
    char *chr_id    = none_default;
    char *vcf_file  = none_default;
    // char *cut_str   = none_default;
    char *bias_file = none_default;
    char *bed_file  = none_default; 
    char *asm_file  = none_default; 
    char *cgmap_file= none_default; 
    char *param_str = none_default; 
    char *methdb_file= none_default; // checked, will not intefere with vcf_file


    //update default from command line
    int c = 0;
    while ((c = getopt(argc, argv, "I:J:K:L:c:C:n:N:d:1:2:E:A:O:v:R:F:X:H:s:a:m:p:W:M:S:P:T:u:B:x:b:D:")) >= 0) {
        switch (c) {
            case 'I': mean_insert= atoi(optarg); break;
            case 'J': sd_insert  = atoi(optarg); break;
            case 'K': min_insert = atoi(optarg); break;
            case 'L': max_insert = atoi(optarg); break;
            case 'c': chr_id     = optarg; break;
            case 'C': chunk_size = atoi(optarg); break;
            case 'n': chr_N      = atoi(optarg); break;
            case 'N': N          = atoi(optarg); break;
            case 'd': depth      = atoi(optarg); break;
            case '1': size_l     = atoi(optarg); break;
            case '2': size_r     = atoi(optarg); break;
            case 'E': ERR_RATE   = atof(optarg); break;
            case 'A': MAX_N_RATIO= atof(optarg); break;
            case 'O': output_fmt = atoi(optarg); break;
            case 'v': vcf_file   = optarg; break;
            case 'R': MUT_RATE   = atof(optarg); break;
            case 'F': INDEL_FRAC = atof(optarg); break;
            case 'X': INDEL_EXTN = atof(optarg); break;
            case 'H': is_hap     = atoi(optarg)!=0; break;
            case 's': seed_snp   = atoi(optarg); break;
            case 'a': asm_file   = optarg; break;
            case 'm': cgmap_file = optarg; break;
            case 'p': cgmap_pool = atoi(optarg); break;
            case 'W': methdb_save= atoi(optarg); break;
            case 'M': methdb_file= optarg; break;
            case 'S': seed_meth  = atoi(optarg); break;
            case 'P': param_str  = optarg; break;
            case 'T': tech_mode  = atoi(optarg); break;
            case 'u': is_uniform = atoi(optarg)!=0; break;
            case 'B': bias_file  = optarg; break;
            // case 'x': cut_str    = optarg; break;
            case 'b': bed_file   = optarg; break;
            case 'D': sd_center  = atoi(optarg); break;
        }
    }
    if (argc - optind < 1) return simu_usage();
    if (seed_snp <= 0) seed_snp  = time(0)&0x7fffffff;
    if (seed_meth<= 0) seed_meth = time(0)&0x7fffffff;

    min_insert = std::max(std::max(size_l, size_r), min_insert); // ensure min_insert >= size_l or size_r

    bool bool_chr_set   = strcmp(chr_id,  "None") && strlen(chr_id);
    bool bool_vcf_set   = strcmp(vcf_file,"None") && strlen(vcf_file);
    // bool bool_site_set  = strcmp(cut_str, "None") && strlen(cut_str);
    bool bool_probe_set = strcmp(bed_file,"None") && strlen(bed_file);
    bool bool_bias_set  = strcmp(bias_file,"None") && strlen(bias_file);
    bool bool_asm_set   = strcmp(asm_file, "None") && strlen(asm_file);
    bool bool_cgmap_set = strcmp(cgmap_file, "None") && strlen(cgmap_file);
    bool bool_methdb_set= strcmp(methdb_file, "None") && strlen(methdb_file);


    // check legal input mode and corresponding files
    // if (tech_mode==2 || bool_probe_set){
    //     fprintf(stderr, "Simulating targeted sequencing reads:\n");
    //     if (!bool_probe_set){fprintf(stderr, "ERROR: Please specify probe bed file path\n");exit(EXIT_FAILURE);}
    //     if (bool_site_set){fprintf(stderr, "WARNING: Cut site is set, ignored in this mode\n");}
    //     tech_mode = 2;
    // } else if (tech_mode==1 || bool_site_set){
    //     fprintf(stderr, "Simulating restricted enzyme cutting reads:\n");
    //     if (!bool_site_set){fprintf(stderr, "ERROR: Please specify enzyme cutting site\n");exit(EXIT_FAILURE);}
    //     if (bool_probe_set){fprintf(stderr, "WARNING: Probe file is set, ignored in this mode\n");}
    //     tech_mode = 1;
    //     parse_cut_rec(cut_str, cut_vec);
    // } else {
    //     fprintf(stderr, "Simulating whole genome reads:\n");
    //     if(!is_uniform && !bool_bias_set){fprintf(stderr, "ERROR: Please specify GC-Bias file when specifying -u as 0\n");exit(EXIT_FAILURE);}
    //     if(bool_bias_set){parse_bias_file(bias_file, eff_vec); is_uniform = false;}
    //     tech_mode = 0;
    // }
    if (tech_mode==2 || bool_probe_set){
        fprintf(stderr, "Simulating targeted sequencing reads:\n");
        if (!bool_probe_set){fprintf(stderr, "ERROR: Please specify probe bed file path\n");exit(EXIT_FAILURE);}
        tech_mode = 2;
    } else if (tech_mode==1 || bool_probe_set){
        fprintf(stderr, "Simulating restricted enzyme cutting reads:\n");
        if (!bool_probe_set){fprintf(stderr, "ERROR: Please specify probe bed file path\n");exit(EXIT_FAILURE);}
        tech_mode = 1;
    } else {
        fprintf(stderr, "Simulating whole genome reads:\n");
        if(!is_uniform && !bool_bias_set){fprintf(stderr, "ERROR: Please specify GC-Bias file when specifying -u as 0\n");exit(EXIT_FAILURE);}
        if(bool_bias_set){parse_bias_file(bias_file, eff_vec); is_uniform = false;}
        tech_mode = 0;
    }


    // check existence of fasta, parse the length
    kseq_t *ks;
    gzFile fp_fa;
    fp_fa = gzopen(argv[optind], "r");
    ks = kseq_init(fp_fa);
    if (!fp_fa) { fprintf (stderr, "ERROR: gzopen of '%s' failed: %s. Exit... \n", argv[optind], strerror (errno)); exit (EXIT_FAILURE);}
    fprintf(stderr, "Reference genome file: %s\n", argv[optind]);

    fprintf(stderr, "[%s] Calculating the length and count of the reference sequences...\n", __func__);
    int l;
    chr_rec tmp_len;
    uint64_t tot_chr_len = 0, tot_eff_len = 0;
    float tot_score = 0;
    while ((l = kseq_read(ks)) >= 0) {
        if (bool_chr_set) {if (strcmp(chr_id, ks->name.s)!=0){continue;}}
        if (l < mean_insert+3*sd_insert){
            fprintf(stderr, "[%s] skip contig '%s' as it is shorter than %d!\n", __func__, ks->name.s, mean_insert+3*sd_insert); 
            continue;
        }
        
        tmp_len = {};
        collect_len_score_chr(ks, &tmp_len, bed_file, probe_vec);
        chr_count[std::string(ks->name.s)] = tmp_len;
        tot_chr_len += tmp_len.chr_len;
        tot_eff_len += tmp_len.eff_len;
        tot_score   += tmp_len.score;
    }
    kseq_destroy(ks);
    gzclose(fp_fa);


    // check if fasta is empty
    if (!chr_count.size()) { fprintf (stderr, "ERROR: Input fasta is empty: %s. Exit... \n", argv[optind]); exit (EXIT_FAILURE);}


    // check input chr_id, calculate the count for contig(s)
    if (bool_chr_set){
        // calculate the count for selected contig
        std::string chr_id_str = std::string(chr_id);
        if (!chr_count.count(chr_id_str)){fprintf(stderr, "ERROR: Contig id '%s' is not found in the fasta file, please check!\n", chr_id); exit(EXIT_FAILURE);}

        uint64_t contig_eff_len = chr_count[chr_id_str].eff_len;
        uint64_t contig_len = chr_count[chr_id_str].chr_len;
        chr_N = chr_N == 0? (contig_eff_len * depth)/(size_l + size_r) : chr_N;
        fprintf(stderr, "[%s] Contig %s specified, total length: %lu, effective length: %lu, #reads: %lu\n", __func__, chr_id, contig_len, contig_eff_len, chr_N);
        chr_count[chr_id_str].count = chr_N;
    } else {
        // calculate the count for all contigs
        if (chr_N > 0) {fprintf(stderr, "ERROR: -n is specified but not -c. Exit... (please note the difference of -n and -N)\n"); exit(EXIT_FAILURE);}
        
        int num_contigs = (int)chr_count.size();
        N = N == 0? (tot_eff_len * depth)/(size_l + size_r) : N;
        fprintf(stderr, "[%s] Found %d contig sequences, total length: %lu, effective length: %lu\n", __func__, num_contigs, tot_chr_len, tot_eff_len);
        fprintf(stderr, "[%s] No contig id specified, will generate %lu reads from all contigs\n", __func__, N);
        
        uint64_t cum_count = 0, tmp_count =0;
        for (auto it = chr_count.begin(); it != chr_count.end(); ++it) {
            tmp_count = tech_mode ==2 ? (uint64_t)(it->second.score * N / tot_score) : (uint64_t)(it->second.eff_len * N / tot_eff_len);
            it->second.count = tmp_count;
            cum_count += tmp_count;
        }

        int rest_count = N - cum_count;
        if(rest_count < 0){fprintf(stderr, "[%s] Read count calculation went wrong \n", __func__); exit(EXIT_FAILURE);} // should never happen
        int step_size  = rest_count > num_contigs ? (int)(rest_count/num_contigs): 1; //hopefully rest_count is small, evenly distributed to contigs
        auto it = chr_count.begin();
        while (rest_count > 0){ 
            int alloc_count  = std::min(step_size, rest_count);
            it->second.count+= alloc_count; 
            rest_count -= alloc_count;
            ++it;
        }
    }


    // check input vcf file
    FILE *vcf;
    if (bool_vcf_set) {
        if(vcf=fopen(vcf_file,"r")){fprintf(stderr, "[%s] VCF file exists, use it to simulate reads\n", __func__); fclose(vcf);
        }else{fprintf(stderr, "ERROR: The specified VCF file does not exist, please check!\n"); exit(EXIT_FAILURE);}
    } else {
        fprintf(stderr, "[%s] No VCF input, will generate SNP randomly if mutation rate is nonzero\n", __func__);
    }


    fprintf(stderr, "[htsim] seed = %d\n", seed_snp);
    srand48(seed_snp);

    sim_core(argv[optind], is_hap, is_uniform, tech_mode, output_fmt, chunk_size, vcf_file, bed_file, chr_id, methdb_save, methdb_file, cgmap_file, cgmap_pool, asm_file, seed_meth);

    return 0;
}

