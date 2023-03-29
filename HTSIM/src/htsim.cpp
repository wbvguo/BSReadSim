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

/* This program is derived from WGSIM(v0.3.1-r13)[https://github.com/lh3/wgsim.git], with very heavy 
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
#include <gsl/gsl_randist.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"

KSEQ_INIT(gzFile, gzread)
#include "methdb.h"
#include "mode.h"
#include "haplo.h"


#define PACKAGE_VERSION "1.0.3"


//  global variables, only changed at program start.
static float ERR_RATE  = 0.005;
static float MUT_RATE  = 0.01;
static float INDEL_FRAC= 0.15;
static float INDEL_EXTN= 0.3;
static float MAX_N_RATIO=0.05;

static int MEAN_INSERT  = 500;
static int SD_INSERT    = 50;
static int MIN_INSERT   = 100;
static int MAX_INSERT   = 1000;
static int SIZE_L       = 100;
static int SIZE_R       = 100;
static int SD_CENTER    = 50;
static int BIN_SIZE     = 100;

static uint8_t MATCH    = 0x00;
static uint8_t SNV      = 0x10;
static uint8_t INSR     = 0x30;
static uint8_t CONVT    = 0x50; // reserved, not used
static uint8_t SEQERR   = 0x90;
const  uint8_t mut_table[16] = {
    0, 1, 0, 2, 
    0, 0, 0, 0, 
    0, 3, 3, 3, 
    3, 3, 3, 3
}; // MXIE


// if the leftmost 4 bit is non-zero, then it must be snp or indel
// enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};

std::vector<float> eff_vec;
std::vector<snp_rec> snp_vec;
std::vector<frag_rec> frag_vec;
std::vector<meth_rec> meth_vec;
std::vector<param_rec> param_vec;
std::vector<probe_rec> probe_vec;
std::map<std::string, chr_rec> chr_count;


// initialize random generator for general usage
std::random_device rd; //Standard mersenne_twister_engine seeded with rd()
std::mt19937 gen(rd());
// initialize random generator for standard normal distribution
std::random_device rn;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen_rn(rn());
std::normal_distribution<float> dis_rn(0.0, 1.0); 
// initialize random generator for uniform distribution between [0,1]
std::random_device ru;
std::mt19937 gen_ru(ru());
std::uniform_real_distribution<float> dis_ru(0.0,1.0);


bool compare_frag_rec(const frag_rec &a, const frag_rec &b)
{
    if (a.haplo != b.haplo) {
        return a.haplo < b.haplo;
    } else {
        return a.pos_l < b.pos_r;
    }
}

void gen_frag_vec(std::uniform_int_distribution<int> *dis_ud, std::discrete_distribution<int> *dis_dd, 
                    uint32_t *posidx_arr, std::vector<frag_rec> &frag_vec, int chunk_size, int tech_mode, bool is_uniform)
{
    frag_rec tmp_frag;
    int pos_l, pos_r, insert_dev, insert_len, frag_idx, probe_center, frag_center, i;

    frag_vec.clear();
    if(tech_mode ==2){
        for(i = 0; i < chunk_size; ++i){
            frag_idx = is_uniform ? (*dis_ud)(gen) : (*dis_dd)(gen);
            probe_rec tmp_probe = probe_vec[frag_idx];
            probe_center = (tmp_probe.pos_l + tmp_probe.pos_r) >> 1;
            frag_center= probe_center + (int)(SD_CENTER * dis_rn(gen_rn));
            insert_dev = (int)(SD_INSERT * dis_rn(gen_rn));
            insert_len = std::max(MIN_INSERT, std::min(MEAN_INSERT + insert_dev, MAX_INSERT));
            
            tmp_frag.pos_l = frag_center - (insert_len>>1);
            tmp_frag.pos_r = frag_center + (insert_len>>1); 
            tmp_frag.strand= tmp_probe.strand;      // denotes the strand
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
            tmp_frag.strand= drand48()<0.5?0:1;     // denotes the strand
            tmp_frag.haplo = drand48()<0.5?0:1;     // can include the haplotype information
            frag_vec.push_back(tmp_frag);
            tmp_frag = {};
        }
    }else{
        if(is_uniform){
            for(i = 0; i < chunk_size; ++i){
                pos_l = (*dis_ud)(gen);
                insert_dev = (int)(SD_INSERT * dis_rn(gen_rn));
                insert_len = std::max(MIN_INSERT, std::min(MEAN_INSERT + insert_dev, MAX_INSERT));
                //pos_r = std::min(pos_l + insert_len, tot_size -2); //will not pass boundary
                tmp_frag.pos_l = pos_l;
                tmp_frag.pos_r = pos_l+insert_len;
                tmp_frag.strand= drand48()<0.5?0:1; // denotes the strand
                tmp_frag.haplo = drand48()<0.5?0:1;
                frag_vec.push_back(tmp_frag);
                tmp_frag = {};
            }
        }else{
            while ((int)frag_vec.size()< chunk_size){
                pos_l = (*dis_ud)(gen);
                insert_dev = (int)(SD_INSERT * dis_rn(gen_rn));
                insert_len = std::max(MIN_INSERT, std::min(MEAN_INSERT + insert_dev, MAX_INSERT));
                pos_r = pos_l + insert_len;
                int gc_count =0;
                for(int kk = pos_l; kk <= pos_r; ++kk){gc_count += (posidx_arr[kk] & 0x2)>>1;}
                float gc_prob = eff_vec[(int)(gc_count*BIN_SIZE/insert_len+0.5)];
                if(dis_ru(gen_ru) > gc_prob){   // when initiate eff_vec, judge the value in case it's too small
                    tmp_frag.pos_l = pos_l;
                    tmp_frag.pos_r = pos_r;
                    tmp_frag.strand= drand48()<0.5?0:1;     // denotes the strand
                    tmp_frag.haplo = drand48()<0.5?0:1;
                    frag_vec.push_back(tmp_frag);
                    tmp_frag= {};
                }
            }
        }
    }
    std::sort(frag_vec.begin(), frag_vec.end(), compare_frag_rec);
}

void check_frag_vec(std::vector<frag_rec> &frag_vec, mutseq_t *hap1, mutseq_t *hap2, int tech_mode)
{
    // find out the start position for read2, TODO: check if the boundaries satisfy for RRBS
    if(tech_mode == 1){
        mutseq_t *ret[2];
        ret[0] = hap1; ret[1] = hap2;
        for(size_t i =0; i < frag_vec.size(); ++i){
            int start2 = frag_vec[i].pos_r;
            int haplo  = frag_vec[i].haplo;
            for(int k=0; k < SIZE_R; k++){
                int c = ret[haplo]->s[start2], mut_type = c & mutmsk;
                if(mutmsk == DELETE){
                    --start2;
                    --k;
                }else if(mutmsk == INSERT){
                    int num_ins = mut_type>>12;
                    if(k + num_ins > SIZE_R){
                        start2 -= SIZE_R - k;
                    }else{
                        start2 -= num_ins;
                        k += num_ins;
                    }
                }else{
                    --start2;
                }
            }
            frag_vec[i].start2 = start2;
        }
    }else{
        for(size_t i =0; i < frag_vec.size(); ++i){
            frag_vec[i].start2 = frag_vec[i].pos_r - SIZE_R;
        }
    }
}

void sim_core(const char *fn, bool is_hap, bool is_uniform, int tech_mode, int output_fmt, int chunk_size,
                char *vcf_file, char *bed_file, char *chr_id,
                bool methdb_save, char *methdb_file, char *cgmap_file, bool cgmap_pool, char *asm_file, int meth_seed)
{
    kseq_t *ks;
    mutseq_t rseq[2];
    gzFile   fp_fa;
    uint64_t tot_sub = 0, tot_indel = 0, tot_err = 0, tot_pairs = 0, n_pairs =0, ii = 0;
    uint8_t *tmp_seq[2];    	// save sequence & mutation (upper half mutation, lower half sequence)
    int8_t  *tmp_offset[2]; 	// save offset per base
    uint8_t *tmp_context[2];	// save context (CG, CHG, CHH)
    mut_t *target;
    int size[2], Q, max_size, l = 0;

    max_size = std::max(SIZE_L, SIZE_R);
    tmp_seq[0] = (uint8_t*)calloc(max_size+2, 1);
    tmp_seq[1] = (uint8_t*)calloc(max_size+2, 1);
    tmp_offset[0]= (int8_t*)calloc(max_size+2, 1);
    tmp_offset[1]= (int8_t*)calloc(max_size+2, 1);
    tmp_context[0] = (uint8_t*)calloc(max_size+2, 1);   // save context and mutation type
    tmp_context[1] = (uint8_t*)calloc(max_size+2, 1);
    size[0] = SIZE_L; size[1] = SIZE_R;

    // TODO: thouroughly check i++ and ++i

    Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;

    bool bool_chr_set   = strcmp(chr_id, "None") && strlen(chr_id);
    bool bool_vcf_set   = strcmp(vcf_file,"None") && strlen(vcf_file);
    bool bool_methdb_set= strcmp(methdb_file,"None") && strlen(methdb_file);
    bool bool_cgmap_set = strcmp(cgmap_file,"None") && strlen(cgmap_file);
    bool bool_asm_set   = strcmp(asm_file,"None") && strlen(asm_file);
    bool bool_cgmap_pool= cgmap_pool;
    bool bool_methdb_save=methdb_save;

    mut_params tmp_params   = {};
    tmp_params.ERR_RATE     = ERR_RATE;
    tmp_params.MUT_RATE     = MUT_RATE;
    tmp_params.INDEL_FRAC   = INDEL_FRAC;
    tmp_params.INDEL_EXTN   = INDEL_EXTN;


    // start simulate
    fp_fa = gzopen(fn, "r");
    ks = kseq_init(fp_fa);

    while ((l = kseq_read(ks)) >= 0) {  //here l is the chromosome length
        if (bool_chr_set) {if (strcmp(chr_id, ks->name.s)!=0){continue;}}
        if (l < MEAN_INSERT + 3 * SD_INSERT) {continue;}
        n_pairs = chr_count[std::string(ks->name.s)].count;
        if(!n_pairs){fprintf(stderr, "[%s] skip contig '%s' for its #reads is 0...\n", __func__, ks->name.s); continue;}
        fprintf(stderr, "[%s] contig '%s': simulate %ld reads...\n", __func__, ks->name.s, n_pairs);

        // Initialize pos_idx array as 0, last 2 bits: whether it's C/G, if it's a SNP position (base can be either REF/ALT)
        // uint8_t posidx_arr[ks->seq.l] = {0}; // only work when length is small, otherwise overflow
        uint32_t* posidx_arr = (uint32_t*)malloc(ks->seq.l * sizeof(uint32_t));
        if (posidx_arr == NULL) {fprintf(stderr, "ERROR: could not allocate memory\n");exit(EXIT_FAILURE);} else {
            memset(posidx_arr, 0, ks->seq.l * sizeof(uint32_t));
        }

        // introduce mutations and print them to stdout
        fprintf(stdout, "Contig Variant Start\n");
        if(bool_vcf_set){
            sim_mut_vcf(ks, vcf_file, rseq, rseq+1, posidx_arr, snp_vec);
            if(snp_vec.size() == 0){fprintf(stdout, "%s\n", ks->name.s);}   //if no variants, print chromosome id
        } else {
            sim_mut_diref(ks, is_hap, rseq, rseq+1, posidx_arr, &tmp_params);
            if(MUT_RATE == 0.0){fprintf(stdout, "%s\n", ks->name.s);}       //if no variants, print chromosome id            
        }
        sim_print_mutref(ks->name.s, ks, rseq, rseq+1, output_fmt);
        fprintf(stdout, "Contig Variant End\n");

        bool bool_update_boundary = false;

        // load or create the methdb
        if(bool_methdb_set){
            // can check methdb_file's filename is the same as ks.name
            load_methdb(posidx_arr, meth_vec, methdb_file); 
        }else{
            create_methdb(ks, posidx_arr, meth_vec);
            if(bool_cgmap_set){fill_cgmap_chr(cgmap_file, ks->name.s, posidx_arr, meth_vec, bool_cgmap_pool, meth_seed);}
            if(bool_asm_set){fill_asm_chr(asm_file, ks->name.s, posidx_arr, meth_vec);}
            fill_beta(meth_vec, param_vec, meth_seed);
            update_methdb(posidx_arr, meth_vec, rseq, rseq+1, bool_asm_set, bool_update_boundary);  //update methdb with genetic variants
            if(bool_methdb_save){save_methdb(meth_vec, methdb_file);}
        }

        // initialize distributions to generate read positions
        uint64_t unif_begin = 2, unif_end = ks->seq.l-MAX_INSERT-2;         //ensure read doesn't pass boundary with 2 base offset
        std::vector<float> weights;

        if(tech_mode){
            parse_bed_chr(bed_file, ks->name.s, probe_vec);
            if(probe_vec.size() == 0){continue;} // parse probe, skip if empty
            if(is_uniform){unif_begin = 0; unif_end = probe_vec.size()-1;}else{
                std::vector<float> weights(probe_vec.size());
                for(int i=0; i < (int)probe_vec.size(); ++i){ weights.push_back(probe_vec[i].score);}
            }
        }
        std::uniform_int_distribution<int> dis_ud(unif_begin, unif_end);
        std::discrete_distribution<int> dis_dd(weights.begin(), weights.end()); //TODO: undefined behavior if weights empty
        // std::vector<float>().swap(weights); //  wights should remain valid for the lifetime of dis_dd

        frag_rec tmp_frag;
        while(ii < n_pairs){// the core loop
            // generate #chunk_size fragments records
            gen_frag_vec(&dis_ud, &dis_dd, posidx_arr, frag_vec, chunk_size, tech_mode, is_uniform);
            check_frag_vec(frag_vec, rseq, rseq+1, tech_mode); // check if the fragments are valid, fill start2
            
            // generate the read sequences
            int tmp_chunk = std::min(chunk_size, (int)(n_pairs - ii));
            ii += tmp_chunk;
            for(int idx=0; idx < tmp_chunk; ++idx){
                tmp_frag = frag_vec[idx];
                if(tmp_frag.haplo == -1){--ii; continue;} // skip for invalid fragments, should not happen
                target =  rseq[tmp_frag.haplo].s; //target = rseq[drand48()<0.5?0:1].s; // haplotype from which the reads are generated
                //cover_pos hold if the read covers a snp *position* (the read don't necessary contain the ALT allele)
                //j hold read1/read2, k hold the length of read, ix hold the cursor that transverses the read
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
                            /* context: 0x00 Match, 0x10 SNP, 0x30 INSERT       \
                                        0x01 CG, 0x03 CHG, 0x07 CHH (>>)        \
                                        0x09 GC, 0x0b GDC, 0x0f GDD (<<) */     \
                            tmp_seq[x][k] = c & 0xf;                            \
                            tmp_offset[x][k] = offset[x];                       \
                            if (mut_type == SUBSTITUTE) {                       \
                                tmp_context[x][k] = SNV;                        \
                                ++n_sub[x];                                     \
                            }                                                   \
                            ++end[x];                                           \
                            ++k;                                                \
                        } else {                                                \
                            tmp_seq[x][k] = c & 0xf; /*The base is ref*/        \
                            tmp_offset[x][k] = offset[x];                       \
                            ++n_indel[x];                                       \
                            ++end[x];                                           \
                            ++k;                                                \
                            int num_ins, ins;                                   \
                            for (num_ins = mut_type>>12, ins = c>>4; num_ins > 0 && k < size[x]; --num_ins, ins >>= 2){ \
                                --offset[x];                                    \
                                tmp_seq[x][k] = ins & 0x3;                      \
                                tmp_offset[x][k] = offset[x];                   \
                                tmp_context[x][k] = INSR;                       \
                                ++k;                                            \
                            }                                                   \
                        }                                                       \
                        if(bool_asm_set){cover_pos[x] |= posidx_arr[i] & 0x1;}  \
                    }                                                           \
                    if (k != size[x]) {ext_coor[x] = -10;}                      \
                } while (0)

                __gen_read(0, start[0], ++i);
                __gen_read(1, start[1], ++i);

                if (ext_coor[0] < 0 || ext_coor[1] < 0) { --ii; continue;}  // failed to generate the read(s)

                for(j = 0; j < 2; ++j){         //check the number of Ns
                    int n_n =0;
                    for (i = 0; i < size[j]; ++i) {
                        int c = tmp_seq[j][i];
                        if (c >= 4) {           // actually c should be never larger than 4 if everything is correct
                            ++n_n; 
                            tmp_seq[j][i] = 4;
                        }
                    }
                    if ((double)n_n / size[j] > MAX_N_RATIO) break;
                }
                if (j < 2) { --ii; continue; }  // too many ambiguous bases on one of the reads
                
                int flag_mut=0;
                for(i=tmp_frag.pos_l; i < tmp_frag.pos_r; ++i){
                    flag_mut |= target[i] & mutmsk;                         //whether frag has mutation
                }
                int flag_indel= n_indel[0] | n_indel[1];                    //whether read1/2 has indel

                // print reads to stdout: mode 0 print string (WGS), else print chars&numbers (WGBS)
                if(output_fmt == 0){
                    // flip and get the reverse complementary
                    int is_flip = drand48() < 0.5? 0 : 1;
                    int tmp_cigar, tmp_base;
                    for (k = 0; k < size[1]; ++k) { 
                        if (k <= int(size[1]/2)) { 
                            tmp_base  = tmp_seq[1][k];                              // reverse
                            tmp_seq[1][k] = tmp_seq[1][size[1]-k];
                            tmp_seq[1][size[1]-k] = tmp_base;

                            tmp_cigar = tmp_context[1][k];
                            tmp_context[1][k] = tmp_context[1][size[1]-k];
                            tmp_context[1][size[1]-k] = tmp_cigar;
                        }
                        tmp_seq[1][k] = tmp_seq[1][k] < 4? 3 - tmp_seq[1][k] : 4;   // complement
                    }
                    for (j = 0; j < 2; ++j) {
                        int jj = j ^ is_flip; // x^0=x; x^1=!x (when x is binary)
                        // header: 1-based coordinates for string output
                        fprintf(stdout, "@%s:%d:%d:%llx:%d:%d/%d\n", ks->name.s, start[0]+1, end[1]+1, (long long)ii, flag_mut, flag_indel, j+1);
                        // sequence (introduce random sequencing error)
                        for (i = 0; i < size[jj]; ++i) {
                            int c = tmp_seq[jj][i];
                            if (drand48() < ERR_RATE){
                                // c = (c + (int)(drand48() * 3.0 + 1)) & 3; // random sequencing errors
                                c = (c + 1) & 3; // recurrent sequencing errors
                                ++n_err[jj];
                                tmp_context[jj][i] |= SEQERR;
                                tmp_seq[jj][i] = c;
                            }
                            fputc("ACGTN"[c], stdout);
                        }
                        fprintf(stdout, "\n");
                        // comment
                        fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:%d:", start[jj]+1, end[jj]+1, cover_pos[jj], n_sub[jj], n_indel[jj], n_err[jj], end[1]-start[0], start[1]-end[0]);
                        for (i = 0; i < size[jj]; ++i) {
                            int c = (tmp_context[jj][i] & 0xf0) >> 4;
                            fputc("MXIE"[mut_table[c]], stdout);
                        }
                        fprintf(stdout, "\n");
                        // quality
                        for (i = 0; i < size[j]; ++i) {
                            if(i != size[j]-1){fputc(Q+33, stdout);}else{fputc(Q+32, stdout);}
                        }
                        fprintf(stdout, "\n");
                    }
                } else {
                    for (j = 0; j < 2; ++j) {
                        // header: ID is 1-based coordinates
                        fprintf(stdout, "@%s:%d:%d:%llx %d %d %d %d ", ks->name.s, start[0]+1, end[1]+1, (long long)ii, j, flag_mut, flag_indel, tmp_frag.strand); 
                        for (i = 0; i < size[j]; ++i) {
                            int pos =  start[j] + i + tmp_offset[j][i];
                            if(posidx_arr[pos] & 0x2){
                                int idx = posidx_arr[pos] >> 2;
                                tmp_context[j][i] |= meth_vec[idx].context;
                                fprintf(stdout, "%.4f,", meth_vec[idx].meth[(int)(bool)flag_mut]);
                            }
                        }
                        fprintf(stdout, "\n");
                        // sequence (no sequencing error, represented by 0-4)
                        for (i = 0; i < size[j]; ++i) {
                            //fprintf(stdout, "%d", tmp_seq[j][i]); //this will output the hex number, below line will output the ascii
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
        free(posidx_arr);
        tot_pairs += n_pairs;
    }

    fprintf(stderr, "[%s] Generated %lu read pairs, with %lu contain SNP, %lu contain INDEL", __func__, tot_pairs, tot_sub, tot_indel);
    if (output_fmt == 0){fprintf(stderr, " and %lu contain sequencing errors\n", tot_err);} else{fprintf(stderr, "\n");}
    kseq_destroy(ks);
    gzclose(fp_fa);
    free(tmp_seq[0]); free(tmp_seq[1]);
    free(tmp_offset[0]); free(tmp_offset[1]);
    free(tmp_context[0]); free(tmp_context[1]);
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
    fprintf(stderr, "         -O INT        output format: 0 for letters; nonzero for ascii numbers (for python module) [0]\n");
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
    char param_defult[] = "0.5|0.5,0.05|0.05,0.05|0.05";
    char *chr_id    = none_default;
    char *vcf_file  = none_default;
    char *bias_file = none_default;
    char *bed_file  = none_default; 
    char *asm_file  = none_default; 
    char *cgmap_file= none_default;
    char *methdb_file= none_default; // checked, will not intefere with vcf_file
    char *param_str = param_defult;


    //update default from command line
    int c = 0;
    while ((c = getopt(argc, argv, "I:J:K:L:c:C:n:N:d:1:2:E:A:O:v:R:F:X:H:s:a:m:p:W:M:S:P:T:u:B:b:D:")) >= 0) {
        switch (c) {
            case 'I': MEAN_INSERT= atoi(optarg); break;
            case 'J': SD_INSERT  = atoi(optarg); break;
            case 'K': MIN_INSERT = atoi(optarg); break;
            case 'L': MAX_INSERT = atoi(optarg); break;
            case 'c': chr_id     = optarg; break;
            case 'C': chunk_size = atoi(optarg); break;
            case 'n': chr_N      = atoi(optarg); break;
            case 'N': N          = atoi(optarg); break;
            case 'd': depth      = atoi(optarg); break;
            case '1': SIZE_L     = atoi(optarg); break;
            case '2': SIZE_R     = atoi(optarg); break;
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
            case 'b': bed_file   = optarg; break;
            case 'D': SD_CENTER  = atoi(optarg); break;
        }
    }
    if (argc - optind < 1) return simu_usage();
    if (seed_snp <= 0) seed_snp  = time(0)&0x7fffffff;
    if (seed_meth<= 0) seed_meth = time(0)&0x7fffffff;
    parse_param(param_str, param_vec);

    MIN_INSERT = std::max(std::max(SIZE_L, SIZE_R), MIN_INSERT); // ensure MIN_INSERT >= SIZE_L or SIZE_R

    bool bool_chr_set   = strcmp(chr_id,  "None") && strlen(chr_id);
    bool bool_vcf_set   = strcmp(vcf_file,"None") && strlen(vcf_file);
    // bool bool_site_set  = strcmp(cut_str, "None") && strlen(cut_str);
    bool bool_probe_set = strcmp(bed_file,"None") && strlen(bed_file);
    bool bool_bias_set  = strcmp(bias_file,"None") && strlen(bias_file);
    // bool bool_asm_set   = strcmp(asm_file, "None") && strlen(asm_file);
    // bool bool_cgmap_set = strcmp(cgmap_file, "None") && strlen(cgmap_file);
    // bool bool_methdb_set= strcmp(methdb_file, "None") && strlen(methdb_file);
    
    
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
        if(bool_bias_set){parse_bias_file(bias_file, eff_vec); is_uniform = false; BIN_SIZE = eff_vec.size();}
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
        if (l < MEAN_INSERT+3*SD_INSERT){
            fprintf(stderr, "[%s] skip contig '%s' as it is shorter than %d!\n", __func__, ks->name.s, MEAN_INSERT+3*SD_INSERT); 
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
        chr_N = chr_N == 0? (contig_eff_len * depth)/(SIZE_L + SIZE_R) : chr_N;
        fprintf(stderr, "[%s] Contig %s specified, total length: %lu, effective length: %lu, #reads: %lu\n", __func__, chr_id, contig_len, contig_eff_len, chr_N);
        chr_count[chr_id_str].count = chr_N;
    } else {
        // calculate the count for all contigs
        if (chr_N > 0) {fprintf(stderr, "ERROR: -n is specified but not -c. Exit... (please note the difference of -n and -N)\n"); exit(EXIT_FAILURE);}
        
        int num_contigs = (int)chr_count.size();
        N = N == 0? (tot_eff_len * depth)/(SIZE_L + SIZE_R) : N;
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
        if((vcf=fopen(vcf_file,"r"))){fprintf(stderr, "[%s] VCF file exists, use it to simulate reads\n", __func__); fclose(vcf);
        }else{fprintf(stderr, "ERROR: The specified VCF file does not exist, please check!\n"); exit(EXIT_FAILURE);}
    } else {
        fprintf(stderr, "[%s] No VCF input, will generate SNP randomly if mutation rate is nonzero\n", __func__);
    }


    fprintf(stderr, "[htsim] seed = %d\n", seed_snp);
    srand48(seed_snp);

    sim_core(argv[optind], is_hap, is_uniform, tech_mode, output_fmt, chunk_size, vcf_file, bed_file, chr_id, methdb_save, methdb_file, cgmap_file, cgmap_pool, asm_file, seed_meth);

    return 0;
}

