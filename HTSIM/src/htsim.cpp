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

KSEQ_INIT(gzFile, gzread)
#include "struct.h"
#include "methdb.h"
#include "mode.h"
#include "haplo.h"
#include "rrcut.h"


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
std::vector<cut_rec> cut_vec;
std::vector<snp_rec> snp_vec;
std::vector<meth_rec> meth_vec;
std::vector<frag_rec> frag_vec;
std::vector<frag_rec> probe_vec;
std::vector<param_rec> param_vec;
std::map<std::string, chr_rec> chr_count;
std::map<int, snpmeth_rec> snpmeth_map;


void output_read(int* tmp_seq[2], int* tmp_context[2], int err_thre, int flag_mut, uint32_t ii, char* chr_id, int half_size, int Q,
                 int(&size)[2], int(&start)[2], int(&end)[2], int(&n_sub)[2], int(&n_indel)[2], int(&n_err)[2], int(&cover_pos)[2])
{
    int j, k;
    int tmp_k, tmp_base, tmp_cigar;
    int is_flip = rand() > (RAND_MAX>>1);
    // Flip and get the reverse complementary sequence
    for (k = 0; k < size[1]; ++k) {
        if (k <= half_size) {
            tmp_k = size[1] - k;

            tmp_base = tmp_seq[1][k];                              
            tmp_seq[1][k] = tmp_seq[1][tmp_k];
            tmp_seq[1][tmp_k] = tmp_base;

            tmp_cigar= tmp_context[1][k];
            tmp_context[1][k] = tmp_context[1][tmp_k];
            tmp_context[1][tmp_k] = tmp_cigar;
        }
        tmp_seq[1][k] = tmp_seq[1][k] < 4 ? 3 - tmp_seq[1][k] : 4;
    }

    // Buffer for storing the output
    char output_buffer[BUFFER_SIZE];
    int buffer_pos = 0;

    // Output
    for (j = 0; j < 2; ++j) {
        int jj = j ^ is_flip;

        // Header: 1-based coordinates in the readID
        buffer_pos += snprintf(output_buffer + buffer_pos, BUFFER_SIZE - buffer_pos, "@%s:%d:%d:%llx:%d:%d/%d\n",
                               chr_id, start[0]+1, end[1]+1, (long long)ii, flag_mut, j+1);
        //fprintf(stdout, "@%s:%d:%d:%llx:%d:%d/%d\n", chr_id, start[0]+1, end[1]+1, (long long)ii, flag_mut, j+1);

        // Sequence (introduce random sequencing error)
        for (k = 0; k < size[jj]; ++k) {
            int c = tmp_seq[jj][k];
            if (rand()< err_thre) {
                c = (c + 1) & 3; // Recurrent sequencing errors
                ++n_err[jj];
                tmp_context[jj][k] |= SEQERR;
                tmp_seq[jj][k] = c;
            }
            //fputc("ACGTN"[c], stdout);
            output_buffer[buffer_pos++] = "ACGTN"[c];
            if (buffer_pos >= BUFFER_SIZE - 1) {
                fwrite(output_buffer, 1, buffer_pos, stdout);
                buffer_pos = 0;
            }
        }
        output_buffer[buffer_pos++] = '\n';
        //fprintf(stdout, "\n");

        // Comment
        buffer_pos += snprintf(output_buffer + buffer_pos, BUFFER_SIZE - buffer_pos, "+:%d:%d:%d:%d:%d:%d:%d:%d:",
                               start[jj]+1+jj, end[jj]+jj, cover_pos[jj], n_sub[jj], n_indel[jj], n_err[jj], end[1]-start[0]+1, start[1]-end[0]+1);
        //fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:%d:", start[jj]+1+jj, end[jj]+jj, cover_pos[jj], n_sub[jj], n_indel[jj], n_err[jj], end[1]-start[0]+1, start[1]-end[0]+1);
        
        for (k = 0; k < size[jj]; ++k) {
            int c = (tmp_context[jj][k] & 0xf0) >> 4;
            output_buffer[buffer_pos++] = "MXIE"[mut_table[c]];
            //fputc("MXIE"[mut_table[c]], stdout);
        }
        output_buffer[buffer_pos++] = '\n';
        //fprintf(stdout, "\n");

        // Quality
        for (k = 0; k < size[j]-1; ++k) {
            output_buffer[buffer_pos++] = Q + 33; // TODO: check if we need to convert to char
            //fputc(Q+33, stdout);
            if (buffer_pos >= BUFFER_SIZE - 1) {
                fwrite(output_buffer, 1, buffer_pos, stdout);
                buffer_pos = 0;
            }
        }
        output_buffer[buffer_pos++] = Q + 32;
        output_buffer[buffer_pos++] = '\n';
        //fputc(Q+32, stdout);
        //fprintf(stdout, "\n");
    }

    // Write any remaining data in the buffer
    if (buffer_pos > 0) {
        fwrite(output_buffer, 1, buffer_pos, stdout);
    }
}


void output_read(int *tmp_seq[2], int *tmp_context[2], int *tmp_pos[2], int flag_mut, uint32_t ii, char *chr_id, int strand,
                 uint32_t *posidx_arr, uint32_t *kmeridx_arr,
                 int (&size)[2], int (&start)[2], int (&end)[2], int (&n_sub)[2], int (&n_indel)[2], int (&n_err)[2], int (&cover_pos)[2]){
    int j, k, idx;
    int pos, pos_prev;
    float kmer_meth;
    snpmeth_rec tmp_snpmeth;

    for (j = 0; j < 2; ++j) {
        // header: ID is 1-based coordinates
        fprintf(stderr, "@%s:%d:%d:%llx %d %d %d ", chr_id, start[0]+1, end[1]+1, (long long)ii, j, flag_mut, strand);
        
        // print meth+kmer
        for (k = 0; k < size[j]; ++k) {
            if(tmp_context[j][k]==0){continue;}
            pos = tmp_pos[j][k];                // pos= start[j] + k + tmp_offset[j][k];
            if(tmp_context[j][k]&0xf0){
                idx = pos == pos_prev? idx+1 : 0;
                kmer_meth= snpmeth_map[pos].kmeridx[idx] + snpmeth_map[pos].meth[idx];
                tmp_context[j][k] |= snpmeth_map[pos].context[idx];
            }else{                              // if it's a normal C/G
                if(kmeridx_arr == NULL){kmer_meth = 0;}else{kmer_meth= flag_mut? kmeridx_arr[pos] >> 16 : (posidx_arr[pos] & 0xffff);}
                idx = posidx_arr[pos] >> 2;
                kmer_meth += meth_vec[idx].meth[flag_mut];
                tmp_context[j][k] |= meth_vec[idx].context[flag_mut];
            }
            fprintf(stdout, "%.4f,", kmer_meth);
        }
        fprintf(stdout, "\n");
        // sequence (no sequencing error, represented by 0-4)
        for (k = 0; k < size[j]; ++k) {
            //fprintf(stdout, "%d", tmp_seq[j][i]); //this will output the hex number, below line will output the ascii
            fputc(tmp_seq[j][k], stdout);
        }
        fprintf(stdout, "\n");
        // comment
        fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:", start[j], end[j], cover_pos[j], n_sub[j], n_indel[j], end[1]-start[0]+1, start[1]-end[0]+1);
        for (k = 0; k < size[j]; ++k) {
            fprintf(stdout, "%d,",tmp_pos[j][k]-start[j]-k);
        }
        // const char *pad = "";
        // for (i = 0; i < size[j]; ++i) {
        //     fprintf(stdout, "%s%d", pad, tmp_offset[j][i]);
        //     pad = ",";
        // }
        fprintf(stdout, "\n");
        // quality
        for (k = 0; k < size[j]; ++k) {
            //fprintf(stdout, "%x", tmp_context[j][i]);
            fputc(tmp_context[j][k], stdout);
        }
        fprintf(stdout, "\n");
    }
}


void sim_core(const char *fn, char *vcf_file, char *bed_file, char *chr_id, char *methdb_file, char *cgmap_file, char *asm_file, 
                expt_param *expt_set, mut_param *mut_set, meth_param *meth_set)
{
    gzFile   fp_fa;
    kseq_t   *ks;
    mut_t    *target;
    mutseq_t rseq[2];
    uint64_t tot_sub = 0, tot_indel = 0, tot_err = 0, tot_pairs = 0;
    uint32_t n_pairs =0;        // max 4,294,967,295 reads per contig
    
    // strange error with int8_t or uint8_t, offset messed up
    int *tmp_seq[2];    	// sequence
    int *tmp_pos[2];    	// position 
    // int *tmp_offset[2]; 	// offset per base 4DEBUG
    int *tmp_context[2];	// cytosine context (CG/CHG/CHH) & mutation (upper half mutation, lower half sequence)
    int size[2], Q, max_length, tmparr_size, half_size, err_thre;

    max_length = std::max(expt_set->size_l, expt_set->size_r);
    tmparr_size= max_length*4;
    tmp_seq[0] = (int*)calloc(max_length, 4);           // sizeof(uint8_t)=1, sizeof(int)=4
    tmp_seq[1] = (int*)calloc(max_length, 4);
    tmp_pos[0] = (int*)calloc(max_length, 4);         
    tmp_pos[1] = (int*)calloc(max_length, 4);
    // tmp_offset[0]= (int*)calloc(max_length, 4);         // 4DEBUG
    // tmp_offset[1]= (int*)calloc(max_length, 4);
    tmp_context[0] = (int*)calloc(max_length, 4);
    tmp_context[1] = (int*)calloc(max_length, 4);
    size[0] = expt_set->size_l; size[1] = expt_set->size_r;
    half_size = size[1] >> 1;
    err_thre  = (int) RAND_MAX * expt_set->err_rate;

    Q = (expt_set->err_rate == 0.0)? 'I' : (int)(-10.0 * log(expt_set->err_rate) / log(10.0) + 0.499) + 33;

    // start simulate
    fp_fa = gzopen(fn, "r");
    ks = kseq_init(fp_fa);
    int l = 0;
    while ((l = kseq_read(ks)) >= 0) {  //here l is the chromosome length
        if (expt_set->is_chr_set) {if (strcmp(chr_id, ks->name.s)!=0){continue;}}
        if (l < expt_set->mean_insert + 3 * expt_set->sd_insert) {continue;}

        n_pairs = chr_count[std::string(ks->name.s)].count;
        if(!n_pairs){fprintf(stderr, "[%s] skip contig '%s' for its #reads is 0...\n", __func__, ks->name.s); continue;}
        fprintf(stderr, "[%s] contig '%s': length %ld, simulate %u reads...\n", __func__, ks->name.s, ks->seq.l, n_pairs);

        // Initialize pos_idx array as 0, last 2 bits: whether it's C/G, if it's a SNP position (base can be either REF/ALT)
        // uint8_t posidx_arr[ks->seq.l] = {0}; // only work when length is small, otherwise stack overflow
        uint32_t* posidx_arr = (uint32_t*) calloc(ks->seq.l, sizeof(uint32_t));
        if (posidx_arr == NULL) { fprintf(stderr, "ERROR: could not allocate memory\n");exit(EXIT_FAILURE);}

        // introduce mutations and print them to stdout
        fprintf(stdout, "Contig Variant Start\n");
        if(mut_set->is_vcf_set){
            sim_mut_vcf(ks, vcf_file, rseq, rseq+1, posidx_arr, snp_vec);
            if(snp_vec.size() == 0){fprintf(stdout, "%s\n", ks->name.s);}       //if no variants, print chromosome id
        } else {
            sim_mut_diref(ks, mut_set, rseq, rseq+1, posidx_arr);
            if(mut_set->mut_rate == 0.0){fprintf(stdout, "%s\n", ks->name.s);}  //if no variants, print chromosome id
        }
        sim_print_mutref(ks->name.s, ks, rseq, rseq+1, expt_set->output_fmt);
        fprintf(stdout, "Contig Variant End\n");


        uint32_t* kmeridx_arr = NULL;
        if(expt_set->is_kmer_set){
            kmeridx_arr = (uint32_t*) calloc(ks->seq.l, sizeof(uint32_t));
            if (kmeridx_arr == NULL) { fprintf(stderr, "ERROR: could not allocate memory for kmer index\n");exit(EXIT_FAILURE);}
            for (int i = 3; i < ks->seq.l-3; ++i){
                if(posidx_arr[i]&0x2){continue;}//skip non-CG positions
                int tmp_kmeridx =0;
                for (int j = -3; j < 4; ++j){
                    int ref_base = nst_nt4_table[(uint8_t)ks->seq.s[i]];
                    ref_base = ref_base < 4 ? ref_base : (rand()&0x3);          // will not interfere with the drand48()
                    tmp_kmeridx = (tmp_kmeridx << 2) | ref_base;
                }
                kmeridx_arr[i] = tmp_kmeridx <<16 | tmp_kmeridx;
            }
        }


        // load or create the methdb
        if(meth_set->is_methdb_set){
            // can check methdb_file's filename is the same as ks.name
            fprintf(stderr, "[%s] contig '%s': load methdb...\n", __func__, ks->name.s);
            load_methdb(posidx_arr, meth_vec, methdb_file);
        }else{
            fprintf(stderr, "[%s] contig '%s': create methdb...\n", __func__, ks->name.s);
            create_methdb(ks, posidx_arr, meth_vec);
            if(meth_set->is_cgmap_set){fill_cgmap_chr(cgmap_file, ks->name.s, posidx_arr, meth_vec, meth_set);}
            if(meth_set->is_asm_set){fill_asm_chr(asm_file, ks->name.s, posidx_arr, meth_vec);}
            fill_beta(meth_vec, param_vec, meth_set->seed_meth);
            //update methdb with genetic variants
            update_variant(ks, rseq, rseq+1, posidx_arr, meth_vec, kmeridx_arr, meth_set, param_vec, snpmeth_map);
            if(meth_set->methdb_save){save_methdb(meth_vec, methdb_file);}
        }

        // TODELETE: print methdb to stderr
        for (const auto& pair : snpmeth_map) {
            fprintf(stderr, "%d:", pair.first);
            for (const auto& element : pair.second.meth) {
                fprintf(stderr, "%f,", element);
            }
            fprintf(stderr, ":");
            for (const auto& element : pair.second.context) {
                fprintf(stderr, "%d,", element);
            }
            fprintf(stderr, ":");
            for (const auto& element : pair.second.kmeridx) {
                fprintf(stderr, "%d,", element);
            }
            fprintf(stderr, "\n");
        }

        // initialize distributions to generate read positions
        uint32_t unif_begin = 2, unif_end = ks->seq.l-expt_set->max_insert-2;       // ensure read doesn't pass boundary with 2 base offset, TOCHECK: if needed
        std::vector<float> weights;

        if(expt_set->tech_mode){
            parse_bed_chr(bed_file, ks->name.s, probe_vec, expt_set->tech_mode);    // parse probe, fill weights
            if(probe_vec.size() == 0){
                fprintf(stderr, "[%s] contig '%s': found 0 region in BED file (tech mode %d). exit ...\n", __func__, ks->name.s, expt_set->tech_mode); 
                exit(EXIT_FAILURE);
            }
            if(expt_set->is_uniform){unif_begin = 0; unif_end = probe_vec.size()-1;}else{
                for(int i=0; i < (int)probe_vec.size(); ++i){weights.push_back(probe_vec[i].score);}
            }
        }
        std::uniform_int_distribution<int> dis_ud(unif_begin, unif_end);
        std::discrete_distribution<int> dis_dd(weights.begin(), weights.end());     //TODO: undefined behavior if weights empty
        // std::vector<float>().swap(weights); //  wights should remain valid for the lifetime of dis_dd


        uint32_t ii = 0;    // record #reads that has been generated for this contig
        frag_rec tmp_frag;  // hold the tmp fragment
        int chunk_size = std::min((uint32_t)std::max(1000, expt_set->chunk_size), n_pairs);
        while(ii < n_pairs){// the core loop
            // generate #chunk_size sorted fragments records (at least 1000)
            gen_frag_vec(&dis_ud, &dis_dd, posidx_arr, (int) ks->seq.l, chunk_size, frag_vec, probe_vec, eff_vec, expt_set);
            // ii += chunk_size; // 4DEBUG
            
            // generate the read sequences
            for(int idx=0; idx < chunk_size && ii < n_pairs; ++idx){
                ++ii;
                tmp_frag= frag_vec[idx];
                target  =  rseq[tmp_frag.haplo].s; // haplotype from which the reads are generated

                //cover_pos hold if the read covers a snp *position* (the read don't necessary contain the ALT allele)
                //i holds cursor on genome, j holds read1/read2, k holds the length of read, ix holds the cursor that transverses the read
                int i, j, k, ix, num_ins, ins, shift_pos;
                int start[2] = {tmp_frag.pos_l, tmp_frag.pos_l};
                int end[2]   = {tmp_frag.pos_r, tmp_frag.pos_r};
                int offset[2]= {0, 0}, ext_coor[2] = {-10, -10}; // the coordinate of the first base of the read
                int n_sub[2] = {0,0}, n_indel[2] = {0,0}, n_err[2] = {0,0}, cover_pos[2] = {0,0};

                // reset
                memset(tmp_seq[0], 0, tmparr_size);
                memset(tmp_seq[1], 0, tmparr_size);
                memset(tmp_pos[0], 0, tmparr_size);
                memset(tmp_pos[1], 0, tmparr_size);
                // memset(tmp_offset[0], 0, tmparr_size);
                // memset(tmp_offset[1], 0, tmparr_size);
                memset(tmp_context[0], 0, tmparr_size);
                memset(tmp_context[1], 0, tmparr_size);

                // x: select read1 or read2; ext_coor: extend corrdinates;
                #define __gen_read(x, start_pos, iter) do {                     \
                    /* generate reads assign mutation flag; */                  \
                    for (i = (start_pos), k = 0, ext_coor[x] = -10; i >= 0 && i < ks->seq.l && k < size[x]; iter) { \
                        int c = target[i], mut_type = c & mutmsk;               \
                        if (ext_coor[x] < 0) {                                  \
                            /* avoid indel as the first base */                 \
                            if ((mut_type!=NOCHANGE) && (mut_type!=SUBSTITUTE)) continue; \
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
                            tmp_pos[x][k] = i;                                  \
                            /*tmp_offset[x][k] = offset[x];*/                   \
                            if (mut_type == SUBSTITUTE) {                       \
                                tmp_context[x][k] = SNV;                        \
                                ++n_sub[x];                                     \
                            }                                                   \
                            ++end[x];                                           \
                            ++k;                                                \
                        } else {                                                \
                            tmp_seq[x][k] = c & 0xf; /*The base is ref*/        \
                            tmp_pos[x][k] = i;                                  \
                            /*tmp_offset[x][k] = offset[x];*/                   \
                            ++n_indel[x];                                       \
                            ++end[x];                                           \
                            ++k;                                                \
                            for (num_ins = mut_type>>12, ins = c>>4; num_ins > 0 && k < size[x]; --num_ins, ins >>= 2){ \
                                --offset[x];                                    \
                                tmp_seq[x][k] = ins & 0x3;                      \
                                tmp_pos[x][k] = i;                              \
                                /*tmp_offset[x][k] = offset[x];*/               \
                                tmp_context[x][k] = INSR;                       \
                                ++k;                                            \
                            }                                                   \
                        }                                                       \
                        if(meth_set->is_asm_set){cover_pos[x] |= posidx_arr[i] & 0x1;}  \
                    }                                                           \
                    if (k != size[x]) {ext_coor[x] = -10;}                      \
                } while (0)

                #define __gen_read2(x, start_pos, iter) do {                    \
                    /* generate reads assign mutation flag; */                  \
                    for (i = (start_pos), k = size[x]-1, ext_coor[x] = -10; i >= 0 && i < ks->seq.l && k >= 0; iter) { \
                        int c = target[i], mut_type = c & mutmsk;               \
                        if (ext_coor[x] < 0) {                                  \
                            /* avoid indel as the first base */                 \
                            if ((mut_type!=NOCHANGE) && (mut_type!=SUBSTITUTE)) continue; \
                            start[x] = i;                                       \
                            end[x] = i;                                         \
                            ext_coor[x] = i;                                    \
                        }                                                       \
                        if (mut_type == DELETE){                                \
                            --offset[x];                                        \
                            --start[x];                                         \
                            ++n_indel[x];                                       \
                        }                                                       \
                        else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
                            /* context: 0x00 Match, 0x10 SNP, 0x30 INSERT       \
                                        0x01 CG, 0x03 CHG, 0x07 CHH (>>)        \
                                        0x09 GC, 0x0b GDC, 0x0f GDD (<<) */     \
                            tmp_seq[x][k] = c & 0xf;                            \
                            tmp_pos[x][k] = i;                                  \
                            /*tmp_offset[x][k] = offset[x];*/                   \
                            if (mut_type == SUBSTITUTE) {                       \
                                tmp_context[x][k] = SNV;                        \
                                ++n_sub[x];                                     \
                            }                                                   \
                            --start[x];                                         \
                            --k;                                                \
                        } else {                                                \
                            tmp_seq[x][k] = c & 0xf; /*The base is ref*/        \
                            tmp_pos[x][k] = i;                                  \
                            /*tmp_offset[x][k] = offset[x];*/                   \
                            ++n_indel[x];                                       \
                            --start[x];                                         \
                            --k;                                                \
                            for (num_ins = mut_type>>12, ins = c>>4; num_ins > 0 && k >=0; --num_ins, ins >>= 2){ \
                                ++offset[x];                                    \
                                tmp_seq[x][k] = ins & 0x3;                      \
                                tmp_pos[x][k] = i;                              \
                                /*tmp_offset[x][k] = offset[x];*/               \
                                tmp_context[x][k] = INSR;                       \
                                --k;                                            \
                            }                                                   \
                        }                                                       \
                        if(meth_set->is_asm_set){cover_pos[x] |= posidx_arr[i] & 0x1;}  \
                    }                                                           \
                    if (k != -1) {ext_coor[x] = -10;}                           \
                    /*shift_pos = end[x]-start[x]-size[x]+1;                    \
                    if(shift_pos != 0){                                         \
                        for(k=0; k < size[x]; ++k){                             \
                            tmp_offset[x][k] = tmp_offset[x][k] + shift_pos;    \
                        }                                                       \
                    }*/                                                         \
                } while (0)

                __gen_read(0, start[0], ++i);
                __gen_read2(1,  end[1], --i);

                if (ext_coor[0] < 0 || ext_coor[1] < 0) { --ii; continue;}  // failed to generate the read(s)
                
                for(j = 0; j < 2; ++j){         //check the number of Ns
                    int n_n =0;
                    // actually should be never larger than 4 if everything is correct
                    for(k = 0; k < size[j]; ++k){if(tmp_seq[j][k] >= 4){++n_n; tmp_seq[j][k]=4;}}
                    if ((double)n_n / size[j] > expt_set->maxN_ratio) break;
                }
                if (j < 2) { --ii; continue; }  // too many ambiguous bases on one of the reads
                
                int flag_mut=0, flag_indel=0;
                // if judgement, might not always need it
                for(i=start[0]; i <= end[1]; ++i){
                    flag_mut |= target[i] & mutmsk;                         //whether frag has mutation
                }
                flag_mut  = (int)(flag_mut != 0);
                flag_indel= (int)((n_indel[0]+n_indel[1])!=0);              //whether read1/2 has indel

                // print reads to stdout: mode 0 print string (WGS), else print chars&numbers (WGBS)
                if(expt_set->output_fmt == 0){

                } else {
                    // if(!flag_indel){continue;} // for bug testing
                }

                tot_sub   += (int)(n_sub[0]  + n_sub[1] > 0);
                tot_indel += (int)(n_indel[0]+ n_indel[1] > 0);
                tot_err   += (int)(n_err[0]  + n_err[1] > 0);
            }
        }
        free(rseq[0].s); free(rseq[1].s);
        free(posidx_arr);
        if(expt_set->is_kmer_set){free(kmeridx_arr);}
        tot_pairs += n_pairs;
    }

    fprintf(stderr, "[%s] Generated %lu read pairs, with %lu contain SNP, %lu contain INDEL", __func__, tot_pairs, tot_sub, tot_indel);
    if (expt_set->output_fmt == 0){fprintf(stderr, " and %lu contain sequencing errors\n", tot_err);} else{fprintf(stderr, "\n");}
    kseq_destroy(ks);
    gzclose(fp_fa);
    free(tmp_seq[0]); free(tmp_seq[1]);
    free(tmp_pos[0]); free(tmp_pos[1]);
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
    fprintf(stderr, "         -I INT        mean insert size (outer distance between 2 ends) [%d]\n", MEAN_INSERT);
    fprintf(stderr, "         -J INT        standard deviation of insert size [%d]\n", SD_INSERT);
    fprintf(stderr, "         -K INT        minimum insert size [%d]\n", MIN_INSERT);
    fprintf(stderr, "         -L INT        maximum insert size [%d]\n", MAX_INSERT);
    fprintf(stderr, "         -c STRING     contig name, only output reads from this contig, default is to output all contigs [None]\n");
    fprintf(stderr, "         -C INT        chunk size for read generation [%d]\n", CHUNK_SIZE);
    fprintf(stderr, "         -n INT        number of read pairs to generate for specified contig, disabled by default [0]\n");
    fprintf(stderr, "         -N INT        total number of read pairs to generate, will caculate from depth when N is 0 [0]\n");
    fprintf(stderr, "         -d INT        average sequencing depth, only used when n/N is not specified [%d]\n", DEPTH);
    fprintf(stderr, "         -1 INT        length of the first read [%d]\n", SIZE_L);
    fprintf(stderr, "         -2 INT        length of the second read [%d]\n", SIZE_R);
    fprintf(stderr, "         -E FLOAT      base error rate (set to be 0 for bisulfite sequencing) [%.3f]\n", ERR_RATE);
    fprintf(stderr, "         -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [%.2f]\n", MAX_N_RATIO);
    fprintf(stderr, "         -O INT        output format: 0 for letters; nonzero for ascii numbers (for python module) [%d]\n", OUTPUT_FMT);
    fprintf(stderr, "         -k INT        whether to output kmer index: 0 for no; nonzero for yes [0]\n");
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
    fprintf(stderr, "         -P STRING     Parameter string for beta distribution [0.5_0.5,0.05_0.05,0.05_0.05]\n");
    fprintf(stderr, "         -U INT        update methdb boundary sites for variants: 0 for no, non-zero for yes [0]\n");
    fprintf(stderr, "technology setting:\n");
    fprintf(stderr, "         -T INT        technology: 0 for Whole genome; 1 for Reduced representation; 2 for Targeted [%d]\n", TECH_MODE);
    fprintf(stderr, "         -u INT        uniform coverage: 0 for diable, nonzero for enable [1]\n");
    fprintf(stderr, "         -B STRING     GC Bias reference for WGS/WGBS, only used when -u set to be 0 [None]\n");
    fprintf(stderr, "         -b STRING     BED file for reduced-representation / targeted sequencing (.bed/.bed.gz) [None]\n");
    fprintf(stderr, "         -D INT        fragment center's deviaiton from the probe center [%d]\n", SD_CENTER);
    fprintf(stderr, "         -x STRING     enzyme cutting site string for reduced representation sequencing, take the format as C_CGG,_AATT [None]\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[])
{
    // default parameters
    uint64_t N  = 0;
    uint32_t chr_N  = 0;
    char none_default[] = "None";
    char param_default[]= "0.5_0.5,0.05_0.05,0.05_0.05";
    char *chr_id    = none_default;
    char *vcf_file  = none_default;
    char *bias_file = none_default;
    char *cut_str   = none_default;
    char *bed_file  = none_default; 
    char *asm_file  = none_default; 
    char *cgmap_file= none_default;
    char *param_str = param_default;
    char *methdb_file= none_default;

    expt_param expt_set;
    meth_param meth_set;
    mut_param  mut_set;

    //update parameters from command line
    int c = 0;
    while ((c = getopt(argc, argv, "I:J:K:L:c:C:n:N:d:1:2:E:A:O:v:R:F:X:H:s:a:m:p:W:M:S:P:T:u:B:b:D:x:")) >= 0) {
        switch (c) {
            case 'I': expt_set.mean_insert= atoi(optarg); break;
            case 'J': expt_set.sd_insert  = atoi(optarg); break;
            case 'K': expt_set.min_insert = atoi(optarg); break;
            case 'L': expt_set.max_insert = atoi(optarg); break;
            case 'd': expt_set.depth      = atof(optarg); break;
            case '1': expt_set.size_l     = atoi(optarg); break;
            case '2': expt_set.size_r     = atoi(optarg); break;
            case 'E': expt_set.err_rate   = atof(optarg); break;
            case 'C': expt_set.chunk_size = atoi(optarg); break;
            case 'A': expt_set.maxN_ratio = atof(optarg); break;
            case 'D': expt_set.sd_center  = atoi(optarg); break;
            case 'T': expt_set.tech_mode  = atoi(optarg); break;
            case 'u': expt_set.is_uniform = atoi(optarg)!=0; break;
            case 'k': expt_set.is_kmer_set= atoi(optarg)!=0; break;
            case 'O': expt_set.output_fmt = atoi(optarg); break;

            case 'R': mut_set.mut_rate    = atof(optarg); break;
            case 'F': mut_set.indel_frac  = atof(optarg); break;
            case 'X': mut_set.indel_extn  = atof(optarg); break;
            case 'H': mut_set.is_hap      = atoi(optarg)!=0; break;
            case 's': mut_set.seed_snp    = atoi(optarg); break;

            case 'S': meth_set.seed_meth  = atoi(optarg); break;
            case 'W': meth_set.methdb_save= atoi(optarg)!=0; break;
            case 'p': meth_set.cgmap_pool = atoi(optarg)!=0; break;
            case 'U': meth_set.update_meth= atoi(optarg)!=0; break;

            case 'n': chr_N      = atoi(optarg); break;
            case 'N': N          = atoi(optarg); break;
            case 'c': chr_id     = optarg; break;
            case 'v': vcf_file   = optarg; break;
            case 'a': asm_file   = optarg; break;
            case 'm': cgmap_file = optarg; break;
            case 'M': methdb_file= optarg; break;
            case 'P': param_str  = optarg; break;
            case 'B': bias_file  = optarg; break;
            case 'b': bed_file   = optarg; break;
            case 'x': cut_str    = optarg; break;
        }
    }
    if (argc - optind < 1) return simu_usage();
    if (mut_set.seed_snp <= 0) mut_set.seed_snp    = time(0)&0x7fffffff;
    if (meth_set.seed_meth<= 0) meth_set.seed_meth = time(0)&0x7fffffff;
    expt_set.min_insert = std::max(std::max(expt_set.size_l, expt_set.size_r), expt_set.min_insert); // ensure MIN_INSERT >= SIZE_L or SIZE_R

    parse_param(param_str, param_vec);
    expt_set.is_chr_set   = strcmp(chr_id, "None") && strlen(chr_id);
    expt_set.is_bias_set  = strcmp(bias_file,"None")&& strlen(bias_file);
    expt_set.is_bed_set   = strcmp(bed_file,"None") && strlen(bed_file);
    expt_set.is_site_set  = strcmp(cut_str, "None") && strlen(cut_str);
    mut_set.is_vcf_set    = strcmp(vcf_file,"None") && strlen(vcf_file);
    meth_set.is_asm_set   = strcmp(asm_file, "None") && strlen(asm_file);
    meth_set.is_cgmap_set = strcmp(cgmap_file, "None")&& strlen(cgmap_file);
    meth_set.is_methdb_set= strcmp(methdb_file, "None")&& strlen(methdb_file);

    if(!expt_set.is_bed_set && expt_set.tech_mode){
        fprintf(stderr, "ERROR: Please specify bed file for your choosen tech mode\n");exit(EXIT_FAILURE);
    }

    if (expt_set.is_bed_set && (expt_set.tech_mode==2 || !expt_set.is_site_set)){
        fprintf(stderr, "Simulating targeted sequencing reads:\n");
        expt_set.tech_mode = 2;
    } else if (expt_set.is_bed_set && (expt_set.tech_mode==1 || expt_set.is_site_set)){
        fprintf(stderr, "Simulating restricted enzyme cutting reads:\n");
        expt_set.tech_mode = 1;
    } else {
        fprintf(stderr, "Simulating whole genome reads:\n");
        if(!expt_set.is_uniform){
            if(expt_set.is_bias_set){
                parse_bias_file(bias_file, eff_vec); expt_set.is_uniform = false; expt_set.bin_size = eff_vec.size();
            }else{
                fprintf(stderr, "ERROR: Please specify GC-Bias file when specifying -u as 0\n");exit(EXIT_FAILURE);
            }
        }
        expt_set.tech_mode = 0;
    }

    // prepare to simulate reads
    fprintf(stderr, "Reference genome file: %s\n", argv[optind]);
    
    // check input vcf file
    if (mut_set.is_vcf_set) {
        FILE *vcf;
        if((vcf=fopen(vcf_file,"r"))){fprintf(stderr, "VCF file: %s, use it to create haplotypes\n", vcf_file); fclose(vcf);
        }else{fprintf(stderr, "ERROR: The specified VCF file does not exist, please check!\n"); exit(EXIT_FAILURE);}
    } else {
        fprintf(stderr, "[%s] No VCF input, will generate SNP randomly if mutation rate is nonzero\n", __func__);
    }
    if (expt_set.is_bed_set) {fprintf(stderr, "BED file: %s\n", bed_file);}

    // check existence of fasta, parse the length and calculate the count for each chr_id
    cal_chr_count(argv[optind], chr_id, bed_file, N, chr_N, &expt_set, probe_vec, chr_count);


    fprintf(stderr, "[htsim] snp seed = %d, meth seed = %d\n", mut_set.seed_snp, meth_set.seed_meth);
    srand48(mut_set.seed_snp);
    srand(time(0)&0x7fffffff);

    sim_core(argv[optind], vcf_file, bed_file, chr_id, methdb_file, cgmap_file, asm_file, &expt_set, &mut_set, &meth_set);

    return 0;
}

