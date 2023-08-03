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


void output_read(int* tmp_seq[2], int* tmp_context[2], int half_size, int flip_thre, int err_thre, int Q, char* chr_id, uint32_t ii, int flag_mut, 
                 int(&size)[2], int(&start)[2], int(&end)[2], int(&n_sub)[2], int(&n_indel)[2], int(&n_err)[2], int(&cover_pos)[2])
{
    int j, k;
    int tmp_k, tmp_base, tmp_cigar;
    int is_flip = rand() > flip_thre;
    
    //output check
    for(j = 0; j <2; ++j){
        for (k = 0; k < size[j]; ++k) {fprintf(stdout, "%d,", tmp_seq[j][k]);}
        fprintf(stdout, "\n");
        for (k = 0; k < size[j]; ++k) {fprintf(stdout, "%d,", tmp_context[j][k]);}
        fprintf(stdout, "\n");
    }

    // Flip and get the reverse complementary sequence
    for (k = 0; k < size[1]; ++k) {
        if (k < half_size) {
            tmp_k = size[1] - k - 1;

            tmp_base = tmp_seq[1][k];                              
            tmp_seq[1][k] = tmp_seq[1][tmp_k];
            tmp_seq[1][tmp_k] = tmp_base;

            tmp_cigar= tmp_context[1][k];
            tmp_context[1][k] = tmp_context[1][tmp_k];
            tmp_context[1][tmp_k] = tmp_cigar;
        }
        tmp_seq[1][k] = tmp_seq[1][k] < 4 ? 3 - tmp_seq[1][k] : 4;
    }

    fprintf(stdout, "\n");

    //output check
    for(j = 0; j <2; ++j){
        for (k = 0; k < size[j]; ++k) {
            fprintf(stdout, "%d,", tmp_seq[j][k]);
        }
        fprintf(stdout, "\n");
        for (k = 0; k < size[j]; ++k) {
            fprintf(stdout, "%d,", tmp_context[j][k]);
        }
        fprintf(stdout, "\n");
    }

    // Output
    for (j = 0; j < 2; ++j) {
        int jj = j ^ is_flip;

        // Header: 1-based coordinates in the readID
        fprintf(stdout, "@%s:%d:%d:%llx:%d:%d/%d\n", chr_id, start[0]+1, end[1]+1, (long long)ii, is_flip, flag_mut, j+1);

        // Sequence (introduce random sequencing error)
        for (k = 0; k < size[jj]; ++k) {
            int c = tmp_seq[jj][k];
            if (rand()< err_thre) {
                c = (c + 1) & 3; // Recurrent sequencing errors
                ++n_err[jj];
                tmp_context[jj][k] |= SEQERR;
                tmp_seq[jj][k] = c;
            }
            fputc("ACGTN"[c], stdout);
        }
        fprintf(stdout, "\n");

        // Comment
        fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:%d:", 
                        start[jj]+1+jj, end[jj]+jj, cover_pos[jj], n_sub[jj], n_indel[jj], n_err[jj], end[1]-start[0]+1, start[1]-end[0]+1);
        
        for (k = 0; k < size[jj]; ++k) {
            int c = (tmp_context[jj][k] & 0xf0) >> 4; // TODO: rethink if we should change the order of the 4-bit flag
            fputc("MXIE"[mut_table[c]], stdout);
        }
        fprintf(stdout, "\n");

        // Quality
        for (k = 0; k < size[j]-1; ++k) {
            fputc(Q+33, stdout);
        }
        fputc(Q+32, stdout);
        fprintf(stdout, "\n");
    }
}

int main()
{
    int *tmp_seq[2];    	// sequence
    int *tmp_pos[2];    	// position 
    int *tmp_context[2];	// cytosine context (CG/CHG/CHH) & mutation (upper half mutation, lower half sequence)
    int size[2], max_length, tmparr_size, err_thre;

    max_length = 100;
    tmparr_size= max_length*4;
    tmp_seq[0] = (int*)calloc(max_length, 4);           // sizeof(uint8_t)=1, sizeof(int)=4
    tmp_seq[1] = (int*)calloc(max_length, 4);
    tmp_pos[0] = (int*)calloc(max_length, 4);         
    tmp_pos[1] = (int*)calloc(max_length, 4);
    tmp_context[0] = (int*)calloc(max_length, 4);
    tmp_context[1] = (int*)calloc(max_length, 4);
    size[0] = 100; size[1] = 100;
    err_thre  = (int) RAND_MAX * 0.001;

    for(int j=0; j < 2; ++j){
        for(int i=0; i < size[j]; ++i){
            tmp_seq[j][i] = rand() & 3;
            tmp_pos[j][i] = i;
            tmp_context[j][i] = 0;
        }
    }

    int flag_mut = 1;
    uint32_t ii = 0;
    char* chr_id = "chr1";
    
    int half_size = max_length / 2;
    int Q = 20;
    int start[2] = {1, 100};
    int end[2] = {51, 150};
    int n_sub[2] = {2, 3};
    int n_indel[2] = {1, 2};
    int n_err[2] = {3, 4};
    int cover_pos[2] = {8, 18};

    output_read(tmp_seq, tmp_context, half_size, 0, err_thre, Q, chr_id, ii, flag_mut, size, start, end, n_sub, n_indel, n_err, cover_pos);
}