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
#include <chrono>

#include "/home/wbguo/iproject/BSReadSim/HTSIM/src/struct.h"

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

const int BUFFER_SIZE = 8192;

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
        buffer_pos += snprintf(output_buffer + buffer_pos, BUFFER_SIZE - buffer_pos, "@%s:%d:%d:%llx:%d/%d\n",
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


void output_read2(int* tmp_seq[2], int* tmp_context[2], int err_thre, int flag_mut, uint32_t ii, char* chr_id, int half_size, int Q,
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

    // Output
    for (j = 0; j < 2; ++j) {
        int jj = j ^ is_flip;

        // Header: 1-based coordinates in the readID
        fprintf(stdout, "@%s:%d:%d:%llx:%d/%d\n", chr_id, start[0]+1, end[1]+1, (long long)ii, flag_mut, j+1);

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
        fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:%d:", start[jj]+1+jj, end[jj]+jj, cover_pos[jj], 
                                                      n_sub[jj], n_indel[jj], n_err[jj], end[1]-start[0]+1, start[1]-end[0]+1);
        
        for (k = 0; k < size[jj]; ++k) {
            int c = (tmp_context[jj][k] & 0xf0) >> 4;
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

int main(){
    // strange error with int8_t or uint8_t, offset messed up
    int *tmp_seq[2];    	// sequence
    int *tmp_pos[2];    	// position 
    // int *tmp_offset[2]; 	// offset per base 4DEBUG
    int *tmp_context[2];	// cytosine context (CG/CHG/CHH) & mutation (upper half mutation, lower half sequence)
    int size[2], max_length, tmparr_size, err_thre;

    max_length = 100;
    tmparr_size= max_length*4;
    tmp_seq[0] = (int*)calloc(max_length, 4);           // sizeof(uint8_t)=1, sizeof(int)=4
    tmp_seq[1] = (int*)calloc(max_length, 4);
    tmp_pos[0] = (int*)calloc(max_length, 4);         
    tmp_pos[1] = (int*)calloc(max_length, 4);
    // tmp_offset[0]= (int*)calloc(max_length, 4);         // 4DEBUG
    // tmp_offset[1]= (int*)calloc(max_length, 4);
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
    int half_size = 10;
    int Q = 20;
    int start[2] = {0, 5};
    int end[2] = {9, 19};
    int n_sub[2] = {2, 3};
    int n_indel[2] = {1, 2};
    int n_err[2] = {3, 4};
    int cover_pos[2] = {8, 18};

    auto start_time = std::chrono::high_resolution_clock::now();
    for(int i=0; i < 1000000; ++i){
        output_read(tmp_seq, tmp_context, err_thre, flag_mut, ii, chr_id, half_size, Q, size, start, end, n_sub, n_indel, n_err, cover_pos);
    }
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    fprintf(stderr, "output_read: %ld us\n", duration.count());

    start_time = std::chrono::high_resolution_clock::now();
    for(int i=0; i < 1000000; ++i){
        output_read2(tmp_seq, tmp_context, err_thre, flag_mut, ii, chr_id, half_size, Q, size, start, end, n_sub, n_indel, n_err, cover_pos);
    }
    end_time = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end_time - start_time);
    fprintf(stderr, "output_read2: %ld us\n", duration.count());

    return 0;
}

