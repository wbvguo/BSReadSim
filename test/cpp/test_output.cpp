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


float err_rate = 0.005;



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
        fprintf(stdout, "@%s:%d:%d:%llx:%d:%d/%d\n", chr_id, start[0]+1, end[1]+1, (long long)ii, flag_mut, j+1);

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

}

