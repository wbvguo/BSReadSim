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

/* This program is based on WGSIM(v0.3.1-r13)[https://github.com/lh3/wgsim.git]
 * with modifications to simulate WGS or WGBS reads in BSReadSim for diploid organism */

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
#include "kseq.h"
#include "vcf.h"
KSEQ_INIT(gzFile, gzread)

#define PACKAGE_VERSION "1.0.2"

const uint8_t nst_nt4_table[256] = {
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

static uint8_t MATCH  = 0x00;
static uint8_t SNV    = 0x01;
static uint8_t INSR   = 0x03;
static uint8_t CONVT  = 0x05;
static uint8_t SEQERR = 0x09;
const uint8_t mut_table[16] = {
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
const uint8_t cg_context_table[64] = {
    0, 0, 0, 0,  0, 0, 0, 0, 
    0, 0, 0, 0,  0, 0, 0, 0, 
    CHH, CHH, CHG, CHH,  CHH, CHH, CHG, CHH, 
    CG,  CG,  CG,  CG,	 CHH, CHH, CHG, CHH, 
    GDD, GDC, GDD, GDD,  GC,  GC,  GC,  GC, 
    GDD, GDC, GDD, GDD,  GDD, GDC, GDD, GDD,
    0, 0, 0, 0,  0, 0, 0, 0, 
    0, 0, 0, 0,  0, 0, 0, 0
};

const uint8_t cg_table[4] = {1, 0, 0, 1}; // for cg check

/* wgsim */
// if the leftmost 4 bit is non-zero, then it must be snp or indel
enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
static mut_t mutmsk = (mut_t)0xf000;

typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
} mutseq_t;

static double ERR_RATE = 0.005;
static double MUT_RATE = 0.01;
static double INDEL_FRAC = 0.15;
static double INDEL_EXTEND= 0.3;
static double MAX_N_RATIO = 0.05;

// to store vcf information
std::vector<int> pos_vec;
std::vector<int> ref_vec;
std::vector<int> alt_vec;
std::vector<int> geno_vec;

void parse_vcf_chr(char *fname, char *chr_id)
{
    //clean the container
    pos_vec.clear();
    ref_vec.clear();
    alt_vec.clear();
    geno_vec.clear();

    //open vcf file
    htsFile *fp    = hts_open(fname,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec    = bcf_init();
    
    //collect control
    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
    int prev_pos= -1;

    bool collect_present = false;
    bool collect_previous= false;

    //check the vcf file
    int nsmpl = bcf_hdr_nsamples(hdr); // number of sample
    if (nsmpl != 1) {
        fprintf(stderr, "ERROR: [%s] Currently only support single-sample simulation, please check the vcf file! Exiting...\n", __func__);
        exit(1);
    }

    while (bcf_read(fp, hdr, rec)>=0) {	
        if (collect_present == false && collect_previous == true) { break; } // finished collecting
        
        // a new round, save last status
        collect_previous = collect_present; 

        //unpack for read REF,ALT,INFO,etc 
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);

        //read CHROM
        std::string chr_current = bcf_hdr_id2name(hdr, rec->rid);

        if (strcmp(chr_current.c_str(), chr_id) != 0){
            collect_present = false;
            continue;
        } else {
            collect_present = true;
            std::string ref = rec->d.allele[0];
            std::string alt = rec->d.allele[1];
            int snp_pos = rec->pos; //it's 0-based coordinates
            int base_change_pos = snp_pos;

            int ref_len = ref.length();
            int alt_len = alt.length();
            int base_offset = alt_len - ref_len;

            // check genotype
            ngt = bcf_get_genotypes(hdr,rec,&gt,&ngt_arr); //The total number of array elements in &gt
            int snp_hap1 = bcf_gt_allele(gt[0]);
            int snp_hap2 = bcf_gt_allele(gt[1]);
            int is_phased = bcf_gt_is_phased(gt[1]); //phased:1, unphased:0

            //skip the following records: 
            // 1. insert/delete offset greater than 4
            // 2. multi-allelic sites
            // 3. multi-ploid organism
            // 4. multi nucleotide polymorphism (MNP); (please represent it with single base SNPs)
            // 5. SNP with the same position (only keep the first occurence)
            // TODO: 
            // 1. contains Ns;
            // 2. missing snps?
            if (abs(base_offset) > 4 || ngt > 2 || prev_pos==snp_pos) {
                fprintf(stderr, "[%s] Skip unusual SNP/INDEL [length/ploidy/pos]: CHROM:%s; POS:%d; REF:%s; ALT:%s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
                continue;
            }
            if (snp_hap1 < 0 || snp_hap2 < 0 || snp_hap1 > 1 || snp_hap2 > 1) {
                fprintf(stderr, "[%s] Skip unusual SNP/INDEL [snp haplotype]: CHROM:%s; POS:%d; HAP1:%d; HAP2:%d\n", __func__, chr_id, snp_pos, snp_hap1, snp_hap2);
                continue;
            }
            prev_pos = snp_pos;

            //pack info: encode SNP info into int
            int ref_int, alt_int, c = 0;
            if (alt_len == 1 && ref_len == 1){ //substitution
                ref_int = (mut_t)nst_nt4_table[(int)ref[0]];
                alt_int = (mut_t)nst_nt4_table[(int)alt[0]];
            } else {
                if ( ref[0] != alt[0] ) { //check if the first base are the same if it's indel, if not skip
                    fprintf(stderr, "[%s] Skip unusual SNP/INDEL [indel & ref]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
                    continue;
                }

                if (alt_len > 1 && ref_len == 1) { //insertion
                    ref_int = (mut_t)nst_nt4_table[(int)ref[0]];
                    base_change_pos = snp_pos + 1; // position +1, because the base change occurs after the first base
                    for (int i = 0; i < alt_len; i++ ){ 
                        c = (mut_t)nst_nt4_table[(int)alt[i]];
                        alt_int = (alt_int << 2) | c;
                    }
                } else if (ref_len > 1 && alt_len == 1){ //deletion
                    alt_int = (mut_t)nst_nt4_table[(int)alt[0]];
                    base_change_pos = snp_pos + 1; // position +1, because the base change occurs after the first base
                    for (int i = 0; i < ref_len; i++ ){ 
                        c = (mut_t)nst_nt4_table[(int)ref[i]];
                        ref_int = (ref_int << 2) | c;
                    }
                } else { // might be MNP, or sth else
                    fprintf(stderr, "[%s] Skip unusual SNP/INDEL [MNP]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
                    continue;
                }
            }

            int geno_int = base_offset << 12 | ref_len << 8 | snp_hap2 << 6 | snp_hap1 << 4 | is_phased;
            pos_vec.push_back(base_change_pos);
            ref_vec.push_back(ref_int);
            alt_vec.push_back(alt_int);
            geno_vec.push_back(geno_int);
        }
    }
    
    fprintf(stderr, "[%s] Finish collecting %lu SNP/INDEL from %s\n", __func__, pos_vec.size(), chr_id);

    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);

    if (int ret=hts_close(fp)){
        fprintf(stderr,"ERROR: [%s] hts_close(%s): non-zero status %d\n", __func__, fname, ret);
        exit(ret);
    }
}

void wgsim_mut_vcf(const kseq_t *ks, char * vcf_file, mutseq_t *hap1, mutseq_t *hap2, uint8_t *site_flag_arr)
{
    // initiate
    mutseq_t *ret[2];

    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    
    // parse VCF
    parse_vcf_chr(vcf_file, ks->name.s);

    int vec_ptr = 0;
    int i, deleting = 0;
    int deletion_count = 0;
    int c;

    for (i = 0; i != ks->seq.l; ++i){
        c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (pos_vec.size() == 0){ continue;} // ignore the rest if there is no SNP

        if (deleting){
            if(deletion_count > 0){
                if (deleting & 1){ ret[0]->s[i] |= DELETE;}
                if (deleting & 2){ ret[1]->s[i] |= DELETE;}
                deletion_count--;
                site_flag_arr[i] = 1;
                continue;
            } else {deleting = 0;}
        }

        if(vec_ptr < pos_vec.size() && i == pos_vec[vec_ptr] && c < 4){
            int geno_int = geno_vec[vec_ptr];
            int is_phased= geno_int & 0x000f;
            int snp_hap1 = (geno_int & 0x003f) >> 4;
            int snp_hap2 = (geno_int & 0x00ff) >> 6;
            int ref_len  = (geno_int & 0x0fff) >> 8;
            int base_offset = geno_int >> 12;
            site_flag_arr[i] = 1;
            //fprintf(stderr, "%d,%d,%d,%d,%d,%d\n", i, is_pahsed,snp_hap1, snp_hap2, ref_len, base_offset);

            if (!is_phased){ // for unphased genotype, randomly swap the haplotype
                if(drand48() < 0.5){
                    int tmp_hap = snp_hap1;
                    snp_hap1 = snp_hap2;
                    snp_hap2 = tmp_hap;
                }
            }

            if(base_offset == 0 && ref_len == 1){ // SNP substitution
                c = alt_vec[vec_ptr];

                if (snp_hap1 == 1 && snp_hap2 == 1){
                    ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                } else if (snp_hap1 == 1 && snp_hap2 == 0){
                    ret[0]->s[i] = SUBSTITUTE|c;
                } else if (snp_hap1 == 0 && snp_hap2 == 1){
                    ret[1]->s[i] = SUBSTITUTE|c;
                } else{continue;}
            } else if (base_offset < 0 ) { // deletion
                c = ref_vec[vec_ptr];
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
                int ins = alt_vec[vec_ptr] & ins_msk;
                //fprintf(stderr, "%d,%d,%d,%d\n", num_ins, ins_msk, alt_vec[vec_ptr], ins);

                if (snp_hap1 == 1 && snp_hap2 == 1){
                    ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else if (snp_hap1 == 1 && snp_hap2 == 0){
                    ret[0]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else if (snp_hap1 == 0 && snp_hap2 == 1){
                    ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else{continue;}
            }
            vec_ptr++;
        }
    }
}

void wgsim_mut_diref(const kseq_t *ks, int is_hap, mutseq_t *hap1, mutseq_t *hap2, uint8_t *site_flag_arr)
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
        if (deleting) {
            if (drand48() < INDEL_EXTEND) {
                if (deleting & 1) ret[0]->s[i] |= DELETE;
                if (deleting & 2) ret[1]->s[i] |= DELETE;
                site_flag_arr[i] = 1;
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
                    } while (num_ins < 4 && drand48() < INDEL_EXTEND);

                    if (is_hap || drand48() < 0.333333) { // hom-ins
                        ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                    } else { // het-ins
                        ret[drand48()<0.5?0:1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                    }
                }
            }
            site_flag_arr[i] = 1;
        }
    }
}

void wgsim_print_mutref(const char *name, const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2)
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

void wgsim_core(const char *fn, int is_hap, uint64_t N, int dist, int std_dev, int size_l, int size_r, int output_mode, char *vcf_file, char *contig_id, int64_t contig_N)
{
    kseq_t *ks;
    mutseq_t rseq[2];
    gzFile fp_fa;
    uint64_t tot_len, contig_len, ii;
    uint64_t tot_sub, tot_indel, tot_err, tot_pairs;
    int i, l, n_ref;
    char *qstr;
    int size[2], Q, max_size;
    uint8_t *tmp_seq[2];    	// save sequence
    uint8_t *tmp_context[2];	// save context (CG, CHG, CHH)
    uint8_t *tmp_mutation[2];	// save mutation status
    int8_t  *tmp_offset[2]; 	// save offset per base
    mut_t *target;
    
    max_size = std::max(size_l, size_r);
    qstr = (char*)calloc(max_size+1, 1);
    tmp_seq[0] = (uint8_t*)calloc(max_size+2, 1);
    tmp_seq[1] = (uint8_t*)calloc(max_size+2, 1);
    tmp_context[0] = (uint8_t*)calloc(max_size+2, 1);
    tmp_context[1] = (uint8_t*)calloc(max_size+2, 1);
    tmp_mutation[0]= (uint8_t*)calloc(max_size+2, 1);
    tmp_mutation[1]= (uint8_t*)calloc(max_size+2, 1);
    tmp_offset[0]  = (int8_t*)calloc(max_size+2, 1);
    tmp_offset[1]  = (int8_t*)calloc(max_size+2, 1);
    size[0] = size_l; size[1] = size_r;
    
    int min_inner, max_inner, inner_dist;
    min_inner = -1 * std::min(size_l, size_r); //when the one read's start is the other read's end
    max_inner = 600 - 2* max_size; // observed from real data, max insert size is ~600

    Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;

    fp_fa = gzopen(fn, "r");
    ks = kseq_init(fp_fa);
    tot_len = n_ref = 0;
    tot_sub = tot_indel = tot_err = tot_pairs = 0;
    bool bool_contig = false; 
    bool bool_contig_N = false;
    fprintf(stderr, "[%s] calculating the total length of the reference sequence...\n", __func__);
    while ((l = kseq_read(ks)) >= 0) {
        tot_len += l;
        ++n_ref;
        if (strcmp(contig_id, ks->name.s)==0){ bool_contig = true; contig_len = l;}
    }
    // if fasta is non-existing or empty
    if (!fp_fa) { fprintf (stderr, "ERROR: gzopen of '%s' failed: %s. Exit... \n", fn, strerror (errno)); exit (EXIT_FAILURE);}
    if (!n_ref) { fprintf (stderr, "ERROR: Input fasta is empty: %s. Exit... \n", fn); exit (EXIT_FAILURE);}
    fprintf(stderr, "[%s] %d contig sequences, total length: %lu\n", __func__, n_ref, tot_len);
    
    // check input contig_id
    if (strcmp(contig_id, "None") == 0 || strlen(contig_id) == 0) {
        if (contig_N > 0 ) {fprintf(stderr, "ERROR: -n is specified but not -c, exit...(please note the difference of -n and -N)\n"); exit(1);}
        fprintf(stderr, "[%s] No contig id specified, will generate %lu reads from all contigs\n", __func__, N);
    } else {
        if (bool_contig && contig_N <= 0) {
            contig_N = int64_t((long double)contig_len / tot_len * N + 0.5);
            fprintf(stderr, "[%s] Simulate %ld reads from contig %s (calculate from -N, as -n is not specified)...\n", __func__, contig_N, contig_id);
        } else if (bool_contig && contig_N > 0) {
            fprintf(stderr, "[%s] Simulate %ld reads from contig %s\n", __func__, contig_N, contig_id);
            bool_contig_N = true;
        } else {
            fprintf(stderr, "ERROR: Contig id '%s' is not found in the fasta, please check!\n", contig_id); exit(1);
        }
    }
    kseq_destroy(ks);
    gzclose(fp_fa);

    // check input vcf file
    FILE *vcf;
    bool bool_vcf = false;
    if (strcmp(vcf_file, "None") == 0 || strlen(vcf_file) == 0) {
        fprintf(stderr, "[%s] No VCF input, will generate SNP randomly if mutation rate is nonzero\n", __func__);
    } else if(vcf=fopen(vcf_file,"r")) {
        fprintf(stderr, "[%s] VCF file exists, use it to simulate reads\n", __func__);
        bool_vcf = true;
        fclose(vcf);
    } else {
        fprintf(stderr, "ERROR: The specified VCF file does not exist, please check!\n");
        exit(1);
    }

    fp_fa = gzopen(fn, "r");
    ks = kseq_init(fp_fa);
    while ((l = kseq_read(ks)) >= 0) {//here l is the chromosome length
        if (bool_contig) {if (strcmp(contig_id, ks->name.s)!=0){continue;}}
        // initialize random number generator to generate read positions
        std::random_device rd;
        std::mt19937 gen(rd());
        // pull from truncated distribution to ensure read doesn't pass boundary 
        std::uniform_int_distribution<int> dis(2, ks->seq.l - max_inner - 2*max_size -2); //add 2 base offset

        uint64_t n_pairs = (uint64_t)((long double)l / tot_len * N + 0.5);
        if (bool_contig_N){n_pairs = contig_N;}

        tot_pairs += n_pairs;
        // initialise random normal for insert size simulation 
        std::random_device rn;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen_rn(rn()); //Standard mersenne_twister_engine seeded with rd()
        std::normal_distribution<float> dis_rn(0.0, 1.0);
        if (l < dist + 3 * std_dev) {
            fprintf(stderr, "[%s] skip sequence '%s' as it is shorter than %d!\n", __func__, ks->name.s, dist + 3 * std_dev);
            continue;
        }

        uint8_t site_flag_arr[ks->seq.l] = {0}; // record if the site is a SNP/INDEL position (the base can be REF)
        // introduce mutations and print them to stdout
        fprintf(stdout, "Contig Variant Start\n");
        if(bool_vcf){
            wgsim_mut_vcf(ks, vcf_file, rseq, rseq+1, site_flag_arr);
            if(pos_vec.size() == 0){fprintf(stdout, "%s\n", ks->name.s);} //if no variants, print contig id
        } else {
            wgsim_mut_diref(ks, is_hap, rseq, rseq+1, site_flag_arr);
            if(MUT_RATE == 0.0){fprintf(stdout, "%s\n", ks->name.s);}
        }
        wgsim_print_mutref(ks->name.s, ks, rseq, rseq+1);
        fprintf(stdout, "Contig Variant End\n");
        
        for (ii = 0; ii != n_pairs; ++ii) { // the core loop
            double ran;
            int pos, s[2];
            int n_sub[2], n_indel[2], n_err[2], ext_coor[2], cover_pos[2], j, k, ix;
            //cover_pos hold if the read covers a snp *position* (the read don't have to contain the ALT allele)
            //j hold read1/read2, k hold the length of read, ix hold the cursor transversing read

            pos = dis(gen);
            int insert_dev = (int)(std_dev * dis_rn(gen_rn));
            inner_dist = dist + insert_dev - size_l - size_r; //dist is the mean insert size
            inner_dist = std::max(min_inner, std::min(inner_dist, max_inner));
            s[0] = size[0]; s[1] = size[1];

            // generate the read sequences
            target = rseq[drand48()<0.5?0:1].s; // haplotype from which the reads are generated
            n_sub[0] = n_sub[1] = n_indel[0] = n_indel[1] = n_err[0] = n_err[1] = cover_pos[0] = cover_pos[1] =0;
            int start[2] = {pos, pos + inner_dist + size_l};
            int end[2] = {start[0], start[1]};
            int offset[2] = {0, 0};

            // x: select read1 or read2; ext_coor: extend corrdinates;
            #define __gen_read(x, start_pos, iter) do {                 \
                /* generate reads assign mutation flag; */              \
                for (i = (start_pos), k = 0, ext_coor[x] = -10; i >= 0 && i < ks->seq.l && k < s[x]; iter) { \
                    int c = target[i], mut_type = c & mutmsk;           \
                    if (ext_coor[x] < 0) {                              \
                        /* avoid indel as the first base */             \
                        if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue; \
                        start[x] = i;                                   \
                        end[x] = i;                                     \
                        ext_coor[x] = i;                                \
                    }                                                   \
                    if (mut_type == DELETE){                            \
                        ++offset[x];                                    \
                        ++end[x];                                       \
                        ++n_indel[x];                                   \
                    }                                                   \
                    else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
                        /* context: 0x00 Match, 0x01 SNP, 0x03 INSERT   \
                                    0x01 CG, 0x03 CHG, 0x07 CHH (>>)    \
                                    0x09 GC, 0x0b GDC, 0x0f GDD (<<) */ \
                        tmp_seq[x][k] = c & 0xf;                        \
                        tmp_offset[x][k] = offset[x];                   \
                        if (mut_type == SUBSTITUTE) {                   \
                            ++n_sub[x];                                 \
                            tmp_mutation[x][k] = SNV;                   \
                        } else {                                        \
                            tmp_mutation[x][k] = MATCH;                 \
                        }                                               \
                        ++end[x];                                       \
                        ++k;                                            \
                    } else {                                            \
                        tmp_seq[x][k] = c & 0xf;                        \
                        tmp_offset[x][k] = offset[x];                   \
                        tmp_mutation[x][k] = MATCH;/*The base is ref*/  \
                        ++n_indel[x];                                   \
                        ++end[x];                                       \
                        ++k;                                            \
                        int num_ins, ins;                               \
                        for (num_ins = mut_type>>12, ins = c>>4; num_ins > 0 && k < s[x]; --num_ins, ins >>= 2){ \
                            --offset[x];                                \
                            tmp_seq[x][k] = ins & 0x3;                  \
                            tmp_offset[x][k] = offset[x];               \
                            tmp_mutation[x][k] = INSR;                  \
                            ++k;                                        \
                        }                                               \
                    }                                                   \
                    cover_pos[x] |= site_flag_arr[i];                   \
                }                                                       \
                /* append CG context flag;                              \
                   currently not handling bounday context by indel*/    \
                for (ix=0; ix < k; ++ix) {                              \
                    int c_d1, c_d2;                                     \
                    int c = tmp_seq[x][ix];                             \
                    if (cg_table[(uint8_t) c]) continue;                \
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
                    tmp_context[x][ix] = cg_context_table[context_idx]; \
                }                                                       \
                if (k != s[x]) ext_coor[x] = -10;                       \
            } while (0)

            __gen_read(0, pos, ++i);
            __gen_read(1, pos + size_l + inner_dist, ++i);

            if (ext_coor[0] < 0 || ext_coor[1] < 0) { // failed to generate the read(s)
                --ii;
                continue;
            }
            for(j = 0; j < 2; ++j){ //check the number of Ns
                int n_n =0;
                for (i = 0; i < s[j]; ++i) {
                    int c = tmp_seq[j][i];
                    if (c >= 4) { // actually c should be never larger than 4 if everything is correct
                        ++n_n; 
                        tmp_seq[j][i] = 4;
                    }
                    qstr[i] = Q; // generate the quality score
                }
                qstr[i] = 0;
                if ((double)n_n / s[j] > MAX_N_RATIO) break;
            }
            if (j < 2) { // too many ambiguous bases on one of the reads
                --ii;
                continue;
            }
            
            int flag_pos = cover_pos[0] | cover_pos[1];
            int flag_mut  = n_sub[0] + n_sub[1] + n_indel[0] + n_indel[1];
            int flag_indel= n_indel[0] + n_indel[1];
            // print reads to stdout: mode 0 print string (WGS), else print numbers (WGBS)
            if(output_mode == 0){
                // flip and get the reverse complementary
                int is_flip = drand48() < 0.5? 0 : 1;
                for (k = 0; k < s[1]; ++k) { 
                    if (k <= int(s[1]/2)) { // reverse
                        int tmp_base  = tmp_seq[1][k];
                        tmp_seq[1][k] = tmp_seq[1][s[1]-k];
                        tmp_seq[1][s[1]-k] = tmp_base;
                        int tmp_cigar = tmp_mutation[1][k];
                        tmp_mutation[1][k] = tmp_mutation[1][s[1]-k];
                        tmp_mutation[1][s[1]-k] = tmp_cigar;
                    }
                    tmp_seq[1][k] = tmp_seq[1][k] < 4? 3 - tmp_seq[1][k] : 4; // complement
                }
                for (j = 0; j < 2; ++j) {
                    int jj = j ^ is_flip; // x^0=x; x^1=!x (when x is binary)
                    // header: 1-based coordinates for string output
                    fprintf(stdout, "@%s:%d:%d:%llx:%d:%d:%d/%d\n", ks->name.s, start[0]+1, end[1]+1, (long long)ii, flag_pos, flag_mut, flag_indel, j+1);
                    // sequence (introduce random sequencing error)
                    for (i = 0; i < s[jj]; ++i) {
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
                    fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:%d:", start[jj]+1, end[jj]+1, cover_pos[jj], n_sub[jj], n_indel[jj], n_err[jj], end[1]-start[0], end[0]-start[1]);
                    for (i = 0; i < s[jj]; ++i) {
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
                    fprintf(stdout, "@%s:%d:%d:%llx %d %d %d %d ", ks->name.s, start[0]+1, end[1]+1, (long long)ii, j+1, flag_pos, flag_mut, flag_indel); 
                    for (i = 0; i < s[j]; ++i) {
                        //fprintf(stdout, "%x", tmp_mutation[j][i]); //this will output the hex number, below line will output the ascii
                        fputc(tmp_mutation[j][i], stdout);
                    }
                    fprintf(stdout, "\n");
                    // sequence (no sequencing error, represented by 0-4)
                    for (i = 0; i < s[j]; ++i) {
                        //fprintf(stdout, "%d", tmp_seq[j][i]);
                        fputc(tmp_seq[j][i], stdout);
                    }
                    fprintf(stdout, "\n");
                    // comment
                    fprintf(stdout, "+:%d:%d:%d:%d:%d:%d:%d:", start[j], end[j], cover_pos[j], n_sub[j], n_indel[j], end[1]-start[0], end[0]-start[1]);
                    const char *pad = "";
                    for (i = 0; i < s[j]; ++i) {
                        fprintf(stdout, "%s%d", pad, tmp_offset[j][i]);
                        pad = ",";
                    }
                    fprintf(stdout, "\n");
                    // quality
                    for (i = 0; i < s[j]; ++i) {
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
        free(rseq[0].s); free(rseq[1].s);
    }

    fprintf(stderr, "[%s] Generated %lu read pairs, with %lu contain SNP, %lu contain INDEL", __func__, tot_pairs, tot_sub, tot_indel);
    if (output_mode == 0){fprintf(stderr, " and %lu contain sequencing errors\n", tot_err);} else{fprintf(stderr, "\n");}
    kseq_destroy(ks);
    gzclose(fp_fa);
    free(qstr);
    free(tmp_seq[0]); free(tmp_seq[1]);
}

static int simu_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Forked wgsim (short read simulator) for simulating WGS or WGBS reads\n");
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Contact: Wenbin Guo <wbguo@ucla.edu>; \n\n");
    fprintf(stderr, "Usage:   wgsim [options] <in.ref.fa> \n\n");
    fprintf(stderr, "Options: -e FLOAT      base error rate [%.3f]\n", ERR_RATE);
    fprintf(stderr, "         -d INT        outer distance between the two ends [500]\n");
    fprintf(stderr, "         -s INT        standard deviation [50]\n");
    fprintf(stderr, "         -N INT        number of read pairs [1000000]\n");
    fprintf(stderr, "         -1 INT        length of the first read [70]\n");
    fprintf(stderr, "         -2 INT        length of the second read [70]\n");
    fprintf(stderr, "         -r FLOAT      rate of mutations [%.4f]\n", MUT_RATE);
    fprintf(stderr, "         -R FLOAT      fraction of indels [%.2f]\n", INDEL_FRAC);
    fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", INDEL_EXTEND);
    fprintf(stderr, "         -S INT        seed for random generator [-1]\n");
    fprintf(stderr, "         -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [%.2f]\n", MAX_N_RATIO);
    fprintf(stderr, "         -h INT        haplotype mode: 0 for No; non-zero for Yes [0]\n");
    fprintf(stderr, "         -m INT        output mode: 0 for letters (sequence: AGCT,CIGAR: MXI); 1 for numbers (sequence: 0123, and flags) [0]\n");
    fprintf(stderr, "         -g STRING     path to the genetic variant file (vcf.gz) [None]\n");
    fprintf(stderr, "         -c STRING     contig name [None]\n");
    fprintf(stderr, "         -n INT        number of read pairs for specified contig [-1]\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[])
{
    uint64_t N;
    int64_t contig_N;
    int dist, std_dev, c, size_l, size_r, is_hap = 0;
    int output_mode= 0;
    char none_default[] = "None";
    char *vcf_file = none_default;
    char *contig_id= none_default; // checked, will not intefere with vcf_file
    int seed = -1;

    N = 1000000; dist = 500; std_dev = 50; contig_N = -1;
    size_l = size_r = 70;
    while ((c = getopt(argc, argv, "e:d:s:N:1:2:r:R:h:X:S:A:m:g:c:n:")) >= 0) {
        switch (c) {
        case 'd': dist = atoi(optarg); break;
        case 's': std_dev = atoi(optarg); break;
        case 'N': N = atoi(optarg); break;
        case '1': size_l = atoi(optarg); break;
        case '2': size_r = atoi(optarg); break;
        case 'e': ERR_RATE = atof(optarg); break;
        case 'r': MUT_RATE = atof(optarg); break;
        case 'R': INDEL_FRAC = atof(optarg); break;
        case 'X': INDEL_EXTEND = atof(optarg); break;
        case 'A': MAX_N_RATIO = atof(optarg); break;
        case 'S': seed = atoi(optarg); break;
        case 'h': is_hap = atoi(optarg); break;
        case 'm': output_mode = atoi(optarg); break;
        case 'g': vcf_file = optarg; break;
        case 'c': contig_id= optarg; break;
        case 'n': contig_N = atoi(optarg); break;
        }
    }
    if (argc - optind < 1) return simu_usage();
    if (seed <= 0) seed = time(0)&0x7fffffff;
    fprintf(stderr, "[wgsim] seed = %d\n", seed);
    srand48(seed);
    wgsim_core(argv[optind], is_hap, N, dist, std_dev, size_l, size_r, output_mode, vcf_file, contig_id, contig_N);

    return 0;
}

