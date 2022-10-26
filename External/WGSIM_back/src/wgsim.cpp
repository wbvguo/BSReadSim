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


/* wgsim */

enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
static mut_t mutmsk = (mut_t)0xf000;

typedef struct {
	int l, m; /* length and maximum buffer size */
	mut_t *s; /* sequence */
} mutseq_t;

static double ERR_RATE = 0.02;
static double MUT_RATE = 0.001;
static double INDEL_FRAC = 0.15;
static double INDEL_EXTEND = 0.3;
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

	bool collect_present = false;
	bool collect_previous= false;

	//check the vcf file
	int nsmpl = bcf_hdr_nsamples(hdr); // number of sample
	if (nsmpl != 1) {
		fprintf(stderr, "[%s] Currently BSReadSim only supports single-sample simulation, please check the vcf file. Exiting...\n", __func__);
		exit(1);
	}

    while (bcf_read(fp, hdr, rec)>=0) {	
		if (collect_present == false && collect_previous == true) { // finished collecting
			break;
		}
		
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
			// 1. insert/delete length greater than 4
			// 2. multi-allelic sites
			// 3. contains no-call values
			// 4. multi-ploid organism
			// 5. multi nucleotide polymorphism (MNP); (consider to include?)
			// TODO: 
			// 1. contains Ns;
			// 2. same POS; 
			// 3. missing snps?
			if (abs(base_offset) > 4 || ngt > 2) {
				fprintf(stderr, "[%s] Skip unusual SNP [length/ploidy]: CHROM:%s; POS:%d; REF:%s; ALT:%s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
				continue;
			}
			if (snp_hap1 < 0 || snp_hap2 < 0 || snp_hap1 > 1 || snp_hap2 > 1 || (snp_hap1 == 0 && snp_hap2 == 0)) {
				fprintf(stderr, "[%s] Skip unusual SNP [snp haplotype]: CHROM:%s; POS:%d; HAP1:%d; HAP2:%d\n", __func__, chr_id, snp_pos, snp_hap1, snp_hap2);
				continue;
			}
			

			//pack info: encode SNP info into int
			int ref_int, alt_int, c = 0;
			if (alt_len == 1 && ref_len == 1){ //substitution
				ref_int = (mut_t)nst_nt4_table[(int)ref[0]];
				alt_int = (mut_t)nst_nt4_table[(int)alt[0]];
			} else {
				if ( ref[0] != alt[0] ) { //check if the first base are the same if it's indel, if not skip
					fprintf(stderr, "[%s] Skip unusual SNP [indel & ref]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
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
					fprintf(stderr, "[%s] Skip unusual SNP [MNP]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
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
	
	fprintf(stderr, "[%s] Finish collecting %lu SNP from %s\n", __func__, pos_vec.size(), chr_id);

    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);

    if (int ret=hts_close(fp)){
        fprintf(stderr,"[%s] hts_close(%s): non-zero status %d\n", __func__, fname, ret);
        exit(ret);
	}
}

void wgsim_mut_vcf(const kseq_t *ks, char * vcf_file, mutseq_t *hap1, mutseq_t *hap2)
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
			//fprintf(stderr, "%d,%d,%d,%d,%d,%d\n", i, is_phased,snp_hap1, snp_hap2, ref_len, base_offset);

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

void wgsim_mut_diref(const kseq_t *ks, int is_hap, mutseq_t *hap1, mutseq_t *hap2)
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
		//fprintf(stderr, "%s,%d,%d,%d,%d\n", ks->name.s,i,c[0],c[1],c[2]);
		if (c[0] >= 4) continue;
		if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
			if (c[1] == c[2]) { // hom
				if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
					printf("%s\t%d\t%c\t%c\t-\n", name, i+1, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]);
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

void wgsim_core(const char *fn, int is_hap, uint64_t N, int dist, int std_dev, int size_l, int size_r, int mean_inner, int output_mode, char *vcf_file)
{
	kseq_t *ks;
    mutseq_t rseq[2];
	gzFile fp_fa;
	uint64_t tot_len, ii;
	int i, l, n_ref;
	char *qstr;
	int size[2], Q, max_size;
	uint8_t *tmp_seq[2];
    mut_t *target;
	
	l = size_l > size_r? size_l : size_r;
	qstr = (char*)calloc(l+1, 1);
	tmp_seq[0] = (uint8_t*)calloc(l+2, 1);
	tmp_seq[1] = (uint8_t*)calloc(l+2, 1);
	size[0] = size_l; size[1] = size_r;
	max_size = size_l > size_r? size_l : size_r;
	
	int min_inner, max_inner, inner_dist;
	min_inner = -1 * size_r; //originally 10, here changed to the -size_r
	max_inner = dist - 2 * size_l; // dist is the insert size

	Q = (ERR_RATE == 0.0)? 'I' : (int)(-10.0 * log(ERR_RATE) / log(10.0) + 0.499) + 33;

	fp_fa = gzopen(fn, "r");
	ks = kseq_init(fp_fa);
	tot_len = n_ref = 0;
	fprintf(stderr, "[%s] calculating the total length of the reference sequence...\n", __func__);
	while ((l = kseq_read(ks)) >= 0) {
		tot_len += l;
		++n_ref;
	}
	fprintf(stderr, "[%s] %d sequences, total length: %llu\n", __func__, n_ref, (long long)tot_len);
	kseq_destroy(ks);
	gzclose(fp_fa);

	// vcf check
	FILE *vcf;
	bool bool_vcf = false;
	if (strcmp(vcf_file, "None") == 0 || strlen(vcf_file) == 0) {
		fprintf(stderr, "No VCF input, will generate SNP randomly if mutation rate is nonzero\n");
	} else if(vcf=fopen(vcf_file,"r")) {
		fprintf(stderr, "VCF file exists, use it to simulate reads\n");
		bool_vcf=true;
		fclose(vcf);
	} else {
		fprintf(stderr, "The specified VCF file does not exist, please check\n");
		exit(1);
	}

	fp_fa = gzopen(fn, "r");
	ks = kseq_init(fp_fa);
	while ((l = kseq_read(ks)) >= 0) {//here l is the chromosome length
		// initialize random number generator to generate read positions
		std::random_device rd;
		std::mt19937 gen(rd());
		// pull from truncated distribution to ensure read doesn't pass boundary 
		std::uniform_int_distribution<int> dis(0, ks->seq.l - dist - 100);

		uint64_t n_pairs = (uint64_t)((long double)l / tot_len * N + 0.5);
		// initialise random normal for insert size simulation 
		std::random_device rn;  //Will be used to obtain a seed for the random number engine
    	std::mt19937 gen_rn(rn()); //Standard mersenne_twister_engine seeded with rd()
    	std::normal_distribution<float> dis_rn(0.0, 1.0);
		if (l < dist + 3 * std_dev) {
			fprintf(stderr, "[%s] skip sequence '%s' as it is shorter than %d!\n", __func__, ks->name.s, dist + 3 * std_dev);
			continue;
		}

		// generate mutations and print them out
		fprintf(stdout, "Contig Variant Start\n");
		if(bool_vcf){
			wgsim_mut_vcf(ks, vcf_file, rseq, rseq+1);
			if(pos_vec.size() == 0){fprintf(stdout, "%s\n", ks->name.s);} //if no variants, print contig id
		} else {
			wgsim_mut_diref(ks, is_hap, rseq, rseq+1);
			if(MUT_RATE == 0.0){fprintf(stdout, "%s\n", ks->name.s);}
		}
		wgsim_print_mutref(ks->name.s, ks, rseq, rseq+1);
		fprintf(stdout, "Contig Variant End\n");

		for (ii = 0; ii != n_pairs; ++ii) { // the core loop
			double ran;
			int d, pos, s[2];
			int n_sub[2], n_indel[2], n_err[2], ext_coor[2], j, k;

			pos = dis(gen);
			int insert_dev = (int)(std_dev * dis_rn(gen_rn));
			inner_dist = mean_inner + insert_dev;
			inner_dist = std::max(min_inner, std::min(inner_dist, max_inner));
			s[0] = size[0]; s[1] = size[1];

			// generate the read sequences
			target = rseq[drand48()<0.5?0:1].s; // haplotype from which the reads are generated
			n_sub[0] = n_sub[1] = n_indel[0] = n_indel[1] = n_err[0] = n_err[1] = 0;
			int start[2] = {pos, pos + inner_dist + size_l};
			int end[2] = {start[0], start[1]};
			std::string cigar[2];
			int offset[2] = {0, 0};
			int c_count[2] = {0, 0};
			int g_count[2] = {0, 0};
			int** c_bases = new int*[2]; // a pointer to a pointer to an int; why new and delete? 
			int** g_bases = new int*[2];
			int** c_offset = new int*[2];
			int** g_offset = new int*[2];
			c_bases[0] = new int[l + 1];
			c_bases[1] = new int[l + 1];
			g_bases[0] = new int[l + 1];
			g_bases[1] = new int[l + 1];
			c_offset[0] = new int[l + 1];
			c_offset[1] = new int[l + 1];
			g_offset[0] = new int[l + 1];
			g_offset[1] = new int[l + 1];

			//Why these \ symbol, in line 595, why it's M for Insertion? (I changed to I)
			// x: read1 or read2; ext_coor: extend corrdinates; 
			#define __gen_read(x, start, iter) do {									\
				for (i = (start), k = 0, ext_coor[x] = -10; i >= 0 && i < ks->seq.l && k < s[x]; iter) {	\
					int c = target[i], mut_type = c & mutmsk;			\
					if (ext_coor[x] < 0) {								\
						/*if the first base is deletion, then move until no deletion*/ \
						if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue; \
						ext_coor[x] = i;								\
					}													\
					if (mut_type == DELETE){                            \
						++n_indel[x];                                   \
						/*cigar[x] += 'D';*/                            \
						++offset[x];                                    \
						++end[x];	                                    \
					}			                                        \
					else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
						char base = c & 0xf;                            \
						int ob = nst_nt4_table[(int)ks->seq.s[i]];      \
						tmp_seq[x][k] = base;						    \
						if (mut_type == SUBSTITUTE) {                   \
							++n_sub[x];	                                \
							cigar[x] += 'X';                            \
							++end[x];	                                \
						} else {                                        \
							cigar[x] += 'M';		                    \
							++end[x];}                                  \
						if ((int)base == 1 || ob == 1) {            \
	                        c_bases[x][c_count[x]] = k;                 \
							c_offset[x][c_count[x]] = offset[x];        \
							++c_count[x];                               \
							}                                           \
						else if ((int)base == 2 || ob == 2){        \
	                        g_bases[x][g_count[x]] = k;                 \
							g_offset[x][g_count[x]] = offset[x];        \
							++g_count[x];                               \
							}                                           \
						k++;                                            \
					} else {											\
						int n, ins;										\
						++n_indel[x];	                                \
						cigar[x] += 'I';								\
						char base = c & 0xf;                            \
						tmp_seq[x][k] = base;						    \
						if ((int)base == 1) {                           \
	                        c_bases[x][c_count[x]] = k;                 \
							c_offset[x][c_count[x]] = offset[x];        \
							++c_count[x];                               \
							}                                           \
						else if ((int)base == 2){                       \
	                        g_bases[x][g_count[x]] = k;                 \
							g_offset[x][g_count[x]] = offset[x];        \
							++g_count[x];                               \
							}                                           \
						++end[x];                                       \
						k++;                                            \
						int i_count = 1;                                \
						for (n = mut_type>>12, ins = c>>4; n > 0 && k < s[x]; --n, ins >>= 2){ \
							cigar[x] += std::to_string(i_count);	    \
							++ i_count;                                 \
							--offset[x];                                \
							tmp_seq[x][k++] = ins & 0x3;				\
						}                                               \
					}													\
				}														\
				if (k != s[x]) ext_coor[x] = -10;						\
			} while (0)

			__gen_read(0, pos, ++i);
			__gen_read(1, pos + inner_dist + size_l, ++i);
			//for (k = 0; k < s[1]; ++k) tmp_seq[1][k] = tmp_seq[1][k] < 4? 3 - tmp_seq[1][k] : 4; // complement
			if (ext_coor[0] < 0 || ext_coor[1] < 0) { // failed to generate the read(s)
				--ii;
				goto CLEAN_CG_PTR;
			}

			// generate sequencing errors
			for (j = 0; j < 2; ++j) {
				int n_n = 0;
				for (i = 0; i < s[j]; ++i) {
					int c = tmp_seq[j][i];
					if (c >= 4) { // actually c should be never larger than 4 if everything is correct
						c = 4;
						++n_n;
					} else if (drand48() < ERR_RATE) {
						// c = (c + (int)(drand48() * 3.0 + 1)) & 3; // random sequencing errors
						c = (c + 1) & 3; // recurrent sequencing errors
						++n_err[j];
					}
					tmp_seq[j][i] = c;
				}
				if ((double)n_n / s[j] > MAX_N_RATIO) break;
			}
			if (j < 2) { // too many ambiguous bases on one of the reads
				--ii;
				goto CLEAN_CG_PTR;
			}

			// print
			if (output_mode == 0) {
				
			}


			for (j = 0; j < 2; ++j) {
				for (i = 0; i < s[j]; ++i) qstr[i] = Q;
				qstr[i] = 0;
				fprintf(stdout, "@%s:%d:%d:%d:%llx:%s:%d:", ks->name.s, start[j], end[j],
						end[1] - start[0], (long long)ii, cigar[j].c_str(), j + 1);
				for (i=0; i != c_count[j]; ++i) {
					fprintf(stdout, "%d_%d,", c_bases[j][i], c_offset[j][i]);
				}
				//if (c_count[j]){
				//	fprintf(stdout, "%d_%d", c_bases[j][c_count[j]], c_offset[j][c_count[j]]);
				//}
				fprintf(stdout, ":");
				for (i=0; i != g_count[j]; ++i) {
					fprintf(stdout, "%d_%d,", g_bases[j][i], g_offset[j][i]);
				}
				//if (g_count[j]){
				//	fprintf(stdout, "%d_%d", g_bases[j][g_count[j]], g_offset[j][g_count[j]]);
				//}
				fprintf(stdout, "\n");
				for (i = 0; i < s[j]; ++i)
					fputc("ACGTN"[(int)tmp_seq[j][i]], stdout);
				fprintf(stdout, "\n+\n%s\n", qstr);

			}
			
			CLEAN_CG_PTR:
			delete c_bases[0];
			delete c_bases[1];
			delete c_bases;
			delete g_bases[0];
			delete g_bases[1];
			delete g_bases;
			delete c_offset[0];
			delete c_offset[1];
			delete c_offset;
			delete g_offset[0];
			delete g_offset[1];
			delete g_offset;
		}
		free(rseq[0].s); free(rseq[1].s);
	}
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
	fprintf(stderr, "Contact: Wenbin Guo <wbguo@ucla.edu>; Hongxiang Fu; Junxi Feng;\n\n");
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
	fprintf(stderr, "         -m INT		output mode: 0 for letters (AGCT,MIX) ; 1 for numbers (0123,012) [0]\n");
	fprintf(stderr, "         -g STRING     path to the genetic variant file (vcf.gz) [None]\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int64_t N;
	int dist, std_dev, c, size_l, size_r, is_hap = 0;
	int output_mode = 0;
	char vcf_default[] = "None";
	char *vcf_file = vcf_default;
	int seed = -1;

	N = 1000000; dist = 500; std_dev = 50;
	size_l = size_r = 70;
	int mean_inner = (dist - size_l * 2) / 2;
	while ((c = getopt(argc, argv, "e:d:s:N:1:2:r:R:h:X:S:A:I:m:g:")) >= 0) {
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
		case 'I': mean_inner = atoi(optarg); break;
		case 'm': output_mode = atoi(optarg); break;
		case 'g': vcf_file = optarg; break;
		}
	}
	if (argc - optind < 1) return simu_usage();
	if (seed <= 0) seed = time(0)&0x7fffffff;
	fprintf(stderr, "[wgsim] seed = %d\n", seed);
	srand48(seed);
	wgsim_core(argv[optind], is_hap, N, dist, std_dev, size_l, size_r, mean_inner, output_mode, vcf_file);

	return 0;
}