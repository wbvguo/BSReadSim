/* The MIT License
   Copyright (c) 2008 Genome Research Ltd (GRL).
                 2011 Heng Li <lh3@live.co.uk>
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
 * with modifications to simulate WGS/WGBS reads */

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
//#include <iostream>
//#include <string>
#include "kseq.h"
#include "vcf.h"

KSEQ_INIT(gzFile, gzread)

#define PACKAGE_VERSION "0.0.2"

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

// to restore vcf information
std::vector<int> positions;
std::vector<std::string>  alt_vec;
std::vector<std::string>  ref_vec;
std::vector<int> h1_genotype;
std::vector<int> h2_genotype;

void parse_vcf_chr(char *fname, char *chr_id){
    //open vcf file
    htsFile *fp    = hts_open(fname,"rb");

    //read header
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t *rec    = bcf_init();

    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
	bool collect_present = false;
	bool collect_previous= false;

    fprintf(stderr, "VCF file contains %i samples\n", bcf_hdr_nsamples(hdr));
    //save for each vcf record
    while ( bcf_read(fp, hdr, rec)>=0)
    {
    	//unpack for read REF,ALT,INFO,etc 
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);

        //read CHROM
        collect_previous = collect_present;
		std::string chr_current = bcf_hdr_id2name(hdr, rec->rid);

		if (strcmp(chr_current.c_str(), chr_id)){
			collect_present = false;
			continue;
		} else {
			collect_present = true;
			//read POS
			//printf("position: %lu\n", (unsigned long)rec->pos);
			positions.push_back(rec->pos);

			std::string ref = rec->d.allele[0];
			std::string alt = rec->d.allele[1];
			ref = ref.c_str();
			alt = alt.c_str();
			
			alt_vec.push_back(alt);
			ref_vec.push_back(ref);
			
			ngt=bcf_get_genotypes(hdr,rec,&gt,&ngt_arr);
			int hap1 = bcf_gt_allele(gt[0]);
			int hap2 = bcf_gt_allele(gt[1]);

			h1_genotype.push_back(hap1);
			h2_genotype.push_back(hap2);

		}
		if (collect_present == false & collect_previous == true) {
			break;
		}
    }

    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
	int ret;
    if ( ret=hts_close(fp)){
        fprintf(stderr, "hts_close(%s): non-zero status %d\n", fname, ret);
        exit(ret);
    }
}
void wgsim_mut_vcf(const kseq_t *ks, int is_hap, mutseq_t *hap1, mutseq_t *hap2, char * vcf_file){
	int vector_position = 0;
	//std::cout << positions.size();
	char *chr_id = ks->name.s;
	parse_vcf_chr(vcf_file, chr_id);

	int i, deleting = 0;
	mutseq_t *ret[2];
	ret[0] = hap1; 
	ret[1] = hap2;
	ret[0]-> l = ks->seq.l;
	ret[1]-> l = ks->seq.l;
	ret[0]-> m = ks->seq.m; 
	ret[1]-> m = ks->seq.m;
	ret[0]-> s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
	ret[1]-> s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));

	bool homo = false;
	bool hap1_set = false;
	bool hap2_set = false;

	int count  =0 ;

	int deletion_count = 0;
	for (i = 0; i != ks->seq.l; ++i) {

		int c;

		c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
			
		if (deleting) {

			if(deletion_count > 0){
            	if (deleting & 1){
					ret[0]->s[i] |= DELETE;
				}
            	if (deleting & 2){
					ret[1]->s[i] |= DELETE;
				}
				deletion_count--;
				continue;
			}
			else deleting = 0;
        }

		if(i == positions[vector_position] && c < 4){

			homo = false;
			hap1_set = false;
			hap2_set = false;
			int hap1 = h1_genotype[vector_position];
			int hap2 = h2_genotype[vector_position];

			if (hap1 == 0 && hap2 == 0){
				vector_position++;
				continue;
			} else if (hap1 == 1 && hap2 == 1){
				homo = true;
			} else if (hap1 == 0){
				hap2_set = true;
			} else{
				hap1_set = true;
			}

			//std::cout << ref_vec[vector_position].length();
			//std::cout << alt_vec[vector_position].length();
			// substitution
			if(ref_vec[vector_position].length() == 1 && alt_vec[vector_position].length() == 1){
				
				count ++;

				//std::cout << alt_vec[vector_position];
				if(alt_vec[vector_position].compare("G")==0){
					c = 2;
				} else if (alt_vec[vector_position].compare("A")==0){
					c = 0;
				} else if (alt_vec[vector_position].compare("C")==0){
					c = 1;
				} else {
					c = 3;
				}

				if(homo){
				ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
				vector_position++;
				} else if(hap1_set){
					//std::cout << "h1";
					ret[0]->s[i] = SUBSTITUTE|c;
					vector_position++;
				} else {
					ret[1]->s[i] = SUBSTITUTE|c;
					//std::cout << c;
					vector_position++;
					//std::cout << vector_position;
				}
			
			}  else{ // indel

				// deletion
				if(ref_vec[vector_position].length() > alt_vec[vector_position].length()){

					if(alt_vec[vector_position].compare("G")==0){
						c = 2;
					} else if (alt_vec[vector_position].compare("A")==0){
						c = 0;
					} else if (alt_vec[vector_position].compare("C")==0){
						c = 1;
					} else {
						c = 3;
					}

					if(homo){
						ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                    	deleting = 3;
						deletion_count = ref_vec[vector_position].length()-1;
						vector_position++;
					} else if (hap1_set){
						ret[0]->s[i] = SUBSTITUTE|c;
						deleting = 1;
						deletion_count = ref_vec[vector_position].length()-1;
						vector_position++;
					} else {
						ret[1]->s[i] = SUBSTITUTE|c;
						deleting = 2;
						deletion_count = ref_vec[vector_position].length()-1;
						vector_position++;
					}
				} else if (alt_vec[vector_position].length() > 4){
				vector_position++;
				continue;
				} // inserstion
				if(ref_vec[vector_position].length() < alt_vec[vector_position].length()){
					int num_ins = 0, ins = 0, r = 0, counter = 1;

					do {
						if(alt_vec[vector_position][counter] == 'G'){
						r = 2;
						} else if (alt_vec[vector_position][counter] == 'A'){
						r = 0;
						} else if (alt_vec[vector_position][counter] == 'C'){
						r = 1;
						} else {
						r = 3;
						}
                        counter++;
						num_ins++;
						ins = (ins << 2) | r;
                    } while (num_ins < alt_vec[vector_position].length()-1);

					if (homo) { // hom-ins
						ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
						vector_position++;

					} else if (hap1_set){
						ret[0]->s[i] = (num_ins << 12) | (ins << 4) | c;
						vector_position++;

					} else {
						ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
						vector_position++;
					}
				}
			}
		}

		if(vector_position > positions.size()){
			break;
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

void wgsim_core(const char *fn, int is_hap, uint64_t N, int dist, int std_dev, int size_l, int size_r, int mean_inner, char *vcf_file)
{
	kseq_t *ks;
    mutseq_t rseq[2];
	gzFile fp_fa;
	uint64_t tot_len, ii;
	int i, l, n_ref;
	int min_inner, max_inner, inner_dist; 																// dist is the insert size
	char *qstr;
	int size[2], Q, max_size;
	uint8_t *tmp_seq[2];
    mut_t *target;
	FILE *vcf;
	bool bool_vcf = false;
	
	l = size_l > size_r? size_l : size_r;
	qstr = (char*)calloc(l+1, 1);
	tmp_seq[0] = (uint8_t*)calloc(l+2, 1);
	tmp_seq[1] = (uint8_t*)calloc(l+2, 1);
	size[0] = size_l; size[1] = size_r;
	max_size = size_l > size_r? size_l : size_r;
	min_inner = -1 * size_l;																			//originally 10, here changed to the -size_l
	max_inner = dist - 2 * size_l;

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

	if (strcmp(vcf_file, "None") == 0 || strlen(vcf_file) == 0) {
		fprintf(stderr, "No vcf input, will generate SNP randomly if needed\n");
	} else if(vcf=fopen(vcf_file,"r")) {
		fprintf(stderr, "VCF file exists, use it in read simulation\n");
		bool_vcf=true;
		fclose(vcf);
	} else{
		fprintf(stderr, "The specified VCF file does not exist, please check\n");
	}

	fp_fa = gzopen(fn, "r");
	ks = kseq_init(fp_fa);
	while ((l = kseq_read(ks)) >= 0) { 																	//here l is the chromosome length
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

		// printf("%s\n", ks->name.s);  ks->name.s is the chrmosome
		// generate mutations and print them out
		fprintf(stdout, "Contig Variant Start\n");
		if(MUT_RATE == 0.0){fprintf(stdout, "%s\n", ks->name.s);}
		if(bool_vcf){wgsim_mut_vcf(ks, is_hap, rseq, rseq+1, vcf_file);}else{wgsim_mut_diref(ks, is_hap, rseq, rseq+1);}
		wgsim_print_mutref(ks->name.s, ks, rseq, rseq+1);
		fprintf(stdout, "Contig Variant End\n");


		for (ii = 0; ii != n_pairs; ++ii) { // the core loop
			double ran;
			int d, pos, s[2];
			int n_sub[2], n_indel[2], n_err[2], ext_coor[2], j, k;
			//FILE *fpo[2];

			pos = dis(gen);
			int insert_dev = (int)(std_dev * dis_rn(gen_rn));
			inner_dist = mean_inner + insert_dev;
			inner_dist = std::max(min_inner, std::min(inner_dist, max_inner));

			// flip or not
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
			int** c_bases = new int*[2]; 																// a pointer to a pointer to an int; why new and delete? 
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

																										//what's these symbol, why need?, in line 345, why it's M for Insertion?
			#define __gen_read(x, start, iter) do {									\
				for (i = (start), k = 0, ext_coor[x] = -10; i >= 0 && i < ks->seq.l && k < s[x]; iter) {	\
					int c = target[i], mut_type = c & mutmsk;			\
					if (ext_coor[x] < 0) {								\
						if (mut_type != NOCHANGE && mut_type != SUBSTITUTE) continue; \
						ext_coor[x] = i;                                \
					}													\
					if (mut_type == DELETE){                            \
						++n_indel[x];                                   \
						/*cigar[x] += 'D';*/                            \
						++offset[x];                                    \
						++end[x];	                                    \
					}			                                        \
					else if (mut_type == NOCHANGE || mut_type == SUBSTITUTE) { \
						char base = c & 0xf;                            \
						int ob = nst_nt4_table[(int)ks->seq.s[i]];                   \
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
						cigar[x] += 'M';								\
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
			if (ext_coor[0] < 0 || ext_coor[1] < 0) { // fail to generate the read(s)
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
	fprintf(stderr, "Forked wgsim (Heng Li) (short read simulator) for simulation of bisulfite treated reads\n");
	fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
	fprintf(stderr, "Contact: Wenbin Guo <wbguo@ucla.edu>;\n\n");
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
	fprintf(stderr, "         -g STRING     path to the genetic variant file (vcf.gz) [None]\n");
	fprintf(stderr, "\n");
	return 1;
}

int main(int argc, char *argv[])
{
	int64_t N;
	int dist, std_dev, c, size_l, size_r, is_hap, mean_inner = 0;
	char vcf_default[] = "None";
	char *vcf_file = vcf_default;
	//FILE *fpout1, *fpout2;
	int seed = -1;

	N = 1000000; dist = 500; std_dev = 50;
	mean_inner = (dist - size_l * 2) / 2; 
	size_l = size_r = 70;
	while ((c = getopt(argc, argv, "e:d:s:N:1:2:r:R:h:X:S:A:I:g:")) >= 0) {
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
		case 'g': vcf_file = optarg; break;
		}
	}
	if (argc - optind < 1) return simu_usage();
	//fpout1 = fopen(argv[optind+1], "w");
	//fpout2 = fopen(argv[optind+2], "w");
	//if (!fpout1 || !fpout2) {
	//	fprintf(stderr, "[wgsim] file open error\n");
	//	return 1;
	//}
	if (seed <= 0) seed = time(0)&0x7fffffff;
	fprintf(stderr, "[wgsim] seed = %d\n", seed);
	srand48(seed);
	wgsim_core(argv[optind], is_hap, N, dist, std_dev, size_l, size_r, mean_inner, vcf_file);

	//fclose(fpout1); fclose(fpout2);
	return 0;
}