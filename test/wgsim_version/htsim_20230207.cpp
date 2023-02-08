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
#include <algorithm>
#include "kseq.h"
#include "vcf.h"
#include "htsim.h"
KSEQ_INIT(gzFile, gzread)

#define PACKAGE_VERSION "1.0.2"
#define MAX_BIN 128 //MAX_BIN=128 to save the enrichment efficiency

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


static double ERR_RATE = 0.005;
static double MUT_RATE = 0.01;
static double INDEL_FRAC = 0.15;
static double INDEL_EXTEND= 0.3;
static double MAX_N_RATIO = 0.05;


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
//encode not as 1,3,5; have problem with print 5 (or 13) when putc
const  uint8_t cg_context_table[64] = {
    0,   0,   0,   0,    0,   0,   0,   0, 
    0,   0,   0,   0,    0,   0,   0,   0, 
    CHH, CHH, CHG, CHH,  CHH, CHH, CHG, CHH, 
    CG,  CG,  CG,  CG,	 CHH, CHH, CHG, CHH, 
    GDD, GDC, GDD, GDD,  GC,  GC,  GC,  GC, 
    GDD, GDC, GDD, GDD,  GDD, GDC, GDD, GDD,
    0,   0,   0,   0,    0,   0,   0,   0, 
    0,   0,   0,   0,    0,   0,   0,   0, 
};


//  global variables, only changed at program start.
const  uint8_t cg_table[4] = {0, 1, 1, 0}; // for cg check
int mean_insert, sd_insert, min_insert, max_insert, size_l, size_r, sd_center;

// if the leftmost 4 bit is non-zero, then it must be snp or indel
enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
static mut_t mutmsk = (mut_t)0xf000;

typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
} mutseq_t;

typedef struct {
    int pos, ref, alt, geno;
} snp_rec;

typedef struct {
    int pos_l, pos_r;
    float score;
    int8_t strand;
} probe_rec;

typedef struct {
    char *name, *chr_id;
} probe_meta;

typedef struct {
    int len = -1;           /* length of cutting site */
    int idx = -1;           /* cutting position on *seq */
    std::vector<int> seq;   /* sequence encoded by numbers*/
} cut_site;

typedef struct {
    int pos;
    int8_t type = -1;
} cut_pos;

typedef struct {
    int pos_l, pos_r;       /*int max value is 2147483647*/
    float cg_ratio= -1;
    int8_t cut_l  = -1;
    int8_t cut_r  = -1;
    int8_t index  = -1;     /*save for capture efficiency usage*/
    int8_t ns     = 0;      /*strand or # of cut sites contained*/
} fragment;                 /*each struct take <= 16 bytes*/

typedef struct {
    uint64_t tot_len;
    uint64_t eff_len;
    std::string name;
} len_rec;


std::vector<len_rec> len_vec;
std::vector<snp_rec> snp_vec;
std::vector<cut_site> site_vec;
std::vector<fragment> frag_vec;
std::vector<cut_pos>  cutpos_vec;
std::vector<probe_rec> probe_vec;
std::vector<float> eff_vec;
std::vector<float> seq_score;


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


void parse_cut_site(char *cut_str)
{
    site_vec.clear();
    int c, idx = 0;
    cut_site tmp_site;
    
    for(int i =0; cut_str[i] !='\0'; ++i){
        if (cut_str[i] == ','){tmp_site.idx = idx; continue;}
        if (cut_str[i] == ';'){
            if (tmp_site.idx >=0) {
                tmp_site.len = idx;
                site_vec.push_back(tmp_site);
                tmp_site = {}; idx = 0; continue; // empty the struct, restart
            } else {
                fprintf (stderr, "[%s] ERROR: Invalid enzyme site format (%s). Exit... \n", __func__, cut_str); exit (EXIT_FAILURE);
            }
        }
        ++idx;
        c = nst_nt4_table[(int)cut_str[i]];
        tmp_site.seq.push_back(c);
    }
    if(tmp_site.idx >=0){tmp_site.len = idx; site_vec.push_back(tmp_site);}
}

void gen_cut_pos(const kseq_t *ks)
{
    cutpos_vec.clear(); // clean the container

    std::vector<mut_t> rseq_ref(ks->seq.l); // TODO: create cut_pos without create rseq_ref

    for (int i = 0; i != ks->seq.l; ++i) {
        rseq_ref.push_back((mut_t)nst_nt4_table[(int)ks->seq.s[i]]);
    }

    auto ptr_begin = rseq_ref.begin();
    auto ptr_end   = rseq_ref.end();
    auto iter_curr = rseq_ref.begin(); // current
    auto iter_temp = rseq_ref.begin(); // temp
    auto iter_save = rseq_ref.begin(); // save
    //printf("%d %d\n", ptr_begin, ptr_end);

    int count = 0; 
    cut_pos tmp_cutpos;
    while (iter_curr < ptr_end)
    {
        // printf("========%i========\n", count);
        // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);
        //initiate
        int min_pos = ptr_end - ptr_begin;
        int site_pos= min_pos;
        int type_idx = -1;

        // printf("[%d %d %d]\n", min_pos, site_pos, type_idx);
        for (int j = 0; j < site_vec.size(); ++j) {
            iter_temp = std::search(iter_curr, ptr_end, site_vec[j].seq.begin(), site_vec[j].seq.end());
            // printf("%d %d %d\n", iter_curr, iter_temp);
            site_pos = iter_temp - ptr_begin;
            if(site_pos < min_pos){
                min_pos = site_pos;
                type_idx= j;
                iter_save= iter_temp;
            }
        }
        // printf("[%d %d %d]\n", min_pos, site_pos, type_idx);
        // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);

        if(type_idx >=0){
            tmp_cutpos = {};
            tmp_cutpos.pos = min_pos + site_vec[type_idx].idx;
            tmp_cutpos.type= type_idx;
            cutpos_vec.push_back(tmp_cutpos);
            iter_curr = iter_save + site_vec[type_idx].len; // need to check for other int types
        }
        // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);
        // ++count;
        // if(count >= 20){
        //     break;
        // }
    }
}

void gen_cut_frag(const kseq_t *ks)
{
    frag_vec.clear();
    // generate potential intervals
    fragment tmp_frag;
    for (int i = -1; i <= int(cutpos_vec.size()); ++i){
        if(i==-1 || i == int(cutpos_vec.size())){ //check the first and last fragments
            tmp_frag = {};
            if(i==-1){
                int frag_len = cutpos_vec[0].pos;
                if (frag_len > min_insert && frag_len < max_insert){
                    tmp_frag.pos_l = 0;
                    tmp_frag.cut_l = -1;
                    tmp_frag.pos_r = cutpos_vec[0].pos;
                    tmp_frag.cut_r = cutpos_vec[0].type;
                    tmp_frag.ns = 0;
                    frag_vec.push_back(tmp_frag);
                }
            }else{
                int frag_len = ks->seq.l - cutpos_vec[i-1].pos;
                if (frag_len > min_insert && frag_len < max_insert){
                    tmp_frag.pos_l = cutpos_vec[i-1].pos;
                    tmp_frag.cut_l = cutpos_vec[i-1].type;
                    tmp_frag.pos_r = ks->seq.l;
                    tmp_frag.cut_r = -1;
                    tmp_frag.ns = 0;
                    frag_vec.push_back(tmp_frag);
                }
            }
            continue;
        }
        for(int j=i+1; j < cutpos_vec.size(); ++j){
            int frag_len = cutpos_vec[j].pos - cutpos_vec[i].pos;
            if (frag_len > min_insert && frag_len < max_insert){
                tmp_frag = {};
                tmp_frag.pos_l = cutpos_vec[i].pos;
                tmp_frag.cut_l = cutpos_vec[i].type;
                tmp_frag.pos_r = cutpos_vec[j].pos;
                tmp_frag.cut_r = cutpos_vec[j].type;
                // calculate the probability and index
                tmp_frag.ns = (int8_t)j-(i+1);
                frag_vec.push_back(tmp_frag);
            }
        }
    }
}

void parse_bed_line(char *line, char *chr_id, probe_rec *tmp_probe, probe_meta *tmp_probe_meta)
{
	char *p, *q, *name = 0;
    int i, start, end, strand;
	float score;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: chr_id = p; break;
            case 1: start  = atoi(p); break;
            case 2: end    = atoi(p); break;
            case 3: name   = strdup(p); break;
            case 4: score  = atof(p); break; // what if score is .
            case 5: strand = int(strcmp(p,"+")==0)-int(strcmp(p,"-")==0); break;
            default: break;}
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

	tmp_probe->pos_l = start;
	tmp_probe->pos_r = end;
    tmp_probe->score = score;
    tmp_probe->strand= strand;
    
    tmp_probe_meta->chr_id = chr_id;
    tmp_probe_meta->name   = name;
	if(i < 4){chr_id = 0;}
}

void parse_bed_chr(char *fname, char *chr_id)
{
    probe_vec.clear();  // clean the container

    htsFile *fp = hts_open(fname,"rb");
    if(fp == 0 ){ fprintf(stderr,"[%s] ERROR: open bed file: %s failed: %s. Exit...\n", __func__, fname, strerror (errno)); exit (EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    probe_rec tmp_probe;  
    probe_meta tmp_probe_meta;

    int ret;
    char *chr_current;
    bool collect_present = false;
    bool collect_previous= false;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (collect_present == false && collect_previous == true) { break; } // finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        parse_bed_line(line.s, chr_current, &tmp_probe, &tmp_probe_meta); //might need to test

        if (strcmp(chr_current, chr_id) != 0){
            collect_present = false;
            continue;
        } else {
            if (!chr_current){fprintf(stderr, "Skip invalid probe: chr %s, name %s...\n", tmp_probe_meta.chr_id, tmp_probe_meta.name);}
            probe_vec.push_back(tmp_probe);
            tmp_probe      = {};
            tmp_probe_meta = {};
        }
    }
    free(line.s);

    if (ret=hts_close(fp)){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
}

void parse_eff_file(char *fname)
{
    FILE* fp = fopen(fname, "r");
    float eff_prob;

    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open capture efficiency file: %s failed. Exit...", __func__, fname); exit (EXIT_FAILURE);}
    while(fscanf(fp, "%f", &eff_prob) == 1) { // this will skip the empty lines
        // printf("%f\n", eff_prob);
        eff_vec.push_back(eff_prob);
    }
    fclose(fp);
}


void parse_vcf_chr(char *fname, char *chr_id)
{
    snp_vec.clear(); //clean the container

    //open vcf file
    htsFile   *fp  = hts_open(fname,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t    *rec = bcf_init();
    
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
        fprintf(stderr, "[%s] ERROR: Currently only support single-sample simulation, please check the vcf file! Exiting...\n", __func__);
        exit(EXIT_FAILURE);
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
                    for (int i = 0; i < alt_len; ++i ){ 
                        c = (mut_t)nst_nt4_table[(int)alt[i]];
                        alt_int = (alt_int << 2) | c;
                    }
                } else if (ref_len > 1 && alt_len == 1){ //deletion
                    alt_int = (mut_t)nst_nt4_table[(int)alt[0]];
                    base_change_pos = snp_pos + 1; // position +1, because the base change occurs after the first base
                    for (int i = 0; i < ref_len; ++i ){ 
                        c = (mut_t)nst_nt4_table[(int)ref[i]];
                        ref_int = (ref_int << 2) | c;
                    }
                } else { // might be MNP, or sth else
                    fprintf(stderr, "[%s] Skip unusual SNP/INDEL [MNP]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
                    continue;
                }
            }

            int geno_int = base_offset << 12 | ref_len << 8 | snp_hap2 << 6 | snp_hap1 << 4 | is_phased;
            snp_rec tmp_snp = {.pos = base_change_pos, .ref = ref_int, .alt = alt_int, .geno = geno_int};
            snp_vec.push_back(tmp_snp);
            tmp_snp = {};
        }
    }
    
    fprintf(stderr, "[%s] Finish collecting %lu SNP/INDEL from %s\n", __func__, snp_vec.size(), chr_id);

    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);

    if (int ret=hts_close(fp)){
        fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret); 
    }
}

void sim_mut_vcf(const kseq_t *ks, char * vcf_file, mutseq_t *hap1, mutseq_t *hap2, uint8_t *site_flag_arr)
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
        if (snp_vec.size() == 0){ continue;} // ignore the rest if there is no SNP

        if (deleting){
            if(deletion_count > 0){
                if (deleting & 1){ ret[0]->s[i] |= DELETE;}
                if (deleting & 2){ ret[1]->s[i] |= DELETE;}
                deletion_count--;
                site_flag_arr[i] = 1;
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

void sim_mut_diref(const kseq_t *ks, bool is_hap, mutseq_t *hap1, mutseq_t *hap2, uint8_t *site_flag_arr)
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

void cal_length_chr(kseq_t *ks, len_rec *tmp_len, int tech_mode, char *bed_file)
{
    uint64_t eff_len;

    if(tech_mode==2){                       // targeted sequencing
        parse_bed_chr(bed_file, ks->name.s);
        int pos_l, pos_r, len;
        int pos_l_prev=0, pos_r_prev=0;
        for (size_t i = 0; i < probe_vec.size(); ++i){
            pos_l = probe_vec[i].pos_l;
            pos_r = probe_vec[i].pos_r;
            len   = pos_r - pos_r;
            // deal with overlap probes
            if ((pos_l < pos_r_prev) && (pos_l > pos_l_prev)){ 
                len = (pos_r > pos_r_prev) ? (pos_r - pos_r_prev) : 0;
            }
            eff_len += len;
            pos_l_prev = pos_l;
            pos_r_prev = pos_r;
        }
    }else if (tech_mode==1){                // reduced representative sequencing 
        gen_cut_pos(ks);
        int len;
        for (int i = -1 ; i <= (int)cutpos_vec.size(); ++i){
            if(i == -1 || i == (int)cutpos_vec.size()){
                if(i==-1){
                    len = cutpos_vec[0].pos - 0;
                }else{
                    len = (int) ks->seq.l - cutpos_vec[i-1].pos;
                }
            }else{
                len = cutpos_vec[i+1].pos - cutpos_vec[i].pos;
            }
            len = len >= min_insert && len <= max_insert? len : 0;
            eff_len += len;
        }
    }else{                                  // whole genome
        eff_len = ks->seq.l;
    }
    tmp_len->tot_len = ks->seq.l;
    tmp_len->eff_len = eff_len;
    tmp_len->name = std::string(ks->name.s);
}

void sim_rand(std::uniform_int_distribution<int> *dis, fragment *tmp_frag, mut_t *target, bool is_uniform, int tot_size)
{
    // for whole genome
    int pos_l, pos_r, insert_dev, insert_len;
    // decide the left and right-most positons of the insert fragment
    if(is_uniform){
        pos_l = (*dis)(gen);
        insert_dev = (int)(sd_insert * dis_rn(gen_rn));
        insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));
        pos_r = std::min(pos_l + insert_len, tot_size -2);
    }else{
        bool flag = true;
        while (flag) {
            pos_l = (*dis)(gen);
            insert_dev = (int)(sd_insert * dis_rn(gen_rn));
            insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));
            pos_r = std::min(pos_l + insert_len, tot_size -2);
            int cg_count =0;
            for(int kk = pos_l; kk <= pos_r; ++kk){
                cg_count += cg_table[(uint8_t)(target[kk]& 0x3)]; // not perfect, but fine enough
            }
            float cg_prob = eff_vec[(int)(cg_count*MAX_BIN/(pos_r - pos_l)+0.5)];
            flag = dis_ru(gen_ru) > cg_prob; // when initiate eff_vec, judge the value in case it's too small
        }
    }
    tmp_frag->pos_l= pos_l;
    tmp_frag->pos_r= pos_r;
}

void sim_rand(std::uniform_int_distribution<int> *dis, fragment *tmp_frag, bool is_uniform)
{
    // for Reduced representative
    int frag_idx;
    if(is_uniform){
        frag_idx = (*dis)(gen);
        *tmp_frag = frag_vec[frag_idx];
    }else{
        bool flag = true;
        while (flag) {
            frag_idx = (*dis)(gen);
            *tmp_frag = frag_vec[frag_idx];
            float cap_prob = eff_vec[tmp_frag->index]; // make sure index is 0-127
            flag = dis_ru(gen_ru) > cap_prob;
        }
    }
}

void sim_rand(std::uniform_int_distribution<int> *dis, fragment *tmp_frag)
{
    // targeted sequencing (uniform)
    int frag_idx = (*dis)(gen);
    probe_rec tmp_probe = probe_vec[frag_idx];
    tmp_frag->ns = tmp_probe.strand;
    int probe_center = (int) (tmp_probe.pos_l + tmp_probe.pos_r)/2;
    int frag_center= probe_center + (int)(sd_center * dis_rn(gen_rn));
    int insert_dev = (int)(sd_insert * dis_rn(gen_rn));
    int insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));

    int pos_l = frag_center - (int)insert_len/2;
    int pos_r = frag_center + (int)insert_len/2; 
    tmp_frag->pos_l = pos_l;
    tmp_frag->pos_r = pos_r;
}

void sim_rand(std::discrete_distribution<int> *dis, fragment *tmp_frag)
{
    // targeted sequencing (non-uniform)
    int frag_idx = (*dis)(gen);
    probe_rec tmp_probe = probe_vec[frag_idx];
    tmp_frag->ns = tmp_probe.strand;
    int probe_center = (int) (tmp_probe.pos_l + tmp_probe.pos_r)/2;
    int frag_center= probe_center + (int)(sd_center * dis_rn(gen_rn));
    int insert_dev = (int)(sd_insert * dis_rn(gen_rn));
    int insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));

    int pos_l = frag_center - (int)insert_len/2;
    int pos_r = frag_center + (int)insert_len/2; 
    tmp_frag->pos_l = pos_l;
    tmp_frag->pos_r = pos_r;
}

void sim_core(const char *fn, bool is_hap, bool is_uniform, uint64_t N, uint64_t chr_N, uint64_t tot_eff_len, int tech_mode, int output_fmt, char *vcf_file, char *probe_bed, char *chr_id)
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

    bool bool_chr_set  = strcmp(chr_id, "None") && strlen(chr_id);
    bool bool_vcf_set  = strcmp(vcf_file,"None") && strlen(vcf_file);
    bool bool_chrN_set = chr_N > 0;

    // start simulate
    fp_fa = gzopen(fn, "r");
    ks = kseq_init(fp_fa);

    while ((l = kseq_read(ks)) >= 0) {  //here l is the chromosome length
        if (bool_chr_set) {if (strcmp(chr_id, ks->name.s)!=0){continue;}}
        if (l < mean_insert + 3 * sd_insert) {
            fprintf(stderr, "[%s] skip sequence '%s' as it is shorter than %d!\n", __func__, ks->name.s, mean_insert + 3 * sd_insert);
            continue;
        }
        if (bool_chrN_set) { n_pairs = chr_N; }else{
            for (size_t i = 0; i < len_vec.size(); ++i){
                if (strcmp(len_vec[i].name.c_str(), ks->name.s)==0){
                    contig_eff_len = len_vec[i].eff_len; break;}
            }
            n_pairs = (uint64_t)(contig_eff_len * N / tot_eff_len + 0.5);
        }
        if(!n_pairs){continue;} // skip if no reads simulated from this contig

        // initialize random number generator to generate read positions
        uint64_t unif_begin, unif_end;
        std::discrete_distribution<int> dis2;
        std::uniform_int_distribution<int> dis;

        if(tech_mode==2){
            parse_bed_chr(probe_bed, ks->name.s);
            if(probe_vec.size() == 0){continue;} // parse probe, skip if empty
            if(is_uniform){
                std::uniform_int_distribution<int> dis(0, probe_vec.size()-1);
            }else{
                std::vector<float> weights(probe_vec.size());
                for(int i=0; i < probe_vec.size(); ++i){ weights.push_back(probe_vec[i].score);}
                std::discrete_distribution<int> dis2(weights.begin(), weights.end());
                // std::vector<float>().swap(weights);
            }
        }else if(tech_mode==1){
            gen_cut_pos(ks); gen_cut_frag(ks); 
            if(frag_vec.size()==0){continue;}   // generate frags, skip if empty
            std::uniform_int_distribution<int> dis(0, frag_vec.size()-1);
        }else{
            // ensure read doesn't pass boundary with 2 base offset
            std::uniform_int_distribution<int> dis(2, ks->seq.l-max_insert-2);
        }

        // print out simulation information
        tot_pairs += n_pairs;
        fprintf(stderr, "[%s] contig '%s': simulate %ld reads...\n", __func__, ks->name.s, n_pairs);

        // introduce mutations and print them to stdout
        uint8_t site_flag_arr[ks->seq.l] = {0}; // record if a site is SNP/INDEL position (base can be either REF/ALT)        
        fprintf(stdout, "Contig Variant Start\n");
        if(bool_vcf_set){
            sim_mut_vcf(ks, vcf_file, rseq, rseq+1, site_flag_arr);
            if(snp_vec.size() == 0){fprintf(stdout, "%s\n", ks->name.s);}   //if no variants, print chromosome id
        } else {
            sim_mut_diref(ks, is_hap, rseq, rseq+1, site_flag_arr);         // TODO: should return an int
            if(MUT_RATE == 0.0){fprintf(stdout, "%s\n", ks->name.s);}       //if no variants, print chromosome id            
        }
        sim_print_mutref(ks->name.s, ks, rseq, rseq+1);
        fprintf(stdout, "Contig Variant End\n");


        fragment tmp_frag;
        for (ii = 0; ii != n_pairs; ++ii) { // the core loop
            tmp_frag ={};
            int start_2 =0;
            int n_sub[2]={0,0}, n_indel[2]={0,0}, n_err[2]={0,0}, cover_pos[2]={0,0}; 
            int ext_coor[2], i, j, k, ix;
            //cover_pos hold if the read covers a snp *position* (the read don't have to contain the ALT allele)
            //j hold read1/read2, k hold the length of read, ix hold the cursor transversing read

            // random position generation
            target = rseq[drand48()<0.5?0:1].s; // haplotype from which the reads are generated
            if(tech_mode==2){
                sim_rand(&dis2, &tmp_frag);
                start_2 = tmp_frag.pos_r - size_l;
            }else if(tech_mode==1){
                sim_rand(&dis, &tmp_frag, is_uniform);
                start_2 = tmp_frag.pos_r;
                // find out the start position for read2, also check if the boundaries satisfy
                for(int k=0; k < size_r; k++){
                    int c = target[start_2], mut_type = c & mutmsk;
                    if(mutmsk == DELETE){
                        --start_2;
                        --k;
                    }else if(mutmsk == INSERT){
                        int num_ins = mut_type>>12;
                        if(k + num_ins > size_r){
                            start_2 -= size_r - k;
                        }else{
                            start_2 -= num_ins;
                            k += num_ins;
                        }
                    }else{
                        --start_2;
                    }
                }
            }else{
                sim_rand(&dis, &tmp_frag, target, is_uniform, ks->seq.l);
                start_2 = tmp_frag.pos_r - size_l;
            }

            // generate the read sequences
            int start[2] = {tmp_frag.pos_l, start_2};
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
                    cover_pos[x] |= site_flag_arr[i];                       \
                }                                                           \
                /* append CG context flag;                                  \
                   currently not handling bounday context by indel*/        \
                for (ix=0; ix < k; ++ix) {                                  \
                    int c_d1, c_d2;                                         \
                    int c = tmp_seq[x][ix];                                 \
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
                        tmp_context[x][ix] = cg_context_table[context_idx]; \
                    }                                                       \
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
            
            int flag_pos, flag_mut;
            for(i=tmp_frag.pos_l; i < tmp_frag.pos_r; ++i){
                flag_pos |= site_flag_arr[i];
                flag_mut |= target[i] & mutmsk;
            }
            int flag_indel= n_indel[0] | n_indel[1];
            
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
                    fprintf(stdout, "@%s:%d:%d:%llx %d %d %d %d %d ", ks->name.s, start[0]+1, end[1]+1, (long long)ii, j, flag_pos, flag_mut, flag_indel, Q); 
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
        free(rseq[0].s); free(rseq[1].s);
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
    fprintf(stderr, "         -i INT        mean insert size (outer distance between 2 ends) [500]\n");
    fprintf(stderr, "         -I INT        standard deviation of insert size [50]\n");
    fprintf(stderr, "         -m INT        minimum insert size [100]\n");
    fprintf(stderr, "         -M INT        maximum insert size [1000]\n");
    fprintf(stderr, "         -c STRING     contig name, only output reads from this contig, default is to output all contigs [None]\n");
    fprintf(stderr, "         -n INT        number of read pairs to generate for specified contig, disabled by default [0]\n");
    fprintf(stderr, "         -N INT        total number of read pairs to generate, will caculate from depth when N is 0 [0]\n");
    fprintf(stderr, "         -d INT        average sequencing depth, only used when n/N is not specified [30]\n");
    fprintf(stderr, "         -1 INT        length of the first read [70]\n");
    fprintf(stderr, "         -2 INT        length of the second read [70]\n");
    fprintf(stderr, "         -E FLOAT      base error rate (set to be 0 for bisulfite sequencing) [%.3f]\n", ERR_RATE);
    fprintf(stderr, "         -A FLOAT      disgard if the fraction of ambiguous bases higher than FLOAT [%.2f]\n", MAX_N_RATIO);
    fprintf(stderr, "         -u INT        uniform coverage: 0 for No, nonzero for Yes [1]\n");
    fprintf(stderr, "         -f INT        output format: 0 for letters; 1 for ascii numbers (for Bisulfite simulaitons) [0]\n");
    fprintf(stderr, "mutation setting:\n");
    fprintf(stderr, "         -g STRING     path to the genetic variant file (vcf/vcf.gz) [None]\n");
    fprintf(stderr, "         -r FLOAT      rate of mutations [%.4f]\n", MUT_RATE);
    fprintf(stderr, "         -R FLOAT      fraction of indels [%.2f]\n", INDEL_FRAC);
    fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", INDEL_EXTEND);
    fprintf(stderr, "         -h INT        haplotype mode: 0 for No, nonzeroz for Yes (all variants are homozygotes) [0]\n");
    fprintf(stderr, "         -s INT        seed for random generator [-1]\n");
    fprintf(stderr, "technology setting:\n");
    fprintf(stderr, "         -T INT        technology: 0 for Whole genome; 1 for Reduced representation; 2 for Targeted [0]\n");
    fprintf(stderr, "         -x STRING     enzyme cutting site string for reduced representation sequencing [None]\n");
    fprintf(stderr, "         -e STRING     enrichment profile lookup file for reduced representation sequencing [None]\n");
    fprintf(stderr, "         -b STRING     probe BED file for targeted sequencing [None]\n");
    fprintf(stderr, "         -D INT        fragment center's deviaiton from the probe center [50]\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[])
{
    uint64_t N=0, chr_N=0;
    int c = 0;
    int tech_mode=0, output_fmt = 0;
    int seed = -1;
    int min_insert, max_insert, depth = 20;
    bool is_hap = false, is_uniform = true;

    char none_default[] = "None";
    char *chr_id    = none_default;
    char *vcf_file  = none_default;
    char *cut_str   = none_default;
    char *enrich_tsv= none_default;
    char *probe_bed = none_default; // checked, will not intefere with vcf_file

    mean_insert = 500; sd_insert = 50; sd_center = 50; size_l = size_r = 70;
    while ((c = getopt(argc, argv, "i:I:m:M:c:n:N:d:1:2:E:A:U:f:g:r:R:X:h:s:T:x:e:b:D:")) >= 0) {
        switch (c) {
            case 'i': mean_insert= atoi(optarg); break;
            case 'I': sd_insert  = atoi(optarg); break;
            case 'm': min_insert = atoi(optarg); break;
            case 'M': max_insert = atoi(optarg); break;
            case 'c': chr_id     = optarg; break;
            case 'n': chr_N   = atoi(optarg); break;
            case 'N': N          = atoi(optarg); break;
            case 'd': depth      = atoi(optarg); break;
            case '1': size_l     = atoi(optarg); break;
            case '2': size_r     = atoi(optarg); break;
            case 'E': ERR_RATE   = atof(optarg); break;
            case 'A': MAX_N_RATIO= atof(optarg); break;
            case 'u': is_uniform = atoi(optarg)!=0; break;
            case 'f': output_fmt = atoi(optarg); break;
            case 'g': vcf_file   = optarg; break;
            case 'r': MUT_RATE   = atof(optarg); break;
            case 'R': INDEL_FRAC = atof(optarg); break;
            case 'X': INDEL_EXTEND=atof(optarg); break;
            case 'h': is_hap     = atoi(optarg)!=0; break;
            case 's': seed       = atoi(optarg); break;
            case 'T': tech_mode  = atoi(optarg); break;
            case 'x': cut_str    = optarg; break;
            case 'e': enrich_tsv = optarg; break;
            case 'b': probe_bed  = optarg; break;
            case 'D': sd_center  = atoi(optarg); break;
        }
    }
    if (argc - optind < 1) return simu_usage();
    if (seed <= 0) seed = time(0)&0x7fffffff;

    bool bool_chr_set   = strcmp(chr_id, "None")    && strlen(chr_id);
    bool bool_vcf_set   = strcmp(vcf_file,"None")   && strlen(vcf_file);
    bool bool_site_set  = strcmp(cut_str, "None")   && strlen(cut_str);
    bool bool_probe_set = strcmp(probe_bed, "None") && strlen(probe_bed);
    bool bool_enrich_set= strcmp(enrich_tsv,"None") && strlen(enrich_tsv); 

    // check legal input mode and corresponding files
    if (tech_mode==2 || bool_probe_set){
        fprintf(stderr, "Simulating targeted sequencing reads:\n");
        if (bool_site_set){fprintf(stderr, "WARN: Cut site is set, ignored in this mode\n");}
        if (!bool_probe_set){fprintf(stderr, "ERROR: Please specify probe bed file path\n");exit(EXIT_FAILURE);}
        tech_mode = 2;
    } else if (tech_mode==1 || bool_site_set){
        fprintf(stderr, "Simulating restricted enzyme cutting reads:\n");
        if (bool_probe_set){fprintf(stderr, "WARN: Probe file is set, ignored in this mode\n");}
        if (!bool_site_set){fprintf(stderr, "ERROR: Please specify enzyme cutting site\n");exit(EXIT_FAILURE);}
        tech_mode = 1;
        parse_cut_site(cut_str);
    } else { 
        fprintf(stderr, "Simulating whole genome reads:\n"); 
        tech_mode = 0;
    }

    // check enrichment file
    if(!is_uniform || bool_enrich_set){
        if (!bool_enrich_set){fprintf(stderr, "ERROR: Please specify enrichment files when specifying -u as 0\n");exit(EXIT_FAILURE);}
        parse_eff_file(enrich_tsv);
        is_uniform = false;
    }

    // check existence of fasta, parse the length
    kseq_t *ks;
    gzFile fp_fa;
    fp_fa = gzopen(argv[optind], "r");
    ks = kseq_init(fp_fa);
    if (!fp_fa) { fprintf (stderr, "ERROR: gzopen of '%s' failed: %s. Exit... \n", argv[optind], strerror (errno)); exit (EXIT_FAILURE);}
    fprintf(stderr, "Reference genome file: %s\n", argv[optind]);

    fprintf(stderr, "[%s] Calculating the total length and effective length of the reference sequences...\n", __func__);
    int l;
    len_rec tmp_len;
    uint64_t tot_len = 0, tot_eff_len = 0;
    while ((l = kseq_read(ks)) >= 0) {
        tmp_len = {};
        cal_length_chr(ks, &tmp_len, tech_mode, probe_bed);
        tot_len += tmp_len.tot_len;
        tot_eff_len += tmp_len.eff_len;
        len_vec.push_back(tmp_len);
    }
    kseq_destroy(ks);
    gzclose(fp_fa);

    // check if fasta is empty
    if (!len_vec.size()) { fprintf (stderr, "ERROR: Input fasta is empty: %s. Exit... \n", argv[optind]); exit (EXIT_FAILURE);}

    // check input chr_id
    if (bool_chr_set){
        bool bool_chr_found = false;
        uint64_t contig_len=0, contig_eff_len=0;
        for (size_t i =0; i < len_vec.size(); ++i){
            if (strcmp(chr_id, len_vec[i].name.c_str())==0){ bool_chr_found = true; contig_len = len_vec[i].tot_len; contig_eff_len = len_vec[i].eff_len;}
        }
        if(!bool_chr_found){fprintf(stderr, "ERROR: Contig id '%s' is not found in the fasta, please check!\n", chr_id); exit(EXIT_FAILURE);}
        N = N == 0? (contig_eff_len * depth)/(size_l + size_r) : N;
        fprintf(stderr, "[%s] Contig %s specified, contig length: %lu, effective length: %lu\n", __func__, chr_id, contig_len, contig_eff_len);
    } else {
        if (chr_N > 0) {fprintf(stderr, "ERROR: -n is specified but not -c. Exit... (please note the difference of -n and -N)\n"); exit(EXIT_FAILURE);}
        fprintf(stderr, "[%s] Found %d contig sequences, total length: %lu, effective length: %lu\n", __func__, (int)len_vec.size(), tot_len, tot_eff_len);
        N = N == 0? (tot_eff_len * depth)/(size_l + size_r) : N;
        fprintf(stderr, "[%s] No contig id specified, will generate %lu reads from all contigs\n", __func__, N);
    }

    // check input vcf file
    FILE *vcf;
    if (bool_vcf_set) {
        if(vcf=fopen(vcf_file,"r")){
            fprintf(stderr, "[%s] VCF file exists, use it to simulate reads\n", __func__); fclose(vcf);
        }else{
            fprintf(stderr, "ERROR: The specified VCF file does not exist, please check!\n"); exit(EXIT_FAILURE);
        }
    } else {
        fprintf(stderr, "[%s] No VCF input, will generate SNP randomly if mutation rate is nonzero\n", __func__);
    }

    fprintf(stderr, "[htsim] seed = %d\n", seed);
    srand48(seed);

    sim_core(argv[optind], is_hap, is_uniform, N, chr_N, tot_eff_len, tech_mode, output_fmt, vcf_file, probe_bed, chr_id);

    return 0;
}

