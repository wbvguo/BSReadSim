/* The MIT License

   Copyright (c) 2022 Wenbin Guo <wbguo@ucla.edu>

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

/* Last Modified: 20MAR2023 */


#ifndef HTSIM_H
#define HTSIM_H

#include <vector>
#include <map>
#include <string>
#include <gsl/gsl_randist.h>
#include "kseq.h"
#include "vcf.h"

/*--------------------------------------------------struct--------------------------------------------------*/
// parse VCF
typedef struct {
    int pos, ref, alt, geno;
} snp_rec;


// parse BED
typedef struct {
    int pos_l, pos_r;
    float score;
    int8_t strand;
} probe_rec;

typedef struct {
    char *name, *chr_id;
} probe_meta;


// parse RRBS
typedef struct {
    int len = -1;           /* length of cutting site */
    int idx = -1;           /* cutting position on *seq */
    std::vector<int> seq;   /* sequence encoded by numbers*/
} cut_rec;

typedef struct {
    int pos;
    int8_t type = -1;
} cut_pos;


// generate read
typedef struct {
    int pos_l,pos_r,start2; /*int max value is 2147483647*/
    int8_t cut_l  = -1;
    int8_t cut_r  = -1;
    int8_t haplo  = -1;     /*haplotype*/
    int8_t strand = -1;     /*strand*/
} frag_rec;                 /*each struct take <= 16 bytes*/


// save for contigs
typedef struct {
    uint32_t chr_len= 0;
    uint32_t eff_len= 0;
    uint32_t count  = 0;
    float   score = 0;
} chr_rec;

// save for methylable bases
typedef struct {
    int pos = -1;
    float ref_meth = -1;    /*save meth if not is_asm*/
    float alt_meth = -1;
    uint8_t context= 0;     /*1,3,7;9,11,15*/
    uint8_t type   = 0;     /*0,1,2,4,8*/
} meth_rec;                 /*each struct take 14 bytes*/

typedef struct {
    float alpha = -1;
    float beta = -1;
} param_rec;                /*take 8 bytes*/


param_rec param_cg = {0.5,0.5};
param_rec param_chg= {0.5,0.5};
param_rec param_chh= {0.5,0.5};
std::map<std::string, int> base_map = {{"C", 0}, {"G",1}};
std::map<std::string, int> context_map = {{"CG",1}, {"CHG",3}, {"CHH",7}};
std::map<int, int> params_map = {{1,0}, {3,1}, {7,2}, {9,0}, {11,1}, {15,2}}; //idx in param_vec

const  uint8_t cg_table[5] = {0, 1, 1, 0, 0}; // for C/G check
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

/*--------------------------------------------------function--------------------------------------------------*/
// for VCF
inline void parse_vcf_chr(char *fname, char *chr_id, std::vector<snp_rec>& snp_vec)
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


// for BED
inline int parse_bed_line(char *line, char *chr_id, probe_rec *tmp_probe, probe_meta *tmp_probe_meta)
{
	char *p, *q, *contig, *name = 0;
    int i, start, end, strand;
	float score;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: contig = p; break;
            case 1: start  = atoi(p); break;
            case 2: end    = atoi(p); break;
            case 3: name   = strdup(p); break;
            case 4: score  = atof(p); break;        // TODO: need to test what if score is "." 
            case 5: strand = int(strcmp(p,"+")==0)-int(strcmp(p,"-")==0); break; //1,0,-1
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return 1;} // termenate early if not equal
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

    if(i < 4){fprintf(stderr, "[%s] Skip invalid probe: chr %s, name %s...\n", __func__, chr_id, name); return 0;}

	tmp_probe->pos_l = start;
	tmp_probe->pos_r = end;
    tmp_probe->score = score;
    tmp_probe->strand= strand;
    
    tmp_probe_meta->chr_id = chr_id;
    tmp_probe_meta->name   = name;
    return 0;
}

inline void parse_bed_chr(char *fname, char *chr_id, std::vector<probe_rec>& probe_vec)
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
        int skip_status = parse_bed_line(line.s, chr_id, &tmp_probe, &tmp_probe_meta); //might need to test
        if(skip_status){collect_present = false; continue;}
        probe_vec.push_back(tmp_probe);
        tmp_probe      = {};
        tmp_probe_meta = {};
    }
    free(line.s);

    if (ret=hts_close(fp)){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
}


// for RRBS
inline void parse_cut_rec(char *cut_str, std::vector<cut_rec>& cut_vec)
{
    cut_vec.clear();
    int c, idx = 0;
    cut_rec tmp_cut;
    
    for(int i =0; cut_str[i] !='\0'; ++i){
        if (cut_str[i] == '|'){tmp_cut.idx = idx; continue;}
        if (cut_str[i] == ','){
            if (tmp_cut.idx >=0) {
                tmp_cut.len = idx;
                cut_vec.push_back(tmp_cut);
                tmp_cut = {}; idx = 0; continue; // empty the struct, restart
            } else {
                fprintf (stderr, "[%s] ERROR: Invalid enzyme site format (%s). Exit... \n", __func__, cut_str); exit (EXIT_FAILURE);
            }
        }
        ++idx;
        c = nst_nt4_table[(int)cut_str[i]];
        tmp_cut.seq.push_back(c);
    }
    if(tmp_cut.idx >=0){tmp_cut.len = idx; cut_vec.push_back(tmp_cut);}
}

inline void gen_cut_pos(const kseq_t *ks, std::vector<cut_pos>& cutpos_vec, std::vector<cut_rec>& cut_vec)
{
    cutpos_vec.clear(); // clean the container

    std::vector<mut_t> rseq_ref(ks->seq.l); // TODO: create cut_pos without create rseq_ref

    for (int i = 0; i != ks->seq.l; ++i) {
        rseq_ref.push_back((mut_t)nst_nt4_table[(int)ks->seq.s[i]]);
    }

    std::vector<mut_t>::iterator ptr_begin = rseq_ref.begin();  // can repalce the type with `auto`
    std::vector<mut_t>::iterator ptr_end   = rseq_ref.end();
    std::vector<mut_t>::iterator iter_curr = rseq_ref.begin();  // current
    std::vector<mut_t>::iterator iter_temp = rseq_ref.begin();  // temp
    std::vector<mut_t>::iterator iter_save = rseq_ref.begin();  // save
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
        for (int j = 0; j < cut_vec.size(); ++j) {
            iter_temp = std::search(iter_curr, ptr_end, cut_vec[j].seq.begin(), cut_vec[j].seq.end());
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
            tmp_cutpos.pos = min_pos + cut_vec[type_idx].idx;
            tmp_cutpos.type= type_idx;
            cutpos_vec.push_back(tmp_cutpos);
            iter_curr = iter_save + cut_vec[type_idx].len; // need to check for other int types
        }
        // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);
        // ++count;
        // if(count >= 20){
        //     break;
        // }
    }
}

//gen cut fragment

// create/stream MethDB
inline void create_methdb(const kseq_t *ks, uint32_t *pos_idx_arr, std::vector<meth_rec>& meth_vec){
    meth_vec.clear();               // clean the container
    meth_rec tmp_meth;
    meth_vec.push_back(tmp_meth);   // put the unitialized meth_rec into as the first element

    int c, c_d1, c_d2;
    uint32_t ix = 1;
    int k = ks->seq.l;
    for (int i = 0; i != ks->seq.l; ++i) {
        c = nst_nt4_table[(int)ks->seq.s[i]];
        if (cg_table[(uint8_t) c]){
            if(c==1){
                if(i > k-3){c_d1 =0; c_d2 =0;}  /*handle the last 2 base*/
                else{
                    c_d1 = nst_nt4_table[(int)ks->seq.s[i + 1]];
                    c_d2 = nst_nt4_table[(int)ks->seq.s[i + 2]];
                }
            }else{
                if(i < 2){c_d1 =0; c_d2 =0;}    /*handle the first 2 base*/
                else{
                    c_d1 = nst_nt4_table[(int)ks->seq.s[i - 1]];
                    c_d2 = nst_nt4_table[(int)ks->seq.s[i - 2]];
                }
            }
            uint8_t context_idx = c << 4 | c_d1 <<2 | c_d2;
            tmp_meth.context = cg_context_table[context_idx];
            tmp_meth.pos = i;
            meth_vec.push_back(tmp_meth);
            pos_idx_arr[i] |= ix << 2;
            ix++;
            tmp_meth = {};
        }
    }
}

inline void save_methdb(std::vector<meth_rec>& meth_vec, char *fname){
    FILE* fp = fopen(fname, "w");
    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open methdb file: %s failed. Exit...", __func__, fname); exit (EXIT_FAILURE);}

    meth_rec tmp_meth;
    for (size_t i=0; i < meth_vec.size(); ++i){
        tmp_meth = meth_vec[i];
        fprintf(fp, "%d\t%f\t%f\t%d\t%d", tmp_meth.pos, tmp_meth.ref_meth, tmp_meth.alt_meth, tmp_meth.context, tmp_meth.type);
    }
    fclose(fp);
}

inline void parse_methdb_line(char *line, meth_rec *tmp_meth)
{
	char *p, *q= 0;
    int i, pos;
	float ref_meth, alt_meth, context, type;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: pos     = atoi(p); break;
            case 1: ref_meth= atof(p); break;
            case 2: alt_meth= atof(p); break;
            case 3: context = atoi(p); break;
            case 4: type    = atoi(p); break;
            default: break;}
			++i, p = q + 1;
			if (i > 5 || c == '\0') break;
		}
	}

    if(i < 4){fprintf(stderr, "[%s] Skip invalid site: position %d...\n", __func__, pos);}

    tmp_meth->pos      = pos;
    tmp_meth->ref_meth = ref_meth;
    tmp_meth->alt_meth = alt_meth;
    tmp_meth->context  = context;
    tmp_meth->type     = type;
}

inline void load_methdb(uint32_t *pos_idx_arr, std::vector<meth_rec>& meth_vec, char *fname){
    htsFile *fp = hts_open(fname,"rb");
    if(fp == 0 ){ fprintf(stderr,"[%s] ERROR: open methdb file: %s failed: %s. Exit...\n", __func__, fname, strerror (errno)); exit (EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    meth_rec tmp_meth;

    int ret;
    int num_404_site= 0; // to record the position conflicting sites
    int num_tot_site= 0;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        parse_methdb_line(line.s, &tmp_meth); //might need to test
        if (tmp_meth.pos == -1) {continue;}
        ++num_tot_site;

        int tmp_idx_rec = pos_idx_arr[tmp_meth.pos];
        int tmp_idx = tmp_idx_rec >> 2;
        int tmp_flag= tmp_idx_rec & 0x1;

        if(meth_vec[tmp_idx].pos!=tmp_meth.pos && meth_vec[tmp_idx].context!=tmp_meth.context){ ++num_404_site; continue;}
        meth_vec[tmp_idx].ref_meth  = tmp_meth.ref_meth;
        meth_vec[tmp_idx].alt_meth  = tmp_meth.alt_meth;
        meth_vec[tmp_idx].type      = tmp_meth.type;
        tmp_meth        = {};
    }
    free(line.s);
    if (ret=hts_close(fp)){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    float ratio_404 = (float)num_404_site/num_tot_site;
    if(ratio_404 > 0.5){fprintf(stderr,"[%s] WARNING: over 50%% sites in CGmap file are not compatible with reference genome, please check!\n", __func__); exit(1);}
    fprintf(stderr,"[%s] Contig %s: %d sites found in cgmap, among them %.2f%% sites not compatiable\n", __func__, num_tot_site, ratio_404);
}


// for CGmap
inline int parse_cgmap_line(char *line, char *chr_id, meth_rec *tmp_meth)
{
	char *p, *q= 0;
    char *contig, *base, *context, *diN = 0;
    int i, pos;
	float meth;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: contig  = p; break;
            case 1: base    = strdup(p); break;
            case 2: pos     = atoi(p); break;
            case 3: context = strdup(p); break;
            case 4: diN     = 0; break; // skip: =strdup(p)
            case 5: meth    = atof(p); break;
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return 1;} // termenate early if not equal
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

    if(i < 5){fprintf(stderr, "[%s] Skip invalid site: chr %s, position %d...\n", __func__, chr_id, pos); return 0;}

    tmp_meth->ref_meth = meth;
    tmp_meth->pos      = pos;
    tmp_meth->context  = base_map[base] << 3 | context_map[std::string(context)];
    return 0;
}

inline void pool_cgmap(std::vector<meth_rec>& meth_vec, int seed)
{
    std::vector<float> cg_vec;
    std::vector<float> chg_vec;
    std::vector<float> chh_vec;
    
    for (size_t i=0; i < meth_vec.size(); ++i){
        if (meth_vec[i].type){
            int tmp_context = meth_vec[i].context &0x7;
            if(tmp_context == 1){cg_vec.push_back(meth_vec[i].ref_meth);}
            else if (tmp_context == 3){chg_vec.push_back(meth_vec[i].ref_meth);}
            else{chh_vec.push_back(meth_vec[i].ref_meth);}
        }
    }

    if (seed <= 0) seed = time(0)&0x7fffffff;
    std::mt19937 eng(seed);
    std::uniform_int_distribution<> dist_cg(0, cg_vec.size()-1);
    std::uniform_int_distribution<> dist_chg(0, chg_vec.size()-1);
    std::uniform_int_distribution<> dist_chh(0, chh_vec.size()-1); // check if the vector can be length of 0

    for (size_t i=0; i < meth_vec.size(); ++i){
        int tmp_context = meth_vec[i].context &0x7;
        if(tmp_context == 1){meth_vec[i].ref_meth = cg_vec[dist_cg(eng)];}
        else if (tmp_context == 3){meth_vec[i].ref_meth = chg_vec[dist_chg(eng)];}
        else{meth_vec[i].ref_meth = chh_vec[dist_chh(eng)];}
    }
}

inline void fill_cgmap_chr(char *fname, char *chr_id, uint32_t *pos_idx_arr, std::vector<meth_rec>& meth_vec, bool cgmap_pool, int seed){
    htsFile *fp = hts_open(fname,"rb");
    if(fp == 0 ){ fprintf(stderr,"[%s] ERROR: open CGmap file: %s failed: %s. Exit...\n", __func__, fname, strerror (errno)); exit (EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    meth_rec tmp_meth;

    int ret;
    int num_404_site= 0; // to record the position conflicting sites
    int num_tot_site= 0;
    bool collect_present = false;
    bool collect_previous= false;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (collect_present == false && collect_previous == true) { break; } // finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        int skip_status = parse_cgmap_line(line.s, chr_id, &tmp_meth); //might need to test
        if(skip_status){collect_present = false; continue;}
        if (tmp_meth.pos == -1) {continue;}
        ++num_tot_site;
        int tmp_idx_rec = pos_idx_arr[tmp_meth.pos];
        int tmp_idx = tmp_idx_rec >> 2;
        int tmp_flag= tmp_idx_rec & 0x1;
        if(meth_vec[tmp_idx].pos!=tmp_meth.pos && meth_vec[tmp_idx].context!=tmp_meth.context){ ++num_404_site; continue;}
        meth_vec[tmp_idx].ref_meth  = tmp_meth.ref_meth;
        meth_vec[tmp_idx].context   = tmp_meth.context;
        meth_vec[tmp_idx].type      = 2;
        tmp_meth        = {};
    }
    free(line.s);
    if (ret=hts_close(fp)){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    float ratio_404 = (float)num_404_site/num_tot_site;
    if(ratio_404 > 0.5){fprintf(stderr,"[%s] WARNING: over 50%% sites in CGmap file are not compatible with reference genome, please check!\n", __func__); exit(1);}
    fprintf(stderr,"[%s] Contig %s: %d sites found in cgmap, among them %.2f%% sites not compatiable\n", __func__, chr_id, num_tot_site, ratio_404);

    if(cgmap_pool){pool_cgmap(meth_vec, seed);}
}


// for ASM
inline int parse_asm_line(char *line, char *chr_id, meth_rec *tmp_meth)
{
	char *p, *q= 0;
    char *contig, *base = 0;
    int i, pos;
	float meth, ref_meth, alt_meth;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: contig  = p; break;
            case 1: base    = strdup(p); break;
            case 2: pos     = atoi(p); break;
            case 3: meth    = atof(p); break; // can skip
            case 4: ref_meth= atof(p); break;
            case 5: alt_meth= atof(p); break;
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return 1;}
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

    if(i < 5){fprintf(stderr, "[%s] Skip invalid site: chr %s, position %d...\n", __func__, chr_id, pos); return 0;}

    tmp_meth->pos      = pos;
    tmp_meth->ref_meth = meth;
    tmp_meth->alt_meth = alt_meth;
    return 0;
}

inline void fill_asm_chr(char *fname, char *chr_id, uint32_t *pos_idx_arr, std::vector<meth_rec>& meth_vec){
    htsFile *fp = hts_open(fname,"rb");
    if(fp == 0 ){ fprintf(stderr,"[%s] ERROR: open ASM file: %s failed: %s. Exit...\n", __func__, fname, strerror (errno)); exit (EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    meth_rec tmp_meth;

    int ret;
    int num_404_site= 0; // to record the position conflicting sites
    int num_tot_site= 0;
    bool collect_present = false;
    bool collect_previous= false;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (collect_present == false && collect_previous == true) { break; } // finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        int skip_status = parse_cgmap_line(line.s, chr_id, &tmp_meth); //might need to test
        if(skip_status){collect_present = false; continue;}
        if (tmp_meth.pos == -1) {continue;}
        ++num_tot_site;
        int tmp_idx_rec = pos_idx_arr[tmp_meth.pos];
        int tmp_idx = tmp_idx_rec >> 2;
        int tmp_flag= tmp_idx_rec & 0x1;
        if(meth_vec[tmp_idx].pos   != tmp_meth.pos){ ++num_404_site; continue;}
        meth_vec[tmp_idx].ref_meth  = tmp_meth.ref_meth;
        meth_vec[tmp_idx].alt_meth  = tmp_meth.alt_meth;
        meth_vec[tmp_idx].type      = 4;
        tmp_meth        = {};
    }
    free(line.s);
    if (ret=hts_close(fp)){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    float ratio_404 = (float)num_404_site/num_tot_site;
    if(ratio_404 > 0.5){fprintf(stderr,"[%s] WARNING: over 50%% sites in ASM file are not compatiable with reference genome, please check!\n", __func__); exit(1);}
    fprintf(stderr,"[%s] Contig %s: %d sites found in cgmap, among them %.2f%% sites not compatiable\n", __func__, chr_id, num_tot_site, ratio_404);    
}


// fill with distribution
inline void parse_param(char *param_str, std::vector<param_rec>& param_vec)
{
    param_vec.clear();
    param_rec tmp_param;
    std::string tmp_str;
    
    for(int i =0; param_str[i] !='\0'; ++i){
        if (param_str[i] == '|'){
            tmp_param.alpha = stof(tmp_str);    // can be illegal negative
            tmp_str.clear();
        } else if (param_str[i] == ',') {
            tmp_param.beta  = stof(tmp_str);    // can be illegal negative
            if (tmp_param.alpha>=0 && tmp_param.beta>=0){param_vec.push_back(tmp_param);}else{
                fprintf (stderr, "[%s] ERROR: Invalid parameter format (%s). Exit... \n", __func__, param_str); exit (EXIT_FAILURE);
            }
            tmp_str.clear();
            tmp_param = {};
            continue;
        } else {
            tmp_str += param_str[i];
        }
    }
    // when finished, check tmp_str (deals with forgetting the last ',')
    if (!tmp_str.empty()){
        tmp_param.beta  = stof(tmp_str);
        if (tmp_param.alpha>=0 && tmp_param.beta>=0){param_vec.push_back(tmp_param);}else{
                fprintf (stderr, "[%s] ERROR: Invalid parameter format (%s). Exit... \n", __func__, param_str); exit (EXIT_FAILURE);
        }
    }
}

inline float gen_beta(gsl_rng *rng, uint8_t context, std::vector<param_rec>& param_vec)
{
    return (float)gsl_ran_beta(rng, param_vec[params_map[context]].alpha, param_vec[params_map[context]].beta);
}

inline void fill_beta(std::vector<meth_rec>& meth_vec, std::vector<param_rec>& param_vec, int seed)
{
    const gsl_rng_type * T;
    gsl_rng *rng;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    rng = gsl_rng_alloc(T);
    if (seed <= 0) seed = time(0)&0x7fffffff;
    gsl_rng_set(rng, seed);

    for (size_t i=0; i < meth_vec.size(); ++i){
        //detect if it's not filled
        if(meth_vec[i].type){continue;}
        meth_vec[i].ref_meth = gen_beta(rng, meth_vec[i].context, param_vec);
        meth_vec[i].alt_meth = meth_vec[i].ref_meth;
        meth_vec[i].type     = 8;
    }
}


// for GC-bias
inline void parse_bias_file(char *fname, std::vector<float>& eff_vec)
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


// for length calculation
inline void cal_length_chr(kseq_t *ks, chr_rec *tmp_len, int tech_mode, char *bed_file)
{
    uint64_t eff_len;

    if(tech_mode==2){                       // targeted sequencing
        parse_bed_chr(bed_file, ks->name.s, probe_vec);
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
        gen_cut_pos(ks, cutpos_vec);
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
    tmp_len->chr_len = ks->seq.l;
    tmp_len->eff_len = eff_len;
}


#endif
