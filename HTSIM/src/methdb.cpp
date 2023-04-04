#include <vector>
#include <map>
#include <zlib.h>
#include <random>
#include <gsl/gsl_randist.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"

KSEQ_INIT(gzFile, gzread)
#include "methdb.h" // must be here otherwise .h will not be recognized


std::map<std::string, int> base_map = {{"C", 0}, {"G",1}};
std::map<std::string, int> context_map = {{"CG",1}, {"CHG",3}, {"CHH",7}};
std::map<int, int> params_map = {{1,0}, {3,1}, {7,2}, {9,0}, {11,1}, {15,2}}; //idx in param_vec


// create/stream MethDB
void create_methdb(const kseq_t *ks, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec)
{
    int k = ks->seq.l;
    meth_vec.clear();               // clean the container and reserve
    meth_vec.reserve(k/2);
    
    meth_rec tmp_meth;
    meth_vec.push_back(tmp_meth);   // put the unitialized meth_rec as the first element

    int c, c_d1, c_d2;
    uint32_t ix = 1;
    
    for (int i = 0; i < k; ++i) {
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
            posidx_arr[i] |= ix << 2;
            ix++;
            tmp_meth = {};
        }
    }
}

void update_methdb(uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec, mutseq_t *hap1, mutseq_t *hap2, bool is_asm_set, bool is_meth_update)
{
    // mutseq_t *ret[2];
    // ret[0] = hap1; ret[1] = hap2;

    int tmp_pos, tmp_cg;
    int k = hap1->l;
    for(int i=0; i < k; ++i){
        if((posidx_arr[i] & 0x1) != 0){                     // if there is a variant
            tmp_pos = posidx_arr[i] >> 2;
            tmp_cg  = (posidx_arr[i]&0x2)>>1;
            if(tmp_cg){                                     // if the record is in meth_vec
                meth_vec[tmp_pos].meth[0] = -1;
            }
            if(is_meth_update){}                          // TODO: update boundary sites
        }
    }

    if(is_asm_set){                                       // use the last bit to store asm signal
        for(int i=0; i < k; ++i){
            posidx_arr[i] &= 0xfffffffe;
            tmp_pos = posidx_arr[i] >> 2;
            if(tmp_pos){posidx_arr[i] |= (uint32_t)meth_vec[tmp_pos].meth[0] != meth_vec[tmp_pos].meth[1];}
        }
    }
}

void save_methdb(std::vector<meth_rec>& meth_vec, const char *fname)
{
    FILE* fp = fopen(fname, "w");
    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open methdb file: %s failed. Exit...\n", __func__, fname); exit (EXIT_FAILURE);}

    meth_rec tmp_meth;
    for (size_t i=1; i < meth_vec.size(); ++i){
        tmp_meth = meth_vec[i];
        fprintf(fp, "%d\t%f\t%f\t%d\t%d\n", tmp_meth.pos, tmp_meth.meth[0], tmp_meth.meth[1], tmp_meth.context, tmp_meth.type);
        if (ferror(fp)) {fprintf(stderr, "[%s] ERROR: failed to write to file %s. Exit...\n", __func__, fname);exit(EXIT_FAILURE);}
    }
    fclose(fp);
}

void parse_methdb_line(char *line, meth_rec *tmp_meth)
{
    // if (line == NULL) {fprintf(stderr, "[%s] ERROR: input line is null. Exit...\n", __func__); exit(EXIT_FAILURE);}
	char *p, *q= 0;
    int i, pos=-1;
	float ref_meth=-1, alt_meth=-1;
    int context=0, type=0;

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

    tmp_meth->pos     = pos;
    tmp_meth->meth[0] = ref_meth;
    tmp_meth->meth[1] = alt_meth;
    tmp_meth->context = context;
    tmp_meth->type    = type;
}

void load_methdb(uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec, char *fname)
{
    htsFile *fp = hts_open(fname,"rb");
    if(fp == 0 ){ fprintf(stderr,"[%s] ERROR: open methdb file: %s failed: %s. Exit...\n", __func__, fname, strerror (errno)); exit (EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    meth_rec tmp_meth;

    int ret;
    int num_404_site= 0; // to record the position conflicting sites
    int num_tot_site= 0;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        parse_methdb_line(line.s, &tmp_meth);
        if (tmp_meth.pos == -1) {continue;}
        ++num_tot_site;

        int tmp_idx_rec = posidx_arr[tmp_meth.pos];
        int tmp_idx = tmp_idx_rec >> 2;

        if(meth_vec[tmp_idx].pos!=tmp_meth.pos && meth_vec[tmp_idx].context!=tmp_meth.context){ ++num_404_site; continue;}
        meth_vec[tmp_idx].meth[0]= tmp_meth.meth[0];
        meth_vec[tmp_idx].meth[1]= tmp_meth.meth[1];
        meth_vec[tmp_idx].type   = tmp_meth.type;
        tmp_meth = {};
    }
    free(line.s);
    if ((ret=hts_close(fp))){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    if(num_tot_site==0){fprintf(stderr,"[%s] no valid site found in methdb file: %s. Skip...\n", __func__, fname); /*exit (EXIT_FAILURE);*/}else{
        float ratio_404 = (float)num_404_site/num_tot_site;
        if(ratio_404 > 0.5){fprintf(stderr,"[%s] WARNING: over 50%% sites in CGmap file are not compatible with reference genome, please check!\n", __func__); exit(1);}
        fprintf(stderr,"[%s] %d sites found in cgmap, among them %.2f%% sites not compatiable\n", __func__, num_tot_site, ratio_404);
    }
}


// for CGmap
int parse_cgmap_line(char *line, char *chr_id, meth_rec *tmp_meth)
{
	char *p, *q= 0;
    char *contig, *base, *context, *diN = 0;
    int i, pos;
	float meth = 0;

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

    if(i < 5){fprintf(stderr, "[%s] Skip invalid site: chr %s, position %d...\n", __func__, chr_id, pos); free(base); free(context); return 0;}

    tmp_meth->meth[0]  = meth;
    tmp_meth->pos      = pos - 1; // convert to 0-based
    tmp_meth->context  = base_map[base] << 3 | context_map[std::string(context)];
    free(base); free(context);
    return 0;
}

void pool_cgmap(std::vector<meth_rec>& meth_vec, int seed)
{
    std::vector<float> cg_vec;
    std::vector<float> chg_vec;
    std::vector<float> chh_vec;
    
    for (size_t i=0; i < meth_vec.size(); ++i){
        if (meth_vec[i].type){
            int tmp_context = meth_vec[i].context &0x7;
            if(tmp_context == 1){cg_vec.push_back(meth_vec[i].meth[0]);}
            else if (tmp_context == 3){chg_vec.push_back(meth_vec[i].meth[0]);}
            else{chh_vec.push_back(meth_vec[i].meth[0]);}
        }
    }

    if (seed <= 0) seed = time(0)&0x7fffffff;
    std::mt19937 eng(seed);
    std::uniform_int_distribution<> dist_cg;
    std::uniform_int_distribution<> dist_chg;
    std::uniform_int_distribution<> dist_chh;
    bool cg_filled  = !cg_vec.empty();
    bool chg_filled = !chg_vec.empty();
    bool chh_filled = !chh_vec.empty();
    if (cg_filled) {dist_cg  = std::uniform_int_distribution<>(0, cg_vec.size()  - 1);}
    if (chg_filled){dist_chg = std::uniform_int_distribution<>(0, chg_vec.size() - 1);}
    if (chh_filled){dist_chh = std::uniform_int_distribution<>(0, chh_vec.size() - 1);}

    for (size_t i = 1; i < meth_vec.size(); ++i) {
        int tmp_context = meth_vec[i].context & 0x7;
        if (tmp_context == 1 && cg_filled) {
            meth_vec[i].meth[0] = cg_vec[dist_cg(eng)];
        } else if (tmp_context == 3 && chg_filled) {
            meth_vec[i].meth[0] = chg_vec[dist_chg(eng)];
        } else if (tmp_context == 7 && chh_filled) {
            meth_vec[i].meth[0] = chh_vec[dist_chh(eng)];
        }
        meth_vec[i].meth[1] = meth_vec[i].meth[0];
        meth_vec[i].type = 10;
    }
}

void fill_cgmap_chr(char *fname, char *chr_id, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec, meth_param *meth_set)
{
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
        int skip_status = parse_cgmap_line(line.s, chr_id, &tmp_meth);
        if(skip_status){collect_present = false; continue;}
        if (tmp_meth.pos == -1) {continue;}
        ++num_tot_site;
        int tmp_idx_rec = posidx_arr[tmp_meth.pos];
        int tmp_idx = tmp_idx_rec >> 2;
        if(meth_vec[tmp_idx].pos!=tmp_meth.pos && meth_vec[tmp_idx].context!=tmp_meth.context){ ++num_404_site; continue;}
        meth_vec[tmp_idx].meth[0]= tmp_meth.meth[0];
        meth_vec[tmp_idx].context= tmp_meth.context;
        meth_vec[tmp_idx].type   = 2;
        tmp_meth    = {};
    }
    free(line.s);
    if ((ret=hts_close(fp))){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    if(num_tot_site == 0){fprintf(stderr,"[%s] no valid sites found in CGmap file: %s. Skip...\n", __func__, fname); /*exit (EXIT_FAILURE);*/}else{
        float ratio_404 = (float)num_404_site/num_tot_site;
        if(ratio_404 > 0.5){fprintf(stderr,"[%s] WARNING: over 50%% sites in CGmap file are not compatible with reference genome, please check!\n", __func__); exit(1);}
        fprintf(stderr,"[%s] Contig %s: %d sites found in cgmap, among them %.2f%% sites not compatiable\n", __func__, chr_id, num_tot_site, ratio_404);
    }
    
    if(meth_set->cgmap_pool){pool_cgmap(meth_vec, meth_set->seed_meth);}
}


// for ASM
int parse_asm_line(char *line, char *chr_id, meth_rec *tmp_meth)
{
	char *p, *q= 0;
    char *contig, *base = 0;
    int i, pos = -1;
	float meth, ref_meth=-1, alt_meth=-1;

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

    if(i < 5){fprintf(stderr, "[%s] Skip invalid site: chr %s, position %d...\n", __func__, chr_id, pos); free(base); return 0;}

    tmp_meth->pos     = pos - 1;     // 0-based
    tmp_meth->meth[0] = ref_meth;
    tmp_meth->meth[1] = alt_meth;
    free(base);
    return 0;
}

void fill_asm_chr(char *fname, char *chr_id, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec)
{
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
        int skip_status = parse_asm_line(line.s, chr_id, &tmp_meth); //might need to test
        if(skip_status){collect_present = false; continue;}
        if (tmp_meth.pos == -1) {continue;}
        ++num_tot_site;
        int tmp_idx_rec = posidx_arr[tmp_meth.pos];
        int tmp_idx = tmp_idx_rec >> 2;
        if(meth_vec[tmp_idx].pos != tmp_meth.pos){ ++num_404_site; continue;}
        meth_vec[tmp_idx].meth[0] = tmp_meth.meth[0];
        meth_vec[tmp_idx].meth[1] = tmp_meth.meth[1];
        meth_vec[tmp_idx].type    = 4;
        tmp_meth    = {};
    }
    free(line.s);
    if ((ret=hts_close(fp))){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    if(num_tot_site == 0){fprintf(stderr,"[%s] no valid sites found in ASM file: %s. Skip...\n", __func__, fname); /*exit (EXIT_FAILURE);*/}else{
        float ratio_404 = (float)num_404_site/num_tot_site;
        if(ratio_404 > 0.5){fprintf(stderr,"[%s] WARNING: over 50%% sites in ASM file are not compatiable with reference genome, please check!\n", __func__); exit(1);}
        fprintf(stderr,"[%s] Contig %s: %d sites found in cgmap, among them %.2f%% sites not compatiable\n", __func__, chr_id, num_tot_site, ratio_404);   
    }
}


// fill with distribution
void parse_param(char *param_str, std::vector<param_rec>& param_vec)
{
    // should detect illegal input
    param_vec.clear();
    param_rec tmp_param;
    std::string tmp_str;
    
    for(int i =0; param_str[i] !='\0'; ++i){
        if (param_str[i] == '_'){               // alpha found
            tmp_param.alpha = stof(tmp_str);
            tmp_str.clear();
        } else if (param_str[i] == ',') {       // beta found
            tmp_param.beta  = stof(tmp_str);
            if (tmp_param.alpha < 0 || tmp_param.beta < 0){
                fprintf (stderr, "[%s] ERROR: parameter cannot be negative (%s). Exit...\n", __func__, param_str); exit (EXIT_FAILURE);
            }
            param_vec.push_back(tmp_param);
            tmp_str.clear();
            tmp_param = {};
            continue;
        } else {
            tmp_str += param_str[i];
        }
    }

    // Check for an incomplete parameter at the end (missing comma)
    if (!tmp_str.empty()){
        tmp_param.beta  = stof(tmp_str);
        if (tmp_param.alpha < 0 || tmp_param.beta < 0){
            fprintf (stderr, "[%s] ERROR: parameter cannot be negative (%s). Exit...\n", __func__, param_str); exit (EXIT_FAILURE);
        }
        param_vec.push_back(tmp_param);
    }
}

float gen_beta(gsl_rng *rng, uint8_t context, std::vector<param_rec>& param_vec)
{
    return (float) gsl_ran_beta(rng, param_vec[params_map[context]].alpha, param_vec[params_map[context]].beta);
}

void fill_beta(std::vector<meth_rec>& meth_vec, std::vector<param_rec>& param_vec, int seed)
{
    if(param_vec.size()==0){fprintf(stderr,"[%s] ERROR: No parameter provided. Exit...\n", __func__); exit(1);}
    
    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    if (seed <= 0) seed = time(0)&0x7fffffff;
    gsl_rng_set(rng, seed);

    for (size_t i=1; i < meth_vec.size(); ++i){
        if(meth_vec[i].type){
            if(meth_vec[i].meth[1] < 0){meth_vec[i].meth[1]  = meth_vec[i].meth[0];}
        }else{
            //if it's not filled
            meth_vec[i].meth[0] = gen_beta(rng, meth_vec[i].context, param_vec);
            meth_vec[i].meth[1] = meth_vec[i].meth[0];
            meth_vec[i].type    = 8;
        }
    }
    gsl_rng_free(rng);
}

