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


std::map<std::string, int> context_map = {{"CG",1}, {"CHG",3}, {"CHH",7}};


// create/stream MethDB
void create_methdb(const kseq_t *ks, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec)
{
    int k = ks->seq.l;
    meth_vec.clear();               // clean the container and reserve
    meth_vec.reserve((int) k*0.48); // chr19 has the highest GC raito (0.4794)
    
    meth_rec tmp_meth;
    meth_vec.push_back(tmp_meth);   // put the unitialized meth_rec as the first element

    int c, c_d1, c_d2, c_kmeridx;
    int ix_base;
    
    for (int i = 0; i < k; ++i) {
        c = nst_nt4_table[(int)ks->seq.s[i]];
        if (cg_table[(uint8_t) c]){
            tmp_meth.pos  = i;
            if(c==1){
                if(i > k-3){c_d1 =0; c_d2 =0;}else{     /*handle the last 2 base*/
                    c_d1 = nst_nt4_table[(int)ks->seq.s[i + 1]];
                    c_d2 = nst_nt4_table[(int)ks->seq.s[i + 2]];
                }
            }else{
                if(i < 2){c_d1 =0; c_d2 =0;}else{       /*handle the first 2 base*/
                    c_d1 = nst_nt4_table[(int)ks->seq.s[i - 1]];
                    c_d2 = nst_nt4_table[(int)ks->seq.s[i - 2]];
                }
            }
            uint8_t context_idx = c << 4 | c_d1 <<2 | c_d2;
            tmp_meth.context[0] = cg_context_table[context_idx];
            tmp_meth.context[1] = cg_context_table[context_idx];

            if(i < 3 || i > k-4){c_kmeridx = 0;}else{   /*handle the first 3 base and last 3 base*/
                c_kmeridx = 0;
                for (int j = -3; j < 4; ++j){
                    ix_base = nst_nt4_table[(int)ks->seq.s[i+j]];
                    ix_base = ix_base < 4 ? ix_base : (rand()&0x3);          // will not interfere with the drand48()
                    c_kmeridx = (c_kmeridx << 2) | ix_base;
                }
            }
            tmp_meth.kmeridx[0] = c_kmeridx;
            tmp_meth.kmeridx[1] = c_kmeridx;
            meth_vec.push_back(tmp_meth);
            posidx_arr[i] |= 2;
            tmp_meth = {};
        }
    }
}

void update_variant(const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec, 
                    meth_param *meth_set, std::map<int, param_rec>& params_map, std::map<int, snpmeth_rec>& snpmeth_map)
{
    const gsl_rng_type * T = gsl_rng_default;
    gsl_rng *rng = gsl_rng_alloc(T);
    gsl_rng_env_setup();
    gsl_rng_set(rng, meth_set->seed_meth);

    mutseq_t *rseq;
    snpmeth_rec tmp_snpmeth;
    std::vector<int> tmp_sites;
    std::vector<uint8_t> sites_u = {0,0,0};
    std::vector<uint8_t> sites_d = {0,0,0};

    int i, c[3];
    int max_ik = ks->seq.l-1;
    int pos_k, pos_base, pos_mask, pos_ins_len;
    int tmp_pos, tmp_base, tmp_mask, tmp_ins_len;
    int tmp_kmeridx, tmp_context, tmp_idx;
    int count, t;

    if(snpmeth_map.size()==0){collect_snpmeth(ks, hap1, hap2, posidx_arr, snpmeth_map);}

    // iterate through the snpmeth_map
    for (std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it) {
        if(it->second.offset < 0){continue;} // deletion is not on the read
        i = it->first;
        rseq = it->second.hap1==1? hap1 : hap2;
        
        // handle the non-mutational neighborhood sites: calculate sites within 3 bases from i
        // for each of them calulate the context and kmeridx(if needed)
        for(int k =-3; k<4; ++k){
            pos_k = std::min(std::max(0, i+k), max_ik); // ensure not passing the boundary
            pos_base = rseq->s[pos_k];
            if(pos_base&mutmsk) continue;               // skip mutation, will skip k=0
            if(!cg_table[pos_base]) continue;           // skip nonCG
            
            // fill in vector
            #define __fill_vec(sites_ptr, updown)                       \
                count=0, t=0;                                           \
                while(count < 3){                                       \
                    ++t;                                                \
                    tmp_pos  = pos_k+t*updown;                          \
                    if(tmp_pos < 0 || tmp_pos >= ks->seq.l) break;      \
                    tmp_base = rseq->s[tmp_pos];                        \
                    tmp_mask = (tmp_base&mutmsk);                       \
                    if(tmp_mask == DELETE) continue;                    \
                    (*sites_ptr)[count]= tmp_base & 0xf;                \
                    ++count;                                            \
                    if(tmp_mask!=SUBSTITUTE && tmp_mask!=NOCHANGE){     \
                        tmp_base >>= 4;                                 \
                        tmp_ins_len  = tmp_mask >> 12;                  \
                        while (count < 3 && tmp_ins_len > 0){           \
                            (*sites_ptr)[count]= tmp_base&0x3;          \
                            ++count;                                    \
                            tmp_base >>= 2;                             \
                            --tmp_ins_len;                              \
                        }                                               \
                    }                                                   \
                }                                                       \

            std::fill(sites_u.begin(), sites_u.end(), 0);
            std::fill(sites_d.begin(), sites_d.end(), 0);
            __fill_vec(&sites_u,-1);
            __fill_vec(&sites_d, 1);

            // compute and update context
            if (pos_base == 1){         // C
                tmp_context = cg_context_table[((pos_base<<4) | (sites_d[0]<<2) | sites_d[1])]; 
            }else if (pos_base == 2){   // G
                tmp_context = cg_context_table[((pos_base<<4) | (sites_u[0]<<2) | sites_u[1])];
            }else{} // should never happen

            tmp_idx = posidx_arr[pos_k] >> 2;
            if(meth_vec[tmp_idx].context[0] != tmp_context){
                meth_vec[tmp_idx].context[1] = tmp_context;
                meth_vec[tmp_idx].meth[1]    = gen_beta(rng, tmp_context, params_map);
            }

            //compute and update kmeridx
            tmp_kmeridx = sites_u[2] << 6 | sites_u[1] << 4 | sites_u[0] << 2 | pos_base;
            tmp_kmeridx = (tmp_kmeridx<<6)| sites_d[0] << 4 | sites_d[1] << 2 | sites_d[2];
            meth_vec[tmp_idx].kmeridx[1] = tmp_kmeridx;
        }

        // handle the snp sites
        pos_k = i;
        pos_base = rseq->s[i];
        pos_mask = (pos_base&mutmsk);
        if(pos_mask==DELETE){continue;}                         // DELETE will not appear on the read
        
        std::fill(sites_u.begin(), sites_u.end(), 0);
        std::fill(sites_d.begin(), sites_d.end(), 0);
        __fill_vec(&sites_u,-1);
        __fill_vec(&sites_d, 1);

        tmp_sites.clear();
        tmp_snpmeth = it->second;
        for(int k = 2; k >= 0; --k){tmp_sites.push_back(sites_u[k]);}
        pos_ins_len = pos_mask==SUBSTITUTE ? 0: pos_mask >> 12; // 0 for SUBSTITUTE & NOCHANGE
        if(pos_ins_len == 0){tmp_sites.push_back(pos_base&0xf);}else{
            pos_base >>= 4;
            for(int k = 0; k < pos_ins_len; ++k){
                tmp_sites.push_back(pos_base&0x3); 
                pos_base>>=2;
            }
        }
        for(int k = 0; k <= 2; ++k){tmp_sites.push_back(sites_d[k]);}

        for(int k = 3; k<4+pos_ins_len; ++k){
            tmp_idx = k-3;
            tmp_kmeridx = (tmp_sites[k-3] << 6) | (tmp_sites[k-2] << 4) | (tmp_sites[k-1] << 2) | tmp_sites[k];
            tmp_kmeridx = (tmp_kmeridx << 6) | (tmp_sites[k+1] << 4) | (tmp_sites[k+2] << 2) | tmp_sites[k+3];

            if(tmp_sites[k] == 1){         // C
                tmp_context = cg_context_table[((tmp_sites[k]<<4) | (tmp_sites[k+1]<<2) | tmp_sites[k+2])];
            }else if (tmp_sites[k] == 2){   // G
                tmp_context = cg_context_table[((tmp_sites[k]<<4) | (tmp_sites[k-1]<<2) | tmp_sites[k-2])];
            }else{                      // else
                tmp_context = 0;
            }
            tmp_snpmeth.context[tmp_idx] = tmp_context;
            tmp_snpmeth.kmeridx[tmp_idx] = tmp_kmeridx;
            tmp_snpmeth.meth[tmp_idx]    = gen_beta(rng, tmp_context, params_map);
        }
        snpmeth_map[it->first] = tmp_snpmeth;
    }

    // if(meth_set->is_asm_set){                                       // use the last bit to store asm signal
    //     for(int i=0; i < ks->seq.l; ++i){
    //         posidx_arr[i] &= 0xfffffffe;
    //         int tmp_pos = posidx_arr[i] >> 2;
    //         if(tmp_pos){posidx_arr[i] |= (uint32_t)meth_vec[tmp_pos].meth[0] != meth_vec[tmp_pos].meth[1];}
    //     }
    // }
    gsl_rng_free(rng);
}

void save_methdb(char *fname, char *chr_id, std::vector<meth_rec>& meth_vec)
{
    FILE* fp = fopen(fname, "a");
    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open methdb file: %s failed. Exit...\n", __func__, fname); exit(EXIT_FAILURE);}

    meth_rec tmp_meth;
    for (size_t i=1; i < meth_vec.size(); ++i){
        tmp_meth = meth_vec[i];
        fprintf(fp, "%s\t%d\t%f\t%f\t%d\t%d\t%d\n", chr_id, tmp_meth.pos, tmp_meth.meth[0], tmp_meth.meth[1], tmp_meth.context[0], tmp_meth.context[1], tmp_meth.type);
        if(ferror(fp)){fprintf(stderr, "[%s] ERROR: failed to write to methdb. Exit...\n", __func__); exit(EXIT_FAILURE);}
    }
    fclose(fp);
}

int parse_methdb_line(char *line, char *chr_id, meth_rec *tmp_meth)
{
    // if(line==NULL){fprintf(stderr, "[%s] ERROR: input line is null. Exit...\n", __func__); exit(EXIT_FAILURE);}
	char *p, *q= 0;
    char *contig;
    int i, pos=-1;
	float ref_meth=-1, alt_meth=-1;
    int context1=0, context2=0, type=0;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: contig  = p; break;
            case 1: pos     = atoi(p); break;
            case 2: ref_meth= atof(p); break;
            case 3: alt_meth= atof(p); break;
            case 4: context1= atoi(p); break;
            case 5: context2= atoi(p); break;
            case 6: type    = atoi(p); break;
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return -1;} // termenate early if not equal
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

    if(i < 5){fprintf(stderr, "[%s] Skip invalid site: chr %s, position %d...\n", __func__, chr_id, pos); return 0;}

    tmp_meth->pos       = pos;
    tmp_meth->meth[0]   = ref_meth;
    tmp_meth->meth[1]   = alt_meth;
    tmp_meth->context[0]= context1;
    tmp_meth->context[1]= context2;
    tmp_meth->type      = type;
    return 0;
}

void load_methdb(char *fname, char *chr_id, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec)
{
    htsFile *fp = hts_open(fname,"rb");
    if(fp==0){fprintf(stderr, "[%s] ERROR: open methdb file: %s failed: %s. Exit...\n", __func__, fname, strerror(errno)); exit(EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    meth_rec tmp_meth;

    int ret;
    int num_404_site= 0; // to record the position conflicting sites
    int num_tot_site= 0;
    float ratio_404 = 0;
    bool collect_present = false;
    bool collect_previous= false;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (!collect_present && collect_previous) { break; } // collect_present false && collect_previous true, finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        int skip_status = parse_methdb_line(line.s, chr_id, &tmp_meth);
        if(skip_status < 0){collect_present = false; continue;}
        
        if (tmp_meth.pos == -1) {continue;}
        ++num_tot_site;

        int tmp_idx_rec = posidx_arr[tmp_meth.pos];
        int tmp_idx = tmp_idx_rec >> 2;

        if(meth_vec[tmp_idx].pos!=tmp_meth.pos && meth_vec[tmp_idx].context[0]!=tmp_meth.context[0]){ ++num_404_site; continue;}
        meth_vec[tmp_idx].meth[0]= tmp_meth.meth[0];
        meth_vec[tmp_idx].meth[1]= tmp_meth.meth[1];
        meth_vec[tmp_idx].type   = tmp_meth.type;
        tmp_meth = {};
    }
    free(line.s);
    if ((ret=hts_close(fp))){fprintf(stderr, "[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    ratio_404 = num_tot_site == 0? 0 : num_404_site/(float) num_tot_site;
    if(ratio_404 > 0.5){fprintf(stderr, "[%s] WARNING: over 50%% sites in CGmap file are not compatible with reference genome, please check!\n", __func__); exit(1);}
    fprintf(stderr, "[%s] %s: found %d sites in methdb file, among them %.2f%% sites not compatiable...\n", __func__, chr_id, num_tot_site, ratio_404);
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
            if(i==0 && strcmp(contig, chr_id)){return -1;} // termenate early if not equal
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

    if(i < 5){fprintf(stderr, "[%s] Skip invalid site: chr %s, position %d...\n", __func__, chr_id, pos); free(base); free(context); return 0;}

    tmp_meth->meth[0]  = meth;
    tmp_meth->pos      = pos - 1; // convert to 0-based
    tmp_meth->context[0] = (base == "G") << 3 | context_map[std::string(context)];
    tmp_meth->context[1] = tmp_meth->context[0];
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
            int tmp_context = meth_vec[i].context[0] &0x7;
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
        int tmp_context = meth_vec[i].context[0] & 0x7;
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
    if(fp == 0 ){ fprintf(stderr,"[%s] ERROR: open CGmap file: %s failed: %s. Exit...\n", __func__, fname, strerror(errno)); exit(EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    meth_rec tmp_meth;

    int ret;
    int num_404_site= 0; // to record the position conflicting sites
    int num_tot_site= 0;
    float ratio_404 = 0;
    bool collect_present = false;
    bool collect_previous= false;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (!collect_present && collect_previous) { break; } // collect_present false && collect_previous true, finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        int skip_status = parse_cgmap_line(line.s, chr_id, &tmp_meth);
        if(skip_status < 0){collect_present = false; continue;}
        if (tmp_meth.pos == -1) {continue;}
        ++num_tot_site;
        int tmp_idx_rec = posidx_arr[tmp_meth.pos];
        int tmp_idx = tmp_idx_rec >> 2;
        if(meth_vec[tmp_idx].pos!=tmp_meth.pos && meth_vec[tmp_idx].context[0]!=tmp_meth.context[0]){ ++num_404_site; continue;}
        meth_vec[tmp_idx].meth[0]= tmp_meth.meth[0];
        meth_vec[tmp_idx].context[0] = tmp_meth.context[0];
        meth_vec[tmp_idx].context[1] = tmp_meth.context[1];
        meth_vec[tmp_idx].type   = 2;
        tmp_meth    = {};
    }
    free(line.s);
    if ((ret=hts_close(fp))){fprintf(stderr, "[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    ratio_404 = num_tot_site == 0? 0 : num_404_site/(float) num_tot_site;
    if(ratio_404 > 0.5){fprintf(stderr, "[%s] WARNING: over 50%% sites in CGmap file are not compatible with reference genome, please check!\n", __func__); exit(1);}
    fprintf(stderr, "[%s] %s: found %d sites in cgmap file, among them %.2f%% sites not compatiable...\n", __func__, chr_id, num_tot_site, ratio_404);
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
            if(i==0 && strcmp(contig, chr_id)){return -1;}
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
    if(fp == 0 ){ fprintf(stderr, "[%s] ERROR: open ASM file: %s failed: %s. Exit...\n", __func__, fname, strerror(errno)); exit(EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    meth_rec tmp_meth;

    int ret;
    int num_404_site= 0; // to record the position conflicting sites
    int num_tot_site= 0;
    float ratio_404 = 0;
    bool collect_present = false;
    bool collect_previous= false;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (!collect_present && collect_previous) { break; } // collect_present false && collect_previous true, finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        int skip_status = parse_asm_line(line.s, chr_id, &tmp_meth); //might need to test
        if(skip_status < 0){collect_present = false; continue;}
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
    if ((ret=hts_close(fp))){fprintf(stderr, "[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    ratio_404 = num_tot_site == 0? 0 : num_404_site/(float) num_tot_site;
    if(ratio_404 > 0.5){fprintf(stderr, "[%s] WARNING: over 50%% sites in ASM file are not compatiable with reference genome, please check!\n", __func__); exit(1);}
    fprintf(stderr, "[%s] %s: found %d sites in asm file, among them %.2f%% sites not compatiable\n", __func__, chr_id, num_tot_site, ratio_404);
}


// fill with distribution
void parse_param(char *param_str, std::map<int, param_rec>& params_map)
{
    // should detect illegal input
    params_map.clear();
    param_rec tmp_param;
    std::string tmp_str;
    int context_idx[6] = {1,3,7,9,11,15};
    
    int idx = 0;
    for(int i =0; param_str[i] !='\0'; ++i){
        if (param_str[i] == '_'){               // alpha found
            tmp_param.alpha = stof(tmp_str);
            tmp_str.clear();
        } else if (param_str[i] == ',') {       // beta found
            tmp_param.beta  = stof(tmp_str);
            if (tmp_param.alpha < 0 || tmp_param.beta < 0){
                fprintf (stderr, "[%s] ERROR: parameter cannot be negative (%s). Exit...\n", __func__, param_str); exit(EXIT_FAILURE);
            }
            params_map[context_idx[idx]] = tmp_param;
            tmp_str.clear();
            tmp_param = {};
            ++idx;
            continue;
        } else {
            tmp_str += param_str[i];
        }
    }

    // Check for an incomplete parameter at the end (missing comma)
    if (!tmp_str.empty()){
        tmp_param.beta  = stof(tmp_str);
        if (tmp_param.alpha < 0 || tmp_param.beta < 0){
            fprintf (stderr, "[%s] ERROR: parameter cannot be negative (%s). Exit...\n", __func__, param_str); exit(EXIT_FAILURE);
        }
        params_map[context_idx[idx]] = tmp_param;
        ++idx;
    }
    if(idx%3!=2){fprintf(stderr, "[%s] ERROR: parameter number should be multiple of 3 (%s). Exit...\n", __func__, param_str); exit(EXIT_FAILURE);}
}

float gen_beta(gsl_rng *rng, uint8_t context, std::map<int, param_rec>& params_map)
{
    return context == 0? 1 : (float) gsl_ran_beta(rng, params_map[context].alpha, params_map[context].beta);
}

void fill_beta(std::vector<meth_rec>& meth_vec, std::map<int, param_rec>& params_map, int seed)
{
    if(params_map.size()==0){fprintf(stderr, "[%s] ERROR: No parameter provided. Exit...\n", __func__); exit(1);}
    
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
            meth_vec[i].meth[0] = gen_beta(rng, meth_vec[i].context[0], params_map);
            meth_vec[i].meth[1] = meth_vec[i].meth[0];
            meth_vec[i].type    = 8;
        }
    }
    gsl_rng_free(rng);
}


// for SNP meth
void collect_snpmeth(const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, std::map<int, snpmeth_rec>& snpmeth_map)
{
    int c[3];
    int i, j, ix = 0; // j keeps the end of the last deletion
    int tmp_c, hap1_mut, hap2_mut, hap1_ins_len, hap2_ins_len;
    snpmeth_rec tmp_snpmeth;
    mutseq_t* hap0 = nullptr;
    for (i = 0; i != (int)ks->seq.l; ++i) {
        if((posidx_arr[i] & 0x1) != 0){                     // if there is a variant
            c[0] = nst_nt4_table[(int)ks->seq.s[i]];
            if (c[0] >= 4) continue;                        // skip N
            c[1] = hap1->s[i]; hap1_mut = (c[1] & mutmsk);
            c[2] = hap2->s[i]; hap2_mut = (c[2] & mutmsk);
            if (hap1_mut != NOCHANGE || hap2_mut != NOCHANGE) {
                if (c[1] == c[2]) { // hom
                    tmp_snpmeth.hap1 = 1; tmp_snpmeth.hap2 = 1;
                    if (hap1_mut == SUBSTITUTE) {   // substitution
                        tmp_snpmeth.ref = c[0];
                        tmp_snpmeth.alt = c[1];
                        tmp_snpmeth.offset = 0;
                    } else if (hap1_mut == DELETE) {// del
                        tmp_c  = 0;
                        for (j = i; j < ks->seq.l && hap1->s[j] == hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j){} //find j
                        for (ix= j-1; ix>= i; --ix){tmp_c = tmp_c << 2 | hap1->s[ix];}
                        tmp_snpmeth.ref = tmp_c;
                        tmp_snpmeth.alt = 0;
                        tmp_snpmeth.offset = i-j;
                    } else {                        // ins
                        tmp_snpmeth.ref = 0;
                        tmp_snpmeth.alt = c[1];
                        tmp_snpmeth.offset = hap1_mut >> 12;// can check if ins_len > 4 if needed
                    }
                } else { // het
                    if (hap1_mut == SUBSTITUTE || hap2_mut == SUBSTITUTE) {  // substitution
                        if((hap1_mut == SUBSTITUTE) && (hap2_mut == NOCHANGE)){
                            tmp_snpmeth.hap1 = 1;   tmp_snpmeth.hap2 = 0;
                            tmp_snpmeth.ref  = c[2];tmp_snpmeth.alt  = c[1];
                            tmp_snpmeth.offset = 0;
                        }else if((hap1_mut == NOCHANGE) && (hap2_mut == SUBSTITUTE)){
                            tmp_snpmeth.hap1 = 0;   tmp_snpmeth.hap2 = 1;
                            tmp_snpmeth.ref  = c[1];tmp_snpmeth.alt  = c[2];
                            tmp_snpmeth.offset = 0;
                        }else{
                            fprintf(stderr, "[%s] Error: hap1_mut = %d, hap2_mut = %d\n", __func__, hap1_mut, hap2_mut); exit(EXIT_FAILURE);
                        }
                    } else if (hap1_mut == DELETE || hap2_mut == DELETE) {  // del
                        tmp_c= 0;
                        if(hap1_mut == DELETE && hap2_mut == NOCHANGE){
                            hap0 = hap1; tmp_snpmeth.hap1 = 1; tmp_snpmeth.hap2 = 0;
                        }else if (hap1_mut == NOCHANGE && hap2_mut == DELETE){
                            hap0 = hap2; tmp_snpmeth.hap1 = 0; tmp_snpmeth.hap2 = 1;
                        }else{
                            fprintf(stderr, "[%s] Error: hap1_mut = %d, hap2_mut = %d\n", __func__, hap1_mut, hap2_mut); exit(EXIT_FAILURE);
                        }
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap0->s[j]&mutmsk) == DELETE; ++j){}
                        for (ix= j-1; ix>= i; --ix){tmp_c = tmp_c << 2 | hap0->s[ix];}
                        tmp_snpmeth.ref = tmp_c;
                        tmp_snpmeth.alt = 0;
                        tmp_snpmeth.offset = i-j;
                    } else {
                        hap1_ins_len = hap1_mut >> 12;
                        hap2_ins_len = hap2_mut >> 12;
                        if(hap2_ins_len == 0 && hap1_ins_len > 0 && hap1_ins_len <=4){
                            tmp_snpmeth.hap1 = 1; tmp_snpmeth.hap2 = 0;
                            tmp_snpmeth.ref = 0;
                            tmp_snpmeth.alt = c[1];
                            tmp_snpmeth.offset = hap1_ins_len;// can check if ins_len > 4 if needed
                        }else if (hap1_ins_len == 0 && hap2_ins_len > 0 && hap2_ins_len <=4){
                            tmp_snpmeth.hap1 = 0; tmp_snpmeth.hap2 = 1;
                            tmp_snpmeth.ref = 0;
                            tmp_snpmeth.alt = c[2];
                            tmp_snpmeth.offset = hap2_ins_len;// can check if ins_len > 4 if needed
                        }else{
                            fprintf(stderr, "[%s] Error: hap1_mut = %d, hap2_mut = %d\n", __func__, hap1_mut, hap2_mut); exit(EXIT_FAILURE);
                        } // else: a position has different mut_type
                    }
                }
                snpmeth_map[i] = tmp_snpmeth;
                if(tmp_snpmeth.offset < 0){i=j-1;}
                tmp_snpmeth = {};
            }
        }
    }
}

void save_snpmeth(char *fname, char *chr_id, std::map<int, snpmeth_rec>& snpmeth_map)
{
    FILE* fp = fopen(fname, "a");
    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open methdb file: %s failed. Exit...\n", __func__, fname); exit(EXIT_FAILURE);}

    int i, indel_int, indel_len;
    snpmeth_rec tmp_snpmeth;
    for (std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it) {
        tmp_snpmeth = it->second;

        fprintf(fp, "%s\t%d\t", chr_id, it->first);
        if(tmp_snpmeth.offset == 0){
            fprintf(fp, "%c\t%c\t%d\t%d\t%d\t%d\t%.4f\n", "ACGTN"[tmp_snpmeth.ref], "ACGTN"[tmp_snpmeth.alt], tmp_snpmeth.hap1, tmp_snpmeth.hap2, 
                                                          tmp_snpmeth.is_phased, tmp_snpmeth.context[0], tmp_snpmeth.meth[0]);
        }else if (tmp_snpmeth.offset < 0){
            indel_len = abs(tmp_snpmeth.offset);
            indel_int = tmp_snpmeth.alt;
            for(i = 0; i < indel_len; ++i){putc("ACGTN"[indel_int & 0x3], fp); indel_int >>= 2;}
            fprintf(fp, "\t-\t%d\t%d\t%d\t", tmp_snpmeth.hap1, tmp_snpmeth.hap2, tmp_snpmeth.is_phased);
            for(i = 0; i < indel_len; ++i){fprintf(fp, "%d,", tmp_snpmeth.context[i]);}
            fprintf(fp, "\t");
            for(i = 0; i < indel_len; ++i){fprintf(fp, "%.4f,", tmp_snpmeth.meth[i]);}
            fprintf(fp, "\n");
        }else{
            indel_len = abs(tmp_snpmeth.offset);
            indel_int = tmp_snpmeth.alt;
            fprintf(fp, "-\t");
            for(i = 0; i < indel_len; ++i){putc("ACGTN"[indel_int & 0x3], fp);indel_int >>= 2;}
            fprintf(fp, "\t%d\t%d\t%d\t", tmp_snpmeth.hap1, tmp_snpmeth.hap2, tmp_snpmeth.is_phased);
            for(i = 0; i < indel_len; ++i){fprintf(fp, "%d,", tmp_snpmeth.context[i]);}
            fprintf(fp, "\t");
            for(i = 0; i < indel_len; ++i){fprintf(fp, "%.4f,", tmp_snpmeth.meth[i]);}
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}

int parse_snpmeth_line(char *line, char *chr_id, snpmeth_rec *tmp_snpmeth)
{
    // if (line == NULL) {fprintf(stderr, "[%s] ERROR: input line is null. Exit...\n", __func__); exit(EXIT_FAILURE);}
    char *p, *q= 0;
    char *contig, *ref_base, *alt_base, *context_str, *kmeridx_str, *meth_str;
	int i, ix, hap1, hap2, is_phased, pos;
    std::string tmp_str;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: contig  = p; break;
            case 1: pos     = atoi(p); break;
            case 2: ref_base= strdup(p); break; 
            case 3: alt_base= strdup(p); break;
            case 4: hap1    = atoi(p); break;
            case 5: hap2    = atoi(p); break;
            case 6: is_phased   = atoi(p); break;
            case 7: context_str = strdup(p); break;
            case 8: kmeridx_str = strdup(p); break;
            case 9: meth_str= strdup(p); break;
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return -1;}
			++i, p = q + 1;
			if (i > 9 || c == '\0') break;
		}
	}

    if(i < 9){
        fprintf(stderr, "[%s] Skip invalid site: chr %s, position %d...\n", __func__, chr_id, pos);
        free(ref_base); free(alt_base); free(context_str); free(kmeridx_str); free(meth_str);
        return 0;
    }

    tmp_snpmeth->hap1 = hap1;
    tmp_snpmeth->hap2 = hap2;
    tmp_snpmeth->is_phased = is_phased;

    if(strcmp(ref_base, "-")==0){       // insertion
        tmp_snpmeth->ref = 0;
        tmp_snpmeth->alt = 0;
        tmp_snpmeth->offset = strlen(alt_base);
        for(ix = 0; ix < tmp_snpmeth->offset; ++ix){
            tmp_snpmeth->alt = tmp_snpmeth->alt << 2 | nst_nt4_table[(int)alt_base[ix]];
        }
        for(i = 0, ix = 0; context_str[i] !='\0' && ix <= tmp_snpmeth->offset; ++i){
            if (context_str[i] == ',') {
                tmp_snpmeth->context[ix] = stoi(tmp_str);
                tmp_str.clear();
                ++ix; continue;
            } else {
                tmp_str += context_str[i];
            }
        }
        for(i = 0, ix = 0; kmeridx_str[i] !='\0' && ix <= tmp_snpmeth->offset; ++i){
            if (kmeridx_str[i] == ',') {
                tmp_snpmeth->kmeridx[ix] = stoi(tmp_str);
                tmp_str.clear();
                ++ix; continue;
            } else {
                tmp_str += kmeridx_str[i];
            }
        }
        for(i = 0, ix = 0; meth_str[i] !='\0' && ix <= tmp_snpmeth->offset; ++i){
            if (meth_str[i] == ',') {
                tmp_snpmeth->meth[ix] = stof(tmp_str);
                tmp_str.clear();
                ++ix; continue;
            } else {
                tmp_str += meth_str[i];
            }
        }
    }else if (strcmp(alt_base, "-")==0){// deletion
        tmp_snpmeth->ref = 0;
        tmp_snpmeth->alt = 0;
        tmp_snpmeth->offset = -strlen(ref_base);
        for(ix = 0; ix < -tmp_snpmeth->offset; ++ix){
            tmp_snpmeth->ref = tmp_snpmeth->ref << 2 | nst_nt4_table[(int)ref_base[ix]];
        }
    }else{                              // substitution
        tmp_snpmeth->ref = nst_nt4_table[(int)ref_base[0]];
        tmp_snpmeth->alt = nst_nt4_table[(int)alt_base[0]];
        tmp_snpmeth->offset = 0;
        for(i = 0, ix = 0; context_str[i] !='\0' && ix <= tmp_snpmeth->offset; ++i){
            if (context_str[i] == ',') {
                tmp_snpmeth->context[ix] = stoi(tmp_str);
                tmp_str.clear();
                ++ix; continue;
            } else {
                tmp_str += context_str[i];
            }
        }
        for(i = 0, ix = 0; kmeridx_str[i] !='\0' && ix <= tmp_snpmeth->offset; ++i){
            if (kmeridx_str[i] == ',') {
                tmp_snpmeth->kmeridx[ix] = stoi(tmp_str);
                tmp_str.clear();
                ++ix; continue;
            } else {
                tmp_str += kmeridx_str[i];
            }
        }
        for(i = 0, ix = 0; meth_str[i] !='\0' && ix <= tmp_snpmeth->offset; ++i){
            if (meth_str[i] == ',') {
                tmp_snpmeth->meth[ix] = stof(tmp_str);
                tmp_str.clear();
                ++ix; continue;
            } else {
                tmp_str += meth_str[i];
            }
        }
    }
    free(ref_base); free(alt_base); free(context_str); free(kmeridx_str); free(meth_str);
    return pos;
}

void load_snpmeth(char *fname, char *chr_id, uint32_t *posidx_arr, std::map<int, snpmeth_rec>& snpmeth_map)
{
    htsFile *fp = hts_open(fname,"rb");
    if(fp==0){fprintf(stderr, "[%s] ERROR: open snpmeth file: %s failed: %s. Exit...\n", __func__, fname, strerror(errno)); exit(EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    snpmeth_rec tmp_snpmeth;

    int ret, pos;
    int num_tot_site= 0;
    bool collect_present = false;
    bool collect_previous= false;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (!collect_present && collect_previous) { break; } // collect_present false && collect_previous true, finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        pos = parse_snpmeth_line(line.s, chr_id, &tmp_snpmeth);
        if(pos < 0){collect_present = false; continue;}
        
        ++num_tot_site;

        posidx_arr[pos] |= 1;
        snpmeth_map[pos] = tmp_snpmeth;
        tmp_snpmeth = {};
    }
    free(line.s);
    if((ret=hts_close(fp))){fprintf(stderr, "[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
    fprintf(stderr, "[%s] %s: found %d sites in snpmeth file...\n", __func__, chr_id, num_tot_site);
}
