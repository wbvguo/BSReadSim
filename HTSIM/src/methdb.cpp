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
            tmp_meth.context[0] = cg_context_table[context_idx];
            tmp_meth.context[1] = cg_context_table[context_idx];
            tmp_meth.pos = i;
            meth_vec.push_back(tmp_meth);
            posidx_arr[i] |= ix << 2;
            ix++;
            tmp_meth = {};
        }
    }
}

void update_variant(const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, std::vector<meth_rec>& meth_vec, 
                    uint32_t *kmeridx_arr, meth_param *meth_set, std::vector<param_rec>& param_vec, std::map<int, snpmeth_rec>& snpmeth_map)
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
    std::vector<uint8_t>* sites_ptr;
    
    int i, j = 0; // j keeps the end of the last deletion
    int c[3];
    int pos_k, pos_base, pos_mask, pos_ins_len;
    int tmp_pos, tmp_base, tmp_mask, ins_len, ins_base, updown;
    int tmp_kmeridx, tmp_context;
    int count, t;
    bool collect_flag;

    //int tmp_idx, tmp_cg, tmp_kmeridx;
    for(int i=0; i < ks->seq.l; ++i){                       // skip contig boundary
        if((posidx_arr[i] & 0x1) != 0){                     // if there is a variant
            c[0] = nst_nt4_table[(int)ks->seq.s[i]];
            if (c[0] >= 4) continue;
            c[1] = hap1->s[i]; c[2] = hap2->s[i];
            if(c[0]==c[1]){rseq=hap2;}else{rseq=hap1;}      // can check MNV, and if they are both same with ref

            // handle the non-mutational boundary sites: calculate from [-3,3] from i
            // for each of them calulate the context and kmeridx(if needed)
            for(int k =-3; k<4 && k!=0; ++k){
                sites_u.clear();
                sites_d.clear();
                pos_k = i+k;
                pos_base = rseq->s[pos_k];
                if(pos_base&mutmsk) continue;               // skip mutation
                if(!cg_table[pos_base]) continue;           // skip nonCG

                // fill in vector
                #define __fill_vec(sites_ptr, updown)                       \
                    count=0, t=1;                                           \
                    collect_flag = true;                                    \
                    while(count < 3){                                       \
                        tmp_pos  = pos_k+t*updown;                          \
                        if(tmp_pos < 0 || tmp_pos >= ks->seq.l) break;      \
                        tmp_base = rseq->s[tmp_pos];                        \
                        tmp_mask = (tmp_base&mutmsk);                       \
                        ins_len  = tmp_mask >> 12;                          \
                        if(tmp_mask == DELETE){                             \
                        }else if (ins_len){                                 \
                            ins_base = tmp_base & 0xf;                      \
                            (*sites_ptr)[count]= ins_base;                  \
                            ++count;                                        \
                            if(count == 3) collect_flag = false;            \
                            ins_base >>= 4;                                 \
                            while (collect_flag && ins_len > 0){            \
                                (*sites_ptr)[count]= ins_base&0x3;          \
                                ++count;                                    \
                                ins_base >>= 2;                             \
                                --ins_len;                                  \
                                if(count == 3) collect_flag = false;        \
                            }                                               \
                            if(!collect_flag) break;                        \
                        }else{                                              \
                            (*sites_ptr)[count]= tmp_base;                  \
                            ++count;                                        \
                        }                                                   \
                        ++t;                                                \
                    }                                                       \
                
                sites_ptr = &sites_u; updown = -1;
                __fill_vec(sites_ptr, updown);

                sites_ptr = &sites_d; updown = 1;
                __fill_vec(sites_ptr, updown);

                if(kmeridx_arr != NULL){
                    //compute and update kmeridx
                    tmp_kmeridx = sites_u[2] << 6 | sites_u[1] << 4 | sites_u[0] << 2 | pos_base;
                    tmp_kmeridx = (tmp_kmeridx<<6)| sites_d[0] << 4 | sites_d[1] << 2 | sites_d[2];
                    //for (int k = 0; k < 3; ++k){tmp_kmeridx = (tmp_kmeridx << 2) | sites_u[3-k-1];}
                    if(tmp_kmeridx != (kmeridx_arr[pos_k] & 0xffff)){
                        kmeridx_arr[pos_k] = (kmeridx_arr[pos_k] &0x0000ffff) | (tmp_kmeridx << 16);
                    }
                }

                // compute and update context
                if (pos_base == 1){         // C
                    tmp_context = cg_context_table[((pos_base<<4) | (sites_d[0]<<2) | sites_d[1])]; 
                }else if (pos_base == 2){   // G
                    tmp_context = cg_context_table[((pos_base<<4) | (sites_u[0]<<2) | sites_u[1])];
                }else{} // should never happen

                int tmp_idx = posidx_arr[pos_k] >> 2;
                if(meth_vec[tmp_idx].context[0] != tmp_context){
                    meth_vec[tmp_idx].context[1] = tmp_context;
                    meth_vec[tmp_idx].meth[1]    = gen_beta(rng, tmp_context, param_vec);
                }
                fprintf(stderr, "ok\n");
            }

            // handle the snp sites
            pos_k = i;
            pos_base = rseq->s[i];
            pos_mask = (pos_base&mutmsk);
            if(pos_mask==DELETE) continue;         // DELETE will not appear on the read
            tmp_snpmeth= {};
            tmp_sites.clear();
            sites_u.clear();
            sites_d.clear();

            sites_ptr = &sites_u; updown = 1;
            __fill_vec(sites_ptr, updown);

            sites_ptr = &sites_d; updown = -1;
            __fill_vec(sites_ptr, updown);


            pos_ins_len= pos_mask >> 12;
            tmp_sites.push_back(sites_u[2]);
            tmp_sites.push_back(sites_u[1]);
            tmp_sites.push_back(sites_u[0]);
            tmp_sites.push_back(pos_base&0x3);
            pos_base >>=4;
            for (size_t i = 0; i < pos_ins_len; ++i){
                tmp_sites.push_back(pos_base&0x3);
                pos_base >> 2;
            }
            tmp_sites.push_back(sites_d[0]);
            tmp_sites.push_back(sites_d[1]);
            tmp_sites.push_back(sites_d[2]);
            
            for(int k =3; k<4+pos_ins_len; ++k){
                tmp_kmeridx = (tmp_sites[k-3] << 6) | (tmp_sites[k-2] << 4) | (tmp_sites[k-1] << 2) | tmp_sites[k];
                tmp_kmeridx = (tmp_kmeridx << 6) | (tmp_sites[k+1] << 4) | (tmp_sites[k+2] << 2) | tmp_sites[k+3];
                
                tmp_base = tmp_sites[k];
                if(tmp_base == 1){         // C
                    tmp_context = cg_context_table[((tmp_base<<4) | (tmp_sites[k+1]<<2) | tmp_sites[k+2])]; 
                }else if (tmp_base == 2){   // G
                    tmp_context = cg_context_table[((tmp_base<<4) | (tmp_sites[k-1]<<2) | tmp_sites[k-2])];
                }else{} // should never happen

                tmp_snpmeth.kmeridx.push_back(tmp_kmeridx);
                tmp_snpmeth.context.push_back(tmp_context);
                tmp_snpmeth.meth.push_back(gen_beta(rng, tmp_context, param_vec));
            }
            
            tmp_snpmeth.ref = c[0];
            tmp_snpmeth.alt = tmp_base;
            tmp_snpmeth.geno= 1 + (int)(c[1] == c[2]);
            snpmeth_map[i] = tmp_snpmeth;
        }
    }

    if(meth_set->is_asm_set){                                       // use the last bit to store asm signal
        for(int i=0; i < ks->seq.l; ++i){
            posidx_arr[i] &= 0xfffffffe;
            int tmp_pos = posidx_arr[i] >> 2;
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
        fprintf(fp, "%d\t%f\t%f\t%d\t%d\t%d\n", tmp_meth.pos, tmp_meth.meth[0], tmp_meth.meth[1], tmp_meth.context[0], tmp_meth.context[1], tmp_meth.type);
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
    int context1=0, context2=0, type=0;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: pos     = atoi(p); break;
            case 1: ref_meth= atof(p); break;
            case 2: alt_meth= atof(p); break;
            case 3: context1= atoi(p); break;
            case 4: context2= atoi(p); break;
            case 5: type    = atoi(p); break;
            default: break;}
			++i, p = q + 1;
			if (i > 5 || c == '\0') break;
		}
	}

    if(i < 4){fprintf(stderr, "[%s] Skip invalid site: position %d...\n", __func__, pos);}

    tmp_meth->pos       = pos;
    tmp_meth->meth[0]   = ref_meth;
    tmp_meth->meth[1]   = alt_meth;
    tmp_meth->context[0]= context1;
    tmp_meth->context[1]= context2;
    tmp_meth->type      = type;
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

        if(meth_vec[tmp_idx].pos!=tmp_meth.pos && meth_vec[tmp_idx].context[0]!=tmp_meth.context[0]){ ++num_404_site; continue;}
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
    tmp_meth->context[0]  = base_map[base] << 3 | context_map[std::string(context)];
    tmp_meth->context[1]  = tmp_meth->context[0];
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
        if(meth_vec[tmp_idx].pos!=tmp_meth.pos && meth_vec[tmp_idx].context[0]!=tmp_meth.context[0]){ ++num_404_site; continue;}
        meth_vec[tmp_idx].meth[0]= tmp_meth.meth[0];
        meth_vec[tmp_idx].context[0] = tmp_meth.context[0];
        meth_vec[tmp_idx].context[1] = tmp_meth.context[1];
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
            meth_vec[i].meth[0] = gen_beta(rng, meth_vec[i].context[0], param_vec);
            meth_vec[i].meth[1] = meth_vec[i].meth[0];
            meth_vec[i].type    = 8;
        }
    }
    gsl_rng_free(rng);
}

