#include <vector>
#include <map>
#include <algorithm>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"


KSEQ_INIT(gzFile, gzread)
#include "rrcut.h"


void parse_cut_site(char *cut_str, std::vector<cut_rec>& cut_vec)
{
    // check if cut_str is valid (seprate by ',' and '|' intercjangeably, 
    // never || without comma in between or ,, without | in between)
    bool vline_flag = false;
    bool comma_flag = false;

    for(int i =0; cut_str[i] !='\0'; ++i){
        if(cut_str[i] == ' '){continue;}
        if(cut_str[i] == '|'){
            if(vline_flag){fprintf (stderr, "[%s] ERROR: Invalid enzyme site format '%s'. Exit... \n", __func__, cut_str); exit (EXIT_FAILURE);}
            vline_flag = true; comma_flag = false;
        }else if(cut_str[i] == ','){
            if(comma_flag){fprintf (stderr, "[%s] ERROR: Invalid enzyme site format '%s'. Exit... \n", __func__, cut_str); exit (EXIT_FAILURE);}
            comma_flag = true; vline_flag = false;
        }else{
            if(nst_nt4_table[(int)cut_str[i]]==4 && (cut_str[i] != 'N' || cut_str[i] != 'n')){
                fprintf (stderr, "[%s] ERROR: Invalid enzyme site format '%s'. Exit... \n", __func__, cut_str); exit (EXIT_FAILURE);
            }
        }
    }


    cut_vec.clear();
    int c, idx = 0;
    cut_rec tmp_site;
    
    for(int i =0; cut_str[i] !='\0'; ++i){
        if(cut_str[i] == ' '){continue;}
        if (cut_str[i] == '|'){tmp_site.idx = idx; continue;}
        if (cut_str[i] == ','){
            if (tmp_site.idx >=0) {
                tmp_site.len = idx;
                cut_vec.push_back(tmp_site);
                tmp_site = {}; idx = 0; continue; // empty the struct, restart
            } else {
                fprintf (stderr, "[%s] ERROR: Invalid enzyme site format (%s). Exit... \n", __func__, cut_str); exit (EXIT_FAILURE);
            }
        }
        ++idx;
        c = nst_nt4_table[(int)cut_str[i]];
        tmp_site.seq.push_back(c);
    }
    if(tmp_site.idx >=0){tmp_site.len = idx; cut_vec.push_back(tmp_site);}
}

void gen_cut_pos(mutseq_t *hap1, mutseq_t *hap2, std::vector<cutpos_rec>& cutpos_vec, std::vector<cut_rec>& cut_vec)
{
    cutpos_vec.clear(); // clean the container
    
    std::map<cutpos_rec, int> cutpos_map;
    cutpos_rec tmp_cutpos;
    cut_rec tmp_cut;
    mut_t *arr_ptr;

    
    for (int j = 0; j < cut_vec.size(); ++j) {
        tmp_cut = cut_vec[j];

        // Search for hap1, currently don't handle cut sites induced by indel, only SNP
        arr_ptr = hap1->s;
        for (int i = 0; i < hap1->l; ++i) {
            auto found = std::search(arr_ptr + i, arr_ptr + hap1->l, tmp_cut.seq.begin(), tmp_cut.seq.end());
            if(found == arr_ptr + hap1->l){ break;}
            tmp_cutpos = {};
            tmp_cutpos.pos  = found - arr_ptr;
            tmp_cutpos.type = j;
            cutpos_map[tmp_cutpos] = 1;
        }

        // Search for hap2
        arr_ptr = hap2->s;
        for (int i = 0; i < hap2->l; ++i) {
            auto found = std::search(arr_ptr + i, arr_ptr + hap2->l, tmp_cut.seq.begin(), tmp_cut.seq.end());
            if(found == arr_ptr + hap2->l){ break;}
            tmp_cutpos = {};
            tmp_cutpos.pos  = found - arr_ptr;
            tmp_cutpos.type = j;
            cutpos_map[tmp_cutpos] |= 0x2; // default is 0 if the elements not in the map, |=0x2 will give 2
        }
    }

    if(cutpos_map.empty()){fprintf (stderr, "[%s] ERROR: No cut site found in hap1 and hap2. Exit... \n", __func__); exit (EXIT_FAILURE);}

    // put into the vector
    for (auto it = cutpos_map.begin(); it != cutpos_map.end(); ++it) {
        tmp_cutpos = it->first;
        tmp_cutpos.haplo = it->second; // 1: hap1, 2: hap2, 3: both
        cutpos_vec.push_back(tmp_cutpos);
    }

    std::sort(cutpos_vec.begin(), cutpos_vec.end());
}

void gen_cut_frag(const kseq_t *ks, expt_param *expt_set, std::vector<frag_rrbs_rec> &frag_vec, std::vector<cutpos_rec>& cutpos_vec, std::vector<cut_rec>& cut_vec)
{
    frag_vec.clear();

    // generate potential intervals
    int frag_len;
    frag_rrbs_rec tmp_frag;

    std::map<int8_t, int8_t> cut_map;
    std::map<int8_t, int8_t> tmp_cutmap;
    for(size_t i=0; i<cut_vec.size(); ++i){cut_map[(int8_t)i] = 0;}
    
    for (int i = -1; i <= int(cutpos_vec.size()); ++i){
        if(i==-1 || i == int(cutpos_vec.size())){ //check the first and last frag_recs
            tmp_frag = {};
            if(i==-1){
                frag_len = cutpos_vec[0].pos;
                if (frag_len > expt_set->min_insert && frag_len < expt_set->max_insert){
                    tmp_frag.pos_l = 0;
                    tmp_frag.pos_r = cutpos_vec[0].pos;
                    tmp_frag.cut_l = 0;
                    tmp_frag.cut_r = cutpos_vec[0].type;
                    tmp_frag.haplo = cutpos_vec[0].haplo;
                    ++tmp_frag.n_cuts[cutpos_vec[0].type];
                    frag_vec.push_back(tmp_frag);
                }
            }else{
                frag_len = ks->seq.l - cutpos_vec[i-1].pos;
                if (frag_len > expt_set->min_insert && frag_len < expt_set->max_insert){
                    tmp_frag.pos_l = cutpos_vec[i-1].pos;
                    tmp_frag.pos_r = ks->seq.l;
                    tmp_frag.cut_l = cutpos_vec[i-1].type;
                    tmp_frag.cut_r = 0;
                    tmp_frag.haplo = cutpos_vec[i-1].haplo;
                    ++tmp_frag.n_cuts[cutpos_vec[i-1].type];
                    frag_vec.push_back(tmp_frag);
                }
            }
            continue;
        }

        // generate fragments that have middle cuts while take haplotype into consideration
        // 3,3,2,1,2,1,3...
        // +,+,+,-,+,-,+...
        // +,+,-,+,-,+,+...
        //   +,-,+,-,+,+...
        //   +,+,-,+,-,+...
        //     +,-,+,-,+...
        //       +,-,+,+...
        //         +,-,+...
        //           +,+...
        // should have better way to achieve it though

        int haplotype   = cutpos_vec[i].haplo;
        int8_t left_cut = cutpos_vec[i].type;
        
        std::map<frag_rrbs_rec, int> frag_map;
        frag_map.clear();

        // first hap1
        tmp_cutmap = cut_map;
        ++tmp_cutmap[cutpos_vec[i].type];
        for(int j=i+1; j < cutpos_vec.size(); ++j){
            if((cutpos_vec[j].haplo & 0x2) == 0){continue;} // skip hap2-only cuts
            frag_len = cutpos_vec[j].pos - cutpos_vec[i].pos;
            if(frag_len < expt_set->min_insert){continue;} // skip the small fragments (too short)
            if(frag_len > expt_set->max_insert){break;}
            haplotype &= cutpos_vec[j].haplo;
            ++tmp_cutmap[cutpos_vec[j].type];

            tmp_frag = {};
            tmp_frag.pos_l = cutpos_vec[i].pos;
            tmp_frag.pos_r = cutpos_vec[j].pos;
            tmp_frag.cut_l = cutpos_vec[i].type;
            tmp_frag.cut_r = cutpos_vec[j].type;
            tmp_frag.haplo = haplotype;
            tmp_frag.n_cuts= tmp_cutmap;
            ++frag_map[tmp_frag];
        }

        // then hap2
        tmp_cutmap = cut_map;
        ++tmp_cutmap[cutpos_vec[i].type];
        for(int j=i+1; j < cutpos_vec.size(); ++j){    
            if((cutpos_vec[j].haplo & 0x1) == 0){continue;} // skip hap1-only cuts
            frag_len = cutpos_vec[j].pos - cutpos_vec[i].pos;
            if(frag_len < expt_set->min_insert){continue;} // skip the small fragments (too short)
            if(frag_len > expt_set->max_insert){break;}
            haplotype &= cutpos_vec[j].haplo;
            ++tmp_cutmap[cutpos_vec[j].type];

            tmp_frag = {};
            tmp_frag.pos_l = cutpos_vec[i].pos;
            tmp_frag.pos_r = cutpos_vec[j].pos;
            tmp_frag.cut_l = cutpos_vec[i].type;
            tmp_frag.cut_r = cutpos_vec[j].type;
            tmp_frag.haplo = haplotype;
            tmp_frag.n_cuts= tmp_cutmap;
            ++frag_map[tmp_frag];
        }
        
        // put into the vector
        for (auto it = frag_map.begin(); it != frag_map.end(); ++it) {
            tmp_frag = it->first;
            frag_vec.push_back(tmp_frag);
        }
    }
}

void output_rrcut_bed(const char *fname, const char *chr_id, std::vector<frag_rrbs_rec> &frag_vec)
{
    FILE* fp = fopen(fname, "w");
    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open rrbs bed file: %s failed. Exit... \n", __func__, fname); exit (EXIT_FAILURE);}

    frag_rrbs_rec tmp_frag;
    int map_size = frag_vec[0].n_cuts.size();
    for (size_t i=1; i < frag_vec.size(); ++i){
        tmp_frag = frag_vec[i];
        fprintf(fp, "%s\t%d\t%d\t.\t1\t%d\t%d\t%d", chr_id, tmp_frag.pos_l, tmp_frag.pos_r, tmp_frag.haplo, tmp_frag.cut_l, tmp_frag.cut_r);
        for(int j=0; j < map_size; ++j){fprintf(fp, "\t%d", tmp_frag.n_cuts[j]);}
        fprintf(fp, "\n");
        if (ferror(fp)) {fprintf(stderr, "[%s] ERROR: failed to write to file %s. Exit... \n", __func__, fname);exit(EXIT_FAILURE);}
    }
    fclose(fp);
}

