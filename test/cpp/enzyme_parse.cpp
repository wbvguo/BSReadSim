#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>

typedef struct {
    int len = -1;         /* length of cutting site */
    int idx = -1;         /* cutting position on *seq */
    std::vector<int> seq; /* sequence encoded by numbers*/
} cut_t;

typedef struct {
    int pos_l;
    int pos_r;
    int cut_l;
    int cut_r;
    int len;
} cut_frag;


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


int main(){
    char a_str[] = "C,CGG;T,CGA;T,GCA";
    char *str = a_str;
    int array[] = {0,0,0,0, 1,1,2,2, 3,2,1,0, 0,0,0,0, 3,1,2,0};
    
    // parse the cutting sites
    int c, idx = 0;
    cut_t tmp_site;
    std::vector<cut_t> cut_site;
    for(int i =0; str[i] !='\0'; ++i){
        if (str[i] == ','){tmp_site.idx = idx; continue;}
        if (str[i] == ';'){
            if (tmp_site.idx >=0) { 
                tmp_site.len = idx;
                cut_site.push_back(tmp_site);
                // restart, empty the struct
                tmp_site = {}; idx = 0;
                continue;
            } else {
                fprintf (stderr, "ERROR: cutting position issue in Input enzyme site (%s). Exit... \n", str); 
                exit (EXIT_FAILURE);
            }
        }
        ++idx;
        c = nst_nt4_table[(int)str[i]];
        tmp_site.seq.push_back(c);
    }
    if(tmp_site.idx >=0){tmp_site.len = idx; cut_site.push_back(tmp_site);}


    // find locations, add to site_pos
    std::vector<uint8_t> cut_type;
    std::vector<int> cut_pos;

    int *ptr_begin = std::begin(array);
    int *ptr_end = std::end(array);
    int *iter = std::begin(array);
    int *iter2= std::begin(array);
    int *iter3= std::begin(array);
    // printf("%d %d\n", ptr_begin, ptr_end);

    int count = 0;
    while (iter < ptr_end)
    {
        // printf("========%i========\n", count);
        // printf("%d %d %d\n", iter, iter2, iter3);
        //initiate
        int min_pos = ptr_end - ptr_begin;
        int site_pos= min_pos;
        int type_idx = -1;

        // printf("[%d %d %d]\n", min_pos, site_pos, type_idx);
        for (int j = 0; j < cut_site.size(); j++) {
            iter2 = std::search(iter, ptr_end, cut_site[j].seq.begin(), cut_site[j].seq.end());
            // printf("[[%d %d]]\n", iter, iter2);
            site_pos = iter2 - ptr_begin;
            if(site_pos < min_pos){
                min_pos = site_pos;
                type_idx= j;
                iter3   = iter2;
            }
        }
        // printf("[%d %d %d]\n", min_pos, site_pos, type_idx);
        // printf("%d %d %d\n", iter, iter2, iter3);

        if(type_idx >=0){
            cut_pos.push_back(min_pos + cut_site[type_idx].idx);
            cut_type.push_back(type_idx);
            iter = iter3 + cut_site[type_idx].len; // need to check for other int types
        }
        // printf("%d %d %d\n", iter, iter2, iter3);
        // ++count;
        // if(count >= 20){
        //     break;
        // }
    }

    // generate potential intervals
    std::vector<cut_frag> frag_vec;
    cut_frag tmp_frag;
    int MIN_FRAG_LEN = 0;
    int MAX_FRAG_LEN = 50;
    int contig_len = sizeof(array)/sizeof(array[0]);
    //printf("%ld\n", cut_pos.size());
    printf("hello\n");
    

    for (int i = -1; i <= int(cut_pos.size()); i++){
        printf("%d\n", i);
        if(i==-1 || i == cut_pos.size()){ //append the first and last fragments
            if(i==-1){
                tmp_frag.pos_l = 0;
                tmp_frag.cut_l = -1;
                tmp_frag.pos_r = cut_pos[0];
                tmp_frag.cut_r = cut_type[0];
                tmp_frag.len   = cut_pos[0] - 0;
            }else{
                tmp_frag.pos_l = cut_pos[i-1];
                tmp_frag.cut_l = cut_type[i-1];
                tmp_frag.pos_r = contig_len;
                tmp_frag.cut_r = -1;
                tmp_frag.len   = contig_len - cut_pos[i-1];
            }
            frag_vec.push_back(tmp_frag);
            tmp_frag = {};
            continue;
        }
        for(int j=i+1; j < cut_pos.size(); j++){
            int frag_len = cut_pos[j] - cut_pos[i];
            if (frag_len > MIN_FRAG_LEN && frag_len < MAX_FRAG_LEN){
                tmp_frag.pos_l = cut_pos[i];
                tmp_frag.cut_l = cut_type[i];
                tmp_frag.pos_r = cut_pos[j];
                tmp_frag.cut_r = cut_type[j];
                tmp_frag.len   = frag_len;

                frag_vec.push_back(tmp_frag);
                tmp_frag = {};
            }
        }
    }
    

    // print the enzyme sites vectors
    // for (int i = 0; i < cut_site.size(); i++) {
    //     printf("%d\n", cut_site[i].idx);
    //     printf("%d\n", cut_site[i].len);
    //     for(int j =0; j < cut_site[i].seq.size(); ++j){
    //         printf("%d ", cut_site[i].seq[j]);
    //     }
    //     printf("\n");
    // }

    // print all the cut sites
    // for (int i = 0; i < cut_pos.size(); i++) {
    //     printf("%d %d\n", cut_pos[i], cut_type[i]);
    // }

    // //print all possible fragments
    for (int i = 0; i < frag_vec.size(); i++) {
        printf("%d %d %d %d %d\n",  frag_vec[i].pos_l, frag_vec[i].pos_r, 
                                    frag_vec[i].cut_l, frag_vec[i].cut_r, frag_vec[i].len);
    }

    return 0;
}


