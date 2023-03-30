#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <vector>
#include <map>

typedef struct {
    int pos_l, pos_r;       /*int max value is 2147483647*/
    float cg_ratio= -1;
    int8_t cut_l  = -1;
    int8_t cut_r  = -1;
    int8_t index  = -1;     /*save for capture efficiency usage*/
    int8_t ns     = 0;      /*strand or # of cut sites contained*/
} fragment;                 /*each struct take <= 16 bytes*/


typedef struct {
    float mut_rate  = 0.01;
    float indel_frac= 0.15;
    float indel_extn= 0.3;
    int  seed_snp   = -1;
    bool is_hap     = false;
    bool is_vcf_set = false;
} mut_param;



typedef struct {
    int pos_l,pos_r,start2;     /*int max value is 2147483647*/
    uint8_t haplo   = 0;        /*haplotype 1: hap1, 2: hap2, 3: both*/
    uint8_t strand  = 0;        /*strand 1: watson, 2: crick, 3: both*/
    int8_t cut_l    = -1;
    int8_t cut_r    = -1;
    std::map<uint8_t, uint8_t> n_cuts;/*number of cuts contained*/
} frag_rec;                     /*each struct take <= 16 bytes*/


int main(){
    fragment tmp_frag;
    printf("%d\n", tmp_frag.index);
    printf("%d\n", sizeof(tmp_frag));
    tmp_frag.cut_l = 0;
    tmp_frag.cut_r = 0;
    tmp_frag.ns = -1;
    tmp_frag.pos_l = 1000000000;
    tmp_frag.pos_r = 2000000000;
    tmp_frag.cg_ratio = 0.88888;
    tmp_frag.index = 1;
    printf("%d\n", sizeof(tmp_frag));
    printf("%d\n", tmp_frag.index);
    tmp_frag = {};
    printf("%d\n", tmp_frag.index);

    mut_param tmp_mut;
    printf("%d\n", sizeof(tmp_mut));

    frag_rec tmp_frag2;
    printf("%d\n", sizeof(tmp_frag2));
}
