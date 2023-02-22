#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    int pos_l, pos_r;       /*int max value is 2147483647*/
    float cg_ratio= -1;
    int8_t cut_l  = -1;
    int8_t cut_r  = -1;
    int8_t index  = -1;     /*save for capture efficiency usage*/
    int8_t ns     = 0;      /*strand or # of cut sites contained*/
} fragment;                 /*each struct take <= 16 bytes*/

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
}
