#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>

typedef struct {
    int pos_l, pos_r;
    int len;
    int8_t cut_l = -1;
    int8_t cut_r = -1;      /* each struct take <= 16 bytes*/
    int8_t strand= 0;
} fragment;

int main(){
    fragment tmp_frag;
    printf("%d\n", sizeof(tmp_frag));
    tmp_frag.cut_l = 0;
    tmp_frag.cut_r = 0;
    tmp_frag.strand= -1;
    tmp_frag.pos_l = 1000000000;
    tmp_frag.pos_r = 2000000000;
    tmp_frag.len   = 1000000000;
    printf("%d\n", sizeof(tmp_frag));
}
