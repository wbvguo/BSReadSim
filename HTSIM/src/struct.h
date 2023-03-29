#ifndef STRUCT_H
#define STRUCT_H
#include <stdint.h>
#include <vector>


/*-------------------------variable-------------------------*/
extern const uint8_t nst_nt4_table[256];
extern const uint8_t cg_table[5];
extern const uint8_t cg_context_table[64];

// if the leftmost 4 bit is non-zero, then it must be snp or indel
enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
const mut_t mutmsk;


/*-------------------------struct-------------------------*/
// parse VCF
typedef struct {
    int pos, ref, alt, geno;
} snp_rec;


// parse BED
typedef struct {
    int pos_l, pos_r;
    float score;
    int8_t cut_l  = -1;
    int8_t cut_r  = -1;
    int8_t haplo  = -1;
    int8_t strand = -1;
} probe_rec;

typedef struct {
    char *name, *chr_id;
} probe_meta;


// parse RRBS
typedef struct {
    int len = -1;           /* length of cutting site */
    int idx = -1;           /* cutting position on *seq */
    std::vector<mut_t> seq; /* sequence encoded by numbers*/
} cut_rec;


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
    float    score  = 0;
} chr_rec;                  /*each struct take <= 16 bytes*/


// methdb
typedef struct {
    int pos = -1;
    float meth[2]  = {-1,-1};
    uint8_t context= 0;     /*1,3,7;9,11,15 for the context*/
    uint8_t type   = 0;     /*0,2,4,8,10 for uninitial, cgmap, asm, beta, pool*/
    // int16_t asm_ofs= 0;     /*0,1 for asm*/
} meth_rec;                 /*each struct take 15 bytes*/

typedef struct {
    float alpha = -1;
    float beta = -1;
} param_rec;                /*take 8 bytes*/


typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
} mutseq_t;


typedef struct {
    float ERR_RATE      = 0.005;
    float MUT_RATE      = 0.01;
    float INDEL_FRAC    = 0.15;
    float INDEL_EXTN    = 0.3;
} mut_params;

#endif

