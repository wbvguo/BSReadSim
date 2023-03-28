#ifndef HTSIM_H
#define HTSIM_H
#include <stdint.h>
#include <string.h>
#include <vector>
#include <map>

typedef unsigned short mut_t;

/*-------------------------variable-------------------------*/
// extern const uint8_t nst_nt4_table[256];
// extern const uint8_t cg_table[5];
// extern const uint8_t cg_context_table[64];
// extern std::map<std::string, int> base_map;
// extern std::map<std::string, int> context_map;
// extern std::map<int, int> params_map;


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
    float   score = 0;
} chr_rec;                  /*each struct take <= 16 bytes*/


// methdb
typedef struct {
    int pos = -1;
    float meth[2]  = {-1,-1};
    uint8_t context= 0;     /*1,3,7;9,11,15 for the context*/
    uint8_t type   = 0;     /*0,1,2,4,8 for cgmap, pool, asm, beta*/
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

#endif

