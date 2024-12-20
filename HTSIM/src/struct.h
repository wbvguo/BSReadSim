#ifndef STRUCT_H
#define STRUCT_H

#include <stdint.h>
#include <vector>
#include <map>


/*-------------------------variable-------------------------*/
extern const char PACKAGE_VERSION[];
extern const char OUTPUT_PREFIX[];
extern const char PARAM_DEFAULT[];

extern const uint8_t nst_nt4_table[256];
extern const uint8_t cg_table[5];
extern const uint8_t cg_context_table[64];

extern const float ERR_RATE;
extern const float MUT_RATE;
extern const float INDEL_FRAC;
extern const float INDEL_EXTN;
extern const float MAX_N_RATIO;

extern const int DEPTH;
extern const int MEAN_INSERT;
extern const int SD_INSERT;
extern const int MIN_INSERT;
extern const int MAX_INSERT;
extern const int SIZE_L;
extern const int SIZE_R;
extern const int SD_CENTER;
extern const int BIN_SIZE;
extern const int TECH_MODE;
extern const int OUTPUT_FMT;
extern const int CHUNK_SIZE;
extern const int BUFFER_SIZE;

// if the leftmost 4 bit is non-zero, then it must be snp or indel
enum muttype_t {NOCHANGE = 0, INSERT = 0x1000, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
extern const mut_t mutmsk;


/*-------------------------struct-------------------------*/
// parse VCF and hold the methylation information
typedef struct {
    float meth[4];
    uint16_t kmeridx[4];
    uint16_t ref, alt;
    uint8_t context[4];
    int8_t hap1, hap2, is_phased, offset;
} snpmeth_rec;                  /*each struct take 36 bytes*/   


// parse BED
typedef struct {
    char *name, *chr_id;
} probe_meta;


// parse RRBS
typedef struct {
    int len = -1;               /* length of cutting site */
    int idx = -1;               /* cutting position on *seq */
    std::vector<mut_t> seq;     /* sequence encoded by numbers*/
} cut_rec;

typedef struct cutpos_rec {
    int pos;
    int8_t haplo= 0;            /*haplotype 1: hap1, 2: hap2, 3: both*/
    int8_t type = -1;
    bool operator<(const cutpos_rec& other) const {
        return pos < other.pos;
    }
} cutpos_rec;

typedef struct frag_rrbs_rec{
    int pos_l,pos_r;            /*int max value is 2147483647*/
    float  score    = 0;        /*will alternatively be used as start2*/
    int8_t haplo    = 0;        /*haplotype 1: hap1, 2: hap2, 3: both*/
    int8_t strand   = 0;        /*strand 1: watson, 2: crick, 3: both*/
    int8_t cut_l    = -1;
    int8_t cut_r    = -1;
    std::map<int8_t, int8_t> n_cuts;/*number of cuts contained*/
    bool operator<(const frag_rrbs_rec& other) const {
        if (pos_l != other.pos_l) {
            return pos_l < other.pos_l;
        } else {
            return pos_r < other.pos_r;
        }
    }
} frag_rrbs_rec;                /*each struct take 64 bytes*/


// fragments
typedef struct frag_rec {
    int pos_l,pos_r;            /*int max value is 2147483647*/
    float   score   = 0;        /*will alternatively be used as start2*/
    int8_t  haplo   = 0;        /*haplotype 1: hap1, 2: hap2, 3: both*/
    int8_t  strand  = 0;        /*strand 1: watson, 2: crick, 3: both*/   
    int8_t  cut_l   = -1;
    int8_t  cut_r   = -1; 
    bool operator<(const frag_rec& other) const {
        if (haplo != other.haplo) {
            return haplo < other.haplo;
        } else {
            return pos_l < other.pos_l;
        }
    }
} frag_rec;                     /*each struct take <= 16 bytes*/


// chr counts
typedef struct {
    uint32_t chr_len= 0;
    uint32_t eff_len= 0;
    uint32_t count  = 0;
    float    score  = 0;
} chr_rec;                      /*each struct take <= 16 bytes*/


// methdb
typedef struct {
    int pos = -1;
    float meth[2]  = {-1,-1};
    uint16_t kmeridx[2]= {0,0};
    uint8_t context[2] = {0,0}; /*1,3,7;9,11,15 for the context*/
    uint8_t type   = 0;         /*0,2,4,8,10 for uninitial, cgmap, asm, beta, pool*/
    uint8_t is_asm = 0;         /*0,1 for asm*/
} meth_rec;                     /*each struct take 15 bytes*/

typedef struct {
    float alpha = -1;
    float beta = -1;
} param_rec;                    /*take 8 bytes*/


// haplotypes
typedef struct {
    int l, m; /* length and maximum buffer size */
    mut_t *s; /* sequence */
} mutseq_t;


// hold all parameters
typedef struct {
    float err_rate      = ERR_RATE;
    float maxN_ratio    = MAX_N_RATIO;
    float depth         = DEPTH;

    int mean_insert     = MEAN_INSERT;
    int sd_insert       = SD_INSERT;
    int min_insert      = MIN_INSERT;
    int max_insert      = MAX_INSERT;
    int size_l          = SIZE_L;
    int size_r          = SIZE_R;
    int max_length      = SIZE_L;
    int half_size       = SIZE_L >> 1;
    int arr_size        = SIZE_L << 2; // because int takes 4 bytes
    int flip_thre       = RAND_MAX >> 1;
    int err_thre        = (int) RAND_MAX * ERR_RATE;

    int sd_center       = SD_CENTER;
    int bin_size        = BIN_SIZE;
    int chunk_size      = CHUNK_SIZE;

    int tech_mode       = TECH_MODE;
    int output_fmt      = OUTPUT_FMT;
    bool is_chr_set     = false;
    bool is_uniform     = true;
    bool is_bias_set    = false;
    bool is_site_set    = false;
    bool is_bed_set     = false;
    bool is_kmer_set    = false;
} expt_param;

typedef struct {
    float mut_rate      = MUT_RATE;
    float indel_frac    = INDEL_FRAC;
    float indel_extn    = INDEL_EXTN;
    int   seed_snp      = -1;
    bool  is_hap        = false;
    bool  is_vcf_set    = false;
    bool  is_snpmeth_set= false;
    bool  snpmeth_save  = false;
} mut_param;

typedef struct{
    int  seed_meth      = -1;
    bool cgmap_pool     = false;
    bool methdb_save    = false;
    bool is_asm_set     = false;
    bool is_cgmap_set   = false;
    bool is_methdb_set  = false;
    bool update_meth    = false;
} meth_param;

#endif
