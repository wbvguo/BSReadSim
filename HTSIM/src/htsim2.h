/* The MIT License

   Copyright (c) 2022 Wenbin Guo <wbguo@ucla.edu>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* Last Modified: 20MAR2023 */


#ifndef HTSIM_H
#define HTSIM_H
#include <stdint.h>
#include <string.h>
#include <vector>
#include <map>



/*-------------------------variable-------------------------*/
extern const uint8_t nst_nt4_table[256];
extern const uint8_t cg_table[5];
extern const uint8_t cg_context_table[64];
extern std::map<std::string, int> base_map;
extern std::map<std::string, int> context_map;
extern std::map<int, int> params_map;

typedef unsigned short mut_t;

/*-------------------------struct-------------------------*/
// parse VCF
typedef struct {
    int pos, ref, alt, geno;
} snp_rec;


// parse BED
typedef struct {
    int pos_l, pos_r;
    float score;
    int8_t strand;
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

typedef struct {
    int pos;
    int8_t type = -1;
} cut_pos;


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
} chr_rec;

// save for methylable bases
typedef struct {
    int pos = -1;
    float ref_meth = -1;    /*save meth if not is_asm*/
    float alt_meth = -1;
    uint8_t context= 0;     /*1,3,7;9,11,15*/
    uint8_t type   = 0;     /*0,1,2,4,8*/
} meth_rec;                 /*each struct take 14 bytes*/

typedef struct {
    float alpha = -1;
    float beta = -1;
} param_rec;                /*take 8 bytes*/


#endif
