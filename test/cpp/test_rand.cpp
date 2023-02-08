#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <ctype.h>
#include <string.h>
#include <zlib.h>
#include <random>
#include <vector>
#include <algorithm>



typedef struct {
    int pos_l, pos_r;       /*int max value is 2147483647*/
    float cg_ratio= -1;
    int8_t cut_l  = -1;
    int8_t cut_r  = -1;
    int8_t index  = -1;     /*save for capture efficiency usage*/
    int8_t ns     = 0;      /*strand or # of cut sites contained*/
} fragment;                 /*each struct take <= 16 bytes*/

int min_insert =0, max_insert = 750;
int mean_insert = 300, sd_insert = 50, sd_center =50;
int n_pairs = 100000000;
int vec_size= 10000000;
int contig_len = 240000000;

// initialize random generator for general distributions
std::random_device rd;
std::mt19937 gen(rd());
// initialize random generator for standard normal distribution 
std::random_device rn;  //Will be used to obtain a seed for the random number engine
std::mt19937 gen_rn(rn()); //Standard mersenne_twister_engine seeded with rd()
std::normal_distribution<float> dis_rn(0.0, 1.0);
std::random_device ru;
std::mt19937 gen_ru(ru());
std::uniform_real_distribution<float> dis_ru(0.0,1.0);


typedef struct {
    int pos_l, pos_r;
    float score;
    int8_t strand;
} probe_rec;


void sim_rand(std::uniform_int_distribution<int> *dis, fragment *tmp_frag)
{
    // targeted sequencing (uniform)
    int frag_idx = (*dis)(gen);
    
    probe_rec tmp_probe;
    tmp_probe.pos_l = frag_idx;
    tmp_probe.pos_r = frag_idx + 1000;

    int probe_center = (int) (tmp_probe.pos_l + tmp_probe.pos_r)/2;
    int frag_center= probe_center + (int)(sd_center * dis_rn(gen_rn));
    int insert_dev = (int)(sd_insert * dis_rn(gen_rn));
    int insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));

    int pos_l = frag_center - (int)insert_len/2;
    int pos_r = frag_center + (int)insert_len/2; 
    tmp_frag->pos_l = pos_l;
    tmp_frag->pos_r = pos_r;
}

void sim_rand(std::discrete_distribution<int> *dis, fragment *tmp_frag)
{
    // targeted sequencing (non-uniform)
    int frag_idx = (*dis)(gen);
    
    probe_rec tmp_probe;
    tmp_probe.pos_l = frag_idx;
    tmp_probe.pos_r = frag_idx + 1000;

    int probe_center = (int) (tmp_probe.pos_l + tmp_probe.pos_r)/2;
    int frag_center= probe_center + (int)(sd_center * dis_rn(gen_rn));
    int insert_dev = (int)(sd_insert * dis_rn(gen_rn));
    int insert_len = std::max(min_insert, std::min(mean_insert + insert_dev, max_insert));

    int pos_l = frag_center - (int)insert_len/2;
    int pos_r = frag_center + (int)insert_len/2; 
    tmp_frag->pos_l = pos_l;
    tmp_frag->pos_r = pos_r;
}

int main(){
    fragment tmp_frag;
    int ii;
    int tech_mode = 0;
    bool is_uniform = false;
    int unif_begin, unif_end;
    std::uniform_int_distribution<int> dis;
    std::uniform_int_distribution<int> dis2;


    if(is_uniform){
        unif_begin = 0;
        unif_end = contig_len-1;     // might be empty
        std::uniform_int_distribution<int> dis(unif_begin, unif_end);
    }else{
        std::vector<float> weights;
        for(int i=0; i < vec_size; ++i){ weights.push_back(dis_ru(gen_ru));}
        std::discrete_distribution<int> dis2(weights.begin(), weights.end());
        std::vector<float>().swap(weights);
    }
    

    for (ii = 0; ii != n_pairs; ++ii) { // the core loop
        sim_rand(&dis2, &tmp_frag);
    }
}