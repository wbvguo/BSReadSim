#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include "kseq.h"
#include "vcf.h"

typedef struct {
    int pos_l, pos_r;
    int strand;
    float score;
} probe_rec;

typedef struct {
    char *name;
    char *contig;
} probe_meta;

std::vector<probe_rec> probe_vec;


//parse_bed: https://github.com/dhspence/tagbam/blob/main/tagbam.c
char *parse_bed(char *s, probe_rec *tmp_probe, probe_meta *tmp_probe_meta)
{
	char *p, *q, *contig_id, *name = 0;
    int i, start, end, strand;
	float score;
    // int32_t i, start = -1, end = -1;

	for (i = 0, p = q = s;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: contig_id = p; break;
            case 1: start= atoi(p); break;
            case 2: end  = atoi(p); break;
            case 3: name = strdup(p); break;
            case 4: score= atof(p); break;
            case 5: strand = int(strcmp(p,"+")==0)-int(strcmp(p,"-")==0); break;
            default: break;}
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

	tmp_probe->pos_l = start;
	tmp_probe->pos_r = end;
    tmp_probe->score = score;
    tmp_probe->strand= strand;
    
    tmp_probe_meta->contig = contig_id;
    tmp_probe_meta->name = name;
	return i >= 3? contig_id : 0;
}


int main() {
    char fname[] = "/home/wbguo/iproject/BSReadSim/temp/data/test_probe.bed";
    char *fn = fname;
    printf("%s", fn);
    char *contig_id;
    char *probe_name;
    

    htsFile *fp    = hts_open(fn,"rb");
    if(fp == 0 ){ fprintf(stderr,"cannot open bed file: %s\n",fn); exit (EXIT_FAILURE);}
    
    int ret;
    kstring_t line = {0,0,0};
    probe_rec tmp_probe; probe_meta tmp_probe_meta;
    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0)
    {
        tmp_probe = {};
        contig_id = parse_bed(line.s, &tmp_probe, &tmp_probe_meta);
        // printf("%s\n", line.s);
        // printf("%ld\n", line.l);
        // printf("%ld\n", line.m);
        printf("%s \n", contig_id);
        printf("%s \n", tmp_probe_meta.name);
        printf("%d %d %f %d\n", tmp_probe.pos_l, tmp_probe.pos_r, tmp_probe.score, tmp_probe.strand);
        probe_vec.push_back(tmp_probe);
    }
    free(line.s);


    for (size_t i = 0; i < probe_vec.size(); i++)
    {
        printf("%d %d %f %d\n", probe_vec[i].pos_l, probe_vec[i].pos_r, probe_vec[i].score, probe_vec[i].strand);
    }
    
    return 0;
}


// compile
// g++ parse_bed.cpp -o a -I ../../HTSLIB/htslib/ -L ../../HTSLIB/ -lhts -Wl,-rpath ../../HTSLIB/