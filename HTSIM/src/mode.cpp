#include <stdint.h>
#include <vector>
#include <string>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"


KSEQ_INIT(gzFile, gzread)
#include "mode.h"


// for BED
int parse_bed_line(char *line, char *chr_id, probe_rec *tmp_probe, probe_meta *tmp_probe_meta)
{
	char *p, *q, *contig, *name = 0;
    int i, start, end, strand;
	float score;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: contig = p; break;
            case 1: start  = atoi(p); break;
            case 2: end    = atoi(p); break;
            case 3: name   = strdup(p); break;
            case 4: score  = atof(p); break;        // TODO: need to test what if score is "." 
            case 5: strand = int(strcmp(p,"+")==0)-int(strcmp(p,"-")==0); break; //1,0,-1
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return 1;} // termenate early if not equal
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

    if(i < 4){fprintf(stderr, "[%s] Skip invalid probe: chr %s, name %s...\n", __func__, chr_id, name); return 0;}

	tmp_probe->pos_l = start;
	tmp_probe->pos_r = end;
    tmp_probe->score = score;
    tmp_probe->strand= (int8_t) strand;
    
    tmp_probe_meta->chr_id = chr_id;
    tmp_probe_meta->name   = name;
    return 0;
}

int parse_bed_line_rrbs(char *line, char *chr_id, probe_rec *tmp_probe, probe_meta *tmp_probe_meta)
{
	char *p, *q, *contig, *name = 0;
    int i, start, end, strand, cut_l, cut_r;
	float score;

	for (i = 0, p = q = line;; ++q) {
		if (*q == '\t' || *q == '\0') {
			int c = *q;
			*q = 0;
            switch (i) {
            case 0: contig = p; break;
            case 1: start  = atoi(p); break;
            case 2: end    = atoi(p); break;
            case 3: name   = strdup(p); break;
            case 4: score  = atof(p); break;        // TODO: need to test what if score is "." 
            case 5: strand = int(strcmp(p,"+")==0)-int(strcmp(p,"-")==0); break; //1,0,-1
            case 7: cut_l  = atoi(p); break;
            case 8: cut_r  = atoi(p); break;
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return 1;} // termenate early if not equal
			++i, p = q + 1;
			if (i > 9 || c == '\0') break;
		}
	}

    if(i < 9){fprintf(stderr, "[%s] Skip invalid probe: chr %s, name %s...\n", __func__, chr_id, name); return 0;}

	tmp_probe->pos_l = start;
	tmp_probe->pos_r = end;
    tmp_probe->score = score;
    tmp_probe->strand= (int8_t) strand;
    tmp_probe->cut_l= (int8_t) cut_l;
    tmp_probe->cut_r= (int8_t) cut_r;

    tmp_probe_meta->chr_id = chr_id;
    tmp_probe_meta->name   = name;
    return 0;
}

void parse_bed_chr(char *fname, char *chr_id, std::vector<probe_rec>& probe_vec)
{
    probe_vec.clear();  // clean the container

    htsFile *fp = hts_open(fname,"rb");
    if(fp == 0 ){ fprintf(stderr,"[%s] ERROR: open bed file: %s failed: %s. Exit...\n", __func__, fname, strerror (errno)); exit (EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    probe_rec tmp_probe;
    probe_meta tmp_probe_meta;

    int ret;
    char *chr_current;
    bool collect_present = false;
    bool collect_previous= false;

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (collect_present == false && collect_previous == true) { break; } // finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        int skip_status = parse_bed_line(line.s, chr_id, &tmp_probe, &tmp_probe_meta); //might need to test
        if(skip_status){collect_present = false; continue;}
        probe_vec.push_back(tmp_probe);
        tmp_probe      = {};
        tmp_probe_meta = {};
    }
    free(line.s);

    if (ret=hts_close(fp)){fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret);}
}


// for GC-bias
void parse_bias_file(char *fname, std::vector<float>& eff_vec)
{
    FILE* fp = fopen(fname, "r");
    float eff_prob;

    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open capture efficiency file: %s failed. Exit...", __func__, fname); exit (EXIT_FAILURE);}
    while(fscanf(fp, "%f", &eff_prob) == 1) { // this will skip the empty lines
        // printf("%f\n", eff_prob);
        eff_vec.push_back(eff_prob);
    }
    fclose(fp);
}


// for length calculation
void collect_len_score_chr(const kseq_t *ks, chr_rec *tmp_len, char *bed_file, std::vector<probe_rec>& probe_vec)
{
    uint64_t eff_len;
    float sum_score;
    bool bool_bed_set = strcmp(bed_file,"None") && strlen(bed_file);

    if(bool_bed_set){                       // targeted sequencing or
        parse_bed_chr(bed_file, ks->name.s, probe_vec);
        int len, pos_l, pos_r, pos_l_prev, pos_r_prev;
        float score;
        for (size_t i = 0; i < probe_vec.size(); ++i){
            pos_l = probe_vec[i].pos_l;
            pos_r = probe_vec[i].pos_r;
            len   = pos_r - pos_r;
            score = probe_vec[i].score;
            // deal with overlap regions
            if ((pos_l < pos_r_prev) && (pos_l > pos_l_prev)){ 
                len = (pos_r > pos_r_prev) ? (pos_r - pos_r_prev) : 0;
            }
            eff_len += len;
            sum_score += score;
            pos_l_prev = pos_l;
            pos_r_prev = pos_r;
        }
    }else{                                  // whole genome
        eff_len = ks->seq.l;
        sum_score = eff_len;
    }
    tmp_len->chr_len = ks->seq.l;
    tmp_len->eff_len = eff_len;
    tmp_len->score   = sum_score;
}

