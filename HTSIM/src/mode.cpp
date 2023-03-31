#include <stdint.h>
#include <vector>
#include <map>
#include <string>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"


KSEQ_INIT(gzFile, gzread)
#include "mode.h"


// for BED
int parse_bed_line(char *line, char *chr_id, frag_rec *tmp_probe)
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
            case 4: score  = atof(p); break;        // will give 0 when not a number
            case 5: strand = ((int)(strcmp(p,"+")!=0)<<1) | (int)(strcmp(p,"-")!=0); break; //1,2,3
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return 1;} // termenate early if not equal
			++i, p = q + 1;
			if (i > 6 || c == '\0') break;
		}
	}

    if(i < 4){fprintf(stderr, "[%s] Skip invalid probe: chr %s, name %s...\n", __func__, chr_id, name); return 0;}
    if(end <= start){fprintf(stderr, "[%s] Skip invalid probe (coordinate conflict): chr %s, name %s...\n", __func__, chr_id, name); return 0;}

	tmp_probe->pos_l = start;
	tmp_probe->pos_r = end;
    tmp_probe->score = score;
    tmp_probe->strand= (int8_t) strand;

    return 0;
}

int parse_bed_line_rrbs(char *line, char *chr_id, frag_rec *tmp_probe)
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
            case 4: score  = atof(p); break;
            case 5: strand = ((int)(strcmp(p,"+")!=0)<<1) | (int)(strcmp(p,"-")!=0); break; //1,2,3
            case 7: cut_l  = atoi(p); break;
            case 8: cut_r  = atoi(p); break;
            default: break;}
            if(i==0 && strcmp(contig, chr_id)){return 1;} // termenate early if not equal
			++i, p = q + 1;
			if (i > 9 || c == '\0') break;
		}
	}

    if(i < 9){fprintf(stderr, "[%s] Skip invalid probe: chr %s, name %s...\n", __func__, chr_id, name); return 0;}
    if(end <= start){fprintf(stderr, "[%s] Skip invalid probe (coordinate conflict): chr %s, name %s...\n", __func__, chr_id, name); return 0;}

	tmp_probe->pos_l = start;
	tmp_probe->pos_r = end;
    tmp_probe->score = score;
    tmp_probe->haplo = (int8_t) strand;
    tmp_probe->cut_l = (int8_t) cut_l;
    tmp_probe->cut_r = (int8_t) cut_r;

    return 0;
}

void parse_bed_chr(char *fname, char *chr_id, std::vector<frag_rec>& probe_vec, int tech_mode)
{
    probe_vec.clear();  // clean the container

    htsFile *fp = hts_open(fname,"rb");
    if(fp == 0 ){ fprintf(stderr,"[%s] ERROR: open bed file: %s failed: %s. Exit...\n", __func__, fname, strerror (errno)); exit (EXIT_FAILURE);}

    kstring_t line = {0,0,0};
    frag_rec tmp_probe;

    int ret;
    char *chr_current;
    bool collect_present = false;
    bool collect_previous= false;
    int (*parse_bed_ptr)(char*, char*, frag_rec*);
    if(tech_mode==1){parse_bed_ptr = &parse_bed_line_rrbs;}else{parse_bed_ptr = &parse_bed_line;}

    while ((ret = hts_getline(fp, KS_SEP_LINE, &line)) >= 0) {
        if (collect_present == false && collect_previous == true) { break; } // finished collecting
        // a new round, save last status
        collect_previous= collect_present;
        int skip_status = (*parse_bed_ptr)(line.s, chr_id, &tmp_probe); //might need to test
        if(skip_status){collect_present = false; continue;}
        probe_vec.push_back(tmp_probe);
        tmp_probe      = {};
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
void collect_len_score_chr(const kseq_t *ks, chr_rec *tmp_len, char *bed_file, int tech_mode, std::vector<frag_rec>& probe_vec)
{
    uint32_t eff_len=0;
    float sum_score=0;
    bool bool_bed_set = strcmp(bed_file,"None") && strlen(bed_file);

    if(bool_bed_set){                       // targeted sequencing or
        parse_bed_chr(bed_file, ks->name.s, probe_vec, tech_mode);
        int len, pos_l, pos_r, pos_l_prev, pos_r_prev;
        float score = 0;
        for (size_t i = 0; i < probe_vec.size(); ++i){
            pos_l = probe_vec[i].pos_l;
            pos_r = probe_vec[i].pos_r;
            len   = pos_r - pos_l;
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


// for count calculation
void cal_chr_count(const char *fn, char *chr_id, char *bed_file, uint64_t N, uint64_t chr_N, 
                    expt_param *expt_set, std::vector<frag_rec>& probe_vec, std::map<std::string, chr_rec> &chr_count)
{
    kseq_t *ks;
    gzFile fp_fa;
    fp_fa = gzopen(fn, "r");
    ks = kseq_init(fp_fa);

    if (!fp_fa) { fprintf (stderr, "ERROR: gzopen of '%s' failed: %s. Exit... \n", fn, strerror (errno)); exit (EXIT_FAILURE);}
    fprintf(stderr, "[%s] Calculating the length and count of the reference sequences...\n", __func__);
    
    chr_rec tmp_len;
    uint64_t tot_chr_len = 0, tot_eff_len = 0;
    float tot_score   = 0;

    int l;
    while ((l = kseq_read(ks)) >= 0) {
        if (expt_set->is_chr_set) {if (strcmp(chr_id, ks->name.s)!=0){continue;}}
        if (l < expt_set->mean_insert+3*expt_set->sd_insert){
            fprintf(stderr, "[%s] skip contig '%s' as it is shorter than %d!\n", __func__, ks->name.s, expt_set->mean_insert+3*expt_set->sd_insert); 
            continue;
        }
        
        tmp_len = {};
        collect_len_score_chr(ks, &tmp_len, bed_file, expt_set->tech_mode, probe_vec);
        chr_count[std::string(ks->name.s)] = tmp_len;
        tot_chr_len += tmp_len.chr_len;
        tot_eff_len += tmp_len.eff_len;
        tot_score   += tmp_len.score;
    }
    kseq_destroy(ks);
    gzclose(fp_fa);


    // check if fasta is empty
    if (chr_count.empty()) {fprintf (stderr, "ERROR: Input fasta is empty: %s. Exit... \n", fn); exit (EXIT_FAILURE);}


    // check input chr_id, calculate the count for contig(s)
    if (expt_set->is_chr_set){
        // calculate the count for selected contig
        std::string chr_id_str = std::string(chr_id);
        if (!chr_count.count(chr_id_str)){fprintf(stderr, "ERROR: Contig id '%s' is not found in the fasta file, please check!\n", chr_id); exit(EXIT_FAILURE);}

        uint32_t contig_eff_len = chr_count[chr_id_str].eff_len;
        uint32_t contig_len = chr_count[chr_id_str].chr_len;
        chr_N = chr_N == 0? (contig_eff_len * expt_set->depth)/(expt_set->size_l + expt_set->size_r) : chr_N;
        fprintf(stderr, "[%s] Contig %s specified, total length: %u, effective length: %u, #reads: %lu\n", __func__, chr_id, contig_len, contig_eff_len, chr_N);
        chr_count[chr_id_str].count = chr_N;
    } else {
        // calculate the count for all contigs
        if (chr_N > 0) {fprintf(stderr, "ERROR: -n is specified but not -c. Exit... (please note the difference of -n and -N)\n"); exit(EXIT_FAILURE);}
        
        int num_contigs = (int)chr_count.size();
        N = N == 0? (tot_eff_len * expt_set->depth)/(expt_set->size_l + expt_set->size_r) : N;
        fprintf(stderr, "[%s] Found %d contig sequences, total length: %lu, effective length: %lu\n", __func__, num_contigs, tot_chr_len, tot_eff_len);
        fprintf(stderr, "[%s] No contig id specified, will generate %lu reads from all contigs\n", __func__, N);
        
        uint64_t cum_count = 0, tmp_count =0;
        for (auto it = chr_count.begin(); it != chr_count.end(); ++it) {
            tmp_count = expt_set->tech_mode ==2 ? (uint64_t)(it->second.score * N / tot_score) : (uint64_t)(it->second.eff_len * N / tot_eff_len);
            it->second.count = tmp_count;
            cum_count += tmp_count;
        }

        int rest_count = N - cum_count;
        if(rest_count < 0){fprintf(stderr, "[%s] Read count calculation went wrong \n", __func__); exit(EXIT_FAILURE);} // should never happen
        int step_size  = rest_count > num_contigs ? (int)(rest_count/num_contigs): 1; //hopefully rest_count is small, evenly distributed to contigs
        auto it = chr_count.begin();
        while (rest_count > 0){ 
            int alloc_count  = std::min(step_size, rest_count);
            it->second.count+= alloc_count; 
            rest_count -= alloc_count;
            ++it;
        }
    }
}

