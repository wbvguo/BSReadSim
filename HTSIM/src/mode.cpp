#include <algorithm>
#include <stdint.h>
#include <vector>
#include "mode.h"
#include "htsim2.h"
#include "kseq.h"
#include "vcf.h"



// for VCF
void parse_vcf_chr(char *fname, char *chr_id, std::vector<snp_rec>& snp_vec)
{
    snp_vec.clear(); //clean the container

    //open vcf file
    htsFile   *fp  = hts_open(fname,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t    *rec = bcf_init();
    
    //collect control
    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
    int prev_pos= -1;

    bool collect_present = false;
    bool collect_previous= false;

    //check the vcf file
    int nsmpl = bcf_hdr_nsamples(hdr); // number of sample
    if (nsmpl != 1) {
        fprintf(stderr, "[%s] ERROR: Currently only support single-sample simulation, please check the vcf file! Exiting...\n", __func__);
        exit(EXIT_FAILURE);
    }

    while (bcf_read(fp, hdr, rec)>=0) {	
        if (collect_present == false && collect_previous == true) { break; } // finished collecting
        
        // a new round, save last status
        collect_previous = collect_present; 

        //unpack for read REF,ALT,INFO,etc 
        bcf_unpack(rec, BCF_UN_STR);
        bcf_unpack(rec, BCF_UN_INFO);

        //read CHROM
        std::string chr_current = bcf_hdr_id2name(hdr, rec->rid);

        if (strcmp(chr_current.c_str(), chr_id) != 0){
            collect_present = false;
            continue;
        } else {
            collect_present = true;
            std::string ref = rec->d.allele[0];
            std::string alt = rec->d.allele[1];
            int snp_pos = rec->pos; //it's 0-based coordinates
            int base_change_pos = snp_pos;

            int ref_len = ref.length();
            int alt_len = alt.length();
            int base_offset = alt_len - ref_len;

            // check genotype
            ngt = bcf_get_genotypes(hdr,rec,&gt,&ngt_arr); //The total number of array elements in &gt
            int snp_hap1 = bcf_gt_allele(gt[0]);
            int snp_hap2 = bcf_gt_allele(gt[1]);
            int is_phased = bcf_gt_is_phased(gt[1]); //phased:1, unphased:0

            //skip the following records: 
            // 1. insert/delete offset greater than 4
            // 2. multi-allelic sites
            // 3. multi-ploid organism
            // 4. multi nucleotide polymorphism (MNP); (please represent it with single base SNPs)
            // 5. SNP with the same position (only keep the first occurence)
            // TODO: 
            // 1. contains Ns;
            // 2. missing snps?
            if (abs(base_offset) > 4 || ngt > 2 || prev_pos==snp_pos) {
                fprintf(stderr, "[%s] Skip unusual SNP/INDEL [length/ploidy/pos]: CHROM:%s; POS:%d; REF:%s; ALT:%s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
                continue;
            }
            if (snp_hap1 < 0 || snp_hap2 < 0 || snp_hap1 > 1 || snp_hap2 > 1) {
                fprintf(stderr, "[%s] Skip unusual SNP/INDEL [snp haplotype]: CHROM:%s; POS:%d; HAP1:%d; HAP2:%d\n", __func__, chr_id, snp_pos, snp_hap1, snp_hap2);
                continue;
            }
            prev_pos = snp_pos;

            //pack info: encode SNP info into int
            int ref_int, alt_int, c = 0;
            if (alt_len == 1 && ref_len == 1){ //substitution
                ref_int = (mut_t)nst_nt4_table[(int)ref[0]];
                alt_int = (mut_t)nst_nt4_table[(int)alt[0]];
            } else {
                if ( ref[0] != alt[0] ) { //check if the first base are the same if it's indel, if not skip
                    fprintf(stderr, "[%s] Skip unusual SNP/INDEL [indel & ref]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
                    continue;
                }

                if (alt_len > 1 && ref_len == 1) { //insertion
                    ref_int = (mut_t)nst_nt4_table[(int)ref[0]];
                    base_change_pos = snp_pos + 1; // position +1, because the base change occurs after the first base
                    for (int i = 0; i < alt_len; ++i ){ 
                        c = (mut_t)nst_nt4_table[(int)alt[i]];
                        alt_int = (alt_int << 2) | c;
                    }
                } else if (ref_len > 1 && alt_len == 1){ //deletion
                    alt_int = (mut_t)nst_nt4_table[(int)alt[0]];
                    base_change_pos = snp_pos + 1; // position +1, because the base change occurs after the first base
                    for (int i = 0; i < ref_len; ++i ){ 
                        c = (mut_t)nst_nt4_table[(int)ref[i]];
                        ref_int = (ref_int << 2) | c;
                    }
                } else { // might be MNP, or sth else
                    fprintf(stderr, "[%s] Skip unusual SNP/INDEL [MNP]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref.c_str(), alt.c_str());
                    continue;
                }
            }

            int geno_int = base_offset << 12 | ref_len << 8 | snp_hap2 << 6 | snp_hap1 << 4 | is_phased;
            snp_rec tmp_snp = {.pos = base_change_pos, .ref = ref_int, .alt = alt_int, .geno = geno_int};
            snp_vec.push_back(tmp_snp);
            tmp_snp = {};
        }
    }
    
    fprintf(stderr, "[%s] Finish collecting %lu SNP/INDEL from %s\n", __func__, snp_vec.size(), chr_id);

    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);

    if (int ret=hts_close(fp)){
        fprintf(stderr,"[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret); 
    }
}


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
    tmp_probe->strand= strand;
    
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


// for RRBS
void parse_cut_rec(char *cut_str, std::vector<cut_rec>& cut_vec)
{
    cut_vec.clear();
    int c, idx = 0;
    cut_rec tmp_cut;
    
    for(int i =0; cut_str[i] !='\0'; ++i){
        if (cut_str[i] == '|'){tmp_cut.idx = idx; continue;}
        if (cut_str[i] == ','){
            if (tmp_cut.idx >=0) {
                tmp_cut.len = idx;
                cut_vec.push_back(tmp_cut);
                tmp_cut = {}; idx = 0; continue; // empty the struct, restart
            } else {
                fprintf (stderr, "[%s] ERROR: Invalid enzyme site format (%s). Exit... \n", __func__, cut_str); exit (EXIT_FAILURE);
            }
        }
        ++idx;
        c = nst_nt4_table[(int)cut_str[i]];
        tmp_cut.seq.push_back(c);
    }
    if(tmp_cut.idx >=0){tmp_cut.len = idx; cut_vec.push_back(tmp_cut);}
}

void gen_cut_pos(const kseq_t *ks, std::vector<cut_pos>& cutpos_vec, std::vector<cut_rec>& cut_vec)
{
    cutpos_vec.clear(); // clean the container

    std::vector<mut_t> rseq_ref(ks->seq.l); // TODO: create cut_pos without create rseq_ref

    for (int i = 0; i != ks->seq.l; ++i) {
        rseq_ref.push_back((mut_t)nst_nt4_table[(int)ks->seq.s[i]]);
    }

    std::vector<mut_t>::iterator ptr_begin = rseq_ref.begin();  // can repalce the type with `auto`
    std::vector<mut_t>::iterator ptr_end   = rseq_ref.end();
    std::vector<mut_t>::iterator iter_curr = rseq_ref.begin();  // current
    std::vector<mut_t>::iterator iter_temp = rseq_ref.begin();  // temp
    std::vector<mut_t>::iterator iter_save = rseq_ref.begin();  // save
    //printf("%d %d\n", ptr_begin, ptr_end);

    int count = 0; 
    cut_pos tmp_cutpos;
    while (iter_curr < ptr_end)
    {
        // printf("========%i========\n", count);
        // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);
        //initiate
        int min_pos = ptr_end - ptr_begin;
        int site_pos= min_pos;
        int type_idx = -1;

        // printf("[%d %d %d]\n", min_pos, site_pos, type_idx);
        for (int j = 0; j < cut_vec.size(); ++j) {
            iter_temp = std::search(iter_curr, ptr_end, cut_vec[j].seq.begin(), cut_vec[j].seq.end());
            // printf("%d %d %d\n", iter_curr, iter_temp);
            site_pos = iter_temp - ptr_begin;
            if(site_pos < min_pos){
                min_pos = site_pos;
                type_idx= j;
                iter_save= iter_temp;
            }
        }
        // printf("[%d %d %d]\n", min_pos, site_pos, type_idx);
        // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);

        if(type_idx >=0){
            tmp_cutpos = {};
            tmp_cutpos.pos = min_pos + cut_vec[type_idx].idx;
            tmp_cutpos.type= type_idx;
            cutpos_vec.push_back(tmp_cutpos);
            iter_curr = iter_save + cut_vec[type_idx].len; // need to check for other int types
        }
        // printf("%d %d %d\n", iter_curr, iter_temp, iter_save);
        // ++count;
        // if(count >= 20){
        //     break;
        // }
    }
}

//gen cut fragment


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
void cal_length_chr(const kseq_t *ks, chr_rec *tmp_len, int tech_mode, char *bed_file, int min_insert, int max_insert, 
                    std::vector<probe_rec>& probe_vec, std::vector<cut_pos>& cutpos_vec, std::vector<cut_rec>& cut_vec)
{
    uint64_t eff_len;

    if(tech_mode==2){                       // targeted sequencing
        parse_bed_chr(bed_file, ks->name.s, probe_vec);
        int pos_l, pos_r, len;
        int pos_l_prev=0, pos_r_prev=0;
        for (size_t i = 0; i < probe_vec.size(); ++i){
            pos_l = probe_vec[i].pos_l;
            pos_r = probe_vec[i].pos_r;
            len   = pos_r - pos_r;
            // deal with overlap probes
            if ((pos_l < pos_r_prev) && (pos_l > pos_l_prev)){ 
                len = (pos_r > pos_r_prev) ? (pos_r - pos_r_prev) : 0;
            }
            eff_len += len;
            pos_l_prev = pos_l;
            pos_r_prev = pos_r;
        }
    }else if (tech_mode==1){                // reduced representative sequencing 
        gen_cut_pos(ks, cutpos_vec, cut_vec);
        int len;
        for (int i = -1 ; i <= (int)cutpos_vec.size(); ++i){
            if(i == -1 || i == (int)cutpos_vec.size()){
                if(i==-1){
                    len = cutpos_vec[0].pos - 0;
                }else{
                    len = (int) ks->seq.l - cutpos_vec[i-1].pos;
                }
            }else{
                len = cutpos_vec[i+1].pos - cutpos_vec[i].pos;
            }
            len = len >= min_insert && len <= max_insert? len : 0;
            eff_len += len;
        }
    }else{                                  // whole genome
        eff_len = ks->seq.l;
    }
    tmp_len->chr_len = ks->seq.l;
    tmp_len->eff_len = eff_len;
}
