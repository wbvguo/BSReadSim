#include <stdint.h>
#include <vector>
#include <string>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"

KSEQ_INIT(gzFile, gzread)
#include "haplo.h"

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

void sim_mut_vcf(const kseq_t *ks, char * vcf_file, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, std::vector<snp_rec> &snp_vec) 
{
    // initiate
    mutseq_t *ret[2];

    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    
    // parse VCF
    parse_vcf_chr(vcf_file, ks->name.s, snp_vec);

    int vec_ptr = 0;
    int i, deleting = 0;
    int deletion_count = 0;
    int c;

    for (i = 0; i != (int) ks->seq.l; ++i){
        c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (cg_table[c]) {posidx_arr[i]= 2;}
        if (snp_vec.size() == 0){ continue;} // ignore the rest if there is no SNP

        if (deleting){
            if(deletion_count > 0){
                if (deleting & 1){ ret[0]->s[i] |= DELETE;}
                if (deleting & 2){ ret[1]->s[i] |= DELETE;}
                --deletion_count;
                posidx_arr[i] |= 1;
                continue;
            } else {deleting = 0;}
        }

        if(vec_ptr < (int) snp_vec.size() && i == snp_vec[vec_ptr].pos && c < 4){
            int geno_int = snp_vec[vec_ptr].geno;
            int is_phased= geno_int & 0x000f;
            int snp_hap1 = (geno_int & 0x003f) >> 4;
            int snp_hap2 = (geno_int & 0x00ff) >> 6;
            int ref_len  = (geno_int & 0x0fff) >> 8;
            int base_offset = geno_int >> 12;
            posidx_arr[i] |= 1;
            //fprintf(stderr, "%d,%d,%d,%d,%d,%d\n", i, is_pahsed,snp_hap1, snp_hap2, ref_len, base_offset);

            if (!is_phased){ // for unphased genotype, randomly swap the haplotype
                if(drand48() < 0.5){
                    int tmp_hap = snp_hap1;
                    snp_hap1 = snp_hap2;
                    snp_hap2 = tmp_hap;
                }
            }

            if(base_offset == 0 && ref_len == 1){ // SNP substitution
                c = snp_vec[vec_ptr].alt;

                if (snp_hap1 == 1 && snp_hap2 == 1){
                    ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                } else if (snp_hap1 == 1 && snp_hap2 == 0){
                    ret[0]->s[i] = SUBSTITUTE|c;
                } else if (snp_hap1 == 0 && snp_hap2 == 1){
                    ret[1]->s[i] = SUBSTITUTE|c;
                } else{continue;}
            } else if (base_offset < 0 ) { // deletion
                c = snp_vec[vec_ptr].ref;
                deletion_count = abs(base_offset) - 1; //minus one because here already delete one base

                if (snp_hap1 == 1 && snp_hap2 == 1){
                    ret[0]->s[i] = ret[1]->s[i] =  DELETE | c;
                    deleting = 3;
                } else if (snp_hap1 == 1 && snp_hap2 == 0){
                    ret[0]->s[i] =  DELETE | c;
                    deleting = 1;
                } else if (snp_hap1 == 0 && snp_hap2 == 1){
                    ret[1]->s[i] =  DELETE | c;
                    deleting = 2;
                } else{continue;}
            } else if (base_offset > 0){ // inserstion
                int num_ins = base_offset;
                int ins_msk = (1 << (num_ins*2)) - 1;
                int ins = snp_vec[vec_ptr].alt & ins_msk;
                //fprintf(stderr, "%d,%d,%d,%d\n", num_ins, ins_msk, alt_vec[vec_ptr], ins);

                if (snp_hap1 == 1 && snp_hap2 == 1){
                    ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else if (snp_hap1 == 1 && snp_hap2 == 0){
                    ret[0]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else if (snp_hap1 == 0 && snp_hap2 == 1){
                    ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                } else{continue;}
            }
            ++vec_ptr;
        }
    }
}

void sim_mut_diref(const kseq_t *ks, mut_param *mut_set, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr)
{
    int i, deleting = 0;
    mutseq_t *ret[2];
    //drand48() is 8 times faster than uniform_distribution, use it exclusively for snp generation

    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    for (i = 0; i != (int)ks->seq.l; ++i) {
        int c;
        c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (cg_table[c]) {posidx_arr[i]= 2;}
        if (deleting) {
            if (drand48() < mut_set->indel_extn) {
                if (deleting & 1) ret[0]->s[i] |= DELETE;
                if (deleting & 2) ret[1]->s[i] |= DELETE;
                posidx_arr[i] |= 1;
                continue;
            } else deleting = 0;
        }
        if (c < 4 && drand48() < mut_set->mut_rate) { // mutation
            if (drand48() >= mut_set->indel_frac) { // substitution
                double r = drand48();
                c = (c + (int)(r * 3.0 + 1)) & 3;
                if (mut_set->is_hap || drand48() < 0.333333) { // hom
                    ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                } else { // het
                    ret[drand48()<0.5?0:1]->s[i] = SUBSTITUTE|c;
                }
            } else { // indel
                if (drand48() < 0.5) { // deletion
                    if (mut_set->is_hap || drand48() < 0.333333) { // hom-del
                        ret[0]->s[i] = ret[1]->s[i] = DELETE | c;
                        deleting = 3;
                    } else { // het-del
                        deleting = drand48()<0.5?1:2;
                        ret[deleting-1]->s[i] = DELETE | c;
                    }
                } else { // insertion
                    int num_ins = 0, ins = 0;
                    do {
                        ++num_ins;
                        ins = (ins << 2) | (int)(drand48() * 4.0);
                    } while (num_ins < 4 && drand48() < mut_set->indel_extn);

                    if (mut_set->is_hap || drand48() < 0.333333) { // hom-ins
                        ret[0]->s[i] = ret[1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                    } else { // het-ins
                        ret[drand48()<0.5?0:1]->s[i] = (num_ins << 12) | (ins << 4) | c;
                    }
                }
            }
            posidx_arr[i] |= 1;
        }
    }
}

void sim_print_mutref(const char *name, const kseq_t *ks, mutseq_t *hap1, mutseq_t *hap2, int output_fmt)
{
    int fmt_offset = output_fmt == 0;
    int i, j = 0; // j keeps the end of the last deletion
    for (i = 0; i != (int)ks->seq.l; ++i) {
        int c[3];
        c[0] = nst_nt4_table[(int)ks->seq.s[i]];
        c[1] = hap1->s[i]; c[2] = hap2->s[i];
        if (c[0] >= 4) continue;
        if ((c[1] & mutmsk) != NOCHANGE || (c[2] & mutmsk) != NOCHANGE) {
            if (c[1] == c[2]) { // hom
                if ((c[1]&mutmsk) == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t-\n", name, i+fmt_offset, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]); // coordinate is 1-based
                } else if ((c[1]&mutmsk) == DELETE) { // del
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+fmt_offset);
                        for (j = i; j < ks->seq.l && hap1->s[j] == hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (((c[1] & mutmsk) >> 12) <= 4) { // ins
                    printf("%s\t%d\t-\t", name, i+fmt_offset);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        --n;
                    }
                    printf("\t-\n");
                } // else: deleted base in a long deletion
            } else { // het
                if ((c[1]&mutmsk) == SUBSTITUTE || (c[2]&mutmsk) == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t+\n", name, i+fmt_offset, "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
                } else if ((c[1]&mutmsk) == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+fmt_offset);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if ((c[2]&mutmsk) == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+fmt_offset);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap2->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (((c[1] & mutmsk) >> 12) <= 4 && ((c[1] & mutmsk) >> 12) > 0) { // ins1
                    printf("%s\t%d\t-\t", name, i+fmt_offset);
                    int n = (c[1]&mutmsk) >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        --n;
                    }
                    printf("\t+\n");
                } else if (((c[2] & mutmsk) >> 12) <= 4 || ((c[2] & mutmsk) >> 12) > 0) { // ins2
                    printf("%s\t%d\t-\t", name, i+fmt_offset);
                    int n = (c[2]&mutmsk) >> 12, ins = c[2] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        --n;
                    }
                    printf("\t+\n");
                } // else: deleted base in a long deletion
            }
        }
    }
}

