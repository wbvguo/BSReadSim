#include <stdint.h>
#include <vector>
#include <string>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"

KSEQ_INIT(gzFile, gzread)
#include "haplo.h"

void parse_vcf_chr(char *fname, char *chr_id, std::map<int, snpmeth_rec>& snpmeth_map)
{
    snpmeth_map.clear(); //clean the container

    //open vcf file
    htsFile   *fp  = hts_open(fname,"rb");
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    bcf1_t    *rec = bcf_init();

    //store snp info
    snpmeth_rec tmp_snp;

    //collect control
    int ngt_arr = 0;
    int ngt     = 0;
    int *gt     = NULL;
    int prev_pos= -1;

    bool collect_present = false;
    bool collect_previous= false;

    //check the vcf file
    int n_samp = bcf_hdr_nsamples(hdr); // number of sample
    if (n_samp != 1) {
        fprintf(stderr, "[%s] ERROR: Currently only support single-sample simulation, please check the vcf file! Exiting...\n", __func__);
        exit(EXIT_FAILURE);
    }

    while (bcf_read(fp, hdr, rec)>=0) {	
        if (!collect_present && collect_previous) { break; } // collect_present false && collect_previous true, finished collecting
        
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
            std::string ref_str = rec->d.allele[0];
            std::string alt_str = rec->d.allele[1];
            int snp_pos = rec->pos; //it's 0-based coordinates
            int ref_len = ref_str.length();
            int alt_len = alt_str.length();
            int base_offset = alt_len - ref_len;
            int base_change_pos = snp_pos;

            // check genotype
            ngt = bcf_get_genotypes(hdr,rec,&gt,&ngt_arr); //The total number of array elements in &gt
            int snp_hap1 = bcf_gt_allele(gt[0]);
            int snp_hap2 = bcf_gt_allele(gt[1]);
            int is_phased= bcf_gt_is_phased(gt[1]); //phased:1, unphased:0

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
                fprintf(stderr, "[%s] Skip unusual SNP/INDEL [length/ploidy/pos]: CHROM:%s; POS:%d; REF:%s; ALT:%s\n", __func__, chr_id, snp_pos, ref_str.c_str(), alt_str.c_str());
                continue;
            }
            if (snp_hap1 < 0 || snp_hap2 < 0 || snp_hap1 > 1 || snp_hap2 > 1) {
                fprintf(stderr, "[%s] Skip unusual SNP/INDEL [snp haplotype]: CHROM:%s; POS:%d; HAP1:%d; HAP2:%d\n", __func__, chr_id, snp_pos, snp_hap1, snp_hap2);
                continue;
            }
            prev_pos = snp_pos;

            //pack info: encode SNP info into int
            int ref_int, alt_int, c = 0;
            if (alt_len == 1 && ref_len == 1){          //substitution
                ref_int = (mut_t)nst_nt4_table[(int)ref_str[0]];
                alt_int = (mut_t)nst_nt4_table[(int)alt_str[0]];
            } else {
                if (ref_str[0] != alt_str[0]) {         //check if the first base are the same: yes -> indel, no -> skip
                    fprintf(stderr, "[%s] Skip unusual SNP/INDEL [indel & ref]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref_str.c_str(), alt_str.c_str());
                    continue;
                }

                if (alt_len > 1 && ref_len == 1) {      //insertion
                    ref_int = (mut_t)nst_nt4_table[(int)ref_str[0]];
                    base_change_pos = snp_pos + 1;                  // position +1, because the base change occurs after the first base
                    for (int i = alt_len-1; i >0; --i){             // index 0 has the same base as ref, the lower bit, the closer to ref
                        c = (mut_t)nst_nt4_table[(int)alt_str[i]];
                        alt_int = (alt_int << 2) | c;
                    }
                } else if (ref_len > 1 && alt_len == 1){ //deletion
                    alt_int = (mut_t)nst_nt4_table[(int)alt_str[0]];
                    base_change_pos = snp_pos + 1;                  // position +1, because the base change occurs after the first base
                    for (int i = ref_len-1; i >0; --i ){            // no use later
                        c = (mut_t)nst_nt4_table[(int)ref_str[i]];
                        ref_int = (ref_int << 2) | c;
                    }
                } else { // might be MNP, or sth else
                    fprintf(stderr, "[%s] Skip unusual SNP/INDEL [MNP]: CHROM:%s; POS:%d; REF:%s; ALT %s\n", __func__, chr_id, snp_pos, ref_str.c_str(), alt_str.c_str());
                    continue;
                }
            }

            tmp_snp = {.ref = (uint16_t) ref_int, .alt = (uint16_t) alt_int, 
                       .hap1= (int8_t) snp_hap1, .hap2 = (int8_t) snp_hap2, .is_phased = (int8_t) is_phased, .offset = (int8_t) base_offset};
            snpmeth_map[base_change_pos]= tmp_snp;
            tmp_snp = {};
        }
    }
    
    fprintf(stderr, "[%s] Finish collecting %lu SNP/INDEL from %s\n", __func__, snpmeth_map.size(), chr_id);

    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);

    if (int ret=hts_close(fp)){
        fprintf(stderr, "[%s] ERROR: hts_close(%s): non-zero status %d\n", __func__, fname, ret); exit(ret); 
    }
}

void sim_mut_vcf(const kseq_t *ks, char * vcf_file, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr, std::map<int, snpmeth_rec>& snpmeth_map)
{
    // initiate
    mutseq_t *ret[2];

    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    
    // parse VCF
    parse_vcf_chr(vcf_file, ks->name.s, snpmeth_map);

    std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin();
    int i, c, tmp_hap;
    int deleting = 0, deletion_count = 0;

    for (i = 0; i != (int) ks->seq.l; ++i){
        c = ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ks->seq.s[i]];
        if (cg_table[c]) {posidx_arr[i] = 2;}
        if (snpmeth_map.size() == 0){ continue; } // ignore the rest if there is no SNP

        if (deleting){
            if(deletion_count > 0){
                if (deleting & 1){ ret[0]->s[i] |= DELETE;}
                if (deleting & 2){ ret[1]->s[i] |= DELETE;}
                --deletion_count;
                posidx_arr[i]|= 1;
                continue;
            } else { deleting = 0;}
        }
        if(i == it->first && c < 4 && it != snpmeth_map.end()){
            int alt_int  = it->second.alt;
            int snp_hap1 = it->second.hap1;
            int snp_hap2 = it->second.hap2;
            int is_phased= it->second.is_phased;
            int base_offset= it->second.offset;
            posidx_arr[i] |= 1;
            //fprintf(stderr, "%d,%d,%d,%d,%d\n", i, is_pahsed,snp_hap1, snp_hap2, base_offset);

            if (!is_phased){ // for unphased genotype, randomly swap the haplotype
                if(drand48() < 0.5){int tmp_hap = snp_hap1; snp_hap1 = snp_hap2; snp_hap2 = tmp_hap;}
            }

            if(base_offset == 0){           // SNP substitution
                if (snp_hap1 == 1){ret[0]->s[i] = SUBSTITUTE|alt_int;}
                if (snp_hap2 == 1){ret[1]->s[i] = SUBSTITUTE|alt_int;}
            } else if (base_offset < 0 ) {  // deletion
                deletion_count = abs(base_offset) - 1; //minus one because here already delete one base here
                if (snp_hap1 == 1){ret[0]->s[i]|= DELETE; deleting+=1;}
                if (snp_hap2 == 1){ret[1]->s[i]|= DELETE; deleting+=2;}
            } else if (base_offset > 0){    // inserstion
                int num_ins = base_offset;
                int ins_msk = (1 << (num_ins*2)) - 1;
                int ins = alt_int & ins_msk;
                //fprintf(stderr, "%d,%d,%d,%d\n", num_ins, ins_msk, alt_vec[vec_ptr], ins);
                if (snp_hap1 == 1){ret[0]->s[i] |= (num_ins << 12) | (ins << 4);}
                if (snp_hap2 == 1){ret[1]->s[i] |= (num_ins << 12) | (ins << 4);}
            }
            ++it;
        }
    }
}

void sim_mut_diref(const kseq_t *ks, mut_param *mut_set, mutseq_t *hap1, mutseq_t *hap2, uint32_t *posidx_arr)
{
    int i, c, deleting = 0;
    mutseq_t *ret[2];
    //drand48() is 8 times faster than uniform_distribution, use it exclusively for snp generation

    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ks->seq.l; ret[1]->l = ks->seq.l;
    ret[0]->m = ks->seq.m; ret[1]->m = ks->seq.m;
    ret[0]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ks->seq.m, sizeof(mut_t));
    for (i = 0; i != (int)ks->seq.l; ++i) {
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
                c = (c + (int)(drand48() * 3.0 + 1)) & 3;   // random mutation
                if (mut_set->is_hap || drand48() < 0.333333) {      // hom-sub
                    ret[0]->s[i] = ret[1]->s[i] = SUBSTITUTE|c;
                } else { // het
                    ret[drand48()<0.5?0:1]->s[i]= SUBSTITUTE|c;
                }
            } else { // indel
                if (drand48() < 0.5) { // deletion
                    if (mut_set->is_hap || drand48() < 0.333333) {  // hom-del
                        ret[0]->s[i] = ret[1]->s[i] |= DELETE;
                        deleting = 3;
                    } else { // het-del
                        deleting = drand48()<0.5?1:2;
                        ret[deleting-1]->s[i] |= DELETE;
                    }
                } else { // insertion
                    int num_ins = 0, ins = 0;
                    do {
                        ++num_ins;
                        ins = (ins << 2) | (int)(drand48() * 4.0);
                    } while (num_ins < 4 && drand48() < mut_set->indel_extn);

                    if (mut_set->is_hap || drand48() < 0.333333) { // hom-ins
                        ret[0]->s[i] = ret[1]->s[i] |= (num_ins << 12) | (ins << 4);
                    } else { // het-ins
                        ret[drand48()<0.5?0:1]->s[i]|= (num_ins << 12) | (ins << 4);
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
    int hap1_mut, hap2_mut;
    for (i = 0; i != (int)ks->seq.l; ++i) {
        int c[3];
        c[0] = nst_nt4_table[(int)ks->seq.s[i]];
        if (c[0] >= 4) continue;
        c[1] = hap1->s[i]; hap1_mut = (c[1] & mutmsk);
        c[2] = hap2->s[i]; hap2_mut = (c[2] & mutmsk);
        if (hap1_mut != NOCHANGE || hap2_mut != NOCHANGE) {
            if (c[1] == c[2]) { // hom
                if (hap1_mut == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t-\n", name, i+fmt_offset, "ACGTN"[c[0]], "ACGTN"[c[1]&0xf]); // coordinate is 1-based
                } else if (hap1_mut == DELETE) { // del
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+fmt_offset);
                        for (j = i; j < ks->seq.l && hap1->s[j] == hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if ((hap1_mut >> 12) <= 4) { // ins
                    printf("%s\t%d\t-\t", name, i+fmt_offset);
                    int n = hap1_mut >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        --n;
                    }
                    printf("\t-\n");
                } // else: deleted base in a long deletion
            } else { // het
                if (hap1_mut == SUBSTITUTE || hap2_mut == SUBSTITUTE) { // substitution
                    printf("%s\t%d\t%c\t%c\t+\n", name, i+fmt_offset, "ACGTN"[c[0]], "XACMGRSVTWYHKDBN"[1<<(c[1]&0x3)|1<<(c[2]&0x3)]);
                } else if (hap1_mut == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+fmt_offset);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap1->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if (hap2_mut == DELETE) {
                    if (i >= j) {
                        printf("%s\t%d\t", name, i+fmt_offset);
                        for (j = i; j < ks->seq.l && hap1->s[j] != hap2->s[j] && (hap2->s[j]&mutmsk) == DELETE; ++j)
                            putchar("ACGTN"[nst_nt4_table[(int)ks->seq.s[j]]]);
                        printf("\t-\t-\n");
                    }
                } else if ((hap1_mut >> 12) <= 4 && (hap1_mut >> 12) > 0) { // ins1
                    printf("%s\t%d\t-\t", name, i+fmt_offset);
                    int n = hap1_mut >> 12, ins = c[1] >> 4;
                    while (n > 0) {
                        putchar("ACGTN"[ins & 0x3]);
                        ins >>= 2;
                        --n;
                    }
                    printf("\t+\n");
                } else if ((hap2_mut >> 12) <= 4 || (hap2_mut >> 12) > 0) { // ins2
                    printf("%s\t%d\t-\t", name, i+fmt_offset);
                    int n = hap2_mut >> 12, ins = c[2] >> 4;
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

