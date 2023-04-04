#include <vector>
#include <map>
#include <algorithm>
#include <zlib.h>
#include "kseq.h"
#include "vcf.h"
#include "struct.h"


KSEQ_INIT(gzFile, gzread)
#include "rrcut.h"
#include "haplo.h"


std::vector<cut_rec> cut_vec;
std::vector<snp_rec> snp_vec;
std::vector<cutpos_rec> cutpos_vec;
std::vector<frag_rrbs_rec> frag_rrbs_vec;


void parse_cut_site(char *cut_str, std::vector<cut_rec>& cut_vec)
{
    // check if cut_str is valid (seprate by ',' and '_' intercjangeably, 
    // never __ without comma in between nor ,, without _ in between)
    bool vline_flag = false;
    bool comma_flag = false;

    for(int i =0; cut_str[i] !='\0'; ++i){
        if(cut_str[i] == ' '){continue;}
        if(cut_str[i] == '_'){
            if(vline_flag){fprintf (stderr, "[%s] ERROR: Invalid enzyme site format '%s'. Exit...\n", __func__, cut_str); exit (EXIT_FAILURE);}
            vline_flag = true; comma_flag = false;
        }else if(cut_str[i] == ','){
            if(comma_flag){fprintf (stderr, "[%s] ERROR: Invalid enzyme site format '%s'. Exit...\n", __func__, cut_str); exit (EXIT_FAILURE);}
            comma_flag = true; vline_flag = false;
        }else{
            if(nst_nt4_table[(int)cut_str[i]]==4 && (cut_str[i] != 'N' || cut_str[i] != 'n')){
                fprintf (stderr, "[%s] ERROR: Invalid enzyme site format '%s'. Exit...\n", __func__, cut_str); exit (EXIT_FAILURE);
            }
        }
    }


    cut_vec.clear();
    int c, idx = 0;
    cut_rec tmp_site;
    
    for(int i =0; cut_str[i] !='\0'; ++i){
        if(cut_str[i] == ' '){continue;}
        if (cut_str[i] == '_'){tmp_site.idx = idx; continue;}
        if (cut_str[i] == ','){
            if (tmp_site.idx >=0) {
                tmp_site.len = idx;
                cut_vec.push_back(tmp_site);
                tmp_site = {}; idx = 0; continue; // empty the struct, restart
            } else {
                fprintf (stderr, "[%s] ERROR: Invalid enzyme site format (%s). Exit...\n", __func__, cut_str); exit (EXIT_FAILURE);
            }
        }
        ++idx;
        c = nst_nt4_table[(int)cut_str[i]];
        tmp_site.seq.push_back(c);
    }
    if(tmp_site.idx >=0){tmp_site.len = idx; cut_vec.push_back(tmp_site);}
}

void gen_cut_pos(mutseq_t *hap1, mutseq_t *hap2, std::vector<cutpos_rec>& cutpos_vec, std::vector<cut_rec>& cut_vec)
{
    cutpos_vec.clear(); // clean the container
    
    std::map<cutpos_rec, int> cutpos_map;
    cutpos_rec tmp_cutpos;
    cut_rec tmp_cut;
    mut_t *arr_ptr;

    
    for (int j = 0; j < cut_vec.size(); ++j) {
        tmp_cut = cut_vec[j];

        // Search for hap1, currently don't handle cut sites induced by indel, only SNP
        arr_ptr = hap1->s;
        for (int i = 0; i < hap1->l; ++i) {
            auto found = std::search(arr_ptr + i, arr_ptr + hap1->l, tmp_cut.seq.begin(), tmp_cut.seq.end());
            if(found == arr_ptr + hap1->l){ break;}
            tmp_cutpos = {};
            tmp_cutpos.pos  = found - arr_ptr;
            tmp_cutpos.type = j;
            cutpos_map[tmp_cutpos] = 1;
        }

        // Search for hap2
        arr_ptr = hap2->s;
        for (int i = 0; i < hap2->l; ++i) {
            auto found = std::search(arr_ptr + i, arr_ptr + hap2->l, tmp_cut.seq.begin(), tmp_cut.seq.end());
            if(found == arr_ptr + hap2->l){ break;}
            tmp_cutpos = {};
            tmp_cutpos.pos  = found - arr_ptr;
            tmp_cutpos.type = j;
            cutpos_map[tmp_cutpos] |= 0x2; // default is 0 if the elements not in the map, |=0x2 will give 2
        }
    }

    if(cutpos_map.empty()){fprintf (stderr, "[%s] ERROR: No cut site found in hap1 and hap2. Exit...\n", __func__); exit (EXIT_FAILURE);}

    // put into the vector
    for (auto it = cutpos_map.begin(); it != cutpos_map.end(); ++it) {
        tmp_cutpos = it->first;
        tmp_cutpos.haplo = it->second; // 1: hap1, 2: hap2, 3: both
        cutpos_vec.push_back(tmp_cutpos);
    }

    std::sort(cutpos_vec.begin(), cutpos_vec.end());
}

void gen_cut_frag(const kseq_t *ks, expt_param *expt_set, std::vector<frag_rrbs_rec> &frag_vec, std::vector<cutpos_rec>& cutpos_vec, std::vector<cut_rec>& cut_vec)
{
    frag_vec.clear();

    // generate potential intervals
    int frag_len;
    frag_rrbs_rec tmp_frag;

    std::map<int8_t, int8_t> cut_map;
    std::map<int8_t, int8_t> tmp_cutmap;
    for(size_t i=0; i<cut_vec.size(); ++i){cut_map[(int8_t)i] = 0;}
    
    for (int i = -1; i <= int(cutpos_vec.size()); ++i){
        if(i==-1 || i == int(cutpos_vec.size())){ //check the first and last frag_recs
            tmp_frag = {};
            if(i==-1){
                frag_len = cutpos_vec[0].pos;
                if (frag_len > expt_set->min_insert && frag_len < expt_set->max_insert){
                    tmp_frag.pos_l = 0;
                    tmp_frag.pos_r = cutpos_vec[0].pos;
                    tmp_frag.cut_l = 0;
                    tmp_frag.cut_r = cutpos_vec[0].type;
                    tmp_frag.haplo = cutpos_vec[0].haplo;
                    ++tmp_frag.n_cuts[cutpos_vec[0].type];
                    frag_vec.push_back(tmp_frag);
                }
            }else{
                frag_len = ks->seq.l - cutpos_vec[i-1].pos;
                if (frag_len > expt_set->min_insert && frag_len < expt_set->max_insert){
                    tmp_frag.pos_l = cutpos_vec[i-1].pos;
                    tmp_frag.pos_r = ks->seq.l;
                    tmp_frag.cut_l = cutpos_vec[i-1].type;
                    tmp_frag.cut_r = 0;
                    tmp_frag.haplo = cutpos_vec[i-1].haplo;
                    ++tmp_frag.n_cuts[cutpos_vec[i-1].type];
                    frag_vec.push_back(tmp_frag);
                }
            }
            continue;
        }

        // generate fragments that have middle cuts while take haplotype into consideration
        // 3,3,2,1,2,1,3...
        // +,+,+,-,+,-,+...
        // +,+,-,+,-,+,+...
        //   +,-,+,-,+,+...
        //   +,+,-,+,-,+...
        //     +,-,+,-,+...
        //       +,-,+,+...
        //         +,-,+...
        //           +,+...
        // should have better way to achieve it though

        int haplotype   = cutpos_vec[i].haplo;
        int8_t left_cut = cutpos_vec[i].type;
        
        std::map<frag_rrbs_rec, int> frag_map;
        frag_map.clear();

        // first hap1
        tmp_cutmap = cut_map;
        ++tmp_cutmap[cutpos_vec[i].type];
        for(int j=i+1; j < cutpos_vec.size(); ++j){
            if((cutpos_vec[j].haplo & 0x2) == 0){continue;} // skip hap2-only cuts
            frag_len = cutpos_vec[j].pos - cutpos_vec[i].pos;
            if(frag_len < expt_set->min_insert){continue;} // skip the small fragments (too short)
            if(frag_len > expt_set->max_insert){break;}
            haplotype &= cutpos_vec[j].haplo;
            ++tmp_cutmap[cutpos_vec[j].type];

            tmp_frag = {};
            tmp_frag.pos_l = cutpos_vec[i].pos;
            tmp_frag.pos_r = cutpos_vec[j].pos;
            tmp_frag.cut_l = cutpos_vec[i].type;
            tmp_frag.cut_r = cutpos_vec[j].type;
            tmp_frag.haplo = haplotype;
            tmp_frag.n_cuts= tmp_cutmap;
            ++frag_map[tmp_frag];
        }

        // then hap2
        tmp_cutmap = cut_map;
        ++tmp_cutmap[cutpos_vec[i].type];
        for(int j=i+1; j < cutpos_vec.size(); ++j){    
            if((cutpos_vec[j].haplo & 0x1) == 0){continue;} // skip hap1-only cuts
            frag_len = cutpos_vec[j].pos - cutpos_vec[i].pos;
            if(frag_len < expt_set->min_insert){continue;} // skip the small fragments (too short)
            if(frag_len > expt_set->max_insert){break;}
            haplotype &= cutpos_vec[j].haplo;
            ++tmp_cutmap[cutpos_vec[j].type];

            tmp_frag = {};
            tmp_frag.pos_l = cutpos_vec[i].pos;
            tmp_frag.pos_r = cutpos_vec[j].pos;
            tmp_frag.cut_l = cutpos_vec[i].type;
            tmp_frag.cut_r = cutpos_vec[j].type;
            tmp_frag.haplo = haplotype;
            tmp_frag.n_cuts= tmp_cutmap;
            ++frag_map[tmp_frag];
        }
        
        // put into the vector
        for (auto it = frag_map.begin(); it != frag_map.end(); ++it) {
            tmp_frag = it->first;
            frag_vec.push_back(tmp_frag);
        }
    }
}

void output_rrcut_bed(const char *fname, const char *chr_id, std::vector<frag_rrbs_rec> &frag_vec, bool to_stdout)
{
    FILE *fp = to_stdout ? stdout : fopen(fname, "a");
    if(fp==NULL){fprintf(stderr, "[%s] ERROR: open rrbs bed file: %s failed. Exit...\n", __func__, fname); exit (EXIT_FAILURE);}

    frag_rrbs_rec tmp_frag;
    int map_size = frag_vec[0].n_cuts.size();
    for (size_t i=1; i < frag_vec.size(); ++i){
        tmp_frag = frag_vec[i];
        fprintf(fp, "%s\t%d\t%d\t.\t1\t%d\t%d\t%d", chr_id, tmp_frag.pos_l, tmp_frag.pos_r, tmp_frag.haplo, tmp_frag.cut_l, tmp_frag.cut_r);
        for(int j=0; j < map_size; ++j){fprintf(fp, "\t%d", tmp_frag.n_cuts[j]);}
        fprintf(fp, "\n");
        //if (ferror(fp)) {fprintf(stderr, "[%s] ERROR: failed to write to file %s. Exit... \n", __func__, fname);exit(EXIT_FAILURE);}
    }
    if(!to_stdout){fclose(fp);}
}

static int simu_usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "rrcut (a module in htsim) for Reduced Representative fragment generation)\n");
    fprintf(stderr, "Version: %s\n", PACKAGE_VERSION);
    fprintf(stderr, "Contact: Wenbin Guo <wbguo@ucla.edu>; \n\n");
    fprintf(stderr, "Usage:   rrcut [options] <ref.fa> \n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "general setting:\n");
    fprintf(stderr, "         -K INT        minimum insert size [%d]\n", MIN_INSERT);
    fprintf(stderr, "         -L INT        maximum insert size [%d]\n", MAX_INSERT);
    fprintf(stderr, "         -c STRING     contig name, only output reads from this contig, default is to output all contigs [None]\n");
    fprintf(stderr, "mutation setting:\n");
    fprintf(stderr, "         -v STRING     path to the genetic variant file (.vcf/vcf.gz) [None]\n");
    fprintf(stderr, "         -R FLOAT      rate of mutations [%.4f]\n", MUT_RATE);
    fprintf(stderr, "         -F FLOAT      fraction of indels [%.2f]\n", INDEL_FRAC);
    fprintf(stderr, "         -X FLOAT      probability an indel is extended [%.2f]\n", INDEL_EXTN);
    fprintf(stderr, "         -H INT        haplotype mode: 0 for disable, nonzero for enable (all variants are homozygotes) [0]\n");
    fprintf(stderr, "         -s INT        seed for random generator [-1]\n");
    fprintf(stderr, "technology setting:\n");
    fprintf(stderr, "         -b STRING     output BED file, output to stdout if not specified [None]\n");
    fprintf(stderr, "         -x STRING     enzyme cutting site string for reduced representation sequencing [None]\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char *argv[])
{
    // default parameters
    char none_default[] = "None";
    char *chr_id    = none_default;
    char *vcf_file  = none_default;
    char *cut_str   = none_default;
    char *bed_file  = none_default;

    expt_param expt_set;
    mut_param  mut_set;

    //update parameters from command line
    int c = 0;
    while ((c = getopt(argc, argv, "K:L:R:F:X:H:s:c:v:b:x:")) >= 0) {
        switch (c) {
            case 'K': expt_set.min_insert = atoi(optarg); break;
            case 'L': expt_set.max_insert = atoi(optarg); break;
            
            case 'R': mut_set.mut_rate    = atof(optarg); break;
            case 'F': mut_set.indel_frac  = atof(optarg); break;
            case 'X': mut_set.indel_extn  = atof(optarg); break;
            case 'H': mut_set.is_hap      = atoi(optarg)!=0; break;
            case 's': mut_set.seed_snp    = atoi(optarg); break;

            case 'c': chr_id     = optarg; break;
            case 'v': vcf_file   = optarg; break;
            case 'b': bed_file   = optarg; break;
            case 'x': cut_str    = optarg; break;
        }
    }
    if (argc - optind < 1) return simu_usage();
    if (mut_set.seed_snp <= 0) mut_set.seed_snp    = time(0)&0x7fffffff;
    
    expt_set.min_insert = std::max(std::max(expt_set.size_l, expt_set.size_r), expt_set.min_insert); // ensure MIN_INSERT >= SIZE_L or SIZE_R
    expt_set.is_chr_set   = strcmp(chr_id, "None") && strlen(chr_id);
    expt_set.is_bed_set   = strcmp(bed_file,"None") && strlen(bed_file);
    expt_set.is_site_set  = strcmp(cut_str, "None") && strlen(cut_str);
    mut_set.is_vcf_set    = strcmp(vcf_file,"None") && strlen(vcf_file);


    // if the bed if designated, delete it, if not, output to stdout
    bool to_stdout = false;
    if(!expt_set.is_bed_set){fprintf(stderr, "No BED file specified, will output to stdout...\n"); to_stdout = true;
    }else{FILE *bed= fopen(bed_file, "r"); if (bed){fclose(bed);remove(bed_file);}}
    
    if (mut_set.is_vcf_set) {
        FILE *vcf;
        if((vcf=fopen(vcf_file,"r"))){fprintf(stderr, "VCF file: %s, use it to simulate reads\n", vcf_file); fclose(vcf);
        }else{fprintf(stderr, "ERROR: The specified VCF file does not exist, please check!\n"); exit(EXIT_FAILURE);}
    } else {
        fprintf(stderr, "[%s] No VCF input, will generate SNP randomly if mutation rate is nonzero\n", __func__);
    }


    // parse the cut site, hold
    parse_cut_site(cut_str, cut_vec);
    

    // parse reference
    kseq_t *ks;
    mutseq_t rseq[2];
    gzFile fp_fa;

    fp_fa = gzopen(argv[optind], "r");
    ks = kseq_init(fp_fa);
    int l;
    while ((l = kseq_read(ks)) >= 0) {  //here l is the chromosome length
        if (expt_set.is_chr_set) {if (strcmp(chr_id, ks->name.s)!=0){continue;}}
        if (l < expt_set.mean_insert + 3 * expt_set.sd_insert) {continue;}

        uint32_t* posidx_arr = (uint32_t*) calloc(ks->seq.l, sizeof(uint32_t));
        if (posidx_arr == NULL) { fprintf(stderr, "ERROR: could not allocate memory\n");exit(EXIT_FAILURE);}

        // introduce mutations and print them to stdout
        //fprintf(stdout, "Contig Variant Start\n");
        if(mut_set.is_vcf_set){
            sim_mut_vcf(ks, vcf_file, rseq, rseq+1, posidx_arr, snp_vec);
            if(snp_vec.size() == 0){fprintf(stdout, "%s\n", ks->name.s);}       //if no variants, print chromosome id
        } else {
            sim_mut_diref(ks, &mut_set, rseq, rseq+1, posidx_arr);
            if(mut_set.mut_rate == 0.0){fprintf(stdout, "%s\n", ks->name.s);}  //if no variants, print chromosome id
        }
        //sim_print_mutref(ks->name.s, ks, rseq, rseq+1, expt_set.output_fmt);
        //fprintf(stdout, "Contig Variant End\n");

        gen_cut_pos(rseq, rseq+1, cutpos_vec, cut_vec);
        gen_cut_frag(ks, &expt_set, frag_rrbs_vec, cutpos_vec, cut_vec);
        output_rrcut_bed(bed_file, ks->name.s, frag_rrbs_vec, to_stdout);
        free(posidx_arr);
    }
    kseq_destroy(ks);
    return 0;
}
