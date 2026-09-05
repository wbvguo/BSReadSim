#ifndef PROFILE_H
#define PROFILE_H

// Methylation inputs, checked against the reference and applied to the methdb:
// profiles (CGmap, bedMethyl, MethBG, MethBED) and allele-specific methylation
// (CGmapTools ASS, ASM BED).

#include <vector>
#include <map>
#include "struct.h"
#include "utility.h"


// methylation profiles, read contig by contig in FASTA order
typedef struct {
    text_in *in;
    int format;
    std::map<std::string, int> chr_idx;
    std::vector<uint32_t> chr_len;
    // the row read ahead
    bool pending;
    int chr, pos, next_base;    /* CGmap: the dinucleotide's second base */
    char strand;                /* BED formats: '+', '-', or '.' */
    int cclass;                 /* MethBED: context class, 0 if not given */
    uint16_t meth;
    bool defined;
    uint8_t context;
    uint64_t rows, defined_rows;
    size_t ncol;
} profile_rec;

profile_rec *open_profile(const char *fname, int format, const std::vector<chr_rec>& chr_vec);
void close_profile(profile_rec *pf);
// fill the methdb (or, when pooling, the context pools) from the rows of contig chr; returns rows used
int fill_profile_chr(profile_rec *pf, const ref_seq *ref, int chr, std::vector<meth_rec>& meth_vec,
                     const meth_param *meth_set, std::vector<uint16_t> (&pools)[4]);


// for ASM: CGmapTools "asm -m ass" rows called TRUE, or ASM BED rows
typedef struct {
    int target, snv;                 // 0-based positions
    uint8_t ref, alt;                // linked SNV alleles
    uint16_t ref_meth, alt_meth;
    uint8_t context;                 // strand hint, then the target context
} asm_rec;

void parse_asm(const char *fname, bool is_bed, const std::vector<chr_rec>& chr_vec,
               std::vector<std::vector<asm_rec>>& asm_vec, uint64_t *rows);
// check rows against the reference, set their context, sort them by target
void check_asm_chr(const ref_seq *ref, std::vector<asm_rec>& asm_list);
// the heterozygous SNVs implied by ASM rows (when there is no VCF)
void asm_variants(int chr_idx, std::vector<asm_rec>& asm_list, uint64_t seed_phase, std::map<int, snpmeth_rec>& snpmeth_map);
// set the allele-specific levels; fails unless each target is a shared site linked to a het SNV
void fill_asm_chr(const ref_seq *ref, std::vector<asm_rec>& asm_list, std::vector<meth_rec>& meth_vec,
                  std::map<int, snpmeth_rec>& snpmeth_map, bool collect_non_cpg);

#endif
