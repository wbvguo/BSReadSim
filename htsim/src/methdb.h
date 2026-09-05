#ifndef METHDB_H
#define METHDB_H

#include <vector>
#include <map>
#include <htslib/hts.h>
#include "struct.h"
#include "utility.h"

// MethDB (version 1): BGZF-compressed text, one section per reference contig
//
//   #methdb   1
//   #contig   name  length  md5                 every reference contig, FASTA order
//   >         name                              starts a contig's section
//   v  pos kind ref alt hap id source ctx1 meth1 type1 ctx2 meth2 type2      one variant
//   s  step ctx meth type                                       a C/G site shared by hap1 and hap2
//   h  step ref_ctx ctx1 meth1 type1 ctx2 meth2 type2           any other C/G site
//
// Variant positions are 0-based snpmeth_map keys. A variant lists the context,
// level (round(p*65535)) and meth_type_t of each of its bases on hap1 and hap2
// (comma-separated, "." when it has none). A site is a meth_vec record: step is
// its position minus that of the previous site of the contig (its position for
// the first), and an s row is a site whose context is its reference context and
// whose level and type are the same on both haplotypes (almost all of them).


// create the methdb: one record per reference site, and per C/G that a variant may turn into one
// (within 2 haplotype bases of it); only within spans (sorted, disjoint intervals) when given
void create_methdb(const ref_seq *ref, const std::map<int, snpmeth_rec>& snpmeth_map, std::vector<meth_rec>& meth_vec,
                   bool collect_non_cpg, const std::vector<std::pair<int, int>> *spans = 0);
// context of reference base i (0 for none; non-CpG dropped unless collect_non_cpg)
uint8_t ref_context(const ref_seq *ref, int i, bool collect_non_cpg);


// fill with distribution: Beta, or the context pools when pooling
uint16_t gen_meth(const meth_param *meth_set, int chr, uint64_t key, uint8_t context, const std::vector<uint16_t> (&pools)[4], uint8_t *type);
// draw the level of every reference site that has none yet
void fill_levels(std::vector<meth_rec>& meth_vec, const meth_param *meth_set, int chr, const std::vector<uint16_t> (&pools)[4], int threads);
// sites created or changed by variants, on each haplotype
void update_variant(const ref_seq *ref, mutseq_t *hap1, mutseq_t *hap2, std::vector<meth_rec>& meth_vec,
                    std::map<int, snpmeth_rec>& snpmeth_map, const meth_param *meth_set, int chr, const std::vector<uint16_t> (&pools)[4]);


// save/load/export MethDB
void save_methdb_header(out_t *out, const std::vector<chr_rec>& chr_vec);
void save_methdb_chr(out_t *out, const ref_seq *ref, std::vector<meth_rec>& meth_vec, std::map<int, snpmeth_rec>& snpmeth_map);

typedef struct {
    text_in *in;
    bool pending;       /* in->line holds the next unread line */
} methdb_reader;
methdb_reader *open_methdb(const char *fname, const std::vector<chr_rec>& chr_vec);
// read the section of a contig; its sites are skipped when meth_vec is null
void load_methdb_chr(methdb_reader *db, const ref_seq *ref, std::map<int, snpmeth_rec>& snpmeth_map,
                     std::vector<meth_rec> *meth_vec);
void close_methdb(methdb_reader *db);
void export_methdb(const char *fname, out_t *out);

#endif
