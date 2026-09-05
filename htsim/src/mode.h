#ifndef MODE_H
#define MODE_H

#include <stdint.h>
#include <vector>
#include "struct.h"

// Where fragments come from. Pass 1 measures what each contig can produce so the
// fragments can be divided among the contigs (for WGBS/WGS while the FASTA is
// indexed; RRBS candidates and usable targets are kept for pass 2); pass 2 places
// them. Fragment `ordinal` of a contig draws from an Rng seeded by the contig and
// the ordinal, so the output does not depend on the threads.


// for BED: capture targets (TBS/WES/TS), per contig, sorted
void parse_bed(const char *fname, const std::vector<chr_rec>& chr_vec, std::vector<std::vector<probe_rec>>& probe_vec, uint64_t *rows);

// for GC-bias: probability of each GC bin
void parse_bias_file(const char *fname, std::vector<double>& eff_vec);


// what the fragments of a contig are drawn from
typedef struct {
    std::vector<frag_rec> cands;        /* RRBS candidates as fragments (haplo -1: either haplotype) */
    std::vector<int> usable;            /* targets that can be sampled */
    std::vector<double> cum;            /* cumulative weights of cands or usable */
    std::vector<std::pair<int, int>> reach;  /* reference intervals its fragments fall in, sorted and disjoint */
} frag_domain;

// pass 1: the weight (score) and effective length of a contig
void wg_domain(const ref_seq *ref, const expt_param *expt_set, int bins, chr_rec *tmp_chr);
void rrbs_domain(const std::vector<frag_rrbs_rec>& cands, frag_domain *dom, chr_rec *tmp_chr);
void probe_domain(const ref_seq *ref, const expt_param *expt_set, std::vector<probe_rec>& probe_vec, frag_domain *dom, chr_rec *tmp_chr);

// turn the GC profile (eff_vec) into acceptance probabilities and re-weight the contigs
void calibrate_gc(std::vector<double>& eff_vec, std::vector<chr_rec>& chr_vec, bool drop_unreachable);

// divide the fragments among the contigs (largest remainder)
void cal_chr_count(const expt_param *expt_set, std::vector<chr_rec>& chr_vec);


// pass 2: fragment `ordinal` of contig chr and its template; false if no sequenceable fragment was found
bool gen_frag(int chr, uint64_t ordinal, const ref_seq *ref, mutseq_t *rseq, const expt_param *expt_set, std::vector<double>& eff_vec,
              std::vector<probe_rec>& probe_vec, frag_domain *dom, frag_rec *frag, frag_seq *fs);

#endif
