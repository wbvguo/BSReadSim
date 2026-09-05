#ifndef READS_H
#define READS_H

// Reads written by htsim itself (--output PREFIX, in each --format), for assays without
// bisulfite conversion. Each fragment's reads are cut from its ends as the Python
// package cuts them: read 1 is the leftmost bases for OT and the reverse complement of
// the rightmost for OB, read 2 comes from the other end. Every base gets --phred, and
// each called base becomes another with --error-rate.
//
// FASTQ goes to PREFIX.R1.fastq.gz and PREFIX.R2.fastq.gz (BGZF), or to PREFIX.R1.fastq
// and PREFIX.R2.fastq, with the package's read names: @contig:start-end:ordinal/mate, the
// one-based inclusive fragment envelope and the hex ordinal of the fragment. Alignments go
// to PREFIX.bam and/or PREFIX.sam: unsorted, mates adjacent, and aligned as the package
// aligns them: at their simulated origin, MAPQ 60, an indel-aware CIGAR, reference-forward
// SEQ and QUAL, and the tags AS (and MQ and MC when paired). The truth tags of the
// package (zt, zr, zf, zx) are not written. FASTQ.gz and BAM are compressed at level 4,
// the default --gzip-level of the package.

#include <stdint.h>
#include <string>
#include <vector>
#include <htslib/sam.h>
#include "struct.h"
#include "stream.h"
#include "utility.h"

// the reads of one fragment
typedef struct {
    std::string seq[2];                 /* read 1, read 2: letters, sequencing orientation */
    bool reverse[2] = {false, false};
    int64_t pos[2] = {0, 0}, end[2] = {0, 0};   /* aligned reference interval (BAM) */
    std::vector<uint32_t> cigar[2];
} read_pair;

typedef struct {
    int mates = 2;
    std::string quality[2];             /* the quality line of each mate */
    std::vector<out_t*> fastq[2];       /* the FASTQ files of read 1 and read 2 */
    std::vector<samFile*> alignments;   /* the BAM and SAM files */
    sam_hdr_t *hdr = 0;
    bam1_t *rec = 0;
    std::vector<std::string> paths;
    uint64_t n_frag = 0;
} reads_out;

void reads_open(reads_out *ro, const expt_param *expt_set, const std::vector<chr_rec>& chr_vec);
// the reads of fragment `ordinal` of contig chr (in parallel; errors are drawn from the
// fragment's own seed). Alignments need f->ref_pos (fill_frag_out with details).
void cut_reads(const expt_param *expt_set, int chr, uint64_t ordinal, const ref_seq *ref,
               const frag_out *f, read_pair *reads);
void reads_add(reads_out *ro, const char *contig, int chr, const frag_out *f, const read_pair *reads);
void reads_close(reads_out *ro, uint64_t skipped);

#endif
