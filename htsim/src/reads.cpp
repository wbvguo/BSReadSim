#include <stdio.h>
#include <algorithm>
#include <string>
#include <vector>
#include <htslib/sam.h>
#include "struct.h"
#include "utility.h"
#include "reads.h"


static const char BASE_LETTERS[] = "ACGTN";
static const uint32_t NO_POS = 0xffffffffu;     /* an inserted base in frag_out.ref_pos */
static const int MAPQ = 60;
static const int GZIP_LEVEL = 4;    /* of FASTQ.gz and BAM (see reads.h) */

void reads_open(reads_out *ro, const expt_param *expt_set, const std::vector<chr_rec>& chr_vec)
{
    ro->mates = expt_set->paired_end ? 2 : 1;
    for (int mate = 0; mate < ro->mates; ++mate) {
        ro->quality[mate].assign(expt_set->read_length[mate], (char)(33 + expt_set->phred));
    }
    const std::string &prefix = expt_set->output;
    for (int k = 0; k < 2; ++k) {       // fastq.gz, fastq
        if (!(expt_set->formats & 1 << k)) continue;
        for (int mate = 0; mate < ro->mates; ++mate) {
            out_t *out = out_open(prefix + (mate ? ".R2." : ".R1.") + READ_FORMATS[k], k == 0, expt_set->threads,
                                  GZIP_LEVEL);
            ro->fastq[mate].push_back(out);
            ro->paths.push_back(out->path);
        }
    }
    if (expt_set->formats & (READ_BAM | READ_SAM)) {
        std::string text = "@HD\tVN:1.6\tSO:unsorted\n";
        for (size_t i = 0; i < chr_vec.size(); ++i) {
            text += "@SQ\tSN:" + chr_vec[i].name + "\tLN:" + std::to_string(chr_vec[i].chr_len) + "\n";
        }
        text += std::string("@PG\tID:htsim\tPN:htsim\tVN:") + PACKAGE_VERSION + "\n";
        text += "@CO\tMAPQ 60 denotes simulated origin, not calibrated mapping confidence\n";
        ro->hdr = sam_hdr_init();
        if (ro->hdr == 0 || sam_hdr_add_lines(ro->hdr, text.c_str(), text.size()) < 0) die("cannot build the SAM header");
        for (int k = 2; k < 4; ++k) {   // bam, sam
            if (!(expt_set->formats & 1 << k)) continue;
            std::string path = prefix + "." + READ_FORMATS[k];
            samFile *fp = sam_open(path.c_str(), k == 2 ? ("wb" + std::to_string(GZIP_LEVEL)).c_str() : "w");
            if (fp == 0) die("cannot create %s", path.c_str());
            if (expt_set->threads > 1 && hts_set_threads(fp, expt_set->threads) < 0) die("cannot start the threads writing %s", path.c_str());
            if (sam_hdr_write(fp, ro->hdr) < 0) die("cannot write the header of %s", path.c_str());
            ro->alignments.push_back(fp);
            ro->paths.push_back(path);
        }
        ro->rec = bam_init1();
    }
}

static void add_op(std::vector<uint32_t> &cigar, int op, uint32_t n)
{
    if (!cigar.empty() && bam_cigar_op(cigar.back()) == (uint32_t)op) cigar.back() += n << BAM_CIGAR_SHIFT;
    else cigar.push_back(bam_cigar_gen(n, op));
}

// where the read of template bases [begin, begin + len) aligns, as the Python package
// aligns it: M for a base with a reference position, I for an inserted one, D for the
// reference bases skipped between; a read inside an insertion anchors at the next base
static void align_read(const frag_out *f, int begin, int len, int contig_len, int64_t *pos, int64_t *end,
                       std::vector<uint32_t> &cigar)
{
    cigar.clear();
    int64_t first = -1, prev = -1;
    for (int k = begin; k < begin + len; ++k) {
        int64_t p = f->ref_pos[k] == NO_POS ? -1 : (int64_t)f->ref_pos[k];
        if (p < 0) {add_op(cigar, BAM_CINS, 1); continue;}
        if (prev >= 0 && p > prev + 1) add_op(cigar, BAM_CDEL, (uint32_t)(p - prev - 1));
        add_op(cigar, BAM_CMATCH, 1);
        if (first < 0) first = p;
        prev = p;
    }
    if (first >= 0) {*pos = first; *end = prev + 1; return;}
    int n = (int)f->ref_pos.size();
    int64_t anchor = -1;
    for (int k = begin + len; k < n && anchor < 0; ++k) if (f->ref_pos[k] != NO_POS) anchor = f->ref_pos[k];
    for (int k = begin - 1; k >= 0 && anchor < 0; --k) if (f->ref_pos[k] != NO_POS) anchor = (int64_t)f->ref_pos[k] + 1;
    *pos = std::min<int64_t>(anchor, contig_len - 1);
    *end = *pos + 1;
}

void cut_reads(const expt_param *expt_set, int chr, uint64_t ordinal, const ref_seq *ref,
               const frag_out *f, read_pair *reads)
{
    Rng rng(seed_of(expt_set->seed, chr, ordinal, 1));     // gen_frag draws from (chr, ordinal, 0)
    int n = (int)f->bases.size();
    bool read1_reverse = f->strand == OB || f->strand == CTOT;
    for (int mate = 0; mate < (expt_set->paired_end ? 2 : 1); ++mate) {
        int len = expt_set->read_length[mate];
        bool reverse = read1_reverse != (mate == 1);
        reads->reverse[mate] = reverse;
        std::string &seq = reads->seq[mate];
        seq.resize(len);
        for (int k = 0; k < len; ++k) {
            int base = reverse ? f->bases[n - 1 - k] : f->bases[k];
            if (base < 4) {
                if (reverse) base = 3 - base;
                if (rng.unif() < expt_set->error_rate) base = (base + 1 + (int)rng.below(3)) & 3;
            }
            seq[k] = BASE_LETTERS[base];
        }
        if (expt_set->formats & (READ_BAM | READ_SAM)) {
            align_read(f, reverse ? n - len : 0, len, ref->len, &reads->pos[mate], &reads->end[mate], reads->cigar[mate]);
        }
    }
}

static char complement(char base)
{
    switch (base) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

static std::string cigar_text(const std::vector<uint32_t> &cigar)
{
    std::string text;
    for (size_t i = 0; i < cigar.size(); ++i) text += std::to_string(bam_cigar_oplen(cigar[i])) + bam_cigar_opchr(cigar[i]);
    return text;
}

static void write_alignments(reads_out *ro, const std::string &name, int chr, const read_pair *reads)
{
    bool paired = ro->mates == 2;
    int64_t span = paired ? std::max(reads->end[0], reads->end[1]) - std::min(reads->pos[0], reads->pos[1]) : 0;
    for (int mate = 0; mate < ro->mates; ++mate) {
        int other = 1 - mate, len = (int)reads->seq[mate].size();
        std::string seq = reads->seq[mate], qual = ro->quality[mate];   // reference-forward
        if (reads->reverse[mate]) {
            std::reverse(seq.begin(), seq.end());
            for (size_t k = 0; k < seq.size(); ++k) seq[k] = complement(seq[k]);
            std::reverse(qual.begin(), qual.end());
        }
        for (size_t k = 0; k < qual.size(); ++k) qual[k] -= 33;
        uint16_t flag = reads->reverse[mate] ? BAM_FREVERSE : 0;
        int32_t mtid = -1;
        int64_t mpos = -1, isize = 0;
        if (paired) {
            flag |= BAM_FPAIRED | BAM_FPROPER_PAIR | (mate ? BAM_FREAD2 : BAM_FREAD1) | (reads->reverse[other] ? BAM_FMREVERSE : 0);
            mtid = chr;
            mpos = reads->pos[other];
            bool leftmost = reads->pos[mate] < reads->pos[other] || (reads->pos[mate] == reads->pos[other] && mate == 0);
            isize = leftmost ? span : -span;
        }
        const std::vector<uint32_t> &cigar = reads->cigar[mate];
        if (bam_set1(ro->rec, name.size(), name.c_str(), flag, chr, reads->pos[mate], MAPQ, cigar.size(), cigar.data(),
                     mtid, mpos, isize, len, seq.c_str(), qual.c_str(), 32) < 0
            || bam_aux_update_int(ro->rec, "AS", len) < 0
            || (paired && (bam_aux_update_int(ro->rec, "MQ", MAPQ) < 0
                           || bam_aux_update_str(ro->rec, "MC", -1, cigar_text(reads->cigar[other]).c_str()) < 0))) {
            die("cannot build the BAM record of %s", name.c_str());
        }
        for (size_t i = 0; i < ro->alignments.size(); ++i) {
            if (sam_write1(ro->alignments[i], ro->hdr, ro->rec) < 0) die("cannot write an alignment to %s", ro->alignments[i]->fn);
        }
    }
}

void reads_add(reads_out *ro, const char *contig, int chr, const frag_out *f, const read_pair *reads)
{
    char ordinal[24];
    snprintf(ordinal, sizeof(ordinal), "%llx", (unsigned long long)ro->n_frag);
    std::string name = std::string(contig) + ":" + std::to_string(f->start + 1) + "-" + std::to_string(f->end) + ":" + ordinal;
    for (int mate = 0; mate < ro->mates; ++mate) {
        for (size_t i = 0; i < ro->fastq[mate].size(); ++i) {
            out_printf(ro->fastq[mate][i], "@%s/%d\n%s\n+\n%s\n", name.c_str(), mate + 1, reads->seq[mate].c_str(),
                       ro->quality[mate].c_str());
        }
    }
    if (!ro->alignments.empty()) write_alignments(ro, name, chr, reads);
    ++ro->n_frag;
}

void reads_close(reads_out *ro, uint64_t skipped)
{
    for (int mate = 0; mate < 2; ++mate) {
        for (size_t i = 0; i < ro->fastq[mate].size(); ++i) out_close(ro->fastq[mate][i]);
    }
    for (size_t i = 0; i < ro->alignments.size(); ++i) {
        std::string path = ro->alignments[i]->fn;
        if (sam_close(ro->alignments[i]) != 0) die("cannot finish %s", path.c_str());
    }
    if (ro->hdr) {
        sam_hdr_destroy(ro->hdr);
        bam_destroy1(ro->rec);
    }
    std::string paths;
    for (size_t i = 0; i < ro->paths.size(); ++i) paths += (i ? ", " : "") + ro->paths[i];
    fprintf(stderr, "htsim: %llu fragments (%llu could not be placed) in %s\n",
            (unsigned long long)ro->n_frag, (unsigned long long)skipped, paths.c_str());
}
