#ifndef STREAM_H
#define STREAM_H

// The fragment stream read by the Python package (bsreadsim.htsim and bsreadsim.batch):
//
//   header   one line of JSON: the version (the package's own, which it checks),
//            master seed, technology, details, read layout, and contigs
//            (name, length, MD5)
//   batch*   u64 payload size, then the payload (little-endian columns)
//   end      u64 zero
//   summary  one line of JSON
//
// A batch payload is the u64 ordinal of its first fragment followed by the
// columns below, each written as a u64 byte count and the raw values.
//
//   per fragment (n):  contig u32, start u32, end u32, haplotype u8, strand u8
//   template_offsets u32[n+1], bases u8[T]
//   site_offsets u32[n+1], site_positions u32[S], site_levels u16[S],
//   site_contexts u8[S], site_kmers u16[S]
//   event_offsets u32[n+1], event_starts u32[E], event_ends u32[E],
//   event_kinds u8[E], event_deleted u32[E]
//
// start/end is the reference envelope; bases are in reference orientation.
// strand: 0 OT, 1 OB, 2 CTOT, 3 CTOB. Sites are the cytosines of both strands:
// a level is the methylation probability times 65535 (rounded); a context is
// 1 CG, 2 CHG, or 3 CHH, plus 0x80 when the level comes from an ASM input; a
// k-mer (only with --site-kmers) is the index of the forward-strand 7-mer
// centred on the site (A=0, C=1, G=2, T=3, first base highest), 65535 if it
// has an N or passes a contig end. Events (only with details) are the
// variants on the template: kind 1 SNV, 2 insertion, 3 deletion, with a
// deletion an empty interval where it removes event_deleted reference bases.
// The reference position of every template base follows from start and the
// events, so it is not sent.

#include <stdint.h>
#include <vector>
#include <map>
#include "struct.h"

// one fragment of the stream
typedef struct {
    uint32_t contig = 0, start = 0, end = 0;
    uint8_t haplo = 0, strand = 0;
    std::vector<uint8_t> bases;
    std::vector<uint32_t> ref_pos;          /* with details; for htsim's own BAM, not streamed */
    std::vector<uint32_t> site_pos;
    std::vector<uint16_t> site_level;
    std::vector<uint8_t> site_ctx;          /* context class, 0x80 set for an ASM level */
    std::vector<uint16_t> site_kmer;
    std::vector<uint32_t> ev_start, ev_end, ev_deleted;
    std::vector<uint8_t> ev_kind;
} frag_out;

// the sites of a contig as fill_frag_out finds them: meth_vec (in position order) and,
// for each block of 64 reference bases, the index of its first record at or after the block
typedef struct {
    const std::vector<meth_rec> *meth_vec = 0;
    std::vector<uint32_t> block;
} site_index;
static const int SITE_BLOCK_SHIFT = 6;
void index_sites(const std::vector<meth_rec>& meth_vec, int len, site_index *sites);

// fill a fragment from its template; without sites it has none, and site k-mers only when asked
void fill_frag_out(int chr, frag_rec *frag, frag_seq *fs, const site_index *sites,
                   std::map<int, snpmeth_rec>& snpmeth_map, bool details, bool kmers, frag_out *out);

typedef struct {
    int batch_size = 1024;
    bool paired_end = true;
    uint64_t n_frag = 0, n_mate = 0, n_base = 0, n_site = 0;
    std::vector<uint64_t> per_contig;
    // the columns of the pending batch
    uint64_t first = 0, n = 0;
    std::vector<uint32_t> contig, start, end;
    std::vector<uint8_t> haplo, strand, bases;
    std::vector<uint32_t> tmpl_offsets;
    std::vector<uint32_t> site_offsets, site_pos;
    std::vector<uint16_t> site_level;
    std::vector<uint8_t> site_ctx;
    std::vector<uint16_t> site_kmer;
    std::vector<uint32_t> ev_offsets, ev_start, ev_end, ev_deleted;
    std::vector<uint8_t> ev_kind;
} stream_rec;

// write the header
void stream_open(stream_rec *st, const expt_param *expt_set, const std::vector<chr_rec>& chr_vec);
void stream_add(stream_rec *st, frag_out *f);
void stream_close(stream_rec *st, uint64_t skipped);

#endif
