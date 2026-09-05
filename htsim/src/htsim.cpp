/* The MIT License
   Copyright (c) 2008 Genome Research Ltd (GRL).
                 2011 Heng Li <lh3@live.co.uk>
                 2022 Wenbin Guo <wbguo@ucla.edu>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:
   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

/* This program is derived from WGSIM(v0.3.1-r13)[https://github.com/lh3/wgsim.git], with very heavy
 * modifications to simulate WGS/RRS/TS or WGBS/RRBS/TBS reads in BSReadSim for diploid organisms.
 * It writes the fragments to stdout (see stream.h) for the Python package to turn into reads,
 * or, on its own and without bisulfite, simple reads to FASTQ and BAM files (see reads.h). */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <functional>
#include <htslib/kseq.h>
#include <htslib/hts.h>
#include <htslib/sam.h>
#include "struct.h"
#include "utility.h"

KSEQ_INIT(gzFile, gzread)
#include "haplo.h"
#include "methdb.h"
#include "profile.h"
#include "mode.h"
#include "rrcut.h"
#include "stream.h"
#include "reads.h"
#include "option.h"


// the inputs read once for the whole reference
typedef struct {
    std::vector<chr_rec> chr_vec;                       /* contigs, in FASTA order */
    std::vector<std::map<int, snpmeth_rec>> vcf_vec;    /* per contig, with --vcf */
    vcf_summary vcf_sum;
    std::vector<std::vector<asm_rec>> asm_vec;          /* per contig, with ASM input */
    uint64_t asm_rows = 0;
    std::vector<std::vector<probe_rec>> probe_vec;      /* per contig, with --targets */
    std::vector<std::map<std::string, rrbs_row>> rrbs_vec;  /* per contig, with --rrbs-candidates */
    std::vector<double> eff_vec;                        /* GC profile, then acceptance per bin */
    std::vector<cut_rec> cut_vec;
} input_set;

// work on one contig of the reference
typedef std::function<void(const ref_seq *ref, int chr)> contig_fn;


/*-------------------------inputs-------------------------*/
static bool gzip_path(const std::string &path) {return path.size() > 3 && !path.compare(path.size() - 3, 3, ".gz");}

// run body(ref, chr) on every contig of the reference, in FASTA order; only here is the FASTA read with kseq
static void for_each_contig(const std::string &fname, const contig_fn &body)
{
    gzFile fp = gzopen(fname.c_str(), "r");
    if (fp == 0) die("cannot open reference FASTA %s", fname.c_str());
    kseq_t *ks = kseq_init(fp);
    for (int chr = 0; kseq_read(ks) >= 0; ++chr) {
        ref_seq ref = {ks->name.s, ks->seq.s, (int)ks->seq.l};
        body(&ref, chr);
    }
    kseq_destroy(ks);
    gzclose(fp);
}

static void finish_contig(chr_rec *tmp_chr, hts_md5_context *md5, uint64_t len)
{
    if (len == 0) die("FASTA contig '%s' has no bases", tmp_chr->name.c_str());
    if (len > (uint64_t)INT32_MAX - 8) die("FASTA contig '%s' is longer than 2^31 - 9 bases", tmp_chr->name.c_str());
    unsigned char digest[16];
    char hex[33];
    hts_md5_final(digest, md5);
    hts_md5_hex(hex, digest);
    hts_md5_destroy(md5);
    tmp_chr->md5 = hex;
    tmp_chr->chr_len = (uint32_t)len;
}

// names, lengths, and MD5 (of the uppercase sequence, as the SAM M5 tag) of every contig;
// letters other than A, C, G, T (IUPAC codes such as M or R) are N for the simulation.
// The FASTA is checked line by line here, so later passes can read it with kseq; when
// on_contig is set it also gets each contig, so this pass can do the work of another.
static void read_fasta_index(const std::string &fname, std::vector<chr_rec>& chr_vec, const contig_fn &on_contig)
{
    text_in *in = in_open(fname.c_str(), "FASTA");
    std::set<std::string> names;
    hts_md5_context *md5 = 0;
    kstring_t seq = {0, 0, 0};      /* the bases of the contig, for on_contig */
    uint64_t len = 0;
    auto finish = [&]() {
        finish_contig(&chr_vec.back(), md5, len);
        if (on_contig) {
            ref_seq ref = {chr_vec.back().name.c_str(), seq.s, (int)len};
            on_contig(&ref, (int)chr_vec.size() - 1);
            seq.l = 0;
        }
    };
    chr_vec.clear();
    while (in_next(in)) {
        kstring_t &line = in->line;
        if (line.l == 0) IN_DIE(in, "empty lines are not allowed");
        if (line.s[0] == '>') {
            if (md5) finish();
            for (size_t i = 0; i < line.l; ++i) {
                unsigned char c = line.s[i];
                if ((c < 0x20 && c != '\t') || c == 0x7f) IN_DIE(in, "header contains a control character");
            }
            chr_rec tmp_chr;
            tmp_chr.name = std::string(line.s + 1, strcspn(line.s + 1, " \t"));
            if (tmp_chr.name.empty()) IN_DIE(in, "contig name is empty");
            if (tmp_chr.name.size() > (1 << 20) || !valid_utf8(tmp_chr.name)) IN_DIE(in, "contig name must be valid UTF-8 of at most 1 MiB");
            if (!names.insert(tmp_chr.name).second) die("FASTA contig name is not unique: %s", tmp_chr.name.c_str());
            chr_vec.push_back(tmp_chr);
            if ((md5 = hts_md5_init()) == 0) die("cannot compute MD5");
            len = 0;
        } else {
            if (!md5) IN_DIE(in, "sequence must start with a '>' header");
            for (size_t i = 0; i < line.l; ++i) {
                if (!isalpha((unsigned char)line.s[i])) IN_DIE(in, "sequence may contain only letters");
                line.s[i] &= ~0x20;     // uppercase, for the MD5
            }
            hts_md5_update(md5, line.s, line.l);
            if (on_contig) kputsn(line.s, line.l, &seq);
            len += line.l;
        }
    }
    in_close(in);
    if (!md5) die("reference FASTA is empty: %s", fname.c_str());
    finish();
    free(seq.s);
}

// the reference index and the inputs checked against it; with sampling, those that place fragments
static void read_inputs(input_set *inputs, const expt_param *expt_set, const mut_param *mut_set, const meth_param *meth_set,
                        bool sampling, const contig_fn &on_contig = contig_fn())
{
    if (sampling && expt_set->sampling == GC_PROFILE) parse_bias_file(expt_set->bias_file.c_str(), inputs->eff_vec);
    read_fasta_index(expt_set->reference, inputs->chr_vec, on_contig);
    std::vector<chr_rec> &chr_vec = inputs->chr_vec;
    if (!mut_set->vcf_file.empty()) parse_vcf(mut_set->vcf_file.c_str(), chr_vec, mut_set->seed_phase, inputs->vcf_vec, &inputs->vcf_sum);
    if (!meth_set->asm_file.empty()) parse_asm(meth_set->asm_file.c_str(), meth_set->asm_is_bed, chr_vec, inputs->asm_vec, &inputs->asm_rows);
    parse_cut_site(expt_set->cut_sites, inputs->cut_vec);
    if (!sampling) return;
    if (!expt_set->bed_file.empty()) parse_bed(expt_set->bed_file.c_str(), chr_vec, inputs->probe_vec, 0);
    if (!expt_set->rrbs_file.empty()) parse_rrcut_bed(expt_set->rrbs_file.c_str(), chr_vec, inputs->rrbs_vec, 0);
}

static profile_rec *open_profile_set(input_set *inputs, const meth_param *meth_set)
{
    if (meth_set->profile_file.empty()) return 0;
    return open_profile(meth_set->profile_file.c_str(), meth_set->profile_format, inputs->chr_vec);
}


/*-------------------------one contig-------------------------*/
// a contig in memory: its variants, haplotypes, and methylome (when built)
typedef struct {
    std::map<int, snpmeth_rec> snpmeth_map;
    mutseq_t rseq[2] = {};
    std::vector<meth_rec> meth_vec;
} contig_rec;

// variants of a contig: from the VCF, else implied by ASM rows, else de novo
static void chr_variants(input_set *inputs, const ref_seq *ref, int chr, const mut_param *mut_set, std::map<int, snpmeth_rec>& snpmeth_map)
{
    if (!inputs->asm_vec.empty()) check_asm_chr(ref, inputs->asm_vec[chr]);
    if (!inputs->vcf_vec.empty()) {
        snpmeth_map = inputs->vcf_vec[chr];
    } else if (!inputs->asm_vec.empty()) {
        asm_variants(chr, inputs->asm_vec[chr], mut_set->seed_phase, snpmeth_map);
    } else {
        sim_mut_diref(ref, mut_set, chr, snpmeth_map);
    }
}

// the variants and haplotypes of a contig, from a MethDB (with its methylome when with_meth) or from the inputs
static void load_contig(input_set *inputs, const ref_seq *ref, int chr, methdb_reader *db, bool with_meth,
                        const mut_param *mut_set, contig_rec *ct)
{
    ct->meth_vec.clear();
    if (db) load_methdb_chr(db, ref, ct->snpmeth_map, with_meth ? &ct->meth_vec : 0);
    else chr_variants(inputs, ref, chr, mut_set, ct->snpmeth_map);
    sim_mut_map(ref, ct->snpmeth_map, ct->rseq, ct->rseq + 1);
}

static void free_contig(contig_rec *ct)
{
    free(ct->rseq[0].s); free(ct->rseq[1].s);
    ct->rseq[0].s = ct->rseq[1].s = 0;
}

// the methylation level of every site of a contig, on each haplotype; with reach, only of the
// sites there and of the ASM targets (which are checked wherever they are)
static void chr_methylome(input_set *inputs, const ref_seq *ref, int chr, contig_rec *ct, profile_rec *pf,
                          const meth_param *meth_set, int threads, const std::vector<std::pair<int, int>> *reach)
{
    std::vector<std::pair<int, int>> spans;
    if (reach) {
        spans = *reach;
        if (!inputs->asm_vec.empty()) {
            for (size_t i = 0; i < inputs->asm_vec[chr].size(); ++i) {
                int target = inputs->asm_vec[chr][i].target;
                spans.push_back(std::make_pair(target, target + 1));
            }
        }
        spans = merge_intervals(spans);
    }
    std::vector<uint16_t> pools[4];
    create_methdb(ref, ct->snpmeth_map, ct->meth_vec, meth_set->collect_non_cpg, reach ? &spans : 0);
    if (pf) fill_profile_chr(pf, ref, chr, ct->meth_vec, meth_set, pools);
    fill_levels(ct->meth_vec, meth_set, chr, pools, threads);
    update_variant(ref, ct->rseq, ct->rseq + 1, ct->meth_vec, ct->snpmeth_map, meth_set, chr, pools);
    if (!inputs->asm_vec.empty()) fill_asm_chr(ref, inputs->asm_vec[chr], ct->meth_vec, ct->snpmeth_map, meth_set->collect_non_cpg);
}

// prepare a contig: what `methdb-build` does for it, whatever its share of the fragments: its
// variants, haplotypes, and (bisulfite assays) methylome, with every input checked against them;
// all from a MethDB when one is given. With reach, the methylome is built only there.
static void prepare_contig(input_set *inputs, const ref_seq *ref, int chr, methdb_reader *db, profile_rec *pf,
                           const expt_param *expt_set, const mut_param *mut_set, const meth_param *meth_set, contig_rec *ct,
                           const std::vector<std::pair<int, int>> *reach = 0)
{
    load_contig(inputs, ref, chr, db, expt_set->bisulfite, mut_set, ct);
    if (!db && expt_set->bisulfite) chr_methylome(inputs, ref, chr, ct, pf, meth_set, expt_set->threads, reach);
}

// the RRBS candidates of a contig, with the scores of --rrbs-candidates
static void rrbs_cands(input_set *inputs, int chr, contig_rec *ct, const expt_param *expt_set, std::vector<frag_rrbs_rec>& cands)
{
    gen_cut_frag(ct->rseq, ct->rseq + 1, !ct->snpmeth_map.empty(), expt_set, inputs->cut_vec, cands);
    if (!inputs->rrbs_vec.empty()) {
        match_rrcut_bed(inputs->chr_vec[chr].name.c_str(), inputs->rrbs_vec[chr], cands, expt_set->sampling == SCORE);
    }
}


/*-------------------------commands-------------------------*/
// write the fragments to stdout, or their reads to the files under --output
static void sim_core(const expt_param *expt_set, const mut_param *mut_set, const meth_param *meth_set)
{
    if (expt_set->fragments == 0 && expt_set->depth == 0) die("one of --fragments and --depth is required");
    input_set inputs;
    std::vector<chr_rec> &chr_vec = inputs.chr_vec;
    int tech = expt_set->tech_mode;
    bool whole_genome = tech == WGBS || tech == WGS;
    bool has_methdb = !meth_set->methdb_file.empty();
    contig_rec ct;

    // pass 1: what each contig can produce. Whole-genome assays measure it while the FASTA is
    // indexed; RRBS candidates and usable targets are kept for pass 2.
    read_inputs(&inputs, expt_set, mut_set, meth_set, true, whole_genome ? contig_fn([&](const ref_seq *ref, int chr) {
        wg_domain(ref, expt_set, (int)inputs.eff_vec.size(), &chr_vec[chr]);
    }) : contig_fn());
    std::vector<frag_domain> doms(chr_vec.size());
    if (!whole_genome) {
        methdb_reader *db = has_methdb && tech == RRBS ? open_methdb(meth_set->methdb_file.c_str(), chr_vec) : 0;
        for_each_contig(expt_set->reference, [&](const ref_seq *ref, int chr) {
            if (tech == RRBS) {     // only RRBS candidates depend on the variants
                std::vector<frag_rrbs_rec> cands;
                load_contig(&inputs, ref, chr, db, false, mut_set, &ct);
                rrbs_cands(&inputs, chr, &ct, expt_set, cands);
                rrbs_domain(cands, &doms[chr], &chr_vec[chr]);
                free_contig(&ct);
            } else {
                probe_domain(ref, expt_set, inputs.probe_vec[chr], &doms[chr], &chr_vec[chr]);
            }
        });
        if (db) close_methdb(db);
    }
    if (!inputs.eff_vec.empty()) calibrate_gc(inputs.eff_vec, chr_vec, expt_set->sd_insert > 0);
    cal_chr_count(expt_set, chr_vec);

    // pass 2: the fragments, contig by contig
    bool own_reads = !expt_set->output.empty();
    stream_rec st;
    reads_out ro;
    if (own_reads) reads_open(&ro, expt_set, chr_vec); else stream_open(&st, expt_set, chr_vec);
    methdb_reader *db = has_methdb ? open_methdb(meth_set->methdb_file.c_str(), chr_vec) : 0;
    out_t *db_out = 0;
    if (!meth_set->methdb_output.empty()) {
        db_out = out_open(meth_set->methdb_output, true, expt_set->threads);
        save_methdb_header(db_out, chr_vec);
    }
    profile_rec *pf = open_profile_set(&inputs, meth_set);
    uint64_t skipped = 0;
    int threads = expt_set->threads;
    uint64_t block = 512 * (uint64_t)threads;
    std::vector<frag_out> built(block);
    std::vector<read_pair> reads(own_reads ? block : 0);
    // alignments (BAM, SAM) need the reference positions
    bool details = expt_set->details || (own_reads && (expt_set->formats & (READ_BAM | READ_SAM)));
    std::vector<char> placed(block);
    std::vector<frag_seq> fs_vec(threads);
    std::vector<probe_rec> no_probes;

    for_each_contig(expt_set->reference, [&](const ref_seq *ref, int chr) {
        uint64_t count = chr_vec[chr].count;
        // the methylome is needed only where fragments fall (everywhere for a MethDB output)
        static const std::vector<std::pair<int, int>> nowhere;
        const std::vector<std::pair<int, int>> *reach = db_out ? 0 : !count ? &nowhere : whole_genome ? 0 : &doms[chr].reach;
        prepare_contig(&inputs, ref, chr, db, pf, expt_set, mut_set, meth_set, &ct, reach);
        if (db_out) save_methdb_chr(db_out, ref, ct.meth_vec, ct.snpmeth_map);

        if (count) {
            frag_domain &dom = doms[chr];
            std::vector<probe_rec> &probes = inputs.probe_vec.empty() ? no_probes : inputs.probe_vec[chr];
            site_index sites;
            if (expt_set->bisulfite) index_sites(ct.meth_vec, ref->len, &sites);
            const site_index *levels = expt_set->bisulfite ? &sites : 0;
            // built in parallel block by block, written in fragment order
            for (uint64_t first = 0; first < count; first += block) {
                uint64_t size = std::min(block, count - first);
                parallel_for(size, threads, [&](int t, uint64_t begin, uint64_t end) {
                    frag_rec frag;
                    for (uint64_t i = begin; i < end; ++i) {
                        placed[i] = gen_frag(chr, first + i, ref, ct.rseq, expt_set, inputs.eff_vec, probes, &dom, &frag, &fs_vec[t]);
                        if (!placed[i]) continue;
                        fill_frag_out(chr, &frag, &fs_vec[t], levels, ct.snpmeth_map, details, expt_set->site_kmers, &built[i]);
                        if (own_reads) cut_reads(expt_set, chr, first + i, ref, &built[i], &reads[i]);
                    }
                });
                for (uint64_t i = 0; i < size; ++i) {
                    if (!placed[i]) ++skipped;
                    else if (own_reads) reads_add(&ro, ref->name, chr, &built[i], &reads[i]);
                    else stream_add(&st, &built[i]);
                }
            }
        }
        doms[chr] = frag_domain();      // release its candidates
        free_contig(&ct);
    });
    if (db) close_methdb(db);
    if (pf) close_profile(pf);
    if (db_out) out_close(db_out);
    if (own_reads) reads_close(&ro, skipped); else stream_close(&st, skipped);
}

static const std::string &output_path(const expt_param *expt_set)
{
    if (expt_set->output.empty()) die("--output is required");
    return expt_set->output;
}

static void rrbs_catalog(const expt_param *expt_set, const mut_param *mut_set, const meth_param *meth_set)
{
    input_set inputs;
    read_inputs(&inputs, expt_set, mut_set, meth_set, false);
    out_t *out = out_open(output_path(expt_set), gzip_path(expt_set->output));
    save_rrcut_header(out);
    contig_rec ct;
    std::vector<frag_rrbs_rec> cands;
    for_each_contig(expt_set->reference, [&](const ref_seq *ref, int chr) {
        load_contig(&inputs, ref, chr, 0, false, mut_set, &ct);
        gen_cut_frag(ct.rseq, ct.rseq + 1, !ct.snpmeth_map.empty(), expt_set, inputs.cut_vec, cands);
        save_rrcut_bed(out, ref->name, cands);
        free_contig(&ct);
    });
    out_close(out);
}

static void methdb_build(const expt_param *expt_set, const mut_param *mut_set, const meth_param *meth_set)
{
    input_set inputs;
    read_inputs(&inputs, expt_set, mut_set, meth_set, false);
    out_t *out = out_open(output_path(expt_set), true, expt_set->threads);
    save_methdb_header(out, inputs.chr_vec);
    profile_rec *pf = open_profile_set(&inputs, meth_set);
    contig_rec ct;
    for_each_contig(expt_set->reference, [&](const ref_seq *ref, int chr) {
        prepare_contig(&inputs, ref, chr, 0, pf, expt_set, mut_set, meth_set, &ct);
        save_methdb_chr(out, ref, ct.meth_vec, ct.snpmeth_map);
        free_contig(&ct);
    });
    if (pf) close_profile(pf);
    out_close(out);
}

static void variant_catalog(const expt_param *expt_set, const mut_param *mut_set, const meth_param *meth_set)
{
    input_set inputs;
    read_inputs(&inputs, expt_set, mut_set, meth_set, false);
    methdb_reader *db = meth_set->methdb_file.empty() ? 0 : open_methdb(meth_set->methdb_file.c_str(), inputs.chr_vec);
    out_t *out = out_open(output_path(expt_set), gzip_path(expt_set->output));
    save_vcf_header(out);
    contig_rec ct;
    for_each_contig(expt_set->reference, [&](const ref_seq *ref, int chr) {
        load_contig(&inputs, ref, chr, db, false, mut_set, &ct);     // checks the variants
        save_vcf_chr(out, ref, ct.snpmeth_map);
        free_contig(&ct);
    });
    if (db) close_methdb(db);
    out_close(out);
}

// a JSON summary of the inputs; their digests are added by the Python package
static void validate_inputs(const expt_param *expt_set, const mut_param *mut_set, const meth_param *meth_set)
{
    if (!meth_set->methdb_file.empty() || !meth_set->methdb_output.empty()) die("input validation accepts text inputs, not MethDB snapshots");
    if (mut_set->mut_rate != 0) die("input validation does not generate de novo mutations");
    input_set inputs;
    read_inputs(&inputs, expt_set, mut_set, meth_set, false);
    profile_rec *pf = open_profile_set(&inputs, meth_set);
    contig_rec ct;
    std::vector<uint16_t> pools[4];
    uint64_t bases = 0, profile_contigs = 0, asm_contigs = 0;
    for_each_contig(expt_set->reference, [&](const ref_seq *ref, int chr) {
        bases += ref->len;
        load_contig(&inputs, ref, chr, 0, false, mut_set, &ct);
        if (!inputs.asm_vec.empty() && !inputs.asm_vec[chr].empty()) {     // ASM targets are checked against the methylome
            ++asm_contigs;
            if (pf) profile_contigs += pf->pending && pf->chr == chr;
            std::vector<std::pair<int, int>> nowhere;   // only the ASM targets
            chr_methylome(&inputs, ref, chr, &ct, pf, meth_set, expt_set->threads, &nowhere);
        } else if (pf) {
            profile_contigs += fill_profile_chr(pf, ref, chr, ct.meth_vec, meth_set, pools) > 0;
        }
        free_contig(&ct);
    });
    if (meth_set->pool_meth && pf->defined_rows == 0) die("methylation-profile pooling requires at least one defined probability");

    std::string json = "{\"status\":\"valid\",\"reference\":{\"contigs\":" + std::to_string(inputs.chr_vec.size()) + ",\"bases\":" + std::to_string(bases) + "}";
    json += ",\"vcf\":";
    if (!mut_set->vcf_file.empty()) {
        const vcf_summary &s = inputs.vcf_sum;
        json += "{\"rows\":" + std::to_string(s.rows) + ",\"contigs\":" + std::to_string(s.contigs)
              + ",\"retained\":" + std::to_string(s.retained) + ",\"reference_genotypes\":" + std::to_string(s.reference_genotypes)
              + ",\"skipped\":{\"total\":" + std::to_string(s.mnp + s.complex_replacement + s.long_indel)
              + ",\"mnp\":" + std::to_string(s.mnp) + ",\"complex_replacement\":" + std::to_string(s.complex_replacement)
              + ",\"long_indel\":" + std::to_string(s.long_indel) + "}}";
    } else {
        json += "null";
    }
    json += ",\"methylation\":";
    if (pf) {
        static const char *const formats[] = {"cgmap", "bedmethyl", "methbg", "methbed"};
        json += std::string("{\"format\":\"") + formats[meth_set->profile_format] + "\",\"rows\":" + std::to_string(pf->rows)
              + ",\"contigs\":" + std::to_string(profile_contigs) + ",\"defined_probabilities\":" + std::to_string(pf->defined_rows) + "}";
        close_profile(pf);
    } else {
        json += "null";
    }
    json += ",\"asm\":";
    if (!inputs.asm_vec.empty()) {
        json += std::string("{\"format\":\"") + (meth_set->asm_is_bed ? "asm-bed" : "ass") + "\",\"rows\":" + std::to_string(inputs.asm_rows)
              + ",\"contigs\":" + std::to_string(asm_contigs) + "}";
    } else {
        json += "null";
    }
    json += "}\n";
    fputs(json.c_str(), stdout);
}

static void methdb_export(const char *input, const char *fname)
{
    out_t *out = out_open(fname, gzip_path(fname));
    export_methdb(input, out);
    out_close(out);
}

// SAM on stdin to BAM on stdout, for the Python BAM writer
static int digit_argument(const char *text, int maximum, const char *what)
{
    char *end = 0;
    long value = strtol(text, &end, 10);
    if (*text == '\0' || *end != '\0' || value < 0 || value > maximum) die("--sam-to-bam %s must be an integer in [0, %d]", what, maximum);
    return (int)value;
}

static void sam_to_bam(int level, int threads)
{
    widen_pipe(fileno(stdin));
    samFile *in = sam_open("-", "r");
    if (in == 0) die("cannot open the SAM input stream");
    char mode[8];
    snprintf(mode, sizeof(mode), "wb%d", level);
    samFile *out = sam_open("-", mode);
    if (out == 0) die("cannot open the BAM output stream");
    if (threads > 0 && hts_set_threads(out, threads) < 0) die("cannot start BAM compression threads");
    sam_hdr_t *hdr = sam_hdr_read(in);
    if (hdr == 0) die("cannot parse the SAM header");
    if (sam_hdr_write(out, hdr) < 0) die("cannot write the BAM header");
    bam1_t *rec = bam_init1();
    int ret;
    while ((ret = sam_read1(in, hdr, rec)) >= 0) {
        if (sam_write1(out, hdr, rec) < 0) die("cannot write a BAM record");
    }
    if (ret < -1) die("cannot parse a SAM record");
    bam_destroy1(rec);
    sam_hdr_destroy(hdr);
    if (sam_close(out) != 0) die("cannot finish the BAM output stream");
    sam_close(in);
}

static int simu_usage(FILE *fp)
{
    const expt_param e;     // the defaults
    const mut_param m;
    const char *const flag[] = {"false", "true"};
    fprintf(fp, "\n");
    fprintf(fp, "htsim (high throughput reads simulator), the native core of BSReadSim\n");
    fprintf(fp, "Version: %s\n\n", PACKAGE_VERSION);
    fprintf(fp, "Usage: htsim [options] --output PREFIX      reads to files PREFIX.* (WGS, WES, TS; see --format)\n");
    fprintf(fp, "       htsim [options]                      the fragment stream to stdout, for bsreadsim\n");
    fprintf(fp, "       htsim rrbs-catalog [options] --output BED\n");
    fprintf(fp, "       htsim methdb-build [options] --output METHDB\n");
    fprintf(fp, "       htsim variant-catalog [options] --output VCF.gz\n");
    fprintf(fp, "       htsim validate-inputs [options]       print a JSON summary\n");
    fprintf(fp, "       htsim methdb-export METHDB OUTPUT\n");
    fprintf(fp, "       htsim --sam-to-bam LEVEL THREADS      convert SAM on stdin to BAM on stdout\n\n");
    fprintf(fp, "Options take a value (booleans: true or false); only --reference is required. [default]\n");
    fprintf(fp, "  --reference FASTA                  the genome\n");
    fprintf(fp, "  --technology NAME                  WGBS, RRBS, TBS, WGS, WES, or TS [%s]\n", TECH_NAMES[e.tech_mode]);
    fprintf(fp, "  --fragments N | --depth D          fragments, or mean depth (simulation runs)\n");
    fprintf(fp, "  --paired-end BOOL                  [%s]\n", flag[e.paired_end]);
    fprintf(fp, "  --read-length LEN[,LEN2]           read 1 (and read 2) length [%d]\n", e.read_length[0]);
    fprintf(fp, "  --insert-min/-mean/-max N          fragment length [%d/%d/%d]\n", e.min_insert, e.mean_insert, e.max_insert);
    fprintf(fp, "  --insert-sd SD                     [%g]\n", e.sd_insert);
    fprintf(fp, "  --max-ambiguous-fraction F         N bases a read may have [%g]\n", e.maxN_ratio);
    fprintf(fp, "  --mutation-rate R                  de novo variants per base [%g]\n", m.mut_rate);
    fprintf(fp, "  --indel-fraction F                 [%g]\n", m.indel_frac);
    fprintf(fp, "  --indel-extension-probability P    [%g]\n", m.indel_extn);
    fprintf(fp, "  --homozygous-only BOOL             [%s]\n", flag[m.homozygous_only]);
    fprintf(fp, "  --vcf VCF                          variants (instead of de novo ones)\n");
    fprintf(fp, "  --targets BED                      capture targets (TBS, WES, TS)\n");
    fprintf(fp, "  --center-sd SD                     fragment centre shift from a target [%g]\n", e.sd_center);
    fprintf(fp, "  --cut-sites SITE[,SITE]            RRBS cut sites [C|CGG]\n");
    fprintf(fp, "  --seed, --seed-mut, --seed-phase, --seed-meth N   [%llu]\n", (unsigned long long)e.seed);
    fprintf(fp, "  --threads N                        [%d]\n", e.threads);
    std::string formats;
    for (int k = 0; k < 4; ++k) if (e.formats & 1 << k) formats += (formats.empty() ? "" : ",") + std::string(READ_FORMATS[k]);
    fprintf(fp, "  --format F[,F]                     the read files: fastq.gz (PREFIX.R1.fastq.gz, PREFIX.R2.fastq.gz),\n");
    fprintf(fp, "                                     fastq, bam (PREFIX.bam, unsorted), sam [%s]\n", formats.c_str());
    fprintf(fp, "  --phred Q, --error-rate E          base quality and substitution rate of the reads [%d, %g]\n",
            e.phred, e.error_rate);
    fprintf(fp, "Bisulfite assays also take the methylation options of bsreadsim (--cgmap, --beta-cg, ...);\n");
    fprintf(fp, "see htsim/README.md. Outputs ending in .gz are BGZF-compressed.\n");
    return 0;
}

int main(int argc, char *argv[])
{
    hts_set_log_level(HTS_LOG_ERROR);
    if (argc == 2 && !strcmp(argv[1], "--help")) return simu_usage(stdout);
    if (argc == 2 && !strcmp(argv[1], "--version")) {printf("htsim %s\n", PACKAGE_VERSION); return 0;}
    if (argc >= 2 && !strcmp(argv[1], "--sam-to-bam")) {
        if (argc != 4) die("usage: htsim --sam-to-bam LEVEL THREADS");
        sam_to_bam(digit_argument(argv[2], 9, "LEVEL"), digit_argument(argv[3], 64, "THREADS"));
        return 0;
    }
    if (argc >= 2 && !strcmp(argv[1], "methdb-export")) {
        if (argc != 4) die("usage: htsim methdb-export METHDB OUTPUT");
        methdb_export(argv[2], argv[3]);
        return 0;
    }

    expt_param expt_set;
    mut_param  mut_set;
    meth_param meth_set;
    static const char *const commands[] = {"rrbs-catalog", "methdb-build", "variant-catalog", "validate-inputs"};
    void (*const runs[])(const expt_param *, const mut_param *, const meth_param *) = {rrbs_catalog, methdb_build, variant_catalog, validate_inputs};
    for (int i = 0; i < 4; ++i) {
        if (argc >= 2 && !strcmp(argv[1], commands[i])) {
            parse_options(argc - 1, argv + 1, commands[i], &expt_set, &mut_set, &meth_set);
            runs[i](&expt_set, &mut_set, &meth_set);
            return 0;
        }
    }
    parse_options(argc, argv, 0, &expt_set, &mut_set, &meth_set);
    sim_core(&expt_set, &mut_set, &meth_set);
    return 0;
}
