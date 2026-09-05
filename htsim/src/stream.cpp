#include <stdio.h>
#include <string>
#include <vector>
#include <map>
#include "struct.h"
#include "utility.h"
#include "haplo.h"
#include "stream.h"




void index_sites(const std::vector<meth_rec>& meth_vec, int len, site_index *sites)
{
    sites->meth_vec = &meth_vec;
    sites->block.resize(((size_t)len >> SITE_BLOCK_SHIFT) + 1);
    size_t i = 0;
    for (size_t b = 0; b < sites->block.size(); ++b) {
        while (i < meth_vec.size() && meth_vec[i].pos < (int)(b << SITE_BLOCK_SHIFT)) ++i;
        sites->block[b] = (uint32_t)i;
    }
}

void fill_frag_out(int chr, frag_rec *frag, frag_seq *fs, const site_index *sites,
                   std::map<int, snpmeth_rec>& snpmeth_map, bool details, bool kmers, frag_out *out)
{
    int h = frag->haplo, n = (int)fs->seq.size();
    out->contig = chr;
    out->haplo  = h;
    out->strand = frag->strand;
    // the reference envelope: from the first to the last base with a reference position
    // (inserted bases take the position of the base they follow)
    int first = 0;
    while (first < n - 1 && fs->context[first] == INSR) ++first;
    out->start  = fs->pos[first] + (fs->context[first] == INSR);
    out->end    = fs->pos[n - 1] + 1;
    out->bases.assign(fs->seq.begin(), fs->seq.end());
    out->ref_pos.clear(); out->ev_start.clear(); out->ev_end.clear(); out->ev_kind.clear(); out->ev_deleted.clear();
    out->site_pos.clear(); out->site_level.clear(); out->site_ctx.clear(); out->site_kmer.clear();

    // sites: the level of the haplotype's base, from the methdb or from its variant
    if (sites) {
        // reference positions only increase along the template: start at the block of its first base, then walk
        const std::vector<meth_rec> &meth_vec = *sites->meth_vec;
        size_t site = sites->block[fs->pos[0] >> SITE_BLOCK_SHIFT];
        for (int k = 0; k < n; ++k) {
            if (!cg_table[fs->seq[k]]) continue;
            uint8_t context, type;
            uint16_t meth;
            if (fs->context[k] == MATCH) {
                while (site < meth_vec.size() && meth_vec[site].pos < fs->pos[k]) ++site;
                if (site == meth_vec.size() || meth_vec[site].pos != fs->pos[k]) continue;
                const meth_rec &rec = meth_vec[site];
                context = rec.context[h]; meth = rec.meth[h]; type = rec.type[h];
            } else {
                const snpmeth_rec &snp = snpmeth_map.find(fs->pos[k])->second;
                int j = fs->offset[k];
                context = snp.context[h][j]; meth = snp.meth[h][j]; type = snp.type[h][j];
            }
            if (!context) continue;
            out->site_pos.push_back(k);
            out->site_level.push_back(meth);
            out->site_ctx.push_back(ctx_class(context) | (type == ASM_REF || type == ASM_ALT ? 0x80 : 0));
            if (kmers) out->site_kmer.push_back(frag_kmer(fs, k));
        }
    }
    if (details) {
        for (int k = 0; k < n; ++k) out->ref_pos.push_back(fs->context[k] == INSR ? 0xffffffffu : (uint32_t)fs->pos[k]);
        for (size_t e = 0; e < fs->event_start.size(); ++e) {
            out->ev_start.push_back(fs->event_start[e]);
            out->ev_end.push_back(fs->event_end[e]);
            out->ev_kind.push_back(fs->event_kind[e]);
            out->ev_deleted.push_back(fs->event_deleted[e]);
        }
    }
}


static void json_string(std::string &out, const std::string &text)
{
    out += '"';
    for (size_t i = 0; i < text.size(); ++i) {
        unsigned char c = text[i];
        if (c == '"' || c == '\\') {out += '\\'; out += c;}
        else if (c < 0x20) {char escaped[8]; snprintf(escaped, sizeof(escaped), "\\u%04x", c); out += escaped;}
        else out += c;
    }
    out += '"';
}

static void write_stdout(const void *data, size_t size)
{
    if (size && fwrite(data, 1, size, stdout) != size) die("cannot write the fragment stream");
}

template <typename T>
static void append(std::string &payload, const T &value) {payload.append((const char *)&value, sizeof(value));}

template <typename T>
static void append_column(std::string &payload, std::vector<T> &values)
{
    append(payload, (uint64_t)(values.size() * sizeof(T)));
    payload.append((const char *)values.data(), values.size() * sizeof(T));
    values.clear();
}

template <typename T>
static void extend(std::vector<T> &to, const std::vector<T> &from) {to.insert(to.end(), from.begin(), from.end());}

void stream_open(stream_rec *st, const expt_param *expt_set, const std::vector<chr_rec>& chr_vec)
{
    widen_pipe(fileno(stdout));
    st->batch_size = expt_set->batch_size;
    st->paired_end = expt_set->paired_end;
    st->per_contig.assign(chr_vec.size(), 0);
    st->tmpl_offsets.assign(1, 0); st->site_offsets.assign(1, 0); st->ev_offsets.assign(1, 0);

    std::string line = "{\"version\": ";  json_string(line, PACKAGE_VERSION);
    line += ", \"master_seed\": " + std::to_string(expt_set->seed);
    line += ", \"technology\": ";      json_string(line, TECH_NAMES[expt_set->tech_mode]);
    line += std::string(", \"details\": ") + (expt_set->details ? "true" : "false");
    line += std::string(", \"paired_end\": ") + (expt_set->paired_end ? "true" : "false");
    line += ", \"read_lengths\": [" + std::to_string(expt_set->read_length[0]);
    if (expt_set->paired_end) line += ", " + std::to_string(expt_set->read_length[1]);
    line += "]";
    line += ", \"contigs\": [";
    for (size_t i = 0; i < chr_vec.size(); ++i) {
        line += i ? ", [" : "[";
        json_string(line, chr_vec[i].name);
        line += ", " + std::to_string(chr_vec[i].chr_len) + ", ";
        json_string(line, chr_vec[i].md5);
        line += "]";
    }
    line += "]}\n";
    write_stdout(line.data(), line.size());
}

static void stream_flush(stream_rec *st)
{
    if (st->n == 0) return;
    std::string payload;
    append(payload, st->first);
    append_column(payload, st->contig);
    append_column(payload, st->start);
    append_column(payload, st->end);
    append_column(payload, st->haplo);
    append_column(payload, st->strand);
    append_column(payload, st->tmpl_offsets);
    append_column(payload, st->bases);
    append_column(payload, st->site_offsets);
    append_column(payload, st->site_pos);
    append_column(payload, st->site_level);
    append_column(payload, st->site_ctx);
    append_column(payload, st->site_kmer);
    append_column(payload, st->ev_offsets);
    append_column(payload, st->ev_start);
    append_column(payload, st->ev_end);
    append_column(payload, st->ev_kind);
    append_column(payload, st->ev_deleted);
    uint64_t size = payload.size();
    write_stdout(&size, sizeof(size));
    write_stdout(payload.data(), payload.size());
    st->first += st->n;
    st->n = 0;
    st->tmpl_offsets.assign(1, 0); st->site_offsets.assign(1, 0); st->ev_offsets.assign(1, 0);
}

void stream_add(stream_rec *st, frag_out *f)
{
    ++st->n;
    ++st->n_frag;
    ++st->per_contig[f->contig];
    st->n_mate += st->paired_end ? 2 : 1;
    st->n_base += f->bases.size();
    st->n_site += f->site_pos.size();
    st->contig.push_back(f->contig);
    st->start.push_back(f->start);
    st->end.push_back(f->end);
    st->haplo.push_back(f->haplo);
    st->strand.push_back(f->strand);
    extend(st->bases, f->bases);
    st->tmpl_offsets.push_back((uint32_t)st->bases.size());
    extend(st->site_pos, f->site_pos);
    extend(st->site_level, f->site_level);
    extend(st->site_ctx, f->site_ctx);
    extend(st->site_kmer, f->site_kmer);
    st->site_offsets.push_back((uint32_t)st->site_pos.size());
    extend(st->ev_start, f->ev_start);
    extend(st->ev_end, f->ev_end);
    extend(st->ev_kind, f->ev_kind);
    extend(st->ev_deleted, f->ev_deleted);
    st->ev_offsets.push_back((uint32_t)st->ev_start.size());
    if (st->n >= (uint64_t)st->batch_size) stream_flush(st);
}

void stream_close(stream_rec *st, uint64_t skipped)
{
    stream_flush(st);
    uint64_t end = 0;
    write_stdout(&end, sizeof(end));
    std::string line = "{\"fragment_count\": " + std::to_string(st->n_frag);
    line += ", \"mate_count\": " + std::to_string(st->n_mate);
    line += ", \"template_base_count\": " + std::to_string(st->n_base);
    line += ", \"methylation_site_count\": " + std::to_string(st->n_site);
    line += ", \"skipped_fragment_count\": " + std::to_string(skipped);
    line += ", \"per_contig_fragment_counts\": [";
    for (size_t i = 0; i < st->per_contig.size(); ++i) line += (i ? ", " : "") + std::to_string(st->per_contig[i]);
    line += "]}\n";
    write_stdout(line.data(), line.size());
    if (fflush(stdout) != 0) die("cannot write the fragment stream");
}
