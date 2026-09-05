#include <vector>
#include <map>
#include <algorithm>
#include <string.h>
#include <limits.h>
#include <htslib/hts.h>
#include "struct.h"
#include "utility.h"
#include "methdb.h"
#include "haplo.h"


uint8_t ref_context(const ref_seq *ref, int i, bool collect_non_cpg)
{
    int l = (int)ref->len;
    int c = nst_nt4_table[(int)ref->seq[i]];
    if (!cg_table[c]) return 0;
    int step = c == 1 ? 1 : -1;
    int d1 = i + step >= 0 && i + step < l ? nst_nt4_table[(int)ref->seq[i + step]] : 4;
    int d2 = i + 2 * step >= 0 && i + 2 * step < l ? nst_nt4_table[(int)ref->seq[i + 2 * step]] : 4;
    uint8_t context = get_context(c, d1, d2);
    if (!collect_non_cpg && ctx_class(context) != 1) return 0;
    return context;
}


// the reference bases whose context a variant at i may change: those within 2 haplotype bases of it
static std::pair<int, int> variant_window(int i, const snpmeth_rec &snp, int l)
{
    return std::make_pair(std::max(0, i - 3), std::min(l, i + var_span(snp) + 3));
}

// create MethDB
void create_methdb(const ref_seq *ref, const std::map<int, snpmeth_rec>& snpmeth_map, std::vector<meth_rec>& meth_vec,
                   bool collect_non_cpg, const std::vector<std::pair<int, int>> *spans)
{
    int l = ref->len;
    std::vector<std::pair<int, int>> windows, whole(1, std::make_pair(0, l));
    for (std::map<int, snpmeth_rec>::const_iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it) {
        windows.push_back(variant_window(it->first, it->second, l));
    }
    windows = merge_intervals(windows);
    if (!spans) spans = &whole;
    uint64_t bases = 0;
    for (size_t s = 0; s < spans->size(); ++s) bases += (*spans)[s].second - (*spans)[s].first;
    std::vector<meth_rec>().swap(meth_vec);         // release the previous contig's records
    meth_vec.reserve((size_t)(bases * 0.48));       // chr19 has the highest GC ratio (0.4794); unused pages stay unmapped

    // one record per reference site, and per C/G that a variant may turn into one
    size_t w = 0;   // the first window that does not end before the base
    for (size_t s = 0; s < spans->size(); ++s) {
        for (int i = (*spans)[s].first; i < (*spans)[s].second; ++i) {
            if (!cg_table[nst_nt4_table[(int)ref->seq[i]]]) continue;
            meth_rec rec;
            rec.pos = i;
            rec.context[0] = rec.context[1] = ref_context(ref, i, collect_non_cpg);
            if (!rec.context[0]) {
                while (w < windows.size() && windows[w].second <= i) ++w;
                if (w == windows.size() || i < windows[w].first) continue;
            }
            meth_vec.push_back(rec);
        }
    }
}


// fill with distribution
uint16_t gen_meth(const meth_param *meth_set, int chr, uint64_t key, uint8_t context, const std::vector<uint16_t> (&pools)[4], uint8_t *type)
{
    // the level depends only on the seed, the contig, and the site
    Rng rng(seed_of(meth_set->seed_meth, chr, key, context));
    const std::vector<uint16_t> &pool = pools[ctx_class(context)];
    if (!pool.empty()) {*type = POOL; return pool[rng.below(pool.size())];}
    const param_rec &param = meth_set->params_map.at(context);
    *type = BETA;
    return prob_u16(rng.beta(param.alpha, param.beta));
}

void fill_levels(std::vector<meth_rec>& meth_vec, const meth_param *meth_set, int chr, const std::vector<uint16_t> (&pools)[4], int threads)
{
    parallel_for(meth_vec.size(), threads, [&](int, uint64_t begin, uint64_t end) {
        for (uint64_t i = begin; i < end; ++i) {
            meth_rec &rec = meth_vec[i];
            if (rec.type[0] || !rec.context[0]) continue;     // already filled, or no site
            uint8_t type;
            rec.meth[0] = rec.meth[1] = gen_meth(meth_set, chr, rec.pos, rec.context[0], pools, &type);
            rec.type[0] = rec.type[1] = type;
        }
    });
}

void update_variant(const ref_seq *ref, mutseq_t *hap1, mutseq_t *hap2, std::vector<meth_rec>& meth_vec,
                    std::map<int, snpmeth_rec>& snpmeth_map, const meth_param *meth_set, int chr, const std::vector<uint16_t> (&pools)[4])
{
    mutseq_t *rseq[2] = {hap1, hap2};
    std::map<int, uint64_t> ordinal;    // variants are numbered in position order
    for (std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it) {
        uint64_t next = ordinal.size();
        ordinal[it->first] = next;
    }

    frag_seq fs;
    std::vector<uint8_t> ctx;
    int l = (int)ref->len;
    uint8_t type;
    for (std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it) {
        // every base within 2 haplotype bases of the variant may have a new context
        std::pair<int, int> window = variant_window(it->first, it->second, l);
        int a = window.first, b = window.second;
        for (int h = 0; h < 2; ++h) {
            if ((h == 0 ? it->second.hap1 : it->second.hap2) != 1) continue;
            // a reference base that is deleted or substituted on this haplotype is no site
            for (int p = a; p < b; ++p) {
                int mut_type = rseq[h]->s[p] & mutmsk;
                meth_rec *rec = find_site(meth_vec, p);
                if (rec && (mut_type == DELETE || mut_type == SUBSTITUTE)) {
                    rec->context[h] = 0; rec->meth[h] = 0; rec->type[h] = 0;
                }
            }
            gen_frag_seq(rseq[h], a, b, 0, INT_MAX, &fs);
            frag_contexts(&fs, ctx, meth_set->collect_non_cpg);
            for (size_t k = 0; k < fs.seq.size(); ++k) {
                int p = fs.pos[k];
                if (fs.context[k] == MATCH) {
                    meth_rec *rec = find_site(meth_vec, p);
                    if (!rec || ctx[k] == ref_context(ref, p, meth_set->collect_non_cpg)) continue;
                    rec->context[h] = ctx[k];
                    rec->meth[h] = ctx[k] ? gen_meth(meth_set, chr, p, ctx[k], pools, &type) : 0;
                    rec->type[h] = ctx[k] ? type : 0;
                } else {
                    snpmeth_rec &snp = snpmeth_map[p];
                    int j = fs.offset[k];
                    uint64_t key = fs.context[k] == SNV ? (uint64_t)p : ((uint64_t)1 << 40 | ordinal[p] << 2 | j);
                    snp.context[h][j] = ctx[k];
                    snp.meth[h][j] = ctx[k] ? gen_meth(meth_set, chr, key, ctx[k], pools, &type) : 0;
                    snp.type[h][j] = ctx[k] ? type : 0;
                }
            }
        }
    }
}


// save/load/export MethDB
static const char *const kind_names[] = {"", "SNV", "insertion", "deletion"};
static const char *const source_names[] = {"", "vcf", "de-novo", "asm"};

static std::string ref_str(const snpmeth_rec &snp) {return snp.offset > 0 ? "." : packed_str(snp.ref, var_span(snp));}
static std::string alt_str(const snpmeth_rec &snp) {return snp.offset < 0 ? "." : packed_str(snp.alt, var_added(snp));}

void save_methdb_header(out_t *out, const std::vector<chr_rec>& chr_vec)
{
    out_printf(out, "#methdb\t1\n");
    for (size_t i = 0; i < chr_vec.size(); ++i) {
        out_printf(out, "#contig\t%s\t%u\t%s\n", chr_vec[i].name.c_str(), chr_vec[i].chr_len, chr_vec[i].md5.c_str());
    }
}

void save_methdb_chr(out_t *out, const ref_seq *ref, std::vector<meth_rec>& meth_vec, std::map<int, snpmeth_rec>& snpmeth_map)
{
    out_printf(out, ">\t%s\n", ref->name);
    for (std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it) {
        snpmeth_rec &snp = it->second;
        int n = var_added(snp);
        out_printf(out, "v\t%d\t%s\t%s\t%s\t%d\t%s\t%s", it->first, kind_names[var_kind(snp)], ref_str(snp).c_str(), alt_str(snp).c_str(),
                   (snp.hap1 == 1) | (snp.hap2 == 1) << 1, snp.id.c_str(), source_names[snp.source]);
        for (int h = 0; h < 2; ++h) {
            for (int field = 0; field < 3; ++field) {
                out_printf(out, "\t");
                if (n == 0) {out_printf(out, "."); continue;}
                for (int j = 0; j < n; ++j) {
                    int value = field == 0 ? snp.context[h][j] : field == 1 ? snp.meth[h][j] : snp.type[h][j];
                    out_printf(out, j ? ",%d" : "%d", value);
                }
            }
        }
        out_printf(out, "\n");
    }
    // sites, by their step from the previous site; almost all are shared rows
    int prev = 0;
    for (size_t i = 0; i < meth_vec.size(); ++i) {
        meth_rec &rec = meth_vec[i];
        if (!rec.context[0] && !rec.context[1]) continue;
        int ref_ctx = ref_context(ref, rec.pos, true);
        if (rec.context[0] == ref_ctx && rec.context[1] == ref_ctx && rec.meth[0] == rec.meth[1] && rec.type[0] == rec.type[1]) {
            out_printf(out, "s\t%d\t%d\t%d\t%d\n", rec.pos - prev, ref_ctx, rec.meth[0], rec.type[0]);
        } else {
            out_printf(out, "h\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", rec.pos - prev, ref_ctx,
                       rec.context[0], rec.meth[0], rec.type[0], rec.context[1], rec.meth[1], rec.type[1]);
        }
        prev = rec.pos;
    }
}

static bool next_methdb_line(methdb_reader *db) {return db->pending = in_next(db->in);}

static void parse_list(text_in *in, const char *text, int n, int *values, int maximum)
{
    if (n == 0) {if (strcmp(text, ".")) IN_DIE(in, "a deletion has no sites"); return;}
    std::string copy(text);
    std::vector<char*> parts;
    split_line(&copy[0], parts, ',');
    if ((int)parts.size() != n) IN_DIE(in, "a variant needs one value per base");
    for (int j = 0; j < n; ++j) {
        uint32_t value = in_u32(in, parts[j], "MethDB value");
        if ((int)value > maximum) IN_DIE(in, "value out of range");
        values[j] = (int)value;
    }
}

static uint16_t parse_bases(text_in *in, const char *text, int *n)
{
    std::vector<uint8_t> bases;
    if (strcmp(text, ".")) {
        for (const char *p = text; *p; ++p) {
            if (!strchr("ACGT", *p)) IN_DIE(in, "variant alleles must contain only A/C/G/T");
            bases.push_back(nst_nt4_table[(int)*p]);
        }
        if (bases.size() > (size_t)MAX_INDEL) IN_DIE(in, "variant alleles hold at most %d bases", MAX_INDEL);
    }
    *n = (int)bases.size();
    return pack_bases(bases);
}

// one 'v' row
static void parse_variant_row(text_in *in, int length, std::map<int, snpmeth_rec>& snpmeth_map, int *key)
{
    std::vector<char*> &f = in->f;
    if (in_split(in) != 14) IN_DIE(in, "variant rows have fourteen fields");
    snpmeth_rec snp;
    *key = (int)in_u32(in, f[1], "variant position");
    if (*key >= length) IN_DIE(in, "variant position is outside its contig");
    int kind = 0, n_ref, n_alt;
    for (int k = 1; k <= 3; ++k) if (!strcmp(f[2], kind_names[k])) kind = k;
    if (!kind) IN_DIE(in, "unknown variant kind");
    snp.ref = parse_bases(in, f[3], &n_ref);
    snp.alt = parse_bases(in, f[4], &n_alt);
    if ((kind == VAR_SNV && (n_ref != 1 || n_alt != 1)) || (kind == VAR_INS && (n_ref || !n_alt)) || (kind == VAR_DEL && (!n_ref || n_alt))) {
        IN_DIE(in, "variant alleles do not fit its kind");
    }
    snp.offset = kind == VAR_SNV ? 0 : kind == VAR_INS ? n_alt : -n_ref;
    uint32_t mask = in_u32(in, f[5], "variant haplotypes");
    if (mask < 1 || mask > 3) IN_DIE(in, "variant haplotypes must be 1, 2, or 3");
    snp.hap1 = mask & 1; snp.hap2 = (mask >> 1) & 1;
    snp.id = f[6];
    snp.source = 0;
    for (int k = 1; k <= 3; ++k) if (!strcmp(f[7], source_names[k])) snp.source = k;
    if (!snp.source) IN_DIE(in, "unknown variant source");
    int n = var_added(snp), values[MAX_INDEL];
    for (int h = 0; h < 2; ++h) {
        parse_list(in, f[8 + 3 * h], n, values, 15);    for (int j = 0; j < n; ++j) snp.context[h][j] = values[j];
        parse_list(in, f[9 + 3 * h], n, values, 65535); for (int j = 0; j < n; ++j) snp.meth[h][j] = values[j];
        parse_list(in, f[10 + 3 * h], n, values, 10);   for (int j = 0; j < n; ++j) snp.type[h][j] = values[j];
    }
    if (snpmeth_map.count(*key)) IN_DIE(in, "two variants share one position");
    snpmeth_map[*key] = snp;
}

// one 's' (shared) or 'h' (per haplotype) row; prev is the position of the previous site of the contig, -1 for none
static meth_rec parse_site_row(text_in *in, int length, int *prev, int *ref_ctx)
{
    std::vector<char*> &f = in->f;
    bool shared = in->line.s[0] == 's';
    if (in_split(in) != (shared ? 5u : 9u)) IN_DIE(in, shared ? "shared site rows have five fields" : "haplotype site rows have nine fields");
    meth_rec rec;
    uint32_t step = in_u32(in, f[1], "site step");
    if (step == 0 && *prev >= 0) IN_DIE(in, "site positions must increase");
    rec.pos = (*prev < 0 ? 0 : *prev) + (int)std::min<uint32_t>(step, length);
    if (rec.pos >= length) IN_DIE(in, "site position is outside its contig");
    *prev = rec.pos;
    *ref_ctx = (int)in_u32(in, f[2], "reference context");
    for (int h = 0; h < 2; ++h) {
        int k = shared ? 2 : 3 + 3 * h;
        uint32_t context = in_u32(in, f[k], "site context"), meth = in_u32(in, f[k + 1], "site level"), type = in_u32(in, f[k + 2], "site source");
        if (context > 15 || meth > 65535 || type > 10) IN_DIE(in, "site value out of range");
        if (context && !type) IN_DIE(in, "a site with a context must have a level source");
        rec.context[h] = context; rec.meth[h] = meth; rec.type[h] = type;
    }
    return rec;
}

// a MethDB opened at its first #contig line
static methdb_reader *open_reader(const char *fname)
{
    methdb_reader *db = new methdb_reader();
    db->in = in_open(fname, "MethDB");
    if (!next_methdb_line(db) || strcmp(db->in->line.s, "#methdb\t1")) die("%s is not a MethDB version 1 snapshot", fname);
    next_methdb_line(db);
    return db;
}

methdb_reader *open_methdb(const char *fname, const std::vector<chr_rec>& chr_vec)
{
    methdb_reader *db = open_reader(fname);
    size_t count = 0;
    std::vector<char*> &f = db->in->f;
    for (; db->pending && !strncmp(db->in->line.s, "#contig\t", 8); next_methdb_line(db), ++count) {
        in_split(db->in);
        if (f.size() != 4 || count >= chr_vec.size() || chr_vec[count].name != f[1]
            || std::to_string(chr_vec[count].chr_len) != f[2] || chr_vec[count].md5 != f[3]) {
            die("MethDB was built on a different reference (at contig %s)", f.size() > 1 ? f[1] : "?");
        }
    }
    if (count != chr_vec.size()) die("MethDB was built on a different reference (contig count)");
    return db;
}

void load_methdb_chr(methdb_reader *db, const ref_seq *ref, std::map<int, snpmeth_rec>& snpmeth_map,
                     std::vector<meth_rec> *meth_vec)
{
    text_in *in = db->in;
    snpmeth_map.clear();
    if (meth_vec) meth_vec->clear();
    std::string section = std::string(">\t") + ref->name;
    if (!db->pending || section != in->line.s) IN_DIE(in, "expected the section of contig %s", ref->name);
    int length = (int)ref->len, key, ref_ctx, prev = -1;
    while (next_methdb_line(db) && in->line.s[0] != '>') {
        char row = in->line.s[0];
        if (row == 'v') {
            parse_variant_row(in, length, snpmeth_map, &key);
        } else if (row == 's' || row == 'h') {
            if (!meth_vec) continue;
            meth_rec rec = parse_site_row(in, length, &prev, &ref_ctx);
            meth_vec->push_back(rec);
        } else {
            IN_DIE(in, "unknown row type");
        }
    }
}

void close_methdb(methdb_reader *db)
{
    in_close(db->in);
    delete db;
}


// methdb-bed-v1: BED6 rows for sites of reference positions, "#" lines for the rest
typedef struct {
    int pos, set, kind, ordinal, offset;    // kind 0: reference position, 1: inserted base
    uint8_t context, type;
    uint16_t meth;
    const char *allele;
} export_site;

static void export_contig(out_t *out, const std::string &name, std::map<int, snpmeth_rec>& snpmeth_map,
                          std::vector<std::pair<meth_rec, int>>& sites)
{
    static const char *const set_names[] = {"shared", "haplotype-1", "haplotype-2"};
    std::vector<export_site> rows;
    // reference bases, shared when both haplotypes have the same site
    for (size_t i = 0; i < sites.size(); ++i) {
        meth_rec &rec = sites[i].first;
        int ref_ctx = sites[i].second;
        if (rec.context[0] && rec.context[0] == rec.context[1] && rec.meth[0] == rec.meth[1] && rec.type[0] == rec.type[1]) {
            rows.push_back({rec.pos, 0, 0, 0, 0, rec.context[0], rec.type[0], rec.meth[0], "shared"});
            continue;
        }
        for (int h = 0; h < 2; ++h) {
            if (!rec.context[h]) continue;
            const char *allele = rec.type[h] == ASM_ALT ? "alternate" : rec.type[h] == ASM_REF ? "reference"
                               : rec.context[h] == ref_ctx ? "reference" : "alternate";
            rows.push_back({rec.pos, h + 1, 0, 0, 0, rec.context[h], rec.type[h], rec.meth[h], allele});
        }
    }
    // bases from variants
    int ordinal = 0;
    for (std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it, ++ordinal) {
        snpmeth_rec &snp = it->second;
        int kind = snp.offset > 0 ? 1 : 0;
        for (int j = 0; j < var_added(snp); ++j) {
            bool both = snp.hap1 == 1 && snp.hap2 == 1;
            if (both && snp.context[0][j] && snp.context[0][j] == snp.context[1][j]
                && snp.meth[0][j] == snp.meth[1][j] && snp.type[0][j] == snp.type[1][j]) {
                rows.push_back({it->first, 0, kind, ordinal, j, snp.context[0][j], snp.type[0][j], snp.meth[0][j], "shared"});
                continue;
            }
            for (int h = 0; h < 2; ++h) {
                if ((h == 0 ? snp.hap1 : snp.hap2) != 1 || !snp.context[h][j]) continue;
                rows.push_back({it->first, h + 1, kind, ordinal, j, snp.context[h][j], snp.type[h][j], snp.meth[h][j], "alternate"});
            }
        }
    }
    std::stable_sort(rows.begin(), rows.end(), [](const export_site &a, const export_site &b) {
        return a.pos != b.pos ? a.pos < b.pos : a.kind < b.kind;
    });
    for (size_t i = 0; i < rows.size(); ++i) {
        export_site &s = rows[i];
        const char *set = set_names[s.set];
        const char *source = s.type == INPUT ? "input" : (s.type == ASM_REF || s.type == ASM_ALT) ? "asm"
                           : s.type == BETA ? "beta" : s.type == POOL ? "pooled-input" : "unset";
        if (s.kind) {
            out_printf(out, "#insertion\t%s\t%s\ti%d.%d\t%d\t%d", name.c_str(), set, s.ordinal, s.offset, s.ordinal, s.offset);
        } else {
            out_printf(out, "%s\t%d\t%d\tmethdb:%s:%d\t%u\t%c\t%s\t%d\treference\t.\t.", name.c_str(), s.pos, s.pos + 1, set, s.pos,
                       (s.meth * 1000u + 32767u) / 65535u, s.context & 0x08 ? '-' : '+', set, s.pos);
        }
        out_printf(out, "\t%s\t%s\t%s\t%u\t%.9g\n", context_name(s.context), source, s.allele, s.meth, (double)(float)(s.meth / 65535.0));
    }
}

void export_methdb(const char *fname, out_t *out)
{
    methdb_reader *db = open_reader(fname);
    text_in *in = db->in;
    out_printf(out, "#format\tmethdb-bed-v1\n#columns\tchrom\tchromStart\tchromEnd\tname\tscore\tstrand\tset\torigin_id\t"
                    "origin_kind\tvariant_event\tinsertion_offset\tcontext\tsource\tallele\tprobability_u16\tprobability\n");

    std::map<std::string, int> lengths;
    std::string name;
    std::map<int, snpmeth_rec> snpmeth_map;
    std::vector<std::pair<meth_rec, int>> sites;
    int contig_count = 0, length = 0, key, ref_ctx, prev = -1;
    for (; db->pending; next_methdb_line(db)) {
        char *line = in->line.s;
        if (!strncmp(line, "#contig\t", 8)) {
            if (in_split(in) != 4) IN_DIE(in, "contig rows have four fields");
            lengths[in->f[1]] = (int)in_u32(in, in->f[2], "contig length");
            out_printf(out, "#contig\t%d\t%s\t%s\t%s\n", contig_count++, in->f[1], in->f[2], in->f[3]);
        } else if (line[0] == '>') {
            if (!name.empty()) export_contig(out, name, snpmeth_map, sites);
            if (in_split(in) != 2 || !lengths.count(in->f[1])) IN_DIE(in, "section of an undeclared contig");
            name = in->f[1];
            length = lengths[name];
            snpmeth_map.clear();
            sites.clear();
            prev = -1;
        } else if (line[0] == 'v') {
            parse_variant_row(in, length, snpmeth_map, &key);
            snpmeth_rec &snp = snpmeth_map[key];
            // reference interval of the event: an insertion is the empty interval after the base it follows
            int start = snp.offset > 0 ? key + 1 : key, end = snp.offset > 0 ? key + 1 : key + var_span(snp);
            out_printf(out, "#variant\t%s\t%d\t%d\t%d\t%s\t%s\t%s\t%d\t%s\t%s\n", name.c_str(), (int)snpmeth_map.size() - 1, start, end,
                       kind_names[var_kind(snp)], ref_str(snp).c_str(), alt_str(snp).c_str(),
                       (snp.hap1 == 1) | (snp.hap2 == 1) << 1, snp.id.c_str(), source_names[snp.source]);
        } else if (line[0] == 's' || line[0] == 'h') {
            meth_rec rec = parse_site_row(in, length, &prev, &ref_ctx);
            sites.push_back(std::make_pair(rec, ref_ctx));
        } else {
            IN_DIE(in, "unknown row type");
        }
    }
    if (!name.empty()) export_contig(out, name, snpmeth_map, sites);
    close_methdb(db);
}
