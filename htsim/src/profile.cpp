#include <vector>
#include <map>
#include <algorithm>
#include <string.h>
#include "struct.h"
#include "utility.h"
#include "profile.h"
#include "methdb.h"


// the C context on the strand a BED row names ('+' C, '-' G), 0 for '.'
static uint8_t strand_ctx(char strand) {return strand == '+' ? 0x01 : strand == '-' ? 0x09 : 0;}


// for methylation profiles
// [start, end) must be one base of the contig
static int one_base(text_in *in, const char *start_text, const char *end_text, uint32_t length)
{
    uint32_t start = in_u32(in, start_text, "chromStart"), end = in_u32(in, end_text, "chromEnd");
    if (start >= length || end != (uint64_t)start + 1) IN_DIE(in, "target must be one in-range half-open base");
    return (int)start;
}

static void check_counts(text_in *in, const char *methylated, const char *total)
{
    if (in_u32(in, methylated, "methylated count") > in_u32(in, total, "total count")) IN_DIE(in, "methylated count exceeds total count");
}

// One data row of each format, whose fields are f: sets pf->pos and the level (meth,
// defined), and the hints checked against the reference later: the context and the
// second base of the dinucleotide for CGmap, the strand and the context class for the
// BED formats.

// chr, C/G, 1-based position, context, dinucleotide, level|na, mC, total
static void cgmap_row(profile_rec *pf, std::vector<char*>& f, uint32_t length)
{
    text_in *in = pf->in;
    uint32_t pos = in_u32(in, f[2], "position");
    if (pos == 0 || pos > length) IN_DIE(in, "position is outside its reference contig");
    pf->pos = (int)pos - 1;
    if (!(pf->context = ctx_of_name(f[3]))) IN_DIE(in, "context must be CG, CHG, or CHH");
    if (!strcmp(f[1], "G")) pf->context |= 0x08; else if (strcmp(f[1], "C")) IN_DIE(in, "nucleotide must be C or G");
    if (strlen(f[4]) != 2 || f[4][0] != 'C' || !strchr("ACGT", f[4][1])) IN_DIE(in, "dinucleotide must be CA, CC, CG, or CT");
    pf->next_base = nst_nt4_table[(int)f[4][1]];
    if ((pf->next_base == 2) != (ctx_class(pf->context) == 1)) IN_DIE(in, "dinucleotide disagrees with its context class");
    pf->defined = strcmp(f[5], "na") != 0;
    if (pf->defined) pf->meth = prob_u16(in_prob(in, f[5], "methylation level"));
    check_counts(in, f[6], f[7]);
}

// BED9 (chr, start, end, name, score, strand, thickStart, thickEnd, itemRgb), valid coverage,
// percent modified[, modified, canonical, other, deleted, failed, diff, nocall counts]
static void bedmethyl_row(profile_rec *pf, std::vector<char*>& f, uint32_t length)
{
    text_in *in = pf->in;
    pf->pos = one_base(in, f[1], f[2], length);
    if (!*f[3]) IN_DIE(in, "name must not be empty");
    in_u32(in, f[4], "score");
    pf->strand = in_strand(in, f[5]);
    if (in_u32(in, f[6], "thickStart") != (uint32_t)pf->pos || in_u32(in, f[7], "thickEnd") != (uint32_t)pf->pos + 1) IN_DIE(in, "thickStart/thickEnd must match its target interval");
    if (strcmp(f[8], "0")) {
        std::vector<char*> rgb;
        split_line(f[8], rgb, ',');
        if (rgb.size() != 3) IN_DIE(in, "itemRgb must be 0 or an R,G,B triple");
        for (size_t i = 0; i < 3; ++i) if (in_u32(in, rgb[i], "itemRgb component") > 255) IN_DIE(in, "itemRgb components must be in [0, 255]");
    }
    uint32_t coverage = in_u32(in, f[9], "valid coverage");
    double percent = in_double(in, f[10], "percent modified");
    if (percent < 0 || percent > 100) IN_DIE(in, "percent modified must be finite and in [0, 100]");
    pf->meth = prob_u16(percent / 100);
    if (f.size() == 18) {
        uint64_t valid = 0;
        for (size_t i = 11; i < 18; ++i) {uint32_t count = in_u32(in, f[i], "extended count"); if (i < 14) valid += count;}
        if (valid != coverage) IN_DIE(in, "modified, canonical, and other counts must sum to valid coverage");
    }
}

// chr, start, end, level
static void methbg_row(profile_rec *pf, std::vector<char*>& f, uint32_t length)
{
    pf->pos = one_base(pf->in, f[1], f[2], length);
    pf->meth = prob_u16(in_prob(pf->in, f[3], "methylation level"));
}

// chr, start, end, name, score (level*1000), strand[, mC, total[, base[, context]]]
static void methbed_row(profile_rec *pf, std::vector<char*>& f, uint32_t length)
{
    text_in *in = pf->in;
    size_t n = f.size();
    pf->pos = one_base(in, f[1], f[2], length);
    if (!*f[3]) IN_DIE(in, "name must not be empty");
    uint32_t score = in_u32(in, f[4], "score");
    if (score > 1000) IN_DIE(in, "score must be in [0, 1000]");
    pf->meth = prob_u16(score / 1000.0);
    pf->strand = in_strand(in, f[5]);
    if (n >= 8 && (!strcmp(f[6], ".")) != (!strcmp(f[7], "."))) IN_DIE(in, "methylated and total counts must both be set or both be .");
    if (n >= 8 && strcmp(f[6], ".")) check_counts(in, f[6], f[7]);
    if (n >= 9 && strcmp(f[8], ".")) {
        char base = !strcmp(f[8], "C") ? '+' : !strcmp(f[8], "G") ? '-' : 0;
        if (!base) IN_DIE(in, "base must be C, G, or .");
        if (pf->strand != '.' && pf->strand != base) IN_DIE(in, "strand and base disagree");
        pf->strand = base;
    }
    if (n == 10 && strcmp(f[9], ".")) {
        if (!(pf->cclass = ctx_class(ctx_of_name(f[9])))) IN_DIE(in, "context must be CG, CHG, CHH, or .");
    }
}

// the formats, indexed by profile_t
static const struct {
    const char *label;          /* names the input in errors */
    bool bed;                   /* BED track and browser lines are skipped */
    size_t widths[4];           /* accepted field counts (0: none) */
    const char *width_text;
    void (*parse_row)(profile_rec *pf, std::vector<char*>& f, uint32_t length);
} formats[] = {
    {"CGmap",     false, {8},           "exactly eight",      cgmap_row},
    {"bedMethyl", true,  {11, 18},      "eleven or eighteen", bedmethyl_row},
    {"MethBG",    true,  {4},           "exactly four",       methbg_row},
    {"MethBED",   true,  {6, 8, 9, 10}, "6, 8, 9, or 10",     methbed_row},
};

// read ahead the next data row; false at the end of the profile
static bool read_profile_row(profile_rec *pf)
{
    const auto &format = formats[pf->format];
    text_in *in = pf->in;
    int prev_chr = pf->pending ? pf->chr : -1, prev_pos = pf->pending ? pf->pos : -1;
    while (in_next(in)) {
        char *line = in->line.s;
        if (format.bed ? bed_header(line) : (!*line || line[0] == '#')) continue;
        size_t n = in_split(in);
        if (std::find(format.widths, format.widths + 4, n) == format.widths + 4) {
            IN_DIE(in, "row must contain %s tab-separated fields", format.width_text);
        }
        if (pf->ncol && pf->ncol != n) IN_DIE(in, "rows must use one consistent field count");
        pf->ncol = n;
        pf->chr = in_chr(in, pf->chr_idx, in->f[0]);
        pf->strand = '.';
        pf->cclass = pf->next_base = 0;
        pf->context = 0;
        pf->defined = true;
        pf->meth = 0;
        format.parse_row(pf, in->f, pf->chr_len[pf->chr]);
        if (prev_chr >= 0 && (pf->chr < prev_chr || (pf->chr == prev_chr && pf->pos <= prev_pos))) {
            IN_DIE(in, "rows must follow FASTA order with unique positions");
        }
        ++pf->rows;
        pf->defined_rows += pf->defined;
        pf->pending = true;
        return true;
    }
    pf->pending = false;
    return false;
}

profile_rec *open_profile(const char *fname, int format, const std::vector<chr_rec>& chr_vec)
{
    profile_rec *pf = new profile_rec();
    pf->in = in_open(fname, formats[format].label);
    pf->format = format;
    pf->chr_idx = chr_index(chr_vec);
    for (size_t i = 0; i < chr_vec.size(); ++i) pf->chr_len.push_back(chr_vec[i].chr_len);
    pf->pending = false;
    pf->rows = pf->defined_rows = 0;
    pf->ncol = 0;
    if (!read_profile_row(pf)) die("%s contains no data rows", formats[format].label);
    return pf;
}

void close_profile(profile_rec *pf)
{
    in_close(pf->in);
    delete pf;
}

int fill_profile_chr(profile_rec *pf, const ref_seq *ref, int chr, std::vector<meth_rec>& meth_vec,
                     const meth_param *meth_set, std::vector<uint16_t> (&pools)[4])
{
    text_in *in = pf->in;
    int num_site = 0;
    std::vector<meth_rec>::iterator site = meth_vec.begin();     // rows come in position order
    while (pf->pending && pf->chr <= chr) {
        if (pf->chr == chr) {
            // check the row against the reference
            uint8_t observed = ref_context(ref, pf->pos, true);
            if (pf->format == CGMAP) {
                if (observed != pf->context) IN_DIE(in, "nucleotide/context disagrees with the reference");
                int next = observed & 0x08 ? 3 - nst_nt4_table[(int)ref->seq[pf->pos - 1]] : nst_nt4_table[(int)ref->seq[pf->pos + 1]];
                if (next != pf->next_base) IN_DIE(in, "dinucleotide disagrees with the reference");
            } else {
                if (!observed) IN_DIE(in, "target is not a resolved reference C/G context");
                if (pf->strand != '.' && strand_ctx(pf->strand) != (observed & 0x09)) IN_DIE(in, "strand/base disagrees with the reference C/G base");
                if (pf->cclass && pf->cclass != ctx_class(observed)) IN_DIE(in, "context disagrees with the reference context");
                pf->context = observed;
            }
            if (pf->defined) {
                if (meth_set->pool_meth) {
                    pools[ctx_class(pf->context)].push_back(pf->meth);
                } else if (meth_vec.size()) {
                    while (site != meth_vec.end() && site->pos < pf->pos) ++site;
                    if (site != meth_vec.end() && site->pos == pf->pos && site->context[0] == pf->context) {   // not when --cpg-only dropped it
                        site->meth[0] = site->meth[1] = pf->meth;
                        site->type[0] = site->type[1] = INPUT;
                    }
                }
            }
            ++num_site;
        }
        read_profile_row(pf);
    }
    return num_site;
}


// for ASM
static uint8_t asm_base(text_in *in, const char *text, const char *field)
{
    if (strlen(text) != 1 || !strchr("ACGT", text[0])) IN_DIE(in, "%s must be one uppercase A/C/G/T base", field);
    return nst_nt4_table[(int)text[0]];
}

static void check_support(text_in *in, const char *text, const char *field)
{
    std::string copy(text);
    std::vector<char*> parts;
    split_line(&copy[0], parts, '-');
    if (parts.size() != 2) IN_DIE(in, "%s must be methylated-unmethylated counts", field);
    if ((uint64_t)in_u32(in, parts[0], field) + in_u32(in, parts[1], field) == 0) IN_DIE(in, "%s must contain at least one supporting read", field);
}

void parse_asm(const char *fname, bool is_bed, const std::vector<chr_rec>& chr_vec,
               std::vector<std::vector<asm_rec>>& asm_vec, uint64_t *rows)
{
    const char *ass_header = "Chr\tSNP_Pos\tRef\tAllele1\tAllele2\tC_Pos\tAllele1_linked_C\tAllele2_linked_C\t"
                             "Allele1_linked_C_met\tAllele2_linked_C_met\tpvalue\tfdr\tASM";
    std::map<std::string, int> chr_idx = chr_index(chr_vec);
    asm_vec.assign(chr_vec.size(), std::vector<asm_rec>());
    text_in *in = in_open(fname, is_bed ? "ASM BED" : "CGmapTools ASS");
    uint64_t n_rows = 0;
    size_t width = 0;
    int prev_chr = -1;
    while (in_next(in)) {
        char *line = in->line.s;
        if (is_bed ? bed_header(line) : (!*line || line[0] == '#' || !strcmp(line, ass_header))) continue;
        size_t n = in_split(in);
        std::vector<char*> &f = in->f;
        if (is_bed ? (n != 12 && n != 16) : n != 13) IN_DIE(in, "row must contain exactly %s tab-separated fields", is_bed ? "twelve or sixteen" : "thirteen");
        if (width && width != n) IN_DIE(in, "rows must use one consistent width");
        width = n;
        int chr = in_chr(in, chr_idx, f[0]);
        if (chr < prev_chr) IN_DIE(in, "rows must follow FASTA contig order");
        prev_chr = chr;
        uint32_t length = chr_vec[chr].chr_len;
        asm_rec tmp_asm = {};
        if (is_bed) {   // target interval, name, score, strand, SNV interval, REF, ALT, ref level, alt level[, evidence]
            uint32_t ts = in_u32(in, f[1], "chromStart"), te = in_u32(in, f[2], "chromEnd");
            uint32_t vs = in_u32(in, f[6], "linkedStart"), ve = in_u32(in, f[7], "linkedEnd");
            if (ts >= length || te != (uint64_t)ts + 1 || vs >= length || ve != (uint64_t)vs + 1) IN_DIE(in, "target and linked SNV must be in-range one-base intervals");
            if (!*f[3]) IN_DIE(in, "name must not be empty");
            if (in_u32(in, f[4], "score") > 1000) IN_DIE(in, "score must be in [0, 1000]");
            tmp_asm.context = strand_ctx(in_strand(in, f[5]));
            tmp_asm.target = (int)ts; tmp_asm.snv = (int)vs;
            tmp_asm.ref = asm_base(in, f[8], "linked REF");
            tmp_asm.alt = asm_base(in, f[9], "linked ALT");
            if (tmp_asm.ref == tmp_asm.alt) IN_DIE(in, "linked REF and ALT bases must differ");
            tmp_asm.ref_meth = prob_u16(in_prob(in, f[10], "reference methylation level"));
            tmp_asm.alt_meth = prob_u16(in_prob(in, f[11], "alternate methylation level"));
            if (n == 16) {
                check_support(in, f[12], "REF_SUPPORT");
                check_support(in, f[13], "ALT_SUPPORT");
                in_prob(in, f[14], "P_VALUE");
                in_prob(in, f[15], "Q_VALUE");
            }
        } else {        // CGmapTools: chr, SNP_Pos, Ref, Allele1, Allele2, C_Pos, linked C x2, levels x2, pvalue, fdr, ASM
            uint32_t snv = in_u32(in, f[1], "linked variant position"), target = in_u32(in, f[5], "target position");
            if (snv == 0 || target == 0 || snv > length || target > length) IN_DIE(in, "position is outside its reference contig");
            tmp_asm.snv = (int)snv - 1; tmp_asm.target = (int)target - 1;
            tmp_asm.ref = asm_base(in, f[2], "Ref");
            uint8_t allele1 = asm_base(in, f[3], "Allele1"), allele2 = asm_base(in, f[4], "Allele2");
            if (allele1 == allele2 || (allele1 != tmp_asm.ref && allele2 != tmp_asm.ref)) IN_DIE(in, "alleles must differ and exactly one must equal Ref");
            check_support(in, f[6], "Allele1_linked_C");
            check_support(in, f[7], "Allele2_linked_C");
            uint16_t meth1 = prob_u16(in_prob(in, f[8], "Allele1_linked_C_met")), meth2 = prob_u16(in_prob(in, f[9], "Allele2_linked_C_met"));
            in_prob(in, f[10], "pvalue");
            in_prob(in, f[11], "fdr");
            if (strcmp(f[12], "TRUE") && strcmp(f[12], "FALSE")) IN_DIE(in, "ASM field must be TRUE or FALSE");
            if (!strcmp(f[12], "FALSE")) continue;
            bool first_is_ref = allele1 == tmp_asm.ref;
            tmp_asm.alt = first_is_ref ? allele2 : allele1;
            tmp_asm.ref_meth = first_is_ref ? meth1 : meth2;
            tmp_asm.alt_meth = first_is_ref ? meth2 : meth1;
        }
        asm_vec[chr].push_back(tmp_asm);
        ++n_rows;
    }
    in_close(in);
    if (n_rows == 0) die(is_bed ? "ASM BED contains no data rows" : "CGmapTools ASS contains no rows called TRUE");
    if (rows) *rows = n_rows;
}

void check_asm_chr(const ref_seq *ref, std::vector<asm_rec>& asm_list)
{
    for (size_t i = 0; i < asm_list.size(); ++i) {
        asm_rec &row = asm_list[i];
        uint8_t observed = ref_context(ref, row.target, true);
        if (!observed) die("ASM target is not a resolved reference C/G context (%s:%d)", ref->name, row.target + 1);
        if (row.context && (row.context & 0x08) != (observed & 0x08)) die("ASM BED strand disagrees with the reference C/G base (%s:%d)", ref->name, row.target + 1);
        row.context = observed;
        if (nst_nt4_table[(int)ref->seq[row.snv]] != row.ref) die("ASM linked REF base disagrees with the reference (%s:%d)", ref->name, row.snv + 1);
    }
    std::sort(asm_list.begin(), asm_list.end(), [](const asm_rec &a, const asm_rec &b) {
        return a.target != b.target ? a.target < b.target : a.snv < b.snv;
    });
    for (size_t i = 1; i < asm_list.size(); ++i) {
        if (asm_list[i].target == asm_list[i - 1].target) die("ASM rows must name each target once (%s:%d)", ref->name, asm_list[i].target + 1);
    }
}

void asm_variants(int chr_idx, std::vector<asm_rec>& asm_list, uint64_t seed_phase, std::map<int, snpmeth_rec>& snpmeth_map)
{
    snpmeth_map.clear();
    std::vector<asm_rec> links(asm_list);
    std::sort(links.begin(), links.end(), [](const asm_rec &a, const asm_rec &b) {
        return a.snv != b.snv ? a.snv < b.snv : (a.ref != b.ref ? a.ref < b.ref : a.alt < b.alt);
    });
    for (size_t i = 0; i < links.size(); ++i) {
        asm_rec &row = links[i];
        std::map<int, snpmeth_rec>::iterator found = snpmeth_map.find(row.snv);
        if (found != snpmeth_map.end()) {
            if (found->second.ref != row.ref || found->second.alt != row.alt) die("ASM rows define conflicting alleles for one linked SNV");
            continue;
        }
        snpmeth_rec tmp_snp;
        tmp_snp.ref = row.ref; tmp_snp.alt = row.alt; tmp_snp.offset = 0;
        if (Rng(seed_of(seed_phase, chr_idx, snpmeth_map.size())).unif() < 0.5) tmp_snp.hap1 = 1; else tmp_snp.hap2 = 1;
        tmp_snp.source = FROM_ASM;
        tmp_snp.id = "asm_" + std::to_string(chr_idx) + "_" + std::to_string(row.snv + 1);
        snpmeth_map[row.snv] = tmp_snp;
    }
}

void fill_asm_chr(const ref_seq *ref, std::vector<asm_rec>& asm_list, std::vector<meth_rec>& meth_vec,
                  std::map<int, snpmeth_rec>& snpmeth_map, bool collect_non_cpg)
{
    for (size_t i = 0; i < asm_list.size(); ++i) {
        asm_rec &row = asm_list[i];
        std::map<int, snpmeth_rec>::iterator found = snpmeth_map.find(row.snv);
        if (found == snpmeth_map.end() || found->second.offset != 0 || found->second.ref != row.ref || found->second.alt != row.alt) {
            die("ASM row does not resolve to its exact typed linked SNV (%s:%d)", ref->name, row.snv + 1);
        }
        snpmeth_rec &snp = found->second;
        if (snp.hap1 == snp.hap2) die("ASM linked SNV must retain one reference haplotype (%s:%d)", ref->name, row.snv + 1);
        meth_rec *rec = find_site(meth_vec, row.target);
        bool collected = collect_non_cpg || ctx_class(row.context) == 1;
        if (!collected || !rec || rec->context[0] != row.context || rec->context[1] != row.context) {
            die("ASM target is not one shared reference-equivalent diploid site (%s:%d)", ref->name, row.target + 1);
        }
        int alt_h = snp.hap1 == 1 ? 0 : 1;
        rec->meth[alt_h] = row.alt_meth;     rec->type[alt_h] = ASM_ALT;
        rec->meth[1 - alt_h] = row.ref_meth; rec->type[1 - alt_h] = ASM_REF;
    }
}
