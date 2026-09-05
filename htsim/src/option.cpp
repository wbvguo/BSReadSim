#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include "struct.h"
#include "utility.h"
#include "option.h"


/*-------------------------options-------------------------*/
enum {
    OPT_SEED = 256, OPT_SEED_MUT, OPT_SEED_PHASE, OPT_SEED_METH,
    OPT_REFERENCE, OPT_VCF, OPT_CGMAP, OPT_BED_METHYL, OPT_METHBG, OPT_METHBED, OPT_METHDB,
    OPT_ASM, OPT_ASM_BED, OPT_TECH, OPT_DIRECTIONAL, OPT_PAIRED, OPT_READ_LEN,
    OPT_INSERT_MIN, OPT_INSERT_MEAN, OPT_INSERT_MAX, OPT_INSERT_SD, OPT_DEPTH, OPT_FRAGMENTS,
    OPT_MAX_N, OPT_MUT_RATE, OPT_INDEL_FRAC, OPT_INDEL_EXTN, OPT_HOMOZYGOUS, OPT_CPG_ONLY,
    OPT_POOL, OPT_BETA_CG, OPT_BETA_CHG, OPT_BETA_CHH, OPT_SAMPLING, OPT_GC_PROFILE,
    OPT_CUT_SITES, OPT_RRBS_CANDS, OPT_TARGETS, OPT_CENTER_SD, OPT_DETAILS, OPT_BATCH,
    OPT_THREADS, OPT_METHDB_OUT, OPT_OUTPUT, OPT_FORMAT, OPT_PHRED, OPT_ERROR_RATE, OPT_SITE_KMERS, OPT_END
};

static const struct option long_options[] = {
    {"seed", required_argument, 0, OPT_SEED},
    {"seed-mut", required_argument, 0, OPT_SEED_MUT},
    {"seed-phase", required_argument, 0, OPT_SEED_PHASE},
    {"seed-meth", required_argument, 0, OPT_SEED_METH},
    {"reference", required_argument, 0, OPT_REFERENCE},
    {"vcf", required_argument, 0, OPT_VCF},
    {"cgmap", required_argument, 0, OPT_CGMAP},
    {"bed-methyl", required_argument, 0, OPT_BED_METHYL},
    {"methbg", required_argument, 0, OPT_METHBG},
    {"methbed", required_argument, 0, OPT_METHBED},
    {"methdb", required_argument, 0, OPT_METHDB},
    {"asm", required_argument, 0, OPT_ASM},
    {"asm-bed", required_argument, 0, OPT_ASM_BED},
    {"technology", required_argument, 0, OPT_TECH},
    {"directional", required_argument, 0, OPT_DIRECTIONAL},
    {"paired-end", required_argument, 0, OPT_PAIRED},
    {"read-length", required_argument, 0, OPT_READ_LEN},
    {"insert-min", required_argument, 0, OPT_INSERT_MIN},
    {"insert-mean", required_argument, 0, OPT_INSERT_MEAN},
    {"insert-max", required_argument, 0, OPT_INSERT_MAX},
    {"insert-sd", required_argument, 0, OPT_INSERT_SD},
    {"depth", required_argument, 0, OPT_DEPTH},
    {"fragments", required_argument, 0, OPT_FRAGMENTS},
    {"max-ambiguous-fraction", required_argument, 0, OPT_MAX_N},
    {"mutation-rate", required_argument, 0, OPT_MUT_RATE},
    {"indel-fraction", required_argument, 0, OPT_INDEL_FRAC},
    {"indel-extension-probability", required_argument, 0, OPT_INDEL_EXTN},
    {"homozygous-only", required_argument, 0, OPT_HOMOZYGOUS},
    {"cpg-only", required_argument, 0, OPT_CPG_ONLY},
    {"pool-meth", required_argument, 0, OPT_POOL},
    {"beta-cg", required_argument, 0, OPT_BETA_CG},
    {"beta-chg", required_argument, 0, OPT_BETA_CHG},
    {"beta-chh", required_argument, 0, OPT_BETA_CHH},
    {"sampling", required_argument, 0, OPT_SAMPLING},
    {"gc-profile", required_argument, 0, OPT_GC_PROFILE},
    {"cut-sites", required_argument, 0, OPT_CUT_SITES},
    {"rrbs-candidates", required_argument, 0, OPT_RRBS_CANDS},
    {"targets", required_argument, 0, OPT_TARGETS},
    {"center-sd", required_argument, 0, OPT_CENTER_SD},
    {"details", required_argument, 0, OPT_DETAILS},
    {"batch-size", required_argument, 0, OPT_BATCH},
    {"threads", required_argument, 0, OPT_THREADS},
    {"methdb-output", required_argument, 0, OPT_METHDB_OUT},
    {"output", required_argument, 0, OPT_OUTPUT},
    {"format", required_argument, 0, OPT_FORMAT},
    {"phred", required_argument, 0, OPT_PHRED},
    {"error-rate", required_argument, 0, OPT_ERROR_RATE},
    {"site-kmers", required_argument, 0, OPT_SITE_KMERS},
    {0, 0, 0, 0}
};

static const char *opt_name(int c) {return long_options[c - OPT_SEED].name;}

static bool to_bool(const char *text, int c)
{
    if (strcmp(text, "true") && strcmp(text, "false")) die("--%s must be true or false", opt_name(c));
    return !strcmp(text, "true");
}

static int to_int(const char *text, int c)
{
    std::string what = std::string("--") + opt_name(c);
    uint32_t value = to_u32(text, what.c_str());
    if (value > INT32_MAX) die("%s is too large", what.c_str());
    return (int)value;
}

static void to_shape(const char *text, int c, param_rec *param)
{
    std::string copy(text), what = std::string("--") + opt_name(c);
    std::vector<char*> parts;
    split_line(&copy[0], parts, ',');
    if (parts.size() != 2) die("%s must be ALPHA,BETA", what.c_str());
    param->alpha = (float)to_double(parts[0], what.c_str());
    param->beta  = (float)to_double(parts[1], what.c_str());
    if (!(param->alpha > 0 && param->beta > 0)) die("%s must be positive", what.c_str());
}

// LENGTH (both mates), or READ1,READ2; true when two were given
static bool to_lengths(const char *text, int c, int *lengths)
{
    std::string copy(text), what = std::string("--") + opt_name(c);
    std::vector<char*> parts;
    split_line(&copy[0], parts, ',');
    if (parts.empty() || parts.size() > 2) die("%s must be LENGTH or READ1,READ2", what.c_str());
    for (size_t m = 0; m < 2; ++m) lengths[m] = to_int(parts[std::min(m, parts.size() - 1)], c);
    return parts.size() == 2;
}

static int to_choice(const char *text, int c, const char *const *names, int n)
{
    for (int i = 0; i < n; ++i) if (!strcmp(text, names[i])) return i;
    std::string allowed;
    for (int i = 0; i < n; ++i) allowed += (i ? ", " : "") + std::string(names[i]);
    die("--%s must be one of %s", opt_name(c), allowed.c_str());
}

// FORMAT[,FORMAT]: read files (read_format_t bits)
static int to_formats(const char *text, int c)
{
    std::string copy(text);
    std::vector<char*> parts;
    split_line(&copy[0], parts, ',');
    int formats = 0;
    for (size_t i = 0; i < parts.size(); ++i) {
        int format = 1 << to_choice(parts[i], c, READ_FORMATS, 4);
        if (formats & format) die("--%s repeats %s", opt_name(c), parts[i]);
        formats |= format;
    }
    if (formats == 0) die("--%s must name at least one format", opt_name(c));
    return formats;
}

static bool is_prob(double p) {return p >= 0 && p <= 1;}

static const int MAX_READ_LENGTH = 10000;  /* the longest read, as documented */

// parse "--name value" pairs, then check their combinations
void parse_options(int argc, char *argv[], const char *command,
                   expt_param *expt_set, mut_param *mut_set, meth_param *meth_set)
{
    static const char *const samplings[] = {"uniform", "gc", "score"};
    std::set<int> seen;
    bool has_depth = false, has_frag = false, pool = false, two_lengths = false;
    int n_profile = 0;
    param_rec beta[3] = {{0.5f, 0.5f}, {0.01f, 0.05f}, {0.01f, 0.05f}};   /* CG, CHG, CHH */
    std::string cut_str;

    int c;
    opterr = 0;
    optind = 1;
    while ((c = getopt_long(argc, argv, ":", long_options, 0)) >= 0) {
        if (c == '?') die("unknown option: %s", argv[optind - 1]);
        if (c == ':') die("missing value for %s", argv[optind - 1]);
        if (!seen.insert(c).second) die("duplicate option: --%s", opt_name(c));
        if (!*optarg) die("--%s must not be empty", opt_name(c));
        std::string what = std::string("--") + opt_name(c);
        switch (c) {
            case OPT_SEED:        expt_set->seed = to_u64(optarg, what.c_str()); break;
            case OPT_SEED_MUT:    mut_set->seed_snp = to_u64(optarg, what.c_str()); break;
            case OPT_SEED_PHASE:  mut_set->seed_phase = to_u64(optarg, what.c_str()); break;
            case OPT_SEED_METH:   meth_set->seed_meth = to_u64(optarg, what.c_str()); break;
            case OPT_REFERENCE:   expt_set->reference = optarg; break;
            case OPT_VCF:         mut_set->vcf_file = optarg; break;
            case OPT_CGMAP:       meth_set->profile_file = optarg; meth_set->profile_format = CGMAP; ++n_profile; break;
            case OPT_BED_METHYL:  meth_set->profile_file = optarg; meth_set->profile_format = BEDMETHYL; ++n_profile; break;
            case OPT_METHBG:      meth_set->profile_file = optarg; meth_set->profile_format = METHBG; ++n_profile; break;
            case OPT_METHBED:     meth_set->profile_file = optarg; meth_set->profile_format = METHBED; ++n_profile; break;
            case OPT_METHDB:      meth_set->methdb_file = optarg; break;
            case OPT_ASM:         meth_set->asm_file = optarg; break;
            case OPT_ASM_BED:     meth_set->asm_file = optarg; meth_set->asm_is_bed = true; break;
            case OPT_TECH:        expt_set->tech_mode = to_choice(optarg, c, TECH_NAMES, 6); break;
            case OPT_DIRECTIONAL: expt_set->directional = to_bool(optarg, c); break;
            case OPT_PAIRED:      expt_set->paired_end = to_bool(optarg, c); break;
            case OPT_READ_LEN:    two_lengths = to_lengths(optarg, c, expt_set->read_length); break;
            case OPT_INSERT_MIN:  expt_set->min_insert = to_int(optarg, c); break;
            case OPT_INSERT_MEAN: expt_set->mean_insert = to_int(optarg, c); break;
            case OPT_INSERT_MAX:  expt_set->max_insert = to_int(optarg, c); break;
            case OPT_INSERT_SD:   expt_set->sd_insert = to_double(optarg, what.c_str()); break;
            case OPT_DEPTH:       expt_set->depth = to_double(optarg, what.c_str()); has_depth = true; break;
            case OPT_FRAGMENTS:   expt_set->fragments = to_u32(optarg, what.c_str()); has_frag = true; break;
            case OPT_MAX_N:       expt_set->maxN_ratio = to_double(optarg, what.c_str()); break;
            case OPT_MUT_RATE:    mut_set->mut_rate = to_double(optarg, what.c_str()); break;
            case OPT_INDEL_FRAC:  mut_set->indel_frac = to_double(optarg, what.c_str()); break;
            case OPT_INDEL_EXTN:  mut_set->indel_extn = to_double(optarg, what.c_str()); break;
            case OPT_HOMOZYGOUS:  mut_set->homozygous_only = to_bool(optarg, c); break;
            case OPT_CPG_ONLY:    meth_set->collect_non_cpg = !to_bool(optarg, c); break;
            case OPT_POOL:        pool = to_bool(optarg, c); break;
            case OPT_BETA_CG:     to_shape(optarg, c, &beta[0]); break;
            case OPT_BETA_CHG:    to_shape(optarg, c, &beta[1]); break;
            case OPT_BETA_CHH:    to_shape(optarg, c, &beta[2]); break;
            case OPT_SAMPLING:    expt_set->sampling = to_choice(optarg, c, samplings, 3); break;
            case OPT_GC_PROFILE:  expt_set->bias_file = optarg; break;
            case OPT_CUT_SITES:   cut_str = optarg; break;
            case OPT_RRBS_CANDS:  expt_set->rrbs_file = optarg; break;
            case OPT_TARGETS:     expt_set->bed_file = optarg; break;
            case OPT_CENTER_SD:   expt_set->sd_center = to_double(optarg, what.c_str()); break;
            case OPT_DETAILS:     expt_set->details = to_bool(optarg, c); break;
            case OPT_BATCH:       expt_set->batch_size = to_int(optarg, c); break;
            case OPT_THREADS:     expt_set->threads = to_int(optarg, c); break;
            case OPT_METHDB_OUT:  meth_set->methdb_output = optarg; break;
            case OPT_OUTPUT:      expt_set->output = optarg; break;
            case OPT_FORMAT:      expt_set->formats = to_formats(optarg, c); break;
            case OPT_PHRED:       expt_set->phred = to_int(optarg, c); break;
            case OPT_ERROR_RATE:  expt_set->error_rate = to_double(optarg, what.c_str()); break;
            case OPT_SITE_KMERS:  expt_set->site_kmers = to_bool(optarg, c); break;
        }
    }
    if (optind < argc) die("unexpected argument: %s", argv[optind]);
    int tech = expt_set->tech_mode;
    expt_set->bisulfite = tech == WGBS || tech == RRBS || tech == TBS;
    meth_set->pool_meth = pool;
    for (int k = 0; k < 2; ++k) {   // C (k = 0) and G contexts
        meth_set->params_map[0x01 | k << 3] = beta[0];
        meth_set->params_map[0x03 | k << 3] = beta[1];
        meth_set->params_map[0x07 | k << 3] = beta[2];
    }
    if (!cut_str.empty()) {
        std::vector<char*> sites;
        split_line(&cut_str[0], sites, ',');
        for (size_t i = 0; i < sites.size(); ++i) expt_set->cut_sites.push_back(sites[i]);
    }

    // combinations
    #define REQUIRE(cond, ...) do {if (!(cond)) die(__VA_ARGS__);} while (0)
    bool whole_genome = tech == WGBS || tech == WGS, targeted = tech == TBS || tech == WES || tech == TS;
    bool has_profile = n_profile > 0, has_asm = !meth_set->asm_file.empty(), has_vcf = !mut_set->vcf_file.empty();
    bool has_methdb = !meth_set->methdb_file.empty(), has_methdb_out = !meth_set->methdb_output.empty();
    REQUIRE(!expt_set->reference.empty(), "--reference is required");
    REQUIRE(n_profile <= 1, "methylation profile inputs are mutually exclusive");
    for (int mate = 0; mate < 2; ++mate) {
        REQUIRE(expt_set->read_length[mate] >= 1 && expt_set->read_length[mate] <= MAX_READ_LENGTH,
                "--read-length must be in [1, %d]", MAX_READ_LENGTH);
    }
    REQUIRE(expt_set->min_insert > 0 && expt_set->batch_size > 0 && expt_set->threads > 0,
            "insert sizes, batch size, and threads must be positive");
    REQUIRE(expt_set->min_insert <= expt_set->mean_insert && expt_set->mean_insert <= expt_set->max_insert,
            "--insert-min <= --insert-mean <= --insert-max must hold");
    REQUIRE(expt_set->sd_insert >= 0, "--insert-sd must be non-negative");
    bool fixed = expt_set->sd_insert == 0 && tech != RRBS;
    REQUIRE(!two_lengths || expt_set->paired_end, "--read-length takes one length for single-end reads");
    REQUIRE(std::max(expt_set->read_length[0], expt_set->read_length[1]) <= (fixed ? expt_set->mean_insert : expt_set->min_insert),
            "--read-length must not exceed %s", fixed ? "--insert-mean when --insert-sd is 0" : "--insert-min");
    REQUIRE(!(has_depth && has_frag), "--depth and --fragments are mutually exclusive");
    REQUIRE(!has_depth || expt_set->depth > 0, "--depth must be positive");
    REQUIRE(!has_frag || expt_set->fragments > 0, "the fragment count must be positive");
    REQUIRE(is_prob(expt_set->maxN_ratio), "--max-ambiguous-fraction must be in [0, 1]");
    REQUIRE(is_prob(mut_set->mut_rate) && is_prob(mut_set->indel_frac) && is_prob(mut_set->indel_extn),
            "--mutation-rate, --indel-fraction, and --indel-extension-probability must be in [0, 1]");
    REQUIRE(!(seen.count(OPT_ASM) && seen.count(OPT_ASM_BED)), "--asm and --asm-bed are mutually exclusive");
    REQUIRE(!(has_vcf && mut_set->mut_rate > 0), "--vcf and --mutation-rate are mutually exclusive");
    REQUIRE(!(has_asm && mut_set->mut_rate > 0), "--asm/--asm-bed and --mutation-rate are mutually exclusive");
    if (has_methdb) {
        REQUIRE(!has_vcf, "--methdb embeds variants and cannot be combined with --vcf");
        REQUIRE(mut_set->mut_rate == 0, "--methdb embeds variants and forbids de novo mutations");
        REQUIRE(!has_profile && !has_asm && !pool,
                "--methdb embeds a methylation profile and cannot be combined with another profile, ASM, or --pool-meth");
        REQUIRE(!has_methdb_out, "--methdb and a MethDB output are mutually exclusive");
    }
    REQUIRE(!pool || has_profile, "--pool-meth requires a text methylation profile");
    if (!expt_set->bisulfite) {
        REQUIRE(!has_profile && !has_methdb && !has_methdb_out && !has_asm && !pool, "standard sequencing forbids methylation inputs");
    }
    std::set<std::string> unique_sites;
    for (size_t i = 0; i < expt_set->cut_sites.size(); ++i) {
        const std::string &site = expt_set->cut_sites[i];
        size_t bar = site.find('|');
        REQUIRE(bar != std::string::npos && site.find('|', bar + 1) == std::string::npos && bar + 1 < site.size(),
                "--cut-site must contain exactly one | before at least one base: %s", site.c_str());
        for (size_t k = 0; k < site.size(); ++k) {
            REQUIRE(k == bar || strchr("ACGTN", site[k]), "--cut-site must use A, C, G, T, N: %s", site.c_str());
        }
        REQUIRE(unique_sites.insert(site).second, "duplicate --cut-site: %s", site.c_str());
    }
    if (tech == RRBS) {
        if (expt_set->cut_sites.empty()) expt_set->cut_sites.push_back("C|CGG");   /* MspI */
    } else {
        REQUIRE(expt_set->cut_sites.empty() && expt_set->rrbs_file.empty(), "only RRBS accepts --cut-site and --rrbs-candidates");
    }
    if (targeted) {
        REQUIRE(!expt_set->bed_file.empty(), "targeted assays require --targets");
        REQUIRE(expt_set->sd_center >= 0, "--center-sd must be non-negative");
    } else {
        REQUIRE(expt_set->bed_file.empty(), "only TBS, WES, and TS accept --targets");
    }
    REQUIRE(expt_set->sampling != GC_PROFILE || whole_genome, "--sampling gc requires WGBS or WGS");
    REQUIRE(expt_set->sampling != SCORE || !whole_genome, "--sampling score requires RRBS, TBS, WES, or TS");
    REQUIRE((expt_set->sampling == GC_PROFILE) == !expt_set->bias_file.empty(), "--sampling gc and --gc-profile require each other");
    REQUIRE(expt_set->sampling != SCORE || tech != RRBS || !expt_set->rrbs_file.empty(), "--sampling score requires --rrbs-candidates");
    if (!command && !expt_set->output.empty()) {     // a run that writes reads
        REQUIRE(!expt_set->bisulfite, "--output writes reads of WGS, WES, and TS; bisulfite reads come from the bsreadsim package");
    } else if (seen.count(OPT_FORMAT) || seen.count(OPT_PHRED) || seen.count(OPT_ERROR_RATE)) {
        if (command) die("%s does not write reads (--format, --phred, --error-rate)", command);
        die("--format, --phred, and --error-rate need --output PREFIX");
    }
    REQUIRE(expt_set->phred >= 0 && expt_set->phred <= 93, "--phred must be in [0, 93]");
    REQUIRE(is_prob(expt_set->error_rate), "--error-rate must be in [0, 1]");
    #undef REQUIRE
}
