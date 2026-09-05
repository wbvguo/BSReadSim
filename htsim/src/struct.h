#ifndef STRUCT_H
#define STRUCT_H

// Constants, lookup tables, and the structs shared by the modules, each with the
// small functions that read or encode its fields. General tools are in utility.h.

#include <stdint.h>
#include <string>
#include <vector>
#include <map>


/*-------------------------variable-------------------------*/
extern const char PACKAGE_VERSION[];

extern const uint8_t nst_nt4_table[256];
extern const uint8_t cg_table[5];
extern const uint8_t cg_context_table[64];

// a haplotype base (mut_t): the base in the low 4 bits; a substitution or a deletion in the top
// 4 bits; or the number of bases inserted after it (1-4) in the top 4 bits and them in bits 4-11
enum muttype_t {NOCHANGE = 0, SUBSTITUTE = 0xe000, DELETE = 0xf000};
typedef unsigned short mut_t;
extern const mut_t mutmsk;

// flags of a template base (upper half of its context byte)
static const uint8_t MATCH = 0x00, SNV = 0x10, INSR = 0x30;

// technologies (WG: WGBS/WGS; RR: RRBS; T: TBS/WES/TS), sampling, library strands
enum tech_t {WGBS = 0, RRBS, TBS, WGS, WES, TS};
extern const char *const TECH_NAMES[6];  /* "WGBS", ..., "TS" */
enum sampling_t {UNIFORM = 0, GC_PROFILE, SCORE};
enum strand_t {OT = 0, OB, CTOT, CTOB};

// the read files a simulation run writes itself under --output PREFIX (see reads.h)
enum read_format_t {READ_FASTQ_GZ = 1, READ_FASTQ = 2, READ_BAM = 4, READ_SAM = 8};
extern const char *const READ_FORMATS[4];   /* "fastq.gz", "fastq", "bam", "sam": format 1 << k */

// source of a methylation level (meth_rec.type)
enum meth_type_t {UNSET = 0, INPUT = 2, ASM_REF = 4, ASM_ALT = 5, BETA = 8, POOL = 10};

// variant kinds (also the event kinds of the fragment stream) and sources
enum var_kind_t {VAR_SNV = 1, VAR_INS = 2, VAR_DEL = 3};
enum var_source_t {FROM_VCF = 1, FROM_DENOVO = 2, FROM_ASM = 3};

static const int MAX_INDEL = 4;         /* longest insertion/deletion */
static const uint16_t NO_KMER = 0xffff; /* 7-mer index of a k-mer with N or past an end */


/*-------------------------contexts-------------------------*/
// context codes: 1 CG, 3 CHG, 7 CHH for a C; plus 8 (9, 11, 15) for a G; 0 for no site
// context of a C (or G) given its next (previous) two bases
uint8_t get_context(int c, int d1, int d2);
const char *context_name(uint8_t context);                  /* "CG-C", ..., "CHH-G" */
inline int ctx_class(uint8_t context) {return (context & 7) == 1 ? 1 : (context & 7) == 3 ? 2 : (context & 7) == 7 ? 3 : 0;}  /* 1 CG, 2 CHG, 3 CHH */
uint8_t ctx_of_name(const char *name);                      /* "CG", "CHG", "CHH" -> C context, 0 otherwise */


/*-------------------------struct-------------------------*/
// one variant of a contig, keyed by position in snpmeth_map: the substituted base,
// the first deleted base, or the base after which bases are inserted
typedef struct {
    uint16_t meth[2][MAX_INDEL] = {};   /* level of the variant's bases on hap1/hap2 */
    uint8_t context[2][MAX_INDEL] = {}; /* their context, 0 for no site */
    uint8_t type[2][MAX_INDEL] = {};    /* their meth_type_t */
    uint16_t ref = 0, alt = 0;          /* bases, packed (see pack_bases) */
    int8_t hap1 = 0, hap2 = 0;          /* 1 if the haplotype carries it */
    int8_t offset = 0;                  /* +inserted, -deleted, 0 for SNV */
    uint8_t source = FROM_VCF;
    std::string id;
} snpmeth_rec;

inline int var_kind(const snpmeth_rec &snp) {return snp.offset == 0 ? VAR_SNV : snp.offset > 0 ? VAR_INS : VAR_DEL;}
// the bases a variant adds to a haplotype (1 for an SNV), and the reference bases it covers
inline int var_added(const snpmeth_rec &snp) {return snp.offset == 0 ? 1 : snp.offset > 0 ? snp.offset : 0;}
inline int var_span(const snpmeth_rec &snp) {return snp.offset < 0 ? -snp.offset : 1;}
// up to MAX_INDEL bases packed 2 bits each, the first in the lowest bits
uint16_t pack_bases(const std::vector<uint8_t>& bases);
inline int packed_base(uint16_t packed, int k) {return (packed >> (2 * k)) & 3;}
std::string packed_str(uint16_t packed, int n);


// parse RRBS
typedef struct {
    int idx = -1;               /* the cut: before base idx of seq */
    std::vector<uint8_t> seq;   /* sequence encoded by numbers, N matches any base */
} cut_rec;

// an RRBS candidate: the haplotype bases between two cuts
typedef struct {
    int haplo = 3;              /* haplotype 1: hap1, 2: hap2, 3: both (no variants) */
    int pos_l = 0, skip = 0;    /* its first base: base skip of reference position pos_l */
    int env_l = 0, env_r = 0;   /* reference envelope: first to last reference base */
    int len = 0, gc = 0;        /* template length and its C/G bases */
    int n_cuts = 0;             /* motif recognitions at the cuts it spans */
    double score = 1;
} frag_rrbs_rec;


// a capture target (TBS/WES/TS): reference interval [pos_l, pos_r)
typedef struct {
    int pos_l = 0, pos_r = 0;
    double score   = 0;
    char strand    = '.';       /* capture strand: '+', '-', or '.' */
} probe_rec;

// a fragment: len haplotype bases from base skip of reference position pos_l (see gen_frag_seq)
typedef struct {
    int pos_l = 0, skip = 0, len = 0;
    int8_t haplo   = 0;         /* haplotype 0: hap1, 1: hap2 */
    int8_t strand  = 0;         /* strand_t */
} frag_rec;


// one fragment read from its haplotype
typedef struct {
    std::vector<uint8_t> seq;       /* template bases */
    std::vector<int> pos;           /* reference position (the anchor for inserted bases) */
    std::vector<uint8_t> context;   /* MATCH/SNV/INSR | site context */
    std::vector<uint8_t> offset;    /* position within inserted bases */
    std::vector<int> event_start, event_end, event_kind;  /* variant intervals on the template */
    std::vector<int> event_deleted; /* reference bases each event deletes (deletions only) */
    uint8_t before[3], after[3];    /* haplotype bases flanking the template (4 if none) */
} frag_seq;
// the N bases among base codes [from, to)
int count_n(const std::vector<uint8_t>& bases, int from, int to);


// chr counts
typedef struct {
    std::string name, md5;
    uint32_t chr_len = 0;
    uint64_t eff_len = 0;       /* bases that the requested depth refers to */
    uint64_t count   = 0;       /* fragments to generate */
    double   score   = 0;       /* share of the fragments */
    std::vector<uint64_t> gc_starts;  /* WG GC sampling: starts per GC bin */
} chr_rec;

// contig name -> index into chr_vec
std::map<std::string, int> chr_index(const std::vector<chr_rec>& chr_vec);

// the bases of the reference contig being processed, as in the FASTA (either case)
typedef struct {
    const char *name;
    const char *seq;
    int len;
} ref_seq;
// the N bases (other letters included) of reference interval [from, to)
int count_n(const ref_seq *ref, int from, int to);


// methdb
typedef struct {
    int pos = -1;
    uint16_t meth[2] = {0, 0};  /* level on hap1/hap2 (see prob_u16) */
    uint8_t context[2] = {0,0}; /*1,3,7;9,11,15 for the context on hap1/hap2, 0 for no site*/
    uint8_t type[2] = {0, 0};   /* meth_type_t */
} meth_rec;                     /*each struct take 12 bytes*/

inline uint16_t prob_u16(double p) {return (uint16_t)(p * 65535.0 + 0.5);}    /* round(p * 65535) */
// the record of reference position pos in meth_vec (one per site, in position order), or 0 if none
meth_rec *find_site(std::vector<meth_rec>& meth_vec, int pos);

typedef struct {
    float alpha = -1;
    float beta = -1;
} param_rec;                    /*take 8 bytes*/


// haplotypes
typedef struct {
    int l;    /* length (reference bases) */
    mut_t *s; /* sequence */
} mutseq_t;

// a haplotype word: the bases inserted after the base, and inserted base j
inline int mut_ins(int c) {int t = c & mutmsk; return t == NOCHANGE || t == SUBSTITUTE || t == DELETE ? 0 : t >> 12;}
inline int mut_ins_base(int c, int j) {return (c >> (4 + 2 * j)) & 3;}


// hold all parameters
typedef struct {
    std::string reference;
    uint64_t seed       = 0;
    int tech_mode       = WGBS;
    int sampling        = UNIFORM;
    bool bisulfite      = true;
    bool directional    = true;
    bool paired_end     = true;

    double maxN_ratio   = 0.05;
    double depth        = 0;        /* used when fragments is 0 */
    uint64_t fragments  = 0;

    int read_length[2]  = {100, 100};   /* read 1, read 2 */
    int mean_insert     = 400;
    double sd_insert    = 25;
    int min_insert      = 100;
    int max_insert      = 1000;
    double sd_center    = 50;

    std::string bias_file, bed_file, rrbs_file;     /* GC profile, targets, RRBS candidates */
    std::vector<std::string> cut_sites;

    bool details        = false;
    bool site_kmers     = false;    /* stream the 7-mer of every site (for a sequence methylation model) */
    int batch_size      = 1024;
    int threads         = 1;
    std::string output;             /* where a command writes; a simulation run writes its reads to
                                       files PREFIX.* here instead of the fragment stream */
    int formats         = READ_FASTQ_GZ;    /* those files (read_format_t bits) */
    int phred           = 40;       /* their base quality */
    double error_rate   = 0.005;    /* and substitution rate */
} expt_param;

// the most N bases read 1 (mate 0) or read 2 (mate 1) may have
inline int max_n(const expt_param *expt_set, int mate) {return (int)(expt_set->maxN_ratio * expt_set->read_length[mate]);}

typedef struct {
    double mut_rate     = 0;
    double indel_frac   = 0.15;
    double indel_extn   = 0.15;
    uint64_t seed_snp   = 0;
    uint64_t seed_phase = 0;
    bool homozygous_only = false;   /* de novo variants are all homozygous */
    std::string vcf_file;
} mut_param;

enum profile_t {CGMAP = 0, BEDMETHYL, METHBG, METHBED};

typedef struct{
    uint64_t seed_meth  = 0;
    bool pool_meth      = false;    /* levels drawn from the profile's, per context class */
    bool collect_non_cpg= true;
    int  profile_format = CGMAP;
    bool asm_is_bed     = false;
    std::string profile_file, asm_file, methdb_file, methdb_output;
    std::map<int, param_rec> params_map;  /* beta shape per context */
} meth_param;

#endif
