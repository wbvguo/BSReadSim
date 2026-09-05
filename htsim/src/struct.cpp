#include <stdint.h>
#include <string.h>
#include <algorithm>
#include "struct.h"
#include "utility.h"

static uint8_t CG = 0x01; //5to3
static uint8_t CHG= 0x03;
static uint8_t CHH= 0x07;
static uint8_t GC = 0x09; //3to5
static uint8_t GDC= 0x0b;
static uint8_t GDD= 0x0f;
//0110**: 24-27; 01**10: 18, 22, 30; 01****: the rest of 16-31
//1001**: 36-39; 10**01: 33, 41, 45; 10****: the rest of 32-47
//encode not as 1,3,5; have problem with print 5 (or 13) when putcns
const uint8_t cg_context_table[64] = {
    0,   0,   0,   0,    0,   0,   0,   0,
    0,   0,   0,   0,    0,   0,   0,   0,
    CHH, CHH, CHG, CHH,  CHH, CHH, CHG, CHH,
    CG,  CG,  CG,  CG,	 CHH, CHH, CHG, CHH,
    GDD, GDC, GDD, GDD,  GC,  GC,  GC,  GC,
    GDD, GDC, GDD, GDD,  GDD, GDC, GDD, GDD,
    0,   0,   0,   0,    0,   0,   0,   0,
    0,   0,   0,   0,    0,   0,   0,   0,
};

const uint8_t cg_table[5] = {0, 1, 1, 0, 0}; // for C/G check

const uint8_t nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

const mut_t mutmsk = (mut_t)0xf000;


//  global variables.
const char PACKAGE_VERSION[] = BSREADSIM_VERSION;
const char *const TECH_NAMES[6] = {"WGBS", "RRBS", "TBS", "WGS", "WES", "TS"};
const char *const READ_FORMATS[4] = {"fastq.gz", "fastq", "bam", "sam"};


// context of a C (G) from its next (previous) two bases on the same haplotype
uint8_t get_context(int c, int d1, int d2)
{
    if (c != 1 && c != 2) return 0;
    if (d1 >= 4) return 0;                          // missing or N neighbor
    if (d1 == 3 - c) return c == 1 ? CG : GC;       // CG needs only the first neighbor
    if (d2 >= 4) return 0;
    return cg_context_table[c << 4 | d1 << 2 | d2];
}

const char *context_name(uint8_t context)
{
    switch (context) {
    case 0x01: return "CG-C";
    case 0x03: return "CHG-C";
    case 0x07: return "CHH-C";
    case 0x09: return "CG-G";
    case 0x0b: return "CHG-G";
    case 0x0f: return "CHH-G";
    default: die("invalid context code %d", context);
    }
}


uint8_t ctx_of_name(const char *name)
{
    return !strcmp(name, "CG") ? CG : !strcmp(name, "CHG") ? CHG : !strcmp(name, "CHH") ? CHH : 0;
}


// variants
uint16_t pack_bases(const std::vector<uint8_t>& bases)
{
    uint16_t value = 0;
    for (int i = (int)bases.size() - 1; i >= 0; --i) value = (value << 2) | bases[i];
    return value;
}

std::string packed_str(uint16_t packed, int n)
{
    std::string text;
    for (int k = 0; k < n; ++k) text += "ACGT"[packed_base(packed, k)];
    return text;
}


// contigs
std::map<std::string, int> chr_index(const std::vector<chr_rec>& chr_vec)
{
    std::map<std::string, int> chr_idx;
    for (size_t i = 0; i < chr_vec.size(); ++i) chr_idx[chr_vec[i].name] = (int)i;
    return chr_idx;
}


// methdb
meth_rec *find_site(std::vector<meth_rec>& meth_vec, int pos)
{
    std::vector<meth_rec>::iterator site = std::lower_bound(meth_vec.begin(), meth_vec.end(), pos,
        [](const meth_rec &rec, int p) {return rec.pos < p;});
    return site != meth_vec.end() && site->pos == pos ? &*site : 0;
}


// N bases
int count_n(const ref_seq *ref, int from, int to)
{
    int n = 0;
    for (int k = from; k < to; ++k) n += nst_nt4_table[(int)ref->seq[k]] > 3;
    return n;
}

int count_n(const std::vector<uint8_t>& bases, int from, int to)
{
    int n = 0;
    for (int k = from; k < to; ++k) n += bases[k] > 3;
    return n;
}
