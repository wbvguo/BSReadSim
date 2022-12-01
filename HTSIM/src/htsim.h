#include <stdint.h>
#include <vector>


typedef struct {
    int pos, ref, alt, geno;
} snp_rec;

typedef struct {
    int pos_l, pos_r, strand;
    float score;
} probe_rec;

typedef struct {
    char *name, *contig;
} probe_meta;

typedef struct {
    int len = -1;         /* length of cutting site */
    int idx = -1;         /* cutting position on *seq */
    std::vector<int> seq; /* sequence encoded by numbers*/
} cut_site;

typedef struct {
    int pos_l, pos_r;
    int len;
    uint8_t cut_l, cut_r;
} cut_frag;

typedef struct {
    int pos;
    uint8_t type;
} cut_pos;

typedef struct {
    int pos_l, pos_r;
} insrt_frag;

