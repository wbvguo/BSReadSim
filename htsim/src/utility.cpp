#include <stdint.h>
#include <stdarg.h>
#include <stdlib.h>
#include <errno.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <random>
#include <htslib/kseq.h>
#include "utility.h"


/*-------------------------errors-------------------------*/
void die(const char *format, ...)
{
    va_list args;
    va_start(args, format);
    fflush(stdout);
    fputs("htsim: ", stderr);
    vfprintf(stderr, format, args);
    fputc('\n', stderr);
    va_end(args);
    exit(2);
}

void die_at(const char *label, uint64_t line_no, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    fflush(stdout);
    fprintf(stderr, "htsim: %s line %llu: ", label, (unsigned long long)line_no);
    vfprintf(stderr, format, args);
    fputc('\n', stderr);
    va_end(args);
    exit(2);
}


/*-------------------------text-------------------------*/
void split_line(char *line, std::vector<char*> &fields, char sep)
{
    fields.clear();
    char *p = line;
    for (char *q = line;; ++q) {
        if (*q == sep || *q == '\0') {
            int c = *q;
            *q = 0;
            fields.push_back(p);
            if (c == '\0') break;
            p = q + 1;
        }
    }
}

static bool parse_u64(const char *text, uint64_t *value)
{
    char *end = 0;
    errno = 0;
    *value = strtoull(text, &end, 10);
    return *text >= '0' && *text <= '9' && *end == '\0' && !errno;
}

static bool parse_double(const char *text, double *value)
{
    char *end = 0;
    *value = strtod(text, &end);
    return *text != '\0' && *end == '\0' && isfinite(*value);
}

uint64_t to_u64(const char *text, const char *what)
{
    uint64_t value;
    if (!parse_u64(text, &value)) die("%s must be an unsigned integer: '%s'", what, text);
    return value;
}

uint32_t to_u32(const char *text, const char *what)
{
    uint64_t value;
    if (!parse_u64(text, &value) || value > UINT32_MAX) die("%s must be an unsigned integer: '%s'", what, text);
    return (uint32_t)value;
}

double to_double(const char *text, const char *what)
{
    double value;
    if (!parse_double(text, &value)) die("%s must be a finite number: '%s'", what, text);
    return value;
}

double to_prob(const char *text, const char *what)
{
    double value = to_double(text, what);
    if (value < 0 || value > 1) die("%s must be in [0, 1]", what);
    return value;
}

bool valid_utf8(const std::string &text)
{
    const unsigned char *s = (const unsigned char *)text.data(), *end = s + text.size();
    while (s < end) {
        unsigned c = *s++;
        if (c < 0x80) continue;
        int n = c >= 0xc2 && c <= 0xdf ? 1 : c >= 0xe0 && c <= 0xef ? 2 : c >= 0xf0 && c <= 0xf4 ? 3 : -1;
        if (n < 0 || end - s < n) return false;
        unsigned code = c & (0x3f >> n);
        for (int k = 0; k < n; ++k, ++s) {
            if ((*s & 0xc0) != 0x80) return false;
            code = code << 6 | (*s & 0x3f);
        }
        if ((n == 2 && (code < 0x800 || (code >= 0xd800 && code <= 0xdfff))) || (n == 3 && (code < 0x10000 || code > 0x10ffff))) return false;
    }
    return true;
}


/*-------------------------random numbers: SplitMix64-------------------------*/
uint64_t mix64(uint64_t value)
{
    value = (value ^ (value >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)) * UINT64_C(0x94d049bb133111eb);
    return value ^ (value >> 31);
}

uint64_t seed_of(uint64_t seed, uint64_t a, uint64_t b, uint64_t c)
{
    const uint64_t golden = UINT64_C(0x9e3779b97f4a7c15);
    uint64_t value = mix64(seed + golden);
    value = mix64(value ^ (a + golden));
    value = mix64(value ^ (b + 2 * golden));
    return mix64(value ^ (c + 3 * golden));
}

Rng::result_type Rng::operator()()
{
    state += UINT64_C(0x9e3779b97f4a7c15);
    return mix64(state);
}

double Rng::unif() {return (double)((*this)() >> 11) * 0x1.0p-53;}

uint64_t Rng::below(uint64_t n)
{
    __extension__ typedef unsigned __int128 wide;
    return (uint64_t)(((wide)(*this)() * n) >> 64);
}

double Rng::norm(double mean, double sd) {return std::normal_distribution<double>(mean, sd)(*this);}

double Rng::beta(double alpha, double beta)
{
    double x = std::gamma_distribution<double>(alpha, 1.0)(*this);
    double y = std::gamma_distribution<double>(beta, 1.0)(*this);
    return x + y > 0 ? x / (x + y) : (unif() < alpha / (alpha + beta) ? 1.0 : 0.0);
}


/*-------------------------intervals-------------------------*/
std::vector<std::pair<int, int>> merge_intervals(std::vector<std::pair<int, int>> intervals)
{
    std::sort(intervals.begin(), intervals.end());
    std::vector<std::pair<int, int>> merged;
    for (size_t i = 0; i < intervals.size(); ++i) {
        if (merged.empty() || intervals[i].first > merged.back().second) merged.push_back(intervals[i]);
        else merged.back().second = std::max(merged.back().second, intervals[i].second);
    }
    return merged;
}


/*-------------------------text input-------------------------*/
text_in *in_open(const char *fname, const char *label)
{
    text_in *in = new text_in();
    in->fp = hts_open(fname, "r");
    if (in->fp == 0) die("cannot open %s %s", label, fname);
    in->line = {0, 0, 0};
    in->label = label;
    in->line_no = 0;
    return in;
}

bool in_next(text_in *in)
{
    int ret = hts_getline(in->fp, KS_SEP_LINE, &in->line);
    if (ret < -1) die("cannot read %s", in->label);
    if (ret < 0) return false;
    ++in->line_no;
    if (in->line.l && in->line.s[in->line.l - 1] == '\r') in->line.s[--in->line.l] = 0;
    return true;
}

size_t in_split(text_in *in, char sep)
{
    split_line(in->line.s, in->f, sep);
    return in->f.size();
}

void in_close(text_in *in)
{
    hts_close(in->fp);
    free(in->line.s);
    delete in;
}

uint32_t in_u32(text_in *in, const char *text, const char *what)
{
    uint64_t value;
    if (!parse_u64(text, &value) || value > UINT32_MAX) IN_DIE(in, "%s must be an unsigned integer: '%s'", what, text);
    return (uint32_t)value;
}

double in_double(text_in *in, const char *text, const char *what)
{
    double value;
    if (!parse_double(text, &value)) IN_DIE(in, "%s must be a finite number: '%s'", what, text);
    return value;
}

double in_prob(text_in *in, const char *text, const char *what)
{
    double value = in_double(in, text, what);
    if (value < 0 || value > 1) IN_DIE(in, "%s must be in [0, 1]", what);
    return value;
}

char in_strand(text_in *in, const char *text)
{
    if (strcmp(text, "+") && strcmp(text, "-") && strcmp(text, ".")) IN_DIE(in, "strand must be +, -, or .");
    return text[0];
}

int in_chr(text_in *in, const std::map<std::string, int>& chr_idx, const char *name)
{
    std::map<std::string, int>::const_iterator found = chr_idx.find(name);
    if (found == chr_idx.end()) IN_DIE(in, "unknown contig %s", name);
    return found->second;
}

bool bed_header(const char *line)
{
    return !*line || line[0] == '#' || !strcmp(line, "track") || !strncmp(line, "track ", 6)
        || !strcmp(line, "browser") || !strncmp(line, "browser ", 8);
}


/*-------------------------text output-------------------------*/
out_t *out_open(const std::string &path, bool compress, int threads, int level)
{
    out_t *out = new out_t();
    std::string mode = !compress ? "wu" : level < 0 ? "w" : "w" + std::to_string(level);
    out->fp = bgzf_open(path.c_str(), mode.c_str());
    if (out->fp == 0) die("cannot create %s", path.c_str());
    if (compress && threads > 1) bgzf_mt(out->fp, threads, 256);
    out->buf = {0, 0, 0};
    out->path = path;
    return out;
}

void out_printf(out_t *out, const char *format, ...)
{
    va_list args;
    va_start(args, format);
    kvsprintf(&out->buf, format, args);
    va_end(args);
    if (out->buf.l >= (1 << 20)) {
        if (bgzf_write(out->fp, out->buf.s, out->buf.l) < 0) die("cannot write %s", out->path.c_str());
        out->buf.l = 0;
    }
}

void out_close(out_t *out)
{
    if (out->buf.l && bgzf_write(out->fp, out->buf.s, out->buf.l) < 0) die("cannot write %s", out->path.c_str());
    if (bgzf_close(out->fp) != 0) die("cannot write %s", out->path.c_str());
    free(out->buf.s);
    delete out;
}


/*-------------------------pipes-------------------------*/
void widen_pipe(int fd)
{
#ifdef F_SETPIPE_SZ
    fcntl(fd, F_SETPIPE_SZ, 1 << 20);   /* fails harmlessly when fd is not a pipe */
#else
    (void)fd;
#endif
}
