#ifndef UTILITY_H
#define UTILITY_H

// General tools that know nothing of BSReadSim's data: errors, number parsing,
// random numbers, threads, text input and output, and pipes.

#include <stdint.h>
#include <limits>
#include <string>
#include <vector>
#include <map>
#include <thread>

#include <htslib/hts.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>


/*-------------------------errors-------------------------*/
// print "htsim: message" and exit
[[noreturn]] void die(const char *format, ...);
// print "htsim: <label> line <line_no>: message" and exit
[[noreturn]] void die_at(const char *label, uint64_t line_no, const char *format, ...);


/*-------------------------text-------------------------*/
// split a line at tabs (in place)
void split_line(char *line, std::vector<char*> &fields, char sep = '\t');
// checked conversions; `what` names the value in error messages
uint32_t to_u32(const char *text, const char *what);
uint64_t to_u64(const char *text, const char *what);
double to_double(const char *text, const char *what);
double to_prob(const char *text, const char *what);
// well-formed UTF-8 (no overlong forms, surrogates, or code points above U+10FFFF)
bool valid_utf8(const std::string &text);


/*-------------------------random numbers-------------------------*/
// Every random decision draws from an Rng seeded by hashing a stage seed with
// what is being decided (a contig and a fragment ordinal, or a site), so the
// output never depends on the order of the work or on the number of threads.
uint64_t mix64(uint64_t value);
uint64_t seed_of(uint64_t seed, uint64_t a, uint64_t b = 0, uint64_t c = 0);

struct Rng {
    typedef uint64_t result_type;
    uint64_t state;
    explicit Rng(uint64_t seed) : state(seed) {}
    static constexpr result_type min() {return 0;}
    static constexpr result_type max() {return std::numeric_limits<uint64_t>::max();}
    result_type operator()();
    double unif();                      /* [0, 1) */
    uint64_t below(uint64_t n);         /* [0, n) */
    double norm(double mean, double sd);
    double beta(double alpha, double beta);
};


/*-------------------------intervals-------------------------*/
// sorted, disjoint half-open intervals covering the same positions
std::vector<std::pair<int, int>> merge_intervals(std::vector<std::pair<int, int>> intervals);


/*-------------------------threads-------------------------*/
// run work(t, begin, end) on `threads` threads that split [0, n) into contiguous ranges
template <typename F>
void parallel_for(uint64_t n, int threads, F work)
{
    if (threads <= 1) {work(0, 0, n); return;}
    std::vector<std::thread> workers;
    for (int t = 0; t < threads; ++t) workers.emplace_back(work, t, n * t / threads, n * (t + 1) / threads);
    for (size_t t = 0; t < workers.size(); ++t) workers[t].join();
}


/*-------------------------text input-------------------------*/
// a text input (plain or gzip) read line by line; errors name its label and the line number
typedef struct {
    htsFile *fp;
    kstring_t line;             /* the current line, without its end of line */
    const char *label;          /* e.g. "target BED" */
    uint64_t line_no;
    std::vector<char*> f;       /* its fields, after in_split */
} text_in;
text_in *in_open(const char *fname, const char *label);
bool in_next(text_in *in);                              /* false after the last line */
size_t in_split(text_in *in, char sep = '\t');
void in_close(text_in *in);
#define IN_DIE(in, ...) die_at((in)->label, (in)->line_no, __VA_ARGS__)
// checked conversions of a field of the current line
uint32_t in_u32(text_in *in, const char *text, const char *what);
double in_double(text_in *in, const char *text, const char *what);
double in_prob(text_in *in, const char *text, const char *what);
char in_strand(text_in *in, const char *text);          /* '+', '-', or '.' */
// the index of the contig a row names (see chr_index)
int in_chr(text_in *in, const std::map<std::string, int>& chr_idx, const char *name);
// BED lines without data: empty, comment, track, and browser lines
bool bed_header(const char *line);


/*-------------------------text output-------------------------*/
// plain, or BGZF when compress (readable by gzip), at a level 0-9 (-1: zlib's default, 6)
typedef struct {
    BGZF *fp;
    kstring_t buf;
    std::string path;
} out_t;
out_t *out_open(const std::string &path, bool compress, int threads = 1, int level = -1);
void out_printf(out_t *out, const char *format, ...);
void out_close(out_t *out);


/*-------------------------pipes-------------------------*/
// let a pipe on fd hold 1 MiB where Linux allows (a default one holds 64 KiB), so a
// large stream passes in fewer turns; nothing happens to a file or a terminal
void widen_pipe(int fd);

#endif
