#include <stdint.h>
#include <string.h>
#include <vector>
#include <string>
#include <set>
#include <htslib/vcf.h>
#include "struct.h"
#include "utility.h"
#include "haplo.h"


static std::string hex_id(uint64_t value)
{
    char text[24];
    snprintf(text, sizeof(text), "%llx", (unsigned long long)value);
    return text;
}

static std::string bases_str(const uint8_t *bases, int n)
{
    std::string text;
    for (int i = 0; i < n; ++i) text += "ACGTN"[bases[i]];
    return text;
}

void parse_vcf(const char *fname, const std::vector<chr_rec>& chr_vec, uint64_t seed_phase,
               std::vector<std::map<int, snpmeth_rec>>& vcf_vec, vcf_summary *summary)
{
    std::map<std::string, int> chr_idx = chr_index(chr_vec);
    vcf_vec.assign(chr_vec.size(), std::map<int, snpmeth_rec>());

    //open vcf file
    htsFile   *fp  = hts_open(fname,"r");
    if (fp == 0) die("cannot open VCF %s", fname);
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (hdr == 0) die("cannot read the VCF header of %s", fname);
    bcf1_t    *rec = bcf_init();

    //check the vcf file
    if (bcf_hdr_nsamples(hdr) != 1) die("VCF must have exactly one sample: %s", fname);

    vcf_summary counts;
    std::set<std::string> ids;
    std::vector<bool> seen(chr_vec.size(), false);
    int ngt_arr = 0, *gt = NULL;
    int prev_chr = -1, prev_pos = -1, ret;
    while ((ret = bcf_read(fp, hdr, rec)) >= 0) {
        uint64_t row = counts.rows + 1;
        //read CHROM, POS (0-based)
        const char *chr_name = bcf_hdr_id2name(hdr, rec->rid);
        std::map<std::string, int>::iterator found = chr_idx.find(chr_name);
        if (found == chr_idx.end()) die("VCF record %llu: unknown contig %s", (unsigned long long)row, chr_name);
        int chr = found->second;
        int pos = (int)rec->pos;
        if (chr < prev_chr || (chr == prev_chr && pos < prev_pos)) die("VCF record %llu: rows must be sorted in reference order", (unsigned long long)row);
        prev_chr = chr; prev_pos = pos;
        ++counts.rows;
        if (!seen[chr]) {seen[chr] = true; ++counts.contigs;}

        //unpack for read REF,ALT,ID
        bcf_unpack(rec, BCF_UN_STR);
        if (rec->n_allele != 2) die("VCF record %llu: only one ALT allele is supported", (unsigned long long)row);

        // check genotype: diploid 0/0, 0/1, 1/0, 1/1
        int ngt = bcf_get_genotypes(hdr, rec, &gt, &ngt_arr);
        if (ngt != 2 || bcf_gt_is_missing(gt[0]) || bcf_gt_is_missing(gt[1])) die("VCF record %llu: GT must be diploid 0/0, 0/1, 1/0, or 1/1", (unsigned long long)row);
        int snp_hap1 = bcf_gt_allele(gt[0]);
        int snp_hap2 = bcf_gt_allele(gt[1]);
        int is_phased= bcf_gt_is_phased(gt[1]);
        if (snp_hap1 < 0 || snp_hap2 < 0 || snp_hap1 > 1 || snp_hap2 > 1) die("VCF record %llu: GT must be diploid 0/0, 0/1, 1/0, or 1/1", (unsigned long long)row);
        if (snp_hap1 + snp_hap2 == 0) {++counts.reference_genotypes; continue;}

        // REF/ALT: uppercase A/C/G/T/N; normalize by trimming a shared prefix and suffix
        std::vector<uint8_t> ref, alt;
        for (int a = 0; a < 2; ++a) {
            const char *allele = rec->d.allele[a];
            std::vector<uint8_t> &bases = a == 0 ? ref : alt;
            for (const char *p = allele; *p; ++p) {
                if (!strchr("ACGTN", *p)) die("VCF record %llu: %s must contain only uppercase A/C/G/T/N", (unsigned long long)row, a == 0 ? "REF" : "ALT");
                bases.push_back(nst_nt4_table[(int)*p]);
            }
        }
        if (ref == alt) die("VCF record %llu: REF and ALT must differ", (unsigned long long)row);
        if (pos + ref.size() > chr_vec[chr].chr_len) die("VCF record %llu: REF extends past the end of its contig", (unsigned long long)row);
        size_t prefix = 0, suffix = 0;
        while (prefix < ref.size() && prefix < alt.size() && ref[prefix] == alt[prefix]) ++prefix;
        while (suffix < ref.size() - prefix && suffix < alt.size() - prefix
               && ref[ref.size() - 1 - suffix] == alt[alt.size() - 1 - suffix]) ++suffix;
        ref = std::vector<uint8_t>(ref.begin() + prefix, ref.end() - suffix);
        alt = std::vector<uint8_t>(alt.begin() + prefix, alt.end() - suffix);
        int start = pos + (int)prefix;
        for (size_t i = 0; i < ref.size(); ++i) if (ref[i] > 3) die("VCF record %llu: N is allowed only in a shared indel anchor", (unsigned long long)row);
        for (size_t i = 0; i < alt.size(); ++i) if (alt[i] > 3) die("VCF record %llu: N is allowed only in a shared indel anchor", (unsigned long long)row);

        //skip the following records:
        // 1. insert/delete longer than 4
        // 2. multi nucleotide polymorphism (MNP) or complex replacement
        snpmeth_rec tmp_snp;
        int base_change_pos = start;
        if (ref.size() == 1 && alt.size() == 1) {           //substitution
            tmp_snp.ref = ref[0]; tmp_snp.alt = alt[0]; tmp_snp.offset = 0;
        } else if (ref.empty() || alt.empty()) {
            if (std::max(ref.size(), alt.size()) > (size_t)MAX_INDEL) {++counts.long_indel; continue;}
            if (alt.empty()) {                                  //deletion: keyed by the first deleted base
                tmp_snp.ref = pack_bases(ref); tmp_snp.offset = -(int8_t)ref.size();
            } else {                                            //insertion: keyed by the base it follows
                if (start == 0) die("VCF record %llu: an insertion before the first base of a contig is not supported", (unsigned long long)row);
                base_change_pos = start - 1;
                tmp_snp.alt = pack_bases(alt); tmp_snp.offset = (int8_t)alt.size();
            }
        } else {
            ++(ref.size() == alt.size() ? counts.mnp : counts.complex_replacement);
            continue;
        }

        std::map<int, snpmeth_rec>& snpmeth_map = vcf_vec[chr];
        // for unphased heterozygotes, randomly assign the haplotype
        if (snp_hap1 != snp_hap2 && !is_phased) {
            Rng rng(seed_of(seed_phase, chr, snpmeth_map.size()));
            if (rng.unif() < 0.5) {int tmp_hap = snp_hap1; snp_hap1 = snp_hap2; snp_hap2 = tmp_hap;}
        }
        tmp_snp.hap1 = (int8_t)snp_hap1; tmp_snp.hap2 = (int8_t)snp_hap2;
        tmp_snp.source = FROM_VCF;

        uint64_t address = ((uint64_t)chr << 32) | snpmeth_map.size();
        std::string id = strcmp(rec->d.id, ".") ? rec->d.id : "vcf_" + hex_id(address);
        if (!ids.insert(id).second) {id += "@" + hex_id(address); ids.insert(id);}
        tmp_snp.id = id;

        if (snpmeth_map.count(base_change_pos)) die("VCF record %llu: normalized variants overlap or share an insertion point", (unsigned long long)row);
        snpmeth_map[base_change_pos] = tmp_snp;
        ++counts.retained;
    }
    if (ret < -1) die("cannot parse VCF %s after record %llu", fname, (unsigned long long)counts.rows);

    free(gt);
    bcf_destroy(rec);
    bcf_hdr_destroy(hdr);
    if ((ret = hts_close(fp))) die("cannot close VCF %s", fname);
    if (summary) *summary = counts;
}


void sim_mut_diref(const ref_seq *ref, const mut_param *mut_set, int chr_idx, std::map<int, snpmeth_rec>& snpmeth_map)
{
    snpmeth_map.clear();
    if (mut_set->mut_rate <= 0) return;
    Rng rng(seed_of(mut_set->seed_snp, chr_idx));
    int i, c, l = (int)ref->len;
    for (i = 0; i < l; ++i) {
        c = nst_nt4_table[(int)ref->seq[i]];
        if (c < 4 && rng.unif() < mut_set->mut_rate) { // mutation
            snpmeth_rec tmp_snp;
            // hom: both haplotypes; het: one of them
            if (mut_set->homozygous_only || rng.unif() < 0.333333) {tmp_snp.hap1 = tmp_snp.hap2 = 1;}
            else if (rng.unif() < 0.5) {tmp_snp.hap1 = 1;} else {tmp_snp.hap2 = 1;}
            tmp_snp.source = FROM_DENOVO;
            tmp_snp.id = "varsim_" + hex_id(((uint64_t)chr_idx << 32) | snpmeth_map.size());
            int key = i;
            if (rng.unif() >= mut_set->indel_frac) { // substitution
                tmp_snp.ref = c;
                tmp_snp.alt = (c + (int)(rng.unif() * 3.0 + 1)) & 3;   // random mutation
            } else if (rng.unif() < 0.5) { // deletion, extended while the next base is not N
                std::vector<uint8_t> deleted(1, (uint8_t)c);
                while ((int)deleted.size() < MAX_INDEL && i + 1 < l
                       && nst_nt4_table[(int)ref->seq[i + 1]] < 4 && rng.unif() < mut_set->indel_extn) {
                    deleted.push_back(nst_nt4_table[(int)ref->seq[++i]]);
                }
                tmp_snp.ref = pack_bases(deleted);
                tmp_snp.offset = -(int8_t)deleted.size();
            } else { // insertion after base i; the next base must not start a deletion
                std::vector<uint8_t> inserted;
                do {
                    inserted.push_back((uint8_t)(rng.unif() * 4.0));
                } while ((int)inserted.size() < MAX_INDEL && rng.unif() < mut_set->indel_extn);
                tmp_snp.alt = pack_bases(inserted);
                tmp_snp.offset = (int8_t)inserted.size();
                ++i;
            }
            snpmeth_map[key] = tmp_snp;
        }
    }
}


void sim_mut_map(const ref_seq *ref, std::map<int, snpmeth_rec>& snpmeth_map, mutseq_t *hap1, mutseq_t *hap2)
{
    // initiate
    mutseq_t *ret[2];
    ret[0] = hap1; ret[1] = hap2;
    ret[0]->l = ref->len; ret[1]->l = ref->len;
    ret[0]->s = (mut_t *)calloc(ref->len, sizeof(mut_t));
    ret[1]->s = (mut_t *)calloc(ref->len, sizeof(mut_t));
    if (ret[0]->s == 0 || ret[1]->s == 0) die("cannot allocate the haplotypes of %s", ref->name);

    int i, l = (int)ref->len;
    for (i = 0; i != l; ++i) {
        ret[0]->s[i] = ret[1]->s[i] = (mut_t)nst_nt4_table[(int)ref->seq[i]];
    }

    int covered = -1;   // last reference position taken by the previous variant
    int ins_anchor = -2;// base followed by the previous insertion
    for (std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it) {
        i = it->first;
        snpmeth_rec &snp = it->second;
        int n = var_span(snp);                      // reference bases it takes
        // a deletion may not start right after an insertion: both would share one anchor
        if (i <= covered || i < 0 || i + n > l || (snp.offset < 0 && i == ins_anchor + 1)) {
            die("variants on %s overlap, share an insertion point, or pass the contig end (position %d)", ref->name, i + 1);
        }
        for (int k = 0; snp.offset <= 0 && k < n; ++k) {         // the REF bases (an insertion has none)
            uint8_t ref_base = snp.offset < 0 ? packed_base(snp.ref, k) : snp.ref;
            if (nst_nt4_table[(int)ref->seq[i + k]] != ref_base) die("variant %s on %s: REF does not match the reference", snp.id.c_str(), ref->name);
        }
        covered = i + n - 1;
        if (snp.offset > 0) ins_anchor = i;

        for (int h = 0; h < 2; ++h) {
            if ((h == 0 ? snp.hap1 : snp.hap2) != 1) continue;
            if (snp.offset == 0) {                  // SNP substitution
                ret[h]->s[i] = SUBSTITUTE | snp.alt;
            } else if (snp.offset < 0) {            // deletion
                for (int k = 0; k < n; ++k) ret[h]->s[i + k] |= DELETE;
            } else {                                // insertion
                int num_ins = snp.offset;
                ret[h]->s[i] |= (num_ins << 12) | ((snp.alt & ((1 << (2 * num_ins)) - 1)) << 4);
            }
        }
    }
}


// the haplotype bases before reference position i, nearest first
static void hap_bases_before(mutseq_t *hap, int i, uint8_t *out, int n)
{
    int got = 0;
    for (--i; i >= 0 && got < n; --i) {
        int c = hap->s[i];
        if ((c & mutmsk) == DELETE) continue;
        for (int k = mut_ins(c) - 1; k >= 0 && got < n; --k) out[got++] = mut_ins_base(c, k);   // inserted bases follow the base
        if (got < n) out[got++] = c & 0xf;
    }
    for (; got < n; ++got) out[got] = 4;
}

int gen_frag_seq(mutseq_t *hap, int pos_l, int pos_r, int skip, int len, frag_seq *fs)
{
    fs->seq.clear(); fs->pos.clear(); fs->context.clear(); fs->offset.clear();
    fs->event_start.clear(); fs->event_end.clear(); fs->event_kind.clear(); fs->event_deleted.clear();
    hap_bases_before(hap, pos_l, fs->before, 3);

    // the skipped bases move into the left flank; the 3 bases after the template are the right flank
    int n_after = 0, deleting = 0, ins_anchor = -1, del_event = -1;
    for (int i = pos_l; i < hap->l && n_after < 3; ++i) {
        int c = hap->s[i], mut_type = c & mutmsk;
        int k = (int)fs->seq.size();
        if (mut_type == DELETE) {       // an empty interval where the deletion starts, inside the template
            if (deleting) {
                if (del_event >= 0) ++fs->event_deleted[del_event];
            } else if (k > 0 && k < len && i < pos_r) {
                del_event = (int)fs->event_start.size();
                fs->event_start.push_back(k); fs->event_end.push_back(k); fs->event_kind.push_back(VAR_DEL); fs->event_deleted.push_back(1);
            } else {
                del_event = -1;
            }
            deleting = 1;
            continue;
        }
        deleting = 0;
        int num_ins = mut_ins(c);
        for (int j = -1; j < num_ins; ++j) {    // the base, then the bases inserted after it
            int base = j < 0 ? c & 0xf : mut_ins_base(c, j);
            k = (int)fs->seq.size();
            if (skip > 0) {
                --skip;
                fs->before[2] = fs->before[1]; fs->before[1] = fs->before[0]; fs->before[0] = base;
            } else if (i < pos_r && k < len) {
                // context: 0x00 Match, 0x10 SNP, 0x30 INSERT (upper half); the site context is added later
                fs->seq.push_back(base);
                fs->pos.push_back(i);
                fs->offset.push_back(j < 0 ? 0 : j);
                if (j >= 0) {
                    fs->context.push_back(INSR);
                    if (ins_anchor == i) {++fs->event_end.back();}
                    else {fs->event_start.push_back(k); fs->event_end.push_back(k + 1); fs->event_kind.push_back(VAR_INS); fs->event_deleted.push_back(0); ins_anchor = i;}
                } else if (mut_type == SUBSTITUTE) {
                    fs->context.push_back(SNV);
                    fs->event_start.push_back(k); fs->event_end.push_back(k + 1); fs->event_kind.push_back(VAR_SNV); fs->event_deleted.push_back(0);
                } else {
                    fs->context.push_back(MATCH);
                }
            } else if (n_after < 3) {
                fs->after[n_after++] = base;
            }
        }
    }
    for (; n_after < 3; ++n_after) fs->after[n_after] = 4;
    return (int)fs->seq.size();
}

// base k of the fragment extended by its flanking bases (4 outside them)
static inline int ext_base(frag_seq *fs, int k)
{
    int n = (int)fs->seq.size();
    if (k >= 0 && k < n) return fs->seq[k];
    if (k < 0 && k >= -3) return fs->before[-k - 1];
    if (k >= n && k < n + 3) return fs->after[k - n];
    return 4;
}

void frag_contexts(frag_seq *fs, std::vector<uint8_t>& ctx, bool collect_non_cpg)
{
    int n = (int)fs->seq.size();
    ctx.assign(n, 0);
    for (int k = 0; k < n; ++k) {
        int c = fs->seq[k];
        if (!cg_table[c]) continue;
        uint8_t context = c == 1 ? get_context(c, ext_base(fs, k + 1), ext_base(fs, k + 2))
                                 : get_context(c, ext_base(fs, k - 1), ext_base(fs, k - 2));
        if (!collect_non_cpg && ctx_class(context) != 1) context = 0;
        ctx[k] = context;
    }
}

uint16_t frag_kmer(frag_seq *fs, int k)
{
    uint16_t kmer = 0;
    for (int j = -3; j <= 3; ++j) {
        int c = ext_base(fs, k + j);
        if (c > 3) return NO_KMER;
        kmer = (kmer << 2) | c;
    }
    return kmer;
}


void save_vcf_header(out_t *out)
{
    out_printf(out, "##fileformat=VCFv4.3\n##source=BSReadSim\n"
                    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSIMULATED\n");
}

void save_vcf_chr(out_t *out, const ref_seq *ref, std::map<int, snpmeth_rec>& snpmeth_map)
{
    static const char *genotypes[] = {"", "1|0", "0|1", "1|1"};
    int l = (int)ref->len;
    for (std::map<int, snpmeth_rec>::iterator it = snpmeth_map.begin(); it != snpmeth_map.end(); ++it) {
        int i = it->first;
        snpmeth_rec &snp = it->second;
        std::vector<uint8_t> ref_bases, alt_bases;
        int pos = i + 1;    // 1-based
        if (snp.offset == 0) {
            ref_bases.push_back(snp.ref); alt_bases.push_back(snp.alt);
        } else if (snp.offset > 0) {    // the base it follows is the anchor
            uint8_t anchor = nst_nt4_table[(int)ref->seq[i]];
            ref_bases.push_back(anchor); alt_bases.push_back(anchor);
            for (int k = 0; k < snp.offset; ++k) alt_bases.push_back(packed_base(snp.alt, k));
        } else {                        // the anchor is the base before, or after at the contig start
            int n = -snp.offset;
            std::vector<uint8_t> deleted;
            for (int k = 0; k < n; ++k) deleted.push_back(packed_base(snp.ref, k));
            if (i > 0) {
                uint8_t anchor = nst_nt4_table[(int)ref->seq[i - 1]];
                ref_bases.push_back(anchor); ref_bases.insert(ref_bases.end(), deleted.begin(), deleted.end());
                alt_bases.push_back(anchor);
                pos = i;
            } else if (i + n < l) {
                uint8_t anchor = nst_nt4_table[(int)ref->seq[i + n]];
                ref_bases = deleted; ref_bases.push_back(anchor);
                alt_bases.push_back(anchor);
            } else {
                die("variant %s has no anchor base", snp.id.c_str());
            }
        }
        out_printf(out, "%s\t%d\t%s\t%s\t%s\t.\tPASS\t.\tGT\t%s\n", ref->name, pos, snp.id.c_str(),
                   bases_str(ref_bases.data(), (int)ref_bases.size()).c_str(), bases_str(alt_bases.data(), (int)alt_bases.size()).c_str(),
                   genotypes[(snp.hap1 == 1) | (snp.hap2 == 1) << 1]);
    }
}
