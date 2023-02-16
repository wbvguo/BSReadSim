import numpy as np
import subprocess

from typing import Dict
from UtilityFunctions import retrieve_iupac


class StreamHTSIM:
    '''
    stream HTSIM output for bisulfite reads generation
    :param str  sim_cmd : HTSIM commands for simulation
    :param bool pair_end: pair_end or not
    :rtype None
    '''
    def __init__(self, sim_cmd: list = None, pair_end: bool = True):
        self.sim_cmd  = sim_cmd
        self.pair_end = pair_end

    def __iter__(self):
        htsim = subprocess.Popen(self.sim_cmd, stdout=subprocess.PIPE, universal_newlines=True)
        sim_iter = iter(htsim.stdout.readline, b'')

        line  = self.get_line(sim_iter) # line is None when EOF
        while line:
            # collect all variant lines on the contig, after that sim_iter points to read lines
            if line == "Contig Variant Start":
                variant_contig, variant_dict = self.collect_variants(sim_iter)
                yield variant_contig, variant_dict

            # collect read pairs
            for collect_flag, read_pair in self.collect_reads(sim_iter):
                if collect_flag: # {1: collect_reads, 0: swith to collect_vars or EOF}
                    yield False, read_pair
                else:
                    line = "Contig Variant Start" if isinstance(read_pair, list) else None
                    break


    def collect_variants(self, sim_iter):
        '''collect variant lines from stdout'''
        variant_dict = {}
        variant_info = {}

        while True:
            line = self.get_line(sim_iter)
            if line == 'Contig Variant End':
                return variant_info['chrom'], variant_dict

            variant_info = self.process_variant_line(line)
            if variant_info['pos']:
                assert variant_info['pos'] not in variant_dict
                variant_dict[variant_info['pos']] = variant_info


    def collect_reads(self, sim_iter):
        '''collect read lines from stdout'''
        skip_flag = not self.pair_end

        while True:
            line  = self.get_line(sim_iter)
            if not line: # EOF
                yield 0, None
            elif line == "Contig Variant Start": # switch to collect variants
                yield 0, []
            else:
                read1 = self.process_read_lines(sim_iter, line = line)
                read2 = self.process_read_lines(sim_iter, skip = skip_flag)
                yield 1, [read1, read2]


    @staticmethod
    def get_line(sim_iter):
        '''receive lines from console'''
        try:
            line = next(sim_iter).strip()
        except StopIteration:
            print("End of output\n")
            return None
        else:
            return line


    @staticmethod
    def process_variant_line(line: str) -> Dict:
        '''parse variant lines'''
        line_split = line.split('\t')

        try:
            chrom, pos, ref, alt, heter_flag = line_split
        except ValueError:
            return dict(chrom=line_split[0], pos = None)
        else:
            heter = heter_flag == '+'
            indel = int(ref == '-') - int(alt == '-') # 1 for ref=='-', -1 for alt=='-', o.w. 0
            offset= indel * max(len(ref), len(alt))
            if indel:
                iupac  = None
            else:
                iupac  = retrieve_iupac(alt)
                alt    = list(set(iupac) - set(ref))[0]
            return dict(chrom=chrom, pos=int(pos), ref=ref, alt=alt,
                        offset=offset, heter=heter, indel=indel, iupac=iupac)


    @staticmethod
    def process_read_lines(sim_iter, line = None, skip = False):
        '''parse read lines'''
        if skip:
            next(sim_iter)
            next(sim_iter)
            next(sim_iter)
            next(sim_iter)
            return None

        if not line:
            line = next(sim_iter).strip()
        # header, seq, comment process
        read_id, pair, flag_pos, flag_mut, flag_indel, qual, cgr = line.split(' ')
        cgr = np.frombuffer(cgr.encode(), dtype=np.int8)
        seq = np.frombuffer(next(sim_iter).strip().encode(), dtype=np.int8)
        _, start, end, cover_pos, n_sub, n_indel, insert_size, inner_dist, ofs= next(sim_iter).strip().split(':')
        ofs = np.fromstring(ofs, dtype=np.int8, sep = ',')
        ctx = np.frombuffer(next(sim_iter).strip().encode(), np.int8)
        return dict(read_id=read_id, pair=int(pair), qual = int(qual),
                    flag_pos=int(flag_pos), flag_mut=int(flag_mut), flag_indel=int(flag_indel),
                    start=int(start), end=int(end), cover_pos=int(cover_pos),
                    n_sub=int(n_sub), n_indel=int(n_indel),
                    insert_size=int(insert_size), inner_dist=int(inner_dist),
                    cgr=cgr, seq=seq, ofs=ofs, ctx=ctx)
