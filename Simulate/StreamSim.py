import os
import pickle
import subprocess

from array import array
from typing import Dict, Union
from Utils.UtilityFunctions import retrieve_iupac


class StreamWGSIM:
    '''
    stream WGSIM output for bisulfite reads generation
    :param str sim_cmd: WGSIM commands for simulation
    :rtype None
    '''
    
    def __init__(self, sim_cmd: str = None, pair_end: bool = True):
        self.sim_cmd  = sim_cmd
        self.pair_end = pair_end
    
    
    def __iter__(self):
        wgsim = subprocess.Popen(self.sim_cmd, stdout=subprocess.PIPE, universal_newlines=True)
        sim_iter = iter(wgsim.stdout.readline, b'')
        
        line  = self.get_line(sim_iter) # line is None when EOF
        while line:
            if line == "Contig Variant Start":
                # collect all the variant lines, after that sim_iter points to read lines
                variant_contig, variant_dict = self.collect_variants(sim_iter)
                yield variant_contig, variant_dict
                collect_reads_flag = True
            
            # collect read pairs
            while collect_reads_flag:
                line = self.get_line(sim_iter)
                if not line or line == "Contig Variant Start":
                    collect_reads_flag = False
                    break
                
                read1 = self.collect_read(sim_iter, line)
                read2 = self.collect_read(sim_iter)
                assert read1["read_id"] == read2["read_id"]
                yield False, [read1, read2]
                
    
    def collect_variants(self, sim_iter):
        variant_dict = {}
        
        while True:
            line = next(sim_iter).strip()
            if line == 'Contig Variant End':
                if isinstance(variant_info, dict):
                    return variant_info['chrom'], variant_dict
                else:
                    return variant_info, variant_dict
       
            variant_info = self.process_variant_line(line)
            if isinstance(variant_info, dict):
                assert variant_info['pos'] not in variant_dict
                variant_dict[variant_info['pos']] = variant_info
    
    
    def collect_read(self, sim_iter, line = None):
        if not line:
            line = next(sim_iter).strip()
        read_dict= self.process_read_name(line)
        read_dict['seq'] = next(sim_iter).strip()
        next(sim_iter)
        next(sim_iter)
        return read_dict
    

    @staticmethod
    def get_line(sim_iter):
        try:
            line = next(sim_iter).strip()
        except StopIteration:
            print("End of output\n")
            return None
        else:
            return line

        
    @staticmethod
    def process_variant_line(line: str) -> Dict:
        line_split = line.split('\t')
        
        try:
            chrom, pos, ref, alt, heter_flag = line_split
        except ValueError:
            return line_split[0]
        else:
            heter = True if heter_flag == '+' else False
            indel = int(ref == '-') - int(alt == '-') # indel=1 for ref=='-', -1 for alt=='-', o.w. 0
            iupac = retrieve_iupac(alt) if indel == 0 else None
            return dict(chrom=chrom, pos=int(pos), ref=ref, alt=alt, 
                        heter=heter, indel=indel, iupac=iupac)

    @staticmethod
    def process_read_name(line: str) -> Dict[str, Union[str, int, array]]:
        read_info = line.split(':')
        chrom, start, end, insert_size, read_id, cigar, pair, c_base_info, g_base_info = read_info
        return dict(chrom=chrom.replace('@', ''), start=int(start), end=int(end),
                    insert_size=insert_size, read_id=read_id, cigar=cigar, pair=int(pair),
                    c_base_info=c_base_info, g_base_info=g_base_info)



class StreamOutput:
    '''
    write the simulation values/variants to disk
    :param str sim_output: path to the simulation folder
    :param str shuffle: whether to shuffle the reads or not, default is reads are 
                        segemented by contig_id
    :rtype None
    '''
    
    def __init__(self, sim_output=None, shuffle = True):
        self.sim_output = sim_output
        self.shuffle = shuffle

    def generate_sim_directory(self):
        if not os.path.isdir(self.sim_output):
            os.makedirs(self.sim_output, exist_ok=False)

    def output_contig(self, contig_id, contig_profile, values=False, variant=False):
        if contig_id:
            contig_label = contig_id
            if values:
                contig_label = f'{contig_id}_values'
            elif variant:
                contig_label = f'{contig_id}_variants'
            with open(f'{self.sim_output}/{contig_label}.pkl', 'wb') as contig_out:
                pickle.dump(contig_profile, contig_out)

    def load_contig(self, contig_id, values=False):
        contig_label = contig_id
        if values:
            contig_label = f'{contig_id}_values'
        try:
            with open(f'{self.sim_output}/{contig_label}.pkl', 'rb') as contig_out:
                contig_profile = pickle.load(contig_out)
        except FileNotFoundError:
            return None
        else:
            return contig_profile
