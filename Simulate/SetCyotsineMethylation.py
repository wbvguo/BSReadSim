import re
import random
import numpy as np
from typing import Dict, List
from tqdm import tqdm

from Bio import SeqIO
from Utils.CGmapIterator import OpenCGmap
from StreamSim import StreamOutput


class SetCytosineMethylation:
    
    def __init__(self, reference_file: str = None, sim_output: str = None, 
                 meth_ref: str = None, cgmap: str = None, beta_param: dict = None,
                 collect_ch_sites: bool = True, overwrite_db: bool = False):
        
        self.reference_file = reference_file
        self.reference_dict = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))
        
        self.sim_output = sim_output
        self.beta_param = beta_param
        self.sim_meth_db: StreamOutput = None
        
        self.overwrite_db = overwrite_db
        self.collect_ch_sites = collect_ch_sites # currently use all of the Cs
        self.initialize_methylation_reference(meth_ref, cgmap)
    
    
    def initialize_methylation_reference(self, meth_ref, cgmap):
        """If profile exists then assume already run"""
        self.profile_dict= dict()
        self.sim_meth_db = StreamOutput(sim_output=self.sim_output)
        self.sim_meth_db.generate_sim_directory()
        
        
        # initialize data object using fasta sequence
        for contig_id in self.reference_dict.keys():
            contig_profile = dict()
            contig_seq = self.reference_dict[contig_id].seq
            
            for c_match in tqdm(re.finditer("C", str(self.reference_dict[contig_id].seq))):
                c_pos = c_match.start()
                context = 1 if str(contig_seq[c_pos:(c_pos+2)]) == 'CG' else 0
                contig_profile[c_pos] = np.array([0, context, np.NaN], dtype=float)
            for g_match in tqdm(re.finditer("G", str(self.reference_dict[contig_id].seq))):
                g_pos = g_match.start()
                context = 1 if str(contig_seq[c_pos:(c_pos+2)]) == 'GC' else 0
                contig_profile[g_pos] = np.array([1, context, np.NaN], dtype=float)
            self.profile_dict[contig_id] = contig_profile
        
        
        # intialize methylation value based on provided profile
        if meth_ref:
            self.sim_meth_db_from_dist(self.beta_param) # TODO: load the existing meth_db.pkl file
        elif cgmap:
            self.sim_meth_db_from_cgmap(cgmap)
        else:
            self.sim_meth_db_from_dist(self.beta_param)

    def sim_meth_db_from_dist(beta_param):
        pass
    
    def sim_meth_db_from_cgmap(self, cgmap: str):
        for line in OpenCGmap(cgmap):
            chrom, base, pos, context, methlevel = line[0], line[1], int(line[2]) - 1, line[3], float(line[5])
            context = 1 if context == 'CG' else 0
            nucleotide = 1 if base == 'G' else 0
            self.profile_dict[chrom][pos] = np.array([nucleotide, context, methlevel], dtype=float)
        
        for contig_id in self.profile_dict.keys():
            for meth_pos, meth_profile in self.profile_dict[contig_id].items():
                if np.isnan(meth_profile[2]):
                    meth_profile[2] = self.pick_cytosine_methylation(context = int(meth_profile[1]))
                    self.profile_dict[contig_id][meth_pos] = meth_profile
        for contig_id, contig_profile in self.profile_dict.items():
            self.sim_meth_db.output_contig(contig_id, contig_profile, values=True)

    #@property
    def pick_cytosine_methylation(self, context):
        """Sample from Cytosine distribution, 1 for CG and 0 for CH"""
        if context:
            return np.random.beta(self.beta_param['CG'][0], self.beta_param['CG'][1], size =1)[0]
        return np.random.beta(self.beta_param['CHG'][0], self.beta_param['CHG'][1], size =1)[0]
        
        
    def get_contig_methylation(self, contig):
        contig_profile = self.sim_meth_db.load_contig(contig, values=True)
        if contig_profile:
            return contig_profile

    def set_variant_methylation(self, sim_data, contig_profile, current_contig):
        pass
    # def set_variant_methylation(self, sim_data, contig_profile, current_contig):
    #     """Variants are always random, so set random methylation"""
    #     ref_seq = self.reference_dict[current_contig]
    #     for pos, variant_info in sim_data.items():
    #         # don't set methylation for last base
    #         if pos + 2 > len(ref_seq):
    #             continue
    #         if variant_info['indel'] == -1:
    #             continue
    #         elif variant_info['indel'] == 1:
    #             insert_context = f'{ref_seq[pos-1]}{variant_info["alt"]}{ref_seq[pos]}'
    #             for insert_pos, base in enumerate(variant_info['alt']):
    #                 if base != 'C' and base != 'G':
    #                     continue
    #                 if base == 'C':
    #                     nucleotide = 0
    #                     context = 1 if insert_context[insert_pos+1] == "G" else 0
    #                 if base == "G":
    #                     nucleotide = 1
    #                     context = 1 if insert_context[insert_pos+1] == "C" else 0
                    
    #                 meth_profile = np.array([nucleotide, context, self.pick_cytosine_methylation(context)], dtype=float)
    #                 if meth_profile.all():
    #                     contig_profile[f'{pos}_+_{insert_pos}'] = meth_profile
    #         else:
    #             for base in variant_info['iupac']:
    #                 if base == variant_info['ref']:
    #                     continue
    #                 else:
    #                     context = f'{ref_seq[pos - 1]}{base}{ref_seq[pos + 1]}'
    #                     if base != "C" and base != "G":
    #                         continue
    #                     if base == "C":
    #                         nucleotide = 0
    #                         context = 1 if ref_seq[pos] == "G" else 0
    #                     if base == "G":
    #                         nucleotide = 1
    #                         context = 1 if ref_seq[pos] == "C" else 0
                        
    #                     meth_profile = np.array([nucleotide, context, self.pick_cytosine_methylation(context)], dtype=float)
    #                     if meth_profile.all():
    #                         contig_profile[f'{pos}_{variant_info["ref"]}_{base}'] = meth_profile
    #                     else:
    #                         print([base, pos, variant_info["ref"]])
