import os
import re
import numpy as np
import pandas as pd
import pickle

from Bio import SeqIO
from typing import Dict, List
from UtilityFunctions import get_htsim_path


class ParseGenome:
    '''
    parse genome fasta into dict, calculate the effective length and score if needed
    # parameters
    :param str   ref_fasta      : path to the reference genome (.fasta/.fa/.fa.gz) [None]
    :param str   outdir         : path to the output directory [None]
    :param str   probe_bed_file : probe bed file for TBS simulation [None]
    :param str   cut_site_str       : enzyme cut site for RRBS, cutting position is denoted by |, multiple sites
                                  are separated by ,, for example MspI and TaqI: "C|CGG,T|CGA" [None]
                                  
    '''

    def __init__(self, ref_fasta: str = None, outdir: str = None, depth: float = None, num_reads: int = None, 
                 read_len: int = None, pair_end: bool = False, insert_min: int = 100, insert_max: int = 1000, 
                 cut_site_str: str = None, rrbs_model: str = None, middle_cuts: bool = False, 
                 probe_bed_file: str = None, is_uniform: bool = True, overwrite_db: bool = False):

        self.ref_dict = SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))              # can be saved
        self.score_dict = {contig_id: 0 for contig_id in self.ref_dict.keys()}      # calculate #reads per contig if not is_uniform
        self.eff_len_dict = {contig_id: 0 for contig_id in self.ref_dict.keys()}    # calculate #reads from depth if not specified
        self.chr_len_dict = {contig_id: len(self.ref_dict[contig_id]) for contig_id in self.ref_dict.keys()} # record contig length
        self.is_uniform   = is_uniform
        
        if probe_bed_file:
            self.tech_mode      = 2
            self.bed_file       = probe_bed_file
            self.parse_bed()
        elif cut_site_str:
            self.tech_mode      = 1
            self.cut_site_str   = cut_site_str
            self.insert_min     = insert_min
            self.insert_max     = insert_max
            self.middle_cuts    = middle_cuts
            self.overwrite_db   = overwrite_db
            self.rrbs_model     = rrbs_model
            self.bed_file       = f'{outdir}/tmp/rrbs.bed'
            self.new_bed()
            self.cut_genome()
        else:
            self.tech_mode      = 0
            self.score_dict     = self.chr_len_dict
            self.eff_len_dict   = self.chr_len_dict
            self.bed_file       = None

        # should check depth, read_len, pair_end
        self.num_reads  = num_reads if num_reads else int(sum(self.eff_len_dict.values())*depth/read_len/(1+int(pair_end)))
        norm_score_dict = {contig_id: self.score_dict[contig_id]/sum(self.score_dict.values()) for contig_id in self.ref_dict.keys()}
        self.count_dict = {contig_id: int(self.num_reads*norm_score_dict[contig_id]/(1+int(pair_end)))for contig_id in self.ref_dict.keys()}
        # the above could have some rounding errors, 
        num_reads_rest  = self.num_reads/(1+int(pair_end)) - sum(self.count_dict.values())
        while num_reads_rest:
            contig_list = list(self.count_dict.keys())
            self.count_dict[contig_list[int(num_reads_rest % len(contig_list))]] += 1
            num_reads_rest = num_reads_rest -1


    def parse_bed(self):
        with open(self.bed_file, 'r') as f: # bed file regions non-overlap
            for line in f:
                columns = line.strip().split('\t')
                length  = int(columns[2]) - int(columns[1])
                self.eff_len_dict[columns[0]]+= length
                self.score_dict[columns[0]]  += length if self.is_uniform else float(columns[4]) # what if score is .


    def cut_genome(self):
        site_dict, spot_dict, site_len_dict = self.parse_site(self.cut_site_str)

        for contig_id in self.ref_dict.keys():
            pos_arr = self.gen_cut_site(contig_id, site_dict)
            frag_arr, contig_eff_length = self.gen_cut_frag(contig_id, pos_arr, site_dict, spot_dict, site_len_dict)
            self.eff_len_dict[contig_id]= contig_eff_length
            
            frag_df = pd.DataFrame(frag_arr, columns= ['start', 'end', 'cut_l', 'cut_r', 'cg_count'] + list(site_dict.values()))
            frag_df.loc[:, 'chr_id'] = contig_id
            frag_df.loc[:, 'name']   = "."
            frag_df.loc[:, 'strand'] = "."
            frag_df.loc[:, 'length'] = frag_df.loc[:,'end'] - frag_df.loc[:,'start']
            frag_df.loc[:, 'ratio']  = frag_df.loc[:,'cg_count']/frag_df.loc[:,'length']
            frag_df.loc[:, 'score']  = 1 if self.is_uniform else self.cal_score(frag_df)
            
            self.score_dict[contig_id] = self.eff_len_dict[contig_id] if self.is_uniform else np.sum(frag_df.loc[:, 'score'])
            
            frag_df.to_csv(self.bed_file, columns=['chr_id', 'start', 'end', 'name', 'score', 'strand'], 
                           sep = '\t', mode='a', header=False, index=False, quoting=None)


    def gen_cut_site(self, contig_id, site_dict):
        # get all the cutting sites and sort
        arr_list = []
        contig_seq = self.ref_dict[contig_id].seq.upper()
        
        for ix, site in site_dict.items():
            arr_list.append(np.array([[match.start(), ix] for match in re.finditer(site, str(contig_seq))], dtype=np.int32))
        arr_list.append(np.array([[0, -1], [len(contig_seq), -1]])) # -1 mark the boundary
        
        pos_arr = np.concatenate(arr_list, axis=0)
        return pos_arr[pos_arr[:, 0].argsort()] # np.array([position, type])


    def gen_cut_frag(self, contig_id, pos_arr, site_dict, spot_dict, site_len_dict):
        # create fragments from the cutting positions
        num_pos = pos_arr.shape[0]
        frag_list  = [] # each fragment: [pos_left, pos_right, left_cut, right_cut, cg_count] + site_list
        eff_length = 0
        contig_seq = self.ref_dict[contig_id].seq.upper()
        
        for i in range(num_pos):
            cut_l = pos_arr[i,1]
            pos_l = pos_arr[i,0] + spot_dict[cut_l]
            
            cut_dict = {ix:0 for ix in site_dict.keys()}
            if cut_l >=0:
                cut_dict[cut_l] += 1
            
            right_bound = num_pos if self.middle_cuts else min([i+2, num_pos])
            
            for j in range(i+1, right_bound):
                cut_r = pos_arr[j, 1]
                pos_r = pos_arr[j, 0] + site_len_dict[cut_r] - spot_dict[cut_r]
                length= pos_r - pos_l
                
                if length > self.insert_max:
                     break
                
                if cut_r >= 0:
                    cut_dict[cut_r] += 1
                
                if length < self.insert_min:
                    continue
                
                if j-i == 1:
                    eff_length += length
                
                # has duplicate computing, can be optimized
                cg_count = contig_seq[pos_l:pos_r].count('C') + contig_seq[pos_l:pos_r].count('G')
                frag_list.append([pos_l, pos_r, cut_l, cut_r, cg_count] + list(cut_dict.values()))

        frag_arr= np.array(frag_list, dtype=np.int32)
        return frag_arr[np.lexsort((frag_arr[:,1], frag_arr[:,0]))], eff_length


    def cal_score(self, df):
        # load model: TODO
        model, model_vars = pickle.load(open(self.rrbs_model, 'rb'))
        
        # check if the model and the columns fit
        if len(set(model_vars).intersection(set(df.columns))):
            raise ValueError('The trained model is not compatible with input data (cut sites maybe different), please check!')
        return model.predict(df.loc[:, model_vars])


    def parse_site(self, cut_site_str):
        cut_site_str_split= cut_site_str.split(",")
        spot_dict = {-1:0}
        site_dict = {}
        site_len_dict = {-1:0}
        
        for ix, site in enumerate(cut_site_str_split):
            spot_dict[ix] = site.rfind("|")
            site_dict[ix] = site.replace("|", "" ).upper()
            site_len_dict[ix] = len(site_dict[ix])
        
        return [site_dict, spot_dict, site_len_dict]


    def new_bed(self):
        if os.path.exists(self.bed_file):
            if not self.overwrite_db:
                raise ValueError('Cannot write to output file (file already exists), please set overwrite_db to be true!')
            else:
                os.remove(self.bed_file)

