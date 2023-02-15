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
    def __init__(self, ref_fasta: str = None, outdir: str = None, 
                 depth: float = None, num_reads: int = None, 
                 read_len: int = None, pair_end: bool = False, is_uniform: bool = True,
                 probe_bed_file: str = None, cut_site_str: str = None, rrbs_model: str = None,
                 insert_min: int = 100, insert_max: int = 1000, overwrite_db: bool = False):
        
        self.ref_dict = SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))              # save
        self.score_dict = {contig_id: 0 for contig_id in self.ref_dict.keys()}      # used to calculate the #reads per chromosome if not is_uniform
        self.eff_len_dict = {contig_id: 0 for contig_id in self.ref_dict.keys()}    # used to calculate the #reads from depth if not specified
        self.chr_len_dict = {contig_id: len(self.ref_dict[contig_id]) for contig_id in self.ref_dict.keys()} # record the genome length


        if probe_bed_file:
            self.tech_mode      = 2
            self.probe_bed_file = probe_bed_file
            self.parse_bed()
        elif cut_site_str:
            self.tech_mode      = 1
            self.cut_site_str   = cut_site_str
            self.insert_min     = insert_min
            self.insert_max     = insert_max
            self.overwrite_db   = overwrite_db
            self.rrbs_model     = rrbs_model
            self.output_bed     = 'f{outdir}/rrbs.bed'
            self.check_file()
            self.cut_genome()
        else:
            self.tech_mode      = 0
            self.score_dict     = self.chr_len_dict
            self.eff_len_dict   = self.chr_len_dict


        if not num_reads:
            num_reads  = int(sum(self.eff_len_dict.values())*depth/read_len/(1+int(pair_end)))
        if not is_uniform:
            tot_score  = sum(self.score_dict.values())
            self.count_dict = {contig_id: num_reads * self.score_dict[contig_id]/tot_score for contig_id in self.ref_dict.keys()}
        else:
            tot_score  = sum(self.eff_len_dict.values())
            self.count_dict = {contig_id: num_reads * self.eff_len_dict[contig_id]/tot_score for contig_id in self.ref_dict.keys()}


    def parse_bed(self):
        with open(self.probe_bed_file, 'r') as f: # bed file regions non-overlap
            for line in f:
                columns = line.strip().split('\t')
                self.score_dict[columns[0]]   += float(columns[4])
                self.eff_len_dict[columns[0]] += int(columns[2]) - int(columns[1])


    def cut_genome(self):
        site_list, spot_list, site_len_list = self.parse_site(self.cut_site_str)
        
        for contig_id in self.ref_dict.keys():
            pos_arr = self.gen_cut_site(contig_id, site_list)
            frag_arr= self.gen_frag_simple(contig_id, pos_arr, site_list, spot_list, site_len_list)
            self.eff_len_dict[contig_id]= np.sum(frag_arr[:,1] - frag_arr[:,0])
            frag_arr2=self.gen_frag_w_cut(frag_arr,site_list)
            
            frag_df = pd.DataFrame(frag_arr2, columns= ['start', 'end', 'cut_l', 'cut_r', 'cg_count'] + site_list)
            frag_df.loc[:, 'chr_id'] = contig_id
            frag_df.loc[:, 'name']   = "."
            frag_df.loc[:, 'strand'] = "."
            frag_df.loc[:, 'length'] = frag_df.loc[:,'end'] - frag_df.loc[:,'start']
            frag_df.loc[:, 'ratio']  = frag_df.loc[:,'cg_count']/frag_df.loc[:,'length']
            frag_df.loc[:, 'score']  = self.cal_score(frag_df)
            
            self.score_dict[contig_id] += np.sum(frag_df.loc[:, 'score'])
            
            frag_df.to_csv(self.output_bed, columns=['chr_id', 'start', 'end', 'name', 'score', 'strand'], 
                           sep = '\t', mode='a', header=False, index=False, quoting=None)


    def gen_cut_site(self, contig_id, site_list):
        # get all the cutting sites and sort
        arr_list = []
        contig_seq = self.ref_dict[contig_id].seq.upper()
        
        for ix, site in enumerate(site_list):
            arr_list.append(np.array([[match.start(), ix] for match in re.finditer(site, str(contig_seq))], dtype=np.int32))
        arr_list.append(np.array([[len(contig_seq), 0]]))
        
        pos_arr = np.concatenate(arr_list, axis=0)
        return pos_arr[pos_arr[:, 0].argsort()]


    def gen_frag_simple(self, contig_id, pos_arr, site_list, spot_list, site_len_list):
        # create fragment without cut sites in the middle
        pos_num  = pos_arr.shape[0]
        frag_arr = np.zeros((pos_num, 5 + len(site_list)), dtype=np.int32) # pos_left, pos_right, left_cut, right_cut, cg_count
        contig_seq = self.ref_dict[contig_id].seq.upper()
        
        ix = 0
        for i in range(pos_num):
            cut_l = 0 if i==0 else pos_arr[i-1,1]
            cut_r = pos_arr[i, 1]
            cut_dict = {ix:0 for ix in range(len(site_list))}
            cut_dict[cut_l] += 1
            cut_dict[cut_r] += 1

            pos_l = 0 if i==0 else pos_arr[i-1,0] + spot_list[cut_l]
            pos_r = pos_arr[i,0] + site_len_list[cut_r] - spot_list[cut_r]

            length= pos_r - pos_l
            if length >= self.insert_min and length <= self.insert_max:
                cg_count= (contig_seq[pos_l:pos_r].count('C') + contig_seq[pos_l:pos_r].count('G'))
                frag_arr[ix,:] = [pos_l, pos_r, cut_l, cut_r, cg_count] + list(cut_dict.values())
                ix += 1
        
        return frag_arr[:ix,:]


    def gen_frag_w_cut(self, frag_arr, site_list):
        # create frag with cut sites in the middle
        frag_list = []
        frag_num  = frag_arr.shape[0]
        for i in range(frag_num):
            pos_l, cut_l, cut_r, cg_count = frag_arr[i,[0,2,3,4]]
            cut_dict = {ix:0 for ix in range(len(site_list))}
            cut_dict[cut_l] += 1
            cut_dict[cut_r] += 1

            for j in range(i+1, frag_num):
                pos_r, cut_r, frag_cg_count = frag_arr[j,[1,3,4]]
                if pos_r - pos_l > self.insert_max:
                    break

                cut_dict[cut_r] += 1
                cg_count += frag_cg_count
                frag_list.append([pos_l, pos_r, cut_l, cut_r, cg_count] + list(cut_dict.values()))
        
        frag_arr2= np.concatenate([frag_arr, np.array(frag_list, dtype=np.int32)], axis=0)
        frag_arr2= frag_arr2[np.lexsort((frag_arr2[:,1], frag_arr2[:,0]))]
        return frag_arr2


    def cal_score(self, model_file, df):
        # load model
        rrbs_model = f'{os.path.dirname(get_htsim_path())}/data/model/rrbs.model' if not rrbs_model else rrbs_model
        
        model = pickle.load(open(model_file, 'rb'))
        model_vars = []
        
        # check if the model and the columns fit
        if set(model_vars) != set(df.columns):
            raise ValueError('The trained model is not compatible with input data (cut sites maybe different), please check!')
        return model.predict(df.loc[:, model_vars])


    def parse_site(self, cut_site_str):
        cut_site_str_split= cut_site_str.split(",")
        spot_list = [site.rfind("|") for site in cut_site_str_split]
        site_list = [site.replace("|", "" ).upper() for site in cut_site_str_split]
        site_len_list = [len(site) for site in site_list]
        
        return [site_list, spot_list, site_len_list]
    
    
    def check_file(self):
        if os.path.exists(self.output_bed):
            if self.overwrite_db:
                os.remove(self.output_bed)
            else:
                raise ValueError('Cannot write to output file, please set overwrite_db to be true!')
            
