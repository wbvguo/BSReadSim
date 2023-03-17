import os
import re
import numpy as np
import pandas as pd
import pickle

from UtilityFunctions import get_htsim_path
from typing import Dict, List

class SetExperiment:
    '''
    determines the tech mode, and hold other parameters
    '''
    
    def __init__(self, ref_dict: Dict = None, outdir: str = None,                                   # required arguments
                 collect_ch: bool = True,
                 overwrite_bed: bool = True, asm_file: str = None,
                 vcf_file: str = None, middle_cuts: bool = False,                                   # genotype reference
                 mut_rate: float = 0.0010, haplo_mode: bool = False, seed: int = -1,                # mutation parameter
                 mut_indel_frac: float = 0.15, indel_ext_prob: float = 0.15,                        # indel parameter
                 insert_mean: int = 400, insert_std: int = 25,                                      # fragment setting
                 insert_min: int = 100, insert_max: int = 1000,
                 read_len: int = 100, depth: int = 20, num_reads: int = None,                       # reads setting
                 propN_cutoff: float = 0.05, conversion_rate: float = 0.998,                        
                 undirectional: bool = False, pair_end: bool = True,                                # sequencing protocol
                 random_error= True, error_rate: float = 0.005, error_model: str = None,                # sequencing error TODO:
                 qual_uniform: bool =True, qual_model: str = None,                                  # quality model, TODO:
                 is_uniform: bool = True,
                 cover_uniform: bool =True, gc_bias_file: str = None,                               # coverage model, gc_bias_file TODO:
                 site_dependency: bool = False, site_model: str = None,                             # site dependency model TODO:
                 cut_site_str: str = None, rrbs_model: str = None, rrbs_uniform: bool = True,       # RRBS technology model TODO:
                 probe_bed_file: str = None, probe_deviation: int =None,                            # TBS technology file
                 collect_stats: bool = False):
        
        
        if cut_site_str and probe_bed_file:
            raise ValueError('Please specify files for only one technology mode!')
        
        # sequencing settings in python module
        self.ref_dict       = ref_dict
        self.outdir         = outdir
        self.pkl_dir        = f'{self.outdir}/pkl/'
        self.tmp_dir        = f'{self.outdir}/tmp/'
        self.read_len       = read_len
        self.pair_end       = pair_end
        self.undirectional  = undirectional
        self.conversion_rate= conversion_rate
        self.collect_ch     = collect_ch
        self.collect_stats  = collect_stats
        self.overwrite_bed  = overwrite_bed
        self.is_uniform     = is_uniform
        self.asm_sim        = True if asm_file else False
        self.middle_cuts    = middle_cuts
        self.cut_site_str   = cut_site_str
        self.bed_file       = probe_bed_file
        
        # #reads to genereate
        self.num_reads      = num_reads
        self.depth          = depth if not self.num_reads else np.nan
        
        # sequencing error
        self.random_error   = False if error_model else random_error
        self.error_model    = error_model
        self.error_rate     = error_rate if self.random_error else np.nan
        
        # site-site dependency
        self.site_dependency= True if site_model else site_dependency
        self.site_model     = site_model
        
        # uniform quality
        self.qual_uniform   = False if qual_model else qual_uniform 
        self.qual_model     = qual_model
        
        # uniform coverage
        self.cover_uniform  = False if gc_bias_file else cover_uniform
        self.gc_bias_file   = gc_bias_file
        
        # rrbs 
        self.rrbs_uniform   = False if rrbs_model else rrbs_uniform
        self.rrbs_model     = rrbs_model
        
        self.create_outdir()# create folder if not exist
        self.check_file()   # use default file if not specified
        self.cal_count()
        
        
        self.htsim_path     = get_htsim_path()
        self.model_path     = f'{os.path.dirname(os.path.dirname(self.htsim_path))}/data/model'
        
        if self.is_uniform:
            self.gc_bias_file=None
            self.error_model= None
            self.qual_model = None
            self.rrbs_model = None
            self.site_model = None
        
        # prepare htsim command
        self.htsim_opts = { '-i': insert_mean, '-I': insert_std, '-m': insert_min, '-M': insert_max,
                            '-1': read_len, '-2': read_len, '-e': 0, '-A': propN_cutoff, '-u': int(is_uniform), '-f': 1,
                            '-g': vcf_file, '-r': mut_rate, '-R': mut_indel_frac, '-X': indel_ext_prob, 
                            '-h': int(haplo_mode), '-s': seed, '-T': self.tech_mode, '-x': cut_site_str, 
                            '-b': self.bed_file, '-B': self.gc_bias_file, '-D': probe_deviation }


    def check_file(self):
        if not self.cover_uniform and not self.gc_bias_file:
            self.gc_bias_file = f'{self.model_path}/gc_bias.txt'
        if self.site_dependency and not self.site_model:
            self.site_model   = f'{self.model_path}/site.pickle'
        if not self.qual_uniform and not self.qual_model:
            self.qual_model   = f'{self.model_path}/qual.pickle'
        if not self.random_error and not  self.error_model:
            self.error_model  = f'{self.model_path}/error.pickle'


    def cal_count(self):
        '''
        parse genome, and calculate the #reads per contig
        '''
        self.count_dict = {contig_id: 0 for contig_id in self.ref_dict.keys()}      # store the #reads per contig
        self.score_dict = {contig_id: 0 for contig_id in self.ref_dict.keys()}      # store contig scores if not is_uniform
        self.eff_len_dict = {contig_id: 0 for contig_id in self.ref_dict.keys()}    # store contig eff_len if using depth
        
        if self.bed_file:
            self.tech_mode      = 2
            self.parse_bed()
        elif self.cut_site_str:
            self.tech_mode      = 1
            self.bed_file       = f'{self.tmp_dir}/rrbs.bed'
            self.cut_genome()
        else:
            self.tech_mode      = 0
            self.score_dict     = {contig_id: len(self.ref_dict[contig_id]) for contig_id in self.ref_dict.keys()} 
            self.eff_len_dict   = {contig_id: len(self.ref_dict[contig_id]) for contig_id in self.ref_dict.keys()}
            self.bed_file       = None
        
        # should check depth, read_len, pair_end
        if not self.num_reads:
            self.num_reads = int(sum(self.eff_len_dict.values())*self.depth/self.read_len/(1+int(self.pair_end)))
        
        contig_list = list(self.count_dict.keys())
        tot_score = sum(self.score_dict.values())
        for contig_id in contig_list:
            self.count_dict[contig_id] = int(self.num_reads * self.score_dict[contig_id]/tot_score)
        
        # the above is rounded, which makes the sum not the same as desired
        num_reads_rest  = self.num_reads - sum(self.count_dict.values())
        while num_reads_rest:    
            self.count_dict[contig_list[int(num_reads_rest % len(contig_list))]] += 1
            num_reads_rest = num_reads_rest -1


    def parse_bed(self):
        with open(self.bed_file, 'r') as FILE: # bed file regions non-overlap
            for line in FILE:
                columns = line.strip().split('\t')
                length  = int(columns[2]) - int(columns[1])
                self.eff_len_dict[columns[0]]+= length
                self.score_dict[columns[0]]  += length if self.is_uniform else float(columns[4]) # what if score is .


    def cut_genome(self):
        if os.path.exists(self.bed_file):
            if not self.overwrite_bed:
                raise ValueError('Cannot write (file already exists), please set overwrite_bed to be true!')
            os.remove(self.bed_file)
        
        
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
        '''
        get all the cutting sites and sort based on the reference genome
        rtype: np.array([position, type])
        '''
        arr_list = []
        contig_seq = self.ref_dict[contig_id].seq.upper()
        
        for ix, site in site_dict.items():
            arr_list.append(np.array([[match.start(), ix] for match in re.finditer(site, str(contig_seq))], dtype=np.int32))
        arr_list.append(np.array([[0, -1], [len(contig_seq), -1]])) # -1 mark the boundary
        
        pos_arr = np.concatenate(arr_list, axis=0)
        return pos_arr[pos_arr[:, 0].argsort()]


    def gen_cut_frag(self, contig_id, pos_arr, site_dict, spot_dict, site_len_dict):
        insert_min = self.htsim_opts['-m']
        insert_max = self.htsim_opts['-M']
        
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
                
                if length > insert_max:
                     break
                if cut_r >= 0:
                    cut_dict[cut_r] += 1
                if length < insert_min:
                    continue
                if j-i == 1:
                    eff_length += length
                
                # has duplicate computing, can be optimized
                cg_count = contig_seq[pos_l:pos_r].count('C') + contig_seq[pos_l:pos_r].count('G')
                frag_list.append([pos_l, pos_r, cut_l, cut_r, cg_count] + list(cut_dict.values()))
        
        frag_arr= np.array(frag_list, dtype=np.int32)
        return frag_arr[np.lexsort((frag_arr[:,1], frag_arr[:,0]))], eff_length


    def cal_score(self, frag_df: pd.DataFrame = None):
        # load model: TODO:
        model, model_vars = pickle.load(open(self.rrbs_model, 'rb'))
        
        # check if the model and the columns fit
        if len(set(model_vars).intersection(set(frag_df.columns))):
            raise ValueError('The RRBS model is not compatible with input data, please check (cut sites maybe different)!')
        return model.predict(frag_df.loc[:, model_vars])


    def create_outdir(self):
        '''create output directory'''
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir, exist_ok=False)
        if not os.path.isdir(self.pkl_dir):
            os.makedirs(self.pkl_dir, exist_ok=False)
        if not os.path.isdir(self.tmp_dir):
            os.makedirs(self.tmp_dir, exist_ok=False)


    @staticmethod
    def parse_site(cut_site_str: str = None) -> List:
        cut_site_str_split= cut_site_str.split(",")
        spot_dict = {-1:0}
        site_dict = {}
        site_len_dict = {-1:0}
        
        for ix, site in enumerate(cut_site_str_split):
            spot_dict[ix] = site.rfind("|")
            site_dict[ix] = site.replace("|", "" ).upper()
            site_len_dict[ix] = len(site_dict[ix])
        
        return [site_dict, spot_dict, site_len_dict]

