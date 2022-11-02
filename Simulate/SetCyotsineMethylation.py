import os
import re
import warnings
import pandas as pd
import numpy as np

from Bio import SeqIO
from scipy.stats import beta
from tqdm import tqdm
from typing import Dict, List

from StreamOutput import StreamOutput
from UtilityFunctions import parseCGmap, parseASM


class SetCytosineMethylation:
    '''
    prepare Cytosine methylation level database for read methylation state simulation
    :param str  ref_fasta    : path to the reference genome [.fasta/.fa/.fa.gz]
    :param str  outdir       : path to the output directory
    :param str  meth_db_path : path to the previous obj [.pkl] (if use meth_db for simulation)
    :param str  cgmap_file   : path to the cgmap file [.CGmap/.CGmap.gz]
    :param bool cgmap_pool   : whether to pool methylaiton levels and random draws from it
    :param bool asm_sim      : whether to conduct allelic specific methylation (ASM)
    :param str  asm_file     : path to the allelic specific methylation (ASM) file
    :param dict beta_param   : dict of beta parameters for CG/CHG/CHH methylation values simulation
    :param bool collect_ch   : whether to collect the CHG/CHH sites, slow speed is expected
    :param bool overwrite_db : whether to overwrite the meth_db if it already exists
    :param bool verbose      : whether to output processing details
    :rtype None

    :var array meth_arr : nx5 numpy array: context, flag, meth_avg, meth_ref, meth_alt
    (flag: 0 from dist or from cgmap pool, 1 from CGmap, 2 from ASM, -1 for uninitialized)
    :var series pos_map : pd.Series of length n: index is genome coordinate, value is row idx of meth_arr
    '''


    def __init__(self, ref_fasta: str = None, outdir: str = None,
                 meth_db_path: str = None, cgmap_file: str = None, cgmap_pool: bool = False,
                 beta_param: dict = {"CG": (0.5, 0.5), "CHG": (0.01, 0.05), "CHH":(0.01, 0.05)},
                 asm_sim: bool = False, asm_file: str = None,
                 collect_ch: bool = True, overwrite_db: bool = False, verbose: bool = False):

        self.ref_fasta   = ref_fasta
        self.outdir      = outdir
        self.meth_db_path= meth_db_path
        self.cgmap_file  = cgmap_file
        self.cgmap_pool  = cgmap_pool
        self.beta_param  = beta_param
        self.asm_sim     = asm_sim
        self.asm_file    = asm_file
        self.collect_ch  = collect_ch
        self.overwrite_db= overwrite_db
        self.verbose     = verbose
        self.meth_arr    = None
        self.pos_map     = None

        self.context_dict      = {1:'CG', 3:'CHG', 7:'CHH', 9:'CG', 11:'CHG', 15:'CHH'}
        self.base_context_dict = {'C': {'CG':1, 'CHG':3, 'CHH':7}, 'G': {'CG':9, 'CHG':11, 'CHH':15}}
        self.base_context_table= {'C': np.array([[7,3],  [1,1]]),  'G': np.array([[15,11],[9,9]])}

        # check existence
        if not os.path.exists(ref_fasta):
            raise ValueError('Cannot find the reference genome, please check!')
        if not outdir:
            raise ValueError('Please specify the output directory!')
        if asm_sim and not os.path.exists(asm_file):
            raise ValueError('Please specify allelic specific methylation file correctly for ASM simulation!')

        self.meth_db = StreamOutput(outdir=self.outdir)
        self.meth_db.check_outdir()
        self.ref_dict= SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))

        # initiate
        if meth_db_path:
            self.check_meth_db()
        else:
            self.create_meth_db()


    def check_meth_db(self):
        '''if use previous meth_db object, check if they align'''
        contig_id_list = self.ref_dict.keys()
        for contig_id in contig_id_list:
            contig_profile = self.meth_db.load_contig(contig_id, values=True)
            profile_sites  = contig_profile[0].index()
            self.init_meth_db(contig_id)
            not_comp_sites = self.pos_map.index().difference(profile_sites)
            if len(not_comp_sites):
                print(f'{contig_id}: sites in the meth_db is not the same as the reference')


    def create_meth_db(self):
        '''parse the fasta, CGmap, ASM file into pickle for later simulation'''
        contig_id_list = self.ref_dict.keys()
        for contig_id in contig_id_list:
            self.init_meth_db(contig_id)
            self.fill_cgmap(contig_id)
            self.fill_asm(contig_id)
            self.fill_dist()
            #self.meth_db.output_contig(contig_id, [self.pos_map, self.meth_arr], values=True)


    def init_meth_db(self, contig_id):
        '''initialize data object using fasta sequence'''
        idx = 0
        seq = self.ref_dict[contig_id].seq.upper()
        seq_len = len(seq)

        if self.collect_ch:
            count_c = seq.count("C")
            count_g = seq.count("G")
            arr_size= count_c + count_g
            self.meth_arr = np.full((arr_size, 5), fill_value=np.NaN, dtype=np.float16)
            self.meth_arr[:,1] = -1 # flag
            self.pos_map  = pd.Series(0, index=range(arr_size), dtype=np.uint32)

            for pos, base in tqdm(enumerate(seq), disable = not self.verbose):
                if base not in {"C", "G"}:
                    continue

                if pos<2 or pos>(seq_len-2):
                    base_d1 = base_d2 = 0
                elif base == "C":
                    base_d1 = int(seq[pos+1] == "G")
                    base_d2 = int(seq[pos+2] == "G")
                else:
                    base_d1 = int(seq[pos-1] == "C")
                    base_d2 = int(seq[pos-2] == "C")
                # C:{10, 11}: 1, {01}: 3, {00}: 7; G: {10, 11}: 9, {01}: 11, {00}: 15
                self.meth_arr[idx, 0] = self.base_context_table[base][base_d1, base_d2] # context
                self.pos_map[idx] = pos
                idx += 1
        else:
            count_c = count_g = seq.count("CG")
            arr_size= count_c + count_g
            self.meth_arr= np.full((arr_size, 5), fill_value=np.NaN, dtype=np.float16)
            self.meth_arr[:,1]= -1
            self.pos_map = pd.Series(0, index=range(arr_size), dtype=np.uint32)

            for cg_match in tqdm(re.finditer("CG", str(seq)), disable = not self.verbose):
                pos = cg_match.start()
                self.meth_arr[idx,  0] = [pos,  1]
                self.meth_arr[idx+1,0] = [pos+1,9]
                self.pos_map[idx]  = pos
                self.pos_map[idx+1]= pos+1
                idx += 2
        self.pos_map = pd.Series(self.pos_map.index.values, index=self.pos_map)


    def fill_cgmap(self, contig_id):
        '''fill sites in the cgmap file'''
        if not self.cgmap_file:
            return None

        if self.cgmap_pool:
            ctx_idx_dict = {'CG':1, 'CHG':2, 'CHH':3}
            meth_level_pool = [[],[],[]]
            for line in parseCGmap(self.cgmap_file, contig_id, collect_ch = self.collect_ch):
                _, base, pos, context, meth_level = line
                meth_level_pool[ctx_idx_dict[context]].append(float(meth_level))

            idx_cg  = np.where(np.bitwise_and(self.meth_arr[:, 1].astype(np.int16), 0x7) == 1)
            idx_chg = np.where(np.bitwise_and(self.meth_arr[:, 1].astype(np.int16), 0x7) == 3)
            idx_chh = np.where(np.bitwise_and(self.meth_arr[:, 1].astype(np.int16), 0x7) == 7)

            self.meth_arr[idx_cg, 2] = np.random.choice(meth_level_pool[0], size = np.sum(idx_cg))
            self.meth_arr[idx_chg,2] = np.random.choice(meth_level_pool[1], size = np.sum(idx_chg))
            self.meth_arr[idx_chh,2] = np.random.choice(meth_level_pool[2], size = np.sum(idx_chh))
            return None

        num_cgmap_pos, num_404_pos = [0,0]  # counting
        for line in parseCGmap(self.cgmap_file, contig_id, collect_ch = self.collect_ch):
            num_cgmap_pos += 1
            _, base, pos, context, meth_level = line
            try:
                arr_idx = self.pos_map[int(pos) - 1] # CGmap file is 1-based coordinate
            except IndexError:
                num_404_pos += 1
            else:
                if (self.meth_arr[arr_idx, 0] != self.base_context_dict[base][context]):
                    num_404_pos +=1
                    continue
                self.meth_arr[arr_idx, 1:2] = [1, float(meth_level)]

        ratio_404 = round(num_404_pos / num_cgmap_pos, 4)
        ### also print what context
        print(f"Contig {contig_id}: {num_cgmap_pos} sites found in CGmap file, among them {num_404_pos} sites ({ratio_404 * 100}%) are not compatible...")
        if (ratio_404 >=0.5):
            warnings.warn("[WARNING]: More than half sites in CGmap file cannot be found in the reference fasta\n" +
                          "Potential reason: the CGmap does not share the same coordinates with fasta, please check!")
        return None


    def fill_asm(self, contig_id):
        '''fill the allelic specific sites'''
        if not self.asm_sim:
            return None

        num_asm_pos, num_404_pos = [0,0] # counting
        for line in parseASM(self.asm_file, contig_id, collect_ch = self.collect_ch):
            _, base, pos, context, tot_meth, ref_meth, alt_meth = line
            num_asm_pos += 1
            try:
                arr_idx = self.pos_map[int(pos) - 1]
            except IndexError:
                num_404_pos += 1
            else:
                if (self.meth_arr[arr_idx, 0] != self.base_context_dict[base][context]):
                    num_404_pos +=1
                    continue
                self.meth_arr[arr_idx, 1:4] = [2, tot_meth, ref_meth, alt_meth]
        ratio_404 = round(num_404_pos / num_asm_pos, 4)
        print(f"Contig {contig_id}: {num_asm_pos} sites found in CGmap file, among them {num_404_pos} sites ({ratio_404 * 100}%) are not compatible...")
        if (ratio_404 >=0.5):
            warnings.warn("[WARNING]: More than half sites in ASM file cannot be found in the reference fasta\n" +
                          "Potential reason: the CGmap does not share the same coordinates with fasta, please check!")
        return None


    def fill_dist(self):
        '''fill in from beta distribution'''
        idx_nan = np.isnan(self.meth_arr[:,2]) # hold for each element
        idx_ctx = np.bitwise_and(self.meth_arr[:,0].astype(np.int16), 0x7) # hold for each element
        idx_nan_cg  = np.where((idx_ctx == 1) & idx_nan)[0] # hold for index

        if len(idx_nan_cg):
            self.meth_arr[idx_nan_cg, 2] = self.simu_beta_dist(context="CG", size=len(idx_nan_cg))

        if self.collect_ch:
            idx_nan_chg = np.where((idx_ctx == 3) & idx_nan)[0]
            idx_nan_chh = np.where((idx_ctx == 7) & idx_nan)[0]
            if len(idx_nan_chg):
                self.meth_arr[idx_nan_chg, 2]= self.simu_beta_dist(context="CHG", size=len(idx_nan_chg))
            if len(idx_nan_chh):
                self.meth_arr[idx_nan_chh, 2]= self.simu_beta_dist(context="CHH", size=len(idx_nan_chh))

        self.meth_arr[idx_nan, 1] = 0
        self.meth_arr[idx_nan, 3] = self.meth_arr[idx_nan, 2]
        self.meth_arr[idx_nan, 4] = self.meth_arr[idx_nan, 2]


    def set_var_meth(self, sim_data, contig_id) -> Dict(str, List):
        '''set random methylation due to variants are random'''
        var_meth_dict = {}

        seq = self.ref_dict[contig_id].seq.upper()
        seq_len = len(seq)
        for pos, variant_info in sim_data.items():
            if pos<2 or pos>(seq_len-2):
                continue
            if variant_info['indel'] == -1:  # deletion
                continue
            if variant_info['indel'] == 1: # insertion
                local_seq = f'{seq[(pos-1):(pos+1)]}{variant_info["alt"]}{seq[(pos+1):(pos+3)]}'
                for insert_pos, base in enumerate(variant_info['alt']):
                    if base not in {'C', 'G'}:
                        continue
                    if base == 'C':
                        base_d1 = int(local_seq[insert_pos+1+2] == "G")
                        base_d2 = int(local_seq[insert_pos+2+2] == "G")
                    else:
                        base_d1 = int(local_seq[insert_pos-1+2] == "C")
                        base_d2 = int(local_seq[insert_pos-2+2] == "C")
                    context     = self.base_context_table[base][base_d1, base_d2]
                    meth_level  = self.simu_beta_dist(context=context)[0]
                    var_meth_dict[f'{pos}+{insert_pos}'] = [base, 3, context, meth_level]
            else:# substitution
                for base in variant_info['iupac']:
                    if base not in {'C', 'G'}:
                        continue
                    if base == variant_info['ref']:
                        continue
                    local_seq = f'{seq[(pos-2):pos]}{variant_info["alt"]}{seq[(pos+1):(pos+3)]}'
                    if base == "C":
                        base_d1 = int(local_seq[1+2] == "G")
                        base_d2 = int(local_seq[2+2] == "G")
                    else:
                        base_d1 = int(local_seq[-1+2] == "C")
                        base_d2 = int(local_seq[-2+2] == "C")
                    context     = self.base_context_table[base][base_d1, base_d2]
                    meth_level  = self.simu_beta_dist(context=context)[0]
                    var_meth_dict[f"{pos}_{variant_info['ref']}_{variant_info['alt']}"] = [base, 1, context, meth_level]

        self.meth_db.output_contig(contig_id, var_meth_dict, is_variant=True)


    def simu_beta_dist(self, context = "CG", size = 1):
        '''output the values accordig to the context using beta distribution'''
        if isinstance(context, int):
            context = self.context_dict[context]
        return beta.rvs(a=self.beta_param[context][0], b=self.beta_param[context][1], size=size).astype(np.float16)


    def estimate_beta_params(self, context = None):
        '''estimate the beta parameters for each context'''
        pass
