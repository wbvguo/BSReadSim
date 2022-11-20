import os
import re
import warnings
import pandas as pd
import numpy as np

from Bio import SeqIO
from scipy.stats import beta
from tqdm import tqdm

from StreamMethDB import StreamMethDB
from UtilityFunctions import parseCGmap, parseASM


class SetMethylation:
    '''
    prepare methylation level database for read methylation state simulation
    :param str  ref_fasta   : path to the reference genome (.fasta/.fa/.fa.gz) [None]
    :param str  outdir      : path to the output directory [None]
    :param str  meth_db_path: path to the previous obj (.pkl) [None]
    :param str  cgmap_file  : path to the cgmap file (.CGmap/.CGmap.gz) [False]
    :param bool cgmap_pool  : whether to pool methylaiton levels and random draws from it [False]
    :param str  asm_file    : path to the allelic specific methylation (ASM) file [None]
    :param dict beta_params : dict of beta parameters for CG/CHG/CHH methylation values simulation
                              [{"CG": (0.5, 0.5), "CHG": (0.01, 0.05), "CHH":(0.01, 0.05)}]
    :param bool collect_ch  : whether to collect/simulate the CHG/CHH sites [False]
    :param bool overwrite_db: whether to overwrite the meth_db if it already exists [True]
    :param int  seed        : seed for random methylation generation
    :param bool verbose     : whether to output processing details [True]
    :rtype None

    :var np.array meth_arr  : nx5 numpy array: context, flag, meth_avg, meth_ref, meth_alt
                              - flag: 0 from dist or pool, 1 from CGmap, 2 from ASM, -1 unintialized
    :var pd.Series pos_map  : pd.Series of length n
                              - index is genome coordinate, value is row idx of meth_arr
    '''


    def __init__(self, ref_fasta: str = None, outdir: str = None,
                 meth_db_path: str = None, cgmap_file: str = None, cgmap_pool: bool = False,
                 beta_params: dict = {"CG": (0.5, 0.5), "CHG": (0.01, 0.05), "CHH":(0.01, 0.05)},
                 asm_file: str = None, update_boundary: bool = False,
                 collect_ch: bool = True, overwrite_db: bool = False,
                 seed: int = None, verbose: bool = False):

        self.ref_fasta   = ref_fasta
        self.outdir      = outdir
        self.meth_db_path= meth_db_path
        self.cgmap_file  = cgmap_file
        self.cgmap_pool  = cgmap_pool
        self.beta_params = beta_params
        self.asm_file    = asm_file
        self.asm_sim     = True if asm_file else False
        self.collect_ch  = collect_ch
        self.overwrite_db= overwrite_db
        self.seed        = seed
        self.verbose     = verbose
        self.meth_arr    = None
        self.pos_map     = None

        self.update_boundary   = update_boundary and not self.asm_sim
        self.context_dict      = {1:'CG', 3:'CHG', 7:'CHH', 9:'CG', 11:'CHG', 15:'CHH'}
        self.base_context_dict = {'C': {'CG':1, 'CHG':3, 'CHH':7}, 'G': {'CG':9, 'CHG':11, 'CHH':15}}
        self.base_context_table= {'C': np.array([[7,3],  [1,1]]),  'G': np.array([[15,11],[9,9]])}

        # check existence
        if not os.path.exists(ref_fasta):
            raise ValueError('Cannot find the reference genome, please check!')
        if not outdir:
            raise ValueError('Please specify the output directory!')
        if self.asm_sim and not os.path.exists(asm_file):
            raise ValueError('Please specify allelic specific methylation file correctly for ASM simulation!')

        self.ref_dict= SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))
        self.meth_db = StreamMethDB(outdir=self.outdir, overwrite_db=self.overwrite_db)
        self.meth_db.check_outdir()
        self.meth_db.save_ref(self.ref_dict)

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
            self.init_meth_db(contig_id)
            not_comp_sites = self.pos_map.index().difference(contig_profile[0].index())
            if len(not_comp_sites):
                print(f'{contig_id}: sites in the meth_db is not the same as the reference')


    def create_meth_db(self):
        '''parse the fasta, CGmap, ASM file into pickle for later simulation'''
        contig_id_list = self.ref_dict.keys()
        for contig_id in contig_id_list:
            self.current_contig = contig_id
            self.init_meth_db(contig_id)
            self.fill_cgmap(contig_id)
            self.fill_asm(contig_id)
            self.fill_dist()
            if self.verbose:
                print(f"Processed {self.pos_map.shape[0]} sites from contig {self.current_contig}")
            # the 3rd item: 1 for boundary updated by variants, 0 for not updated
            self.meth_db.output_contig(contig_id, [self.pos_map, self.meth_arr, 0], is_variant=False)


    def init_meth_db(self, contig_id):
        '''initialize data object using fasta sequence'''
        seq = self.ref_dict[contig_id].seq.upper()
        seq_len = len(seq)
        idx = 0

        if self.verbose:
            print(f"\n[Initiating the methylaiton database] for {self.current_contig}...")
        if self.collect_ch:
            count_c = seq.count("C")
            count_g = seq.count("G")
            arr_size= count_c + count_g
            self.meth_arr = np.full((arr_size, 5), fill_value=np.NaN, dtype=np.float16)
            self.meth_arr[:,1] = -1 # flag, record uninitialized sites
            self.pos_map  = pd.Series(0, index=range(arr_size), dtype=np.uint32)

            for pos, base in tqdm(enumerate(seq), disable = not self.verbose):
                if base not in {"C", "G"}:
                    continue
                if pos<2 or pos>=(seq_len-3):
                    base_d1 = base_d2 = 0
                else:
                    updown  = 1 if base == "C" else -1
                    base_d1 = seq[pos+1*updown]
                    base_d2 = seq[pos+2*updown]
                # C:{10, 11}: 1, {01}: 3, {00}: 7; G: {10, 11}: 9, {01}: 11, {00}: 15
                self.meth_arr[idx, 0] = self.get_cg_context(base, base_d1, base_d2) # context
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

        if self.verbose:
            print(f"Filling CGmap for {self.current_contig}... ", end="")
        if self.cgmap_pool:
            meth_level_pool = self.get_cgmap_pool(contig_id=None)
            idx_cg  = np.where(np.bitwise_and(self.meth_arr[:,0].astype(np.int16), 0x7)==1)
            np.random.seed(self.seed)
            self.meth_arr[idx_cg, 2] = np.random.choice(meth_level_pool[0], size=np.sum(idx_cg))

            if self.collect_ch:
                idx_chg = np.where(np.bitwise_and(self.meth_arr[:,0].astype(np.int16), 0x7)==3)
                idx_chh = np.where(np.bitwise_and(self.meth_arr[:,0].astype(np.int16), 0x7)==7)
                np.random.seed(self.seed)
                self.meth_arr[idx_chg,2] = np.random.choice(meth_level_pool[1], size=np.sum(idx_chg))
                np.random.seed(self.seed)
                self.meth_arr[idx_chh,2] = np.random.choice(meth_level_pool[2], size=np.sum(idx_chh))
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
                self.meth_arr[arr_idx, 1:3] = [1, float(meth_level)]

        if self.verbose:
            if num_cgmap_pos:
                ratio_404 = round(num_404_pos / num_cgmap_pos, 4)
                ### also print what context
                print(f"{num_cgmap_pos} sites found in CGmap file," +
                    f"among them {num_404_pos} sites ({ratio_404 * 100}%) are not compatible...")
                if ratio_404 >=0.5:
                    warnings.warn("[WARNING]: More than half sites in CGmap file cannot be found in the reference fasta\n" +
                                  "Potential reason: the CGmap does not share the same coordinates with fasta, please check!")
            else:
                print(f"No sites found in CGmap file for {self.current_contig}, skip...")
        return None


    def fill_asm(self, contig_id):
        '''fill the allelic specific sites'''
        if not self.asm_sim:
            return None

        if self.verbose:
            print(f"Filling ASM for {self.current_contig}... ", end="")
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

        if self.verbose:
            if num_asm_pos:
                ratio_404 = round(num_404_pos / num_asm_pos, 4)
                print(f"{num_asm_pos} sites found in CGmap file, " +
                    f"among them {num_404_pos} sites ({ratio_404 * 100}%) are not compatible...")
                if ratio_404 >=0.5:
                    warnings.warn("[WARNING]: More than half sites in ASM file cannot be found in the reference fasta\n" +
                                  "Potential reason: the CGmap does not share the same coordinates with fasta, please check!")
            else:
                print(f"No sites found in ASM file for {self.current_contig}, skip...")
        return None


    def fill_dist(self):
        '''fill in from beta distribution'''

        if self.verbose:
            print(f"Filling with beta distribution for {self.current_contig}...")
        idx_nan = np.isnan(self.meth_arr[:,2]) # hold for each element
        idx_ctx = np.bitwise_and(self.meth_arr[:,0].astype(np.int16), 0x7) # hold for each element
        idx_nan_cg  = np.where((idx_ctx == 1) & idx_nan)[0] # hold for index

        num_idx_nan_cg = len(idx_nan_cg)

        if num_idx_nan_cg:
            self.meth_arr[idx_nan_cg, 2] = self.simu_beta_dist(context="CG", size=num_idx_nan_cg)

        if self.collect_ch:
            idx_nan_chg = np.where((idx_ctx == 3) & idx_nan)[0]
            idx_nan_chh = np.where((idx_ctx == 7) & idx_nan)[0]

            num_idx_nan_chg = len(idx_nan_chg)
            num_idx_nan_chh = len(idx_nan_chh)
            if num_idx_nan_chg:
                self.meth_arr[idx_nan_chg, 2]= self.simu_beta_dist(context="CHG", size=num_idx_nan_chg)
            if num_idx_nan_chh:
                self.meth_arr[idx_nan_chh, 2]= self.simu_beta_dist(context="CHH", size=num_idx_nan_chh)

        self.meth_arr[idx_nan, 1] = 0
        self.meth_arr[idx_nan, 3] = self.meth_arr[idx_nan, 2]
        self.meth_arr[idx_nan, 4] = self.meth_arr[idx_nan, 2]


    def get_cgmap_pool(self, contig_id):
        '''get the pool of cgmaps, return a list of size 3'''
        ctx_idx_dict = {'CG':1, 'CHG':2, 'CHH':3}
        meth_level_pool = [[],[],[]]
        for line in parseCGmap(self.cgmap_file, contig_id, collect_ch = self.collect_ch):
            _, _,  _, context, meth_level = line
            meth_level_pool[ctx_idx_dict[context]].append(float(meth_level))
        return meth_level_pool


    def estimate_beta_params(self, context = None):
        '''estimate the beta parameters for each context'''
        ctx_idx_dict = {'CG':1, 'CHG':2, 'CHH':3}
        idx_ctx_dict = {v: k for k, v in ctx_idx_dict.items()}
        meth_level_pool = self.get_cgmap_pool(contig_id=None)
        beta_param_estimate = {}

        if context:
            idx = ctx_idx_dict[context]
            meth_list = meth_level_pool[idx]
            a, b, _, _ = beta.fit(meth_list)
            beta_param_estimate[context] = (a, b)
        else:
            for idx, meth_list in enumerate(meth_level_pool):
                context = idx_ctx_dict[idx]
                a, b, _, _ = beta.fit(meth_list)
                beta_param_estimate[context] = (a, b)
        return beta_param_estimate


    @classmethod
    def simu_beta_dist(self, context = "CG", size = 1):
        '''output the values accordig to the context using beta distribution'''
        if isinstance(context, int):
            context = self.context_dict[context]
        return beta.rvs(a=self.beta_params[context][0],
                        b=self.beta_params[context][1],
                        size=size, random_state=self.seed).astype(np.float16)


    @classmethod
    def get_cg_context(self, base, base_d1, base_d2):
        '''input the base and surrounding, output context'''
        if base == "C":
            flag_d1 = int(base_d1 == "G")
            flag_d2 = int(base_d2 == "G")
        elif base == "G":
            flag_d1 = int(base_d1 == "C")
            flag_d2 = int(base_d2 == "C")
        else:
            return None
        return self.base_context_table[base][flag_d1, flag_d2]

