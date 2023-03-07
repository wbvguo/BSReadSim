import os
import random
import numpy as np
import subprocess

from Bio import SeqIO
from tqdm import tqdm
from scipy.stats import bernoulli
from typing import Dict, Union, Tuple
from threading import Lock
from concurrent.futures import ThreadPoolExecutor, as_completed


class ReadProcessor():
    def __init__(self, undirectional, collect_ch) -> Dict:
        self.undirectional = undirectional
        self.collect_ch = collect_ch


    def process_read_pair(self, read_pair):
        """
        This processing step works as follows:
        1. check if strand capture (targeted sequencing), if yes assign strandness accordingly
        2. if not, randomly assign reads to Watson or Crick strand, with corresponding base change pattern
        3. set methylation status according to methylation profile
        4. bisulfite converted and introduce sequencing error
        5. output reads
        """

        # for directional library, read1 always C2T, strand can be either watson or crick unless strand captured
        read1_idx   = random.choice([0, 1]) 
        pattern_idx = random.choice([0, 1]) if self.undirectional else read1_idx
        strand_idx  = random.choice([0, 1]) if read_pair[0]['strand']<0 else read_pair[0]['strand']
        read_pair[1-read1_idx]['read2'] = 1
        read_pair[0]['conv'] = pattern_idx
        read_pair[1]['conv'] = pattern_idx
        read_pair[0]['strand'] = strand_idx
        read_pair[1]['strand'] = strand_idx
        
        # mask the context
        self.mask_context(read_pair[0])
        self.mask_context(read_pair[1])
        
        # retrive methy profile
        self.retrive_meth_db(read_pair[0])
        self.retrive_meth_db(read_pair[1])

        # set methylation states
        self.set_context_state(read_pair)

        # bisulfite converted
        self.treat_bisulfite(read_pair[0])
        self.treat_bisulfite(read_pair[1])

        # rev complementary
        self.rev_complement(read_pair)

        # introduce quality scores
        self.add_qual_score(read_pair[0])
        self.add_qual_score(read_pair[1])

        # introduce seq errors
        self.add_seq_err(read_pair[0])
        self.add_seq_err(read_pair[1])

        return read_pair


    def mask_context(self, read_rec):
        '''mask the context based on read substitution pattern (0 for C2T, 1 for G2A)'''
        if read_rec['conv'] :
            if self.collect_ch:
                read_rec['ctx'] = np.ma.masked_less(read_rec['ctx'], 8)
            else:
                read_rec['ctx'] = np.ma.not_equal(read_rec['ctx'], 1)
        else:
            if self.collect_ch:
                read_rec['ctx'] = np.ma.masked_greater(read_rec['ctx'], 8)
            else:
                read_rec['ctx'] = np.ma.not_equal(read_rec['ctx'], 9)


    def retrive_meth_db(self, read_rec):
        '''retrive methylation levels from meth_db, append meth and pos to read_rec'''
        read_meth = np.zeros(self.read_len)
        read_pos  = read_rec['start'] + np.arange(self.read_len)
        site_flag = np.logical_not(read_rec['ctx'].mask)            # unmasked sites

        if np.any(site_flag):                                       # contain methylable bases
            arr_idx  = 2

            if read_rec['flag_pos']:                                # covers mutation position
                if self.asm_sim:
                    arr_idx  = 4 if read_rec['flag_mut'] else 3

                if read_rec['n_indel']:                             # handle indel first (offset)
                    read_pos += read_rec['ofs']
                    indel_site= read_rec['cgr'] == 3
                    read_meth[indel_site]= self.fetch_meth_val(read_pos[indel_site], arr_idx, 3)

                if read_rec['n_sub']:
                    snp_site = site_flag & (read_rec['cgr'] == 1)   # snp methylable site
                    read_meth[snp_site]  = self.fetch_meth_val(read_pos[snp_site],  arr_idx, 1)

                match_site = site_flag & (read_rec['cgr'] == 0)     # match methylable site
                read_meth[match_site] = self.fetch_meth_val(read_pos[match_site],arr_idx, 0)
            else:
                match_site = site_flag                              # SNP/INDEL free region
                read_meth[match_site] = self.fetch_meth_val(read_pos[match_site],arr_idx, 0)
        read_rec['meth']   = read_meth
        read_rec['pos']    = read_pos


    def fetch_meth_val(self, pos_arr, arr_idx, var_type):
        '''fetch the methylation value from meth_db'''
        if var_type == 1:           # SNP
            meth_val = [self.variant_profile[pos][0] for pos in pos_arr]
        elif var_type== 3:          # insertion sites will have the same coordinate
            insert_pos, pos_count = np.unique(pos_arr, return_counts=True) # uniq_pos, pos_count
            meth_val  = []
            for ix, pos in enumerate(insert_pos):
                try:
                    insert_meth = list(self.variant_profile[pos][0][:pos_count[ix]])
                except:
                    insert_meth = [0]*pos_count[ix]
                meth_val += insert_meth
        else:                       # match
            meth_val  = list(self.meth_arr[self.pos_map[pos_arr], arr_idx])
        return meth_val


    def set_context_state(self, read_rec):
        '''set the context methylation state'''
        if self.pair_end:
            if read_rec[0]['inner_dist'] <= 0: # overlapped read pair: merge meth, get state, split
                overlap_idx = np.where(read_rec[0]['pos'] == read_rec[1]['pos'][0])
                if np.any(overlap_idx):
                    for idx in overlap_idx:
                        overlap_len = self.read_len - idx
                        if read_rec[0]['pos'][idx:] == read_rec[1]['pos'][:overlap_len]:
                            break
                else:
                    idx = self.read_len
                # it's okay to have idx=0, or np.where returns null
                comb_meth = np.concatenate((read_rec[0]['meth'][:idx],read_rec[1]['meth']), axis=0)
                comb_state= self.fetch_meth_state(comb_meth)
                read_rec[0]['ctx'][np.where(comb_state[:self.read_len])[0]] +=1
                read_rec[1]['ctx'][np.where(comb_state[-self.read_len:])[0]]+=1
            else:
                read1_state = self.fetch_meth_state(read_rec[0]['meth'])
                read_rec[0]['ctx'][np.where(read1_state)[0]] += 1
                read2_state = self.fetch_meth_state(read_rec[1]['meth'])
                read_rec[1]['ctx'][np.where(read2_state)[0]] += 1
        else:
            read_state = self.fetch_meth_state(read_rec['meth'])
            read_rec['ctx'][np.where(read_state)[0]] += 1


    def fetch_meth_state(self, read_meth):
        '''set methylation states based on the meth_arr'''
        if self.site_dependency:
            # generate methylation pattern according to read_meth and distance
            pass
        else:
            meth_states = np.zeros(read_meth.size)
            nonzero_idx = np.where(read_meth)[0]    # this returns a tuple
            meth_states[nonzero_idx] = bernoulli.rvs(read_meth[nonzero_idx], size = nonzero_idx.size)
        return meth_states


    def treat_bisulfite(self, read_rec):
        """ bisulfite conversion """
        unmeth_idx = np.where(np.bitwise_and(read_rec['ctx'], 0x1)==1)[0] # behave strange without ==1
        conv_states= bernoulli.rvs(self.conversion_rate, size=unmeth_idx.size)
        conv_idx   = unmeth_idx[conv_states == 1]   # successfully converted base index
        read_rec['seq'][conv_idx] = np.bitwise_and(read_rec['seq'][conv_idx] + 2, 0x3) # C2T, G2A
        unconv_idx = unmeth_idx[conv_states == 0]   # remains the same base index
        read_rec['cgr'][unconv_idx]= 4


    def rev_complement(self, read_pair):
        '''reverse complementary for pair_end'''
        for read_rec in read_pair:
            if read_rec['pair']:
                read_rec['cgr'] = np.flip(read_rec['cgr'])
                read_rec['ctx'] = np.flip(read_rec['ctx'])
                read_rec['seq'] = 3 - np.flip(read_rec['seq'])
                read_rec['conv']= 1 - read_rec['conv']


    def add_seq_err(self, read_rec):
        ''' introduce sequencing error'''
        if self.random_err:
            err_idx = np.where(bernoulli.rvs(self.err_rate, size = self.read_len))[0] #cannot np.squeeze
            if np.any(err_idx):
                for idx in err_idx:
                    base_ori = read_rec['seq'][idx]
                    base_err = np.random.choice(np.setdiff1d(np.array([0,1,2,3]), base_ori), size = 1)[0]
                    read_rec['seq'][idx] = base_err
                    read_rec['cgr'][idx] = 2
                    read_rec['ctx'][idx] = 5 if (base_ori, base_err)==[(1,3), (2,0)][read_rec['conv']] else 6
        else:
            # generate sequencing error based on a profile TODO:
            pass


    def add_qual_score(self, read_rec, qual_num=56):
        '''add quality scores'''
        if self.qual_uniform:
            qual_arr = np.full(self.read_len, qual_num)
            qual_arr[[0,-1][read_rec['strand']^read_rec['conv']]] = qual_num -1 # denote the 5' end
            read_rec['qual'] = qual_arr
        else:
            # generate quality score from a profile TODO:
            pass
