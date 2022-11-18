import os
import random
import sys
import subprocess
import numpy as np

from Bio import SeqIO
from tqdm import tqdm
from scipy.stats import bernoulli
from typing import Dict, Union, Tuple
from threading import Lock
from concurrent.futures import ThreadPoolExecutor

from SetMethylation import SetMethylation
from StreamWGSIM import StreamWGSIM
from StreamOutput import StreamOutput
from UtilityFunctions import get_wgsim_path, reverse_complement


class BSReadSim:
    """
    The bisulfite sequencing simulation works as follows:

    1. modified WGSIM module is invoked to simulate pair-end Illumina sequencing reads
        - if run in single end mode, the second read is not processed nor output,
          correspondly the number of reads simulated is doubled to get desired depth
        - the output of 2 reads from WGSIM is always the watson 5' to 3' and error-free
          the strandness/orientation is determined in the python module
    2. methylation values are set for methylable bases (Cytosine and Guanine), from 3 ways
        - reference CGmap, either faithfully or from a pool
        - from a beta distribution (whose parameters can be estimated from the reference CGmap)
        - from the ASM file
    3. Allelic methylation and site-site dependency pattern is permitted if needed [TODO]
    4. Reads are bisulfite converted (with a conversion rate) and sequencing error can be introduced

    # parameters
    :param str   ref_fasta      : path to the reference genome (.fasta/.fa/.fa.gz) [None]
    :param str   outdir         : path to the output directory [None]
    :param str   prefix         : prefix of output reads ['sim']
    :param str   meth_db_path   : path to the previous obj (.pkl) [None]
    :param str   cgmap_file     : path to the cgmap file (.CGmap/.CGmap.gz) [False]
    :param bool  cgmap_pool     : whether to pool methylaiton levels and random draws from it [False]
    :param bool  asm_sim        : whether to conduct allelic specific methylation (ASM) [False]
    :param str   asm_file       : path to the allelic specific methylation (ASM) file [None]
    :param dict  beta_params    : dict of beta parameters for CG/CHG/CHH methylation values simulation
                                  [{"CG": (0.5, 0.5), "CHG": (0.01, 0.05), "CHH":(0.01, 0.05)}]
    :param bool  collect_ch     : whether to collect/simulate the CHG/CHH sites [False]
    :param bool  overwrite_db   : whether to overwrite the meth_db if it already exists [True]
    :param str   vcf_file       : path to the reference vcf file
    :param float mut_rate       : mutation error rate for random mutation generation [0.0010]
    :param float mut_indel_frac : INDEL fraction in the random generated mutation [0.15]
    :param float indel_ext_prob : probability the INDEL length will be extended [0.15]
    :param bool  haplo_mode     : simulate only homozygous variants [False]
    :param int   seed           : seed for methylation, mutation and sequencing error generation,
                                  seed is none if sign is negative [-1]
    :param int   insert_mean    : mean of insert size
    :param int   insert_std     : standard deviation of insert size
    :param int   read_len       : length of simulated reads [100]
    :param float depth          : average sequencing depth for each contigs [20]
    :param float seq_err        : sequencing error rate [0.005]
    :param float conversion_rate: bisulfite conversion rate [0.995]
    :param bool  pair_end       : whether to simulate pair-end data or not [False]
    :param bool  undirectional  : simulate undirectional (PCR product of Watson and Crick strands) [False]
    :param float propN_cutoff   : reference segments where the proportion of ambiguous bases, - or N,
                                  greater than threshold will be skipped [0.05]
    :param bool  verbose        : whether to output processing details [True]
    :param bool  shuffle        : by default the reads are segmented by contigs, shuffle it [True]

    Usage:
    ```python
    simulation = BSReadSim(ref_fasta, outdir, **kwargs)
    simulation.run()
    ```

    """

    def __init__(self, ref_fasta: str = None, outdir: str = None, prefix: str = "sim",          # required arguments
                 meth_db_path: str = None, cgmap_file: str = None, cgmap_pool: bool = False,    # methylation ref
                 asm_file: str = None,                                                          # ASM simulation
                 beta_params: Dict = {"CG":(0.5, 0.5), "CHG":(0.01, 0.05), "CHH":(0.01, 0.05)}, # beta distribution
                 collect_ch: bool = True, overwrite_db: bool = False,                           # methylation output
                 vcf_file: str = None,                                                          # genetic variant reference
                 mut_rate: float = 0.0010, haplo_mode: bool = False, seed: int = -1,            # mutation parameter
                 update_boundary: bool = False, site_dependency: bool = False, random_err= True, collect_sim_stats: bool = False,
                 mut_indel_frac: float = 0.15, indel_ext_prob: float = 0.15,                    # indel parameter
                 insert_mean: int = 400, insert_std: int = 25,                                  # fragment setting
                 read_len: int = 100, depth: int = 20, num_reads: int = None,                   # reads setting
                 seq_err: float = 0.005, conversion_rate: float = 0.998,                        # rates setting
                 undirectional: bool = False, pair_end: bool = True, propN_cutoff: float = 0.05,# sequencing protocol
                 n_threads: int = 4, gzip: bool = True, verbose: bool = True, shuffle: bool = True):

        # check required auguments are provided
        if not os.path.exists(ref_fasta):
            raise ValueError('Cannot find the reference file, please check!')
        if not outdir:
            raise ValueError('Please specify the output directory!')

        self.ref_fasta= ref_fasta
        self.outdir   = outdir
        self.prefix   = prefix
        self.seed     = seed
        self.n_threads= n_threads
        self.lock     = Lock()

        # prepare the methylation reference
        print('Initiating methylation profile:\n')
        self.meth_db  = SetMethylation(ref_fasta=ref_fasta, outdir=outdir,
                                       meth_db_path=meth_db_path,
                                       cgmap_file=cgmap_file, cgmap_pool=cgmap_pool,
                                       asm_file=asm_file,
                                       beta_params=beta_params,
                                       collect_ch=collect_ch,
                                       update_boundary=update_boundary,
                                       overwrite_db=overwrite_db,verbose=verbose,
                                       seed = None if seed < 0 else seed)

        # prepare wgsim command
        wgsim_args    = [get_wgsim_path(), ref_fasta]
        if not num_reads:
            num_reads = int(self.meth_db.genome_len*depth/read_len/(1+int(pair_end)))
        wgsim_options = {'-d': insert_mean, '-s': insert_std,
                         '-1': read_len, '-2': read_len, '-N': num_reads,
                         '-g': vcf_file,
                         '-r': mut_rate, '-h': int(haplo_mode), '-S': seed,
                         '-R': mut_indel_frac, '-X': indel_ext_prob,
                         '-A': propN_cutoff, '-e': 0, '-m': 1}

        self.wgsim_args     = wgsim_args
        self.wgsim_options  = wgsim_options
        self.sim_cmd_part   = wgsim_args + [str(item) for key_val in wgsim_options.items() for item in key_val]

        # sequencing settings
        self.num_reads      = num_reads
        self.read_len       = read_len
        self.pair_end       = pair_end
        self.undirectional  = undirectional
        self.conversion_rate= conversion_rate
        self.random_err     = random_err
        self.seq_err        = seq_err
        self.collect_ch     = collect_ch
        self.asm_sim        = True if asm_file else False
        self.site_dependency= site_dependency


        # to hold intermediate data
        self.current_contig = None
        self.pos_map        = None
        self.meth_arr       = None
        self.variant_profile= None

        # prepare output
        self.gzip           = gzip
        self.shuffle        = shuffle

        # progress
        self.verbose        = verbose
        if verbose:
            self.tqdm_count = [0, 0, num_reads * 0.01] # current, previous, step size
            self.tqdm_pbar  = tqdm(range(num_reads), bar_format="{desc:<5.5}{percentage:3.0f}%|{bar:20}{r_bar}")


    def run(self):
        """Simulating Bisulfite sequencing reads with parallelization:
        # 1. recieve from WGSIM (in chunk)
        # 2. process and write to disk (multiple threads)
        """
        print('Simulating methylated Reads:\n')
        print(f'[CMD]: {" ".join(self.sim_cmd_part)} for each contigs\n')

        with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
            for contig_id in self.meth_db.ref_dict.keys():
                sim_cmd  = self.sim_cmd_part + ['-c', contig_id]
                read_gen = StreamWGSIM(sim_cmd=sim_cmd, pair_end=self.pair_end)         # should we make it thread-safe?
                var_contig, sim_data = next(read_gen)                                   # the first element is the variants
                self.current_contig = var_contig                                        # update the profiles
                self.pos_map, self.meth_arr, _ = self.meth_db.load_contig(var_contig)   # 3 items list, pos_map, meth_arr, status
                self.variant_profile= self.meth_db.set_var_meth(var_contig, sim_data)   # a dict
                if self.pair_end:
                    job_arr = [executor.submit(self.process_read_pair, read_pair) for read_pair in read_gen] #what if read_gen is empty at very beginning
                else:
                    job_arr = [executor.submit(self.process_read, read_pair) for read_pair in read_gen]

                if self.verbose:
                    for job in job_arr:
                        job.add_done_callback(self.progress_bar)


        for output in self.output_obj:      # close the fastq object
            output.close()

        if self.verbose:                    # close the progress bar
            self.tqdm_pbar.close()

        for fastq in self.fastq_list:       # shuffle or gzip reads
            if self.shuffle:
                self.shuffle_read(fastq, self.ref_fasta)
            if self.gzip:
                self.gzip_read(fastq)

        print('Simulation Finished!\n')


    def process_read(self, read_pair):
        '''process single end reads'''
        read_flip = random.choice([0, 1])
        read1_sub = random.choice([0, 1]) if self.undirectional else read_flip

        self.mask_context(read_pair[0], read1_sub)
        self.retrive_meth_db(read_pair[0])
        self.set_context_state(read_pair[0])
        self.treat_bisulfite(read_pair[0])
        self.add_seq_err(read_pair[0])
        self.add_qual_score(read_pair[0])
        self.output_reads(read_pair, read_flip, read1_sub)


    def process_read_pair(self, read_pair):
        """
        This processing step works as follows:
        1. randomly assign reads to Watson or Crick strand, with corresponding base change pattern
        2. set methylation status according to methylation profile
        3. bisulfite converted and introduce sequencing error
        4. output reads
        """

        # for directional library, read1 will be G2A if read_flip else C2A, strand will be W else C
        # if read_flip: read_pair[0], read_pair[1] = read_pair[1], read_pair[0]
        read_flip = random.choice([0, 1])
        read1_sub = random.choice([0, 1]) if self.undirectional else read_flip

        # mask the context
        self.mask_context(read_pair[0], read1_sub)
        self.mask_context(read_pair[1], 1-read1_sub)
        # retrive methy profile
        self.retrive_meth_db(read_pair[0])
        self.retrive_meth_db(read_pair[1])

        self.set_context_state(read_pair)

        # bisulfite converted
        self.treat_bisulfite(read_pair[0])
        self.treat_bisulfite(read_pair[1])

        # introduce seq errors
        self.add_seq_err(read_pair[0])
        self.add_seq_err(read_pair[1])

        # introduce quality scores
        self.add_qual_score(read_pair[0])
        self.add_qual_score(read_pair[1])
        # output
        self.output_reads(read_pair, read_flip, read1_sub)


    def mask_context(self, read_rec, read_sub):
        '''mask the context based on read substitution pattern (0 for C2T, 1 for G2A)'''
        if read_sub:
            if self.collect_ch:
                read_rec['ctx'] = np.ma.masked_greater(read_rec['ctx'], 8)
            else:
                read_rec['ctx'] = np.ma.not_equal(read_rec['ctx'], 9)
        else:
            if self.collect_ch:
                read_rec['ctx'] = np.ma.masked_less(read_rec['ctx'], 8)
            else:
                read_rec['ctx'] = np.ma.not_equal(read_rec['ctx'], 1)


    def retrive_meth_db(self, read_rec):
        '''retrive methylation levels from meth_db, append meth and pos to read_rec'''
        read_meth = np.zeros(self.read_len)
        site_flag = np.logical_not(read_rec['ctx'].mask)            # not masked sites

        if np.any(site_flag):                                       # contain methylable bases
            arr_idx  = 2
            read_pos = read_rec['start'] + np.arange(self.read_len)

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
            insert_pos= [x for x in pos_arr if pos_arr.count(x) > 1]
            meth_val  = []
            for pos in insert_pos:
                meth_val += list(self.variant_profile[pos][0])
        else:                       # match
            meth_val  = list(self.meth_arr[self.pos_map[pos_arr], arr_idx])
        return meth_val


    def set_context_state(self, read_rec):
        '''set the context methylation state'''
        if self.pair_end: #overlappd
            if read_rec[0]['inner_dist'] <= 0: # overlapped read pair: merge meth, get state, split
                overlap_idx = np.where(read_rec[0]['pos'] == read_rec[1]['pos'][0])
                if np.any(overlap_idx):
                    for idx in overlap_idx:
                        overlap_len = self.read_len - idx
                        if read_rec[0]['pos'][idx:] == read_rec[1]['pos'][:overlap_len]:
                            break
                else:
                    idx = self.read_len
                comb_meth = np.concatenate(read_rec[0]['meth'][:idx],read_rec[1]['meth'])
                comb_state= self.fetch_meth_state(comb_meth)
                read_rec[0]['ctx'][np.where(comb_state[:self.read_len])] +=1
                read_rec[1]['ctx'][np.where(comb_state[-self.read_len:])]+=1
            else:
                read1_state = self.fetch_meth_state(read_rec[0]['meth'])
                read_rec[0]['ctx'][np.where(read1_state)] += 1
                read2_state = self.fetch_meth_state(read_rec[1]['meth'])
                read_rec[1]['ctx'][np.where(read2_state)] += 1
        else:
            read_state = self.fetch_meth_state(read_rec['meth'])
            read_rec['ctx'][np.where(read_state)] += 1


    def fetch_meth_state(self, read_meth):
        '''set methylation states based on the meth_arr'''
        if self.site_dependency:
            # generate methylation pattern according to read_meth and distance
            pass
        else:
            meth_states = bernoulli.rvs(read_meth, size = len(read_meth))
            #read_ctx[np.where(meth_states)] += 1 # context increased by 1 if methylated
        return meth_states


    def treat_bisulfite(self, read_rec):
        """ bisulfite conversion """
        unmeth_idx = np.where(np.bitwise_and(read_rec['ctx'], 0x1)==1) # behave strange without ==1
        conv_states= bernoulli.rvs(self.conversion_rate, size=len(unmeth_idx))
        conv_idx   = unmeth_idx[conv_states]        # successfully converted base index
        read_rec['seq'][conv_idx]  = np.bitwise_and(read_rec['seq'][conv_idx] + 2, 0x3) # C2T, G2A
        unconv_idx = unmeth_idx[conv_states == 0]   # remains the same base index
        read_rec['cgr'][unconv_idx]= 4


    def add_seq_err(self, read_rec):
        ''' introduce sequencing error'''
        if self.random_err:
            err_states = bernoulli.rvs(self.seq_err, size = self.read_len)
            err_idx    = np.where(err_states)
            if err_idx:
                base_set = {0,1,2,3}
                for idx in err_idx:
                    base_ori = read_rec[idx]
                    base_alt = np.ramdon.choice(base_set.difference(base_ori), size = 1)
                    read_rec['cgr_str'][idx] = "e" if (base_ori, base_alt) in {(1,3), (2,0)} else "E"
        else:
            # generate sequencing error based on a profile
            pass


    def add_qual_score(self, read_rec, qual_uniform = True):
        '''add quality scores for the read'''
        if qual_uniform:
            qual_arr = np.full(self.read_len, chr(read_rec['qual']))
            read_rec['qual'] = qual_arr
        else:
            # generate quality score from a profile
            pass


    def output_reads(self, read_pair, read_flip, read1_sub):
        """write bisulfite reads to disk"""
        if self.pair_end:
            self.rev_complement(read_pair[1])

            if read_flip:
                read_pair[0], read_pair[1] = read_pair[1], read_pair[0]
            conv_tag = ('G2A', 'C2T') if read1_sub else ('C2T', 'G2A')
            for idx, read_rec in enumerate(read_pair):
                conv_pattern = conv_tag[idx]
                
            
        else:
            
        
        # format reads
        conversion_1, conversion_2 = ('C2T', 'G2A') if sim_data[0][0]['sub_base'] == sub_base else ('G2A', 'C2T')
        reverse_read = 1 if sub_base == 'C' else 0
        sim_data[reverse_read][1] = reverse_complement(sim_data[reverse_read][1]) #seq
        sim_data[reverse_read][3] = sim_data[reverse_read][3][::-1] #qual
        sim_data[reverse_read][0]['cigar'] = sim_data[reverse_read][0]['cigar'][::-1] #cigar
        read_label = f'@{sim_data[0][0]["read_id"]}_{sim_data[0][0]["chrom"]}/1'
        read_comment = f'+{sim_data[0][0]["chrom"]}:{sim_data[0][0]["start"]+1}:' \
                       f'{sim_data[0][0]["end"]}:{sim_data[0][0]["cigar"]}:{ref_strand}{conversion_1}'
        read = f'{read_label}\n{sim_data[0][1].upper()}\n{read_comment}\n{sim_data[0][3]}\n'
        self.output_obj[0].write(read)
        if self.pair_end:
            read_label = f'@{sim_data[1][0]["read_id"]}_{sim_data[1][0]["chrom"]}/2'
            read_comment = f'+{sim_data[1][0]["chrom"]}:{sim_data[1][0]["start"]+1}:' \
                           f'{sim_data[1][0]["end"]}:{sim_data[1][0]["cigar"]}:{ref_strand}{conversion_2}'
            read = f'{read_label}\n{sim_data[1][1].upper()}\n{read_comment}\n{sim_data[1][3]}\n'
            self.output_obj[1].write(read)


    @property
    def progress_bar(self):
        '''show the progress of read simulaiton'''
        with self.lock:
            self.tqdm_count[0] += 2 if self.pair_end else 1

        incre_amount = self.tqdm_count[0] - self.tqdm_count[1]
        if incre_amount > self.tqdm_count[2]:
            self.tqdm_pbar.update(incre_amount)
            self.tqdm_count[1] = self.tqdm_count[0]


# column context, row cigar
cigar_table = np.array([['M', 'c', 'C', 'b',  'B', '-', '-', 'a',  'A', 'c', 'C', 'b',  'B', '-', '-', 'a',  'A'], # match
                        ['x', 'x', 'X', 'x',  'X', '-', '-', 'x',  'X', 'x', 'X', 'x',  'X', '-', '-', 'x',  'X'], # snp
                        ['-', '-', '-', '-',  '-', '-', '-', '-',  '-', '-', '-', '-',  '-', '-', '-', '-',  '-'], # empty
                        ['M', 'i', 'I', 'i',  'I', '-', '-', 'i',  'I', 'i', 'I', 'i',  'I', '-', '-', 'i',  'I'], # insert
                        ['-', '#', '-', '#',  '-', '-', '-', '#',  '-', '#', '-', '#',  '-', '-', '-', '#',  '-']])# convert failed
