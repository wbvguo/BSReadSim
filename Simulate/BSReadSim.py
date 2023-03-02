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

from LockedIterator import LockedIterator
from SetMethylation import SetMethylation
from StreamReads import StreamReads
from StreamHTSIM import StreamHTSIM
from UtilityFunctions import get_htsim_path
from ParseGenome import ParseGenome


class BSReadSim:
    """
    BSReadSim works as follows:
    1. construct meth_db:
        a. parse the reference fasta and initiate a methylation database (a np.array which holds the 
           methylable base (C or G)'s position and methylation value). The database will be updated 
           when receiving mutations (SNP + INDEL) during simulation.
        b. fill in site's methmethylaiton value, which can be obtained from the following 3 ways: 
            - reference CGmap file, either pool randomly or assign to correponding sites 
            - reference allelic specific methylation profiles 
            - beta distribution (prarameters can be estimated from CGmap file or set by users) 
    2. invoke htsim module to simulate variants and pair-end Illumina sequencing reads 
        - genetic variants is introduced by either random mutations or from a reference VCF file 
        - meth_db is updated according to variants: the variant's methylable base (also the boundary 
          methylable bases whose context (CG/CHG/CHH) is changed by variants) is set using beta 
          distribution according to their new context 
        - the htsim module generates error-free reads whose directionality is watson 5' to 3', the 
          read strandness/orientation is determined in the following python module 
        - if run in single end mode, the second read is not processed nor output, correspondingly 
          the #reads to generate is doubled to get desired depth 
    3. For each read, methylation values are retrived from meth_db for methylable bases (C or G),
       methylation states is simulated using either the naive independent bernoulli distribution or
       using pretrained RNN model to accomedate site-site dependency 
    4. Allelic methylation is enabled if the asm_file specified 
    5. Reads are bisulfite converted (with a conversion rate) and sequencing error is introduced 

    # parameters
    :param str   ref_fasta      : path to the reference genome (.fasta/.fa/.fa.gz) [None]
    :param str   outdir         : path to the output directory [None]
    :param str   prefix         : prefix of output reads ["sim"]
    :param str   meth_db_path   : path to the previous obj (.pkl) [None]
    :param str   cgmap_file     : path to the cgmap file (.CGmap/.CGmap.gz) [False]
    :param bool  cgmap_pool     : whether to pool methylaiton levels and random draws from it [False]
    :param str   asm_file       : path to the allelic specific methylation (ASM) file [None]
    :param dict  beta_params    : dict of beta parameters for CG/CHG/CHH methylation values simulation
                                  [{"CG": (0.5, 0.5), "CHG": (0.01, 0.05), "CHH":(0.01, 0.05)}]
    :param bool  collect_ch     : whether to collect/simulate the CHG/CHH sites [False]
    :param bool  update_boundary: whether to change the boundary sites' mehtylation due to mutation [False]
    :param bool  overwrite_db   : whether to overwrite the meth_db if it already exists [True]
    :param str   vcf_file       : path to the reference vcf file
    :param float mut_rate       : mutation error rate for random mutation generation [0.0010]
    :param bool  haplo_mode     : simulate only homozygous variants [False]
    :param int   seed           : seed for methylation, mutation and sequencing error generation,
                                  seed is none if sign is negative [-1]
    :param float mut_indel_frac : INDEL fraction in the random generated mutation [0.15]
    :param float indel_ext_prob : probability the INDEL length will be extended [0.15]
    :param int   insert_mean    : mean of insert size
    :param int   insert_std     : standard deviation of insert size
    :param int   insert_min     : min insert size
    :param int   insert_max     : max insert size
    :param int   read_len       : length of simulated reads [100]
    :param float depth          : average sequencing depth for each contigs [20]
    :param float num_reads      : total number of reads to simulate [None]
    :param bool  random_err     : whether to generate sequencing error randomly [True]
    :param float err_rate       : sequencing error rate for random error generation[0.005]
    :param float propN_cutoff   : reference segments where the proportion of ambiguous bases, - or N,
                                  greater than threshold will be skipped [0.05]
    :param bool  site_dependency: whether to consider site-site dependency on the read-level methylation [False]
    :param float conversion_rate: bisulfite conversion rate [0.995]
    :param bool  undirectional  : simulate undirectional (PCR product of Watson and Crick strands) [False]
    :param bool  pair_end       : whether to simulate pair-end data or not [False]
    :param str   cut_site_str   : enzyme cut site for RRBS, cutting position is denoted by |, multiple sites
                                  are separated by comma, for example MspI and TaqI: "C|CGG,T|CGA" [None]
    :param str   probe_bed_file : probe bed file for TBS simulation [None]
    :param bool  is_uniform     : Whether the coverage should be uniform or not [True]
    :param int   n_threads      : number of threads for simulation [4]
    :param bool  shuffle        : by default the reads are segmented by contigs, shuffle it [True]
    :param bool  gzip           : whether to compressed the output fastq files [True]
    :param bool  verbose        : whether to output processing details [True]
    :param bool  collect_stats  : whether to collect simulation statistics [False]


    Usage:
    ```python
    simulation = BSReadSim(ref_fasta, outdir, **kwargs)
    simulation.run()
    ```

    """

    def __init__(self, ref_fasta: str = None, outdir: str = None, prefix: str = "sim",              # required arguments
                 meth_db_path: str = None, cgmap_file: str = None, cgmap_pool: bool = False,        # methylation ref
                 asm_file: str = None,                                                              # ASM simulation
                 beta_params: Dict = {"CG":(0.5, 0.5), "CHG":(0.01, 0.05), "CHH":(0.01, 0.05)},     # beta distribution
                 collect_ch: bool = True, update_boundary: bool = False, overwrite_db: bool = False,# meth_db output
                 vcf_file: str = None,                                                              # genotype reference
                 mut_rate: float = 0.0010, haplo_mode: bool = False, seed: int = -1,                # mutation parameter
                 mut_indel_frac: float = 0.15, indel_ext_prob: float = 0.15,                        # indel parameter
                 insert_mean: int = 400, insert_std: int = 25,                                      # fragment setting
                 insert_min: int = 100, insert_max: int = 1000,
                 read_len: int = 100, depth: int = 20, num_reads: int = None,                       # reads setting
                 propN_cutoff: float = 0.05, conversion_rate: float = 0.998,                        
                 undirectional: bool = False, pair_end: bool = True,                                # sequencing protocol
                 random_err= True, err_rate: float = 0.005, error_model: str = None,                # sequencing error TODO:
                 qual_uniform: bool =True, qual_model: str = None,                                  # quality model, TODO:
                 is_uniform: bool = True, gc_bias_file: str = None,                                 # coverage model, gc_bias_file TODO:
                 site_dependency: bool = False, site_model: str = None,                             # site dependency model TODO:
                 cut_site_str: str = None, rrbs_model: str = None,                                  # RRBS technology model TODO:
                 probe_bed_file: str = None, probe_deviation: int =None,                            # TBS technology file
                 n_threads: int = 4, shuffle: bool = True, gzip: bool = True,                       # output setting
                 verbose: bool = True, collect_stats: bool = False):
 
        # check required auguments are provided
        if not os.path.exists(ref_fasta):
            raise ValueError('Cannot find the reference file, please check!')
        if not outdir:
            raise ValueError('Please specify the output directory!')
        if cut_site_str and probe_bed_file:
            raise ValueError('Please specify files for only one technology mode!')

        self.ref_fasta      = ref_fasta
        self.outdir         = outdir
        self.pkl_dir        = f'{self.outdir}/pkl/'
        self.tmp_dir        = f'{self.outdir}/tmp/'
        self.prefix         = prefix
        self.seed           = None if seed < 0 else seed
        self.n_threads      = n_threads
        #self.countLock      = Lock()
        self.writeLock      = Lock()
        self.verbose        = verbose

        # sequencing settings
        self.read_len       = read_len
        self.pair_end       = pair_end
        self.undirectional  = undirectional
        self.conversion_rate= conversion_rate
        self.random_err     = False if error_model else random_err
        self.err_rate       = err_rate
        self.collect_ch     = collect_ch
        self.asm_sim        = True if asm_file else False
        self.site_dependency= True if site_model else site_dependency
        self.is_uniform     = is_uniform
        self.collect_stats  = collect_stats
        self.qual_uniform   = qual_uniform

        # prepare folder, advanced files
        self.create_outdir()
        self.tool_path      = os.path.dirname(os.path.dirname(get_htsim_path()))

        if self.is_uniform:
            self.gc_bias_file=None
            self.error_model= None
            self.qual_model = None
            self.rrbs_model = None
            self.site_model = None
        else:
            self.gc_bias_file=gc_bias_file if gc_bias_file else f'{self.tool_path}/data/model/gc_bias.txt'
            self.error_model= error_model if error_model else f'{self.tool_path}/data/model/error.pickle'
            self.qual_model = qual_model if qual_model else f'{self.tool_path}/data/model/qual.pickle'
            self.rrbs_model = rrbs_model if rrbs_model else f'{self.tool_path}/data/model/rrbs.pickle'
            self.site_model = site_model if site_model else f'{self.tool_path}/data/model/site.pickle'

        # parse genome
        print('Initiating genome...')
        self.genome_db  = ParseGenome(ref_fasta=self.ref_fasta, outdir=self.tmp_dir, num_reads=num_reads,
                                      depth=depth, read_len=read_len, pair_end=pair_end, is_uniform=is_uniform,
                                      probe_bed_file=probe_bed_file, cut_site_str=cut_site_str, rrbs_model=rrbs_model,
                                      insert_min=insert_min, insert_max=insert_max, overwrite_db=overwrite_db)
        self.ref_dict   = self.genome_db.ref_dict
        self.bed_file   = self.genome_db.bed_file
        self.tech_mode  = self.genome_db.tech_mode
        self.num_reads  = self.genome_db.num_reads
        self.count_dict = self.genome_db.count_dict

        # prepare the methylation reference
        print('Initiating methylation profile...')
        self.meth_set   = SetMethylation(ref_dict=self.ref_dict, outdir=self.outdir,
                                         meth_db_path=meth_db_path,
                                         cgmap_file=cgmap_file, cgmap_pool=cgmap_pool,
                                         asm_file=asm_file,
                                         beta_params=beta_params,
                                         collect_ch=collect_ch,
                                         update_boundary=update_boundary,
                                         overwrite_db=overwrite_db,
                                         verbose=self.verbose,
                                         seed=self.seed)
        self.meth_db    = self.meth_set.meth_db

        # prepare output
        self.fastq_out  = StreamReads(outdir=outdir, prefix=prefix, pair_end=pair_end, 
                                      write_lock = self.writeLock,
                                      shuffle=shuffle, ref_fasta=self.ref_fasta, gzip=gzip)

        # prepare htsim command
        htsim_args      = [get_htsim_path(), ref_fasta]
        htsim_options   = {'-i': insert_mean, '-I': insert_std, '-m': insert_min, '-M': insert_max,
                           '-1': read_len, '-2': read_len, '-e': 0, '-A': propN_cutoff, '-u': int(is_uniform), '-f': 1,
                           '-g': vcf_file, '-r': mut_rate, '-R': mut_indel_frac, '-X': indel_ext_prob, 
                           '-h': int(haplo_mode), '-s': seed, '-T': self.tech_mode, '-x': cut_site_str, 
                           '-b': self.bed_file, '-B': self.gc_bias_file, '-D': probe_deviation}
        self.cmd_part   = htsim_args + [str(item) for key_val in htsim_options.items() for item in key_val]
        self.htsim_args     = htsim_args
        self.htsim_options  = htsim_options

        # to hold intermediate data
        self.current_contig = None
        self.pos_map        = None
        self.meth_arr       = None
        self.variant_profile= None

        # # progress # has some problem
        # if self.verbose:
        #     self.tqdm_count = [0, 0, self.num_reads * 0.01] # current, previous, step size
        #     self.tqdm_pbar  = tqdm(range(self.num_reads), bar_format="{desc:<5.5}{percentage:3.0f}%|{bar:20}{r_bar}")


    def run(self):
        """Simulating Bisulfite sequencing reads with parallelization:
        # 1. recieve from htsim
        # 2. process and write to disk (multiple threads)
        """
        print('Simulating methylated Reads...')
        print(f'[CMD]: {" ".join(self.cmd_part)}')
        if self.verbose:
            print(f'[INFO]: #reads/#read pairs for each contig', end=" ")
            print(self.count_dict)

        # with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
        for contig_id in self.count_dict.keys():
            if self.count_dict[contig_id] == 0: # need to test if reads < self.n_threads TODO:
                continue
            sim_cmd  = self.cmd_part + ['-c', contig_id] + ['-n', str(self.count_dict[contig_id])]
            read_gen = LockedIterator(StreamHTSIM(sim_cmd=sim_cmd, pair_end=self.pair_end)) # only output 1 header for -c TODO:
            var_contig, sim_data= next(read_gen)                                    # the first element of generator is variants
            self.current_contig = var_contig                                        # update the profiles
            self.pos_map, self.meth_arr, _ = self.meth_db.load_contig(var_contig)   # [pos_map, meth_arr, status]
            self.variant_profile= self.meth_set.set_var_meth(var_contig, sim_data)  # a dict, can be empty
            
            with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
                if self.pair_end:                                                       # what if read_gen is empty at very beginning
                    job_arr = [executor.submit(self.process_read_pair, read_pair) for _, read_pair in read_gen]
                else:
                    job_arr = [executor.submit(self.process_read, read_pair) for _, read_pair in read_gen]

                for job in as_completed(job_arr):
                    try:
                        data = job.result()
                    except Exception as exc:
                        print('generated an exception: %s' % (exc))
                    else:
                        read_pair, read1_idx, pattern_idx = data
                        self.fastq_out.output_reads(read_pair, read1_idx, pattern_idx)
                
                # if self.verbose:
                #     for job in job_arr:
                #         job.add_done_callback(self.progress_bar)

        # if self.verbose:                    # close the progress bar, has problem
        #     self.tqdm_pbar.close()

        self.fastq_out.close()
        print('Simulation Finished!\n')


    def process_read(self, read_pair):
        '''process single end reads'''
        read1_idx   = random.choice([0, 1]) if read_pair[0]['strand']<0 else read_pair[0]['strand']
        pattern_idx = random.choice([0, 1]) if self.undirectional else read1_idx
        
        self.mask_context(read_pair[read1_idx], pattern_idx)
        self.retrive_meth_db(read_pair[read1_idx])
        self.set_context_state(read_pair[read1_idx])
        self.treat_bisulfite(read_pair[read1_idx])
        self.rev_complement(read_pair)
        self.add_seq_err(read_pair[read1_idx], pattern_idx)
        self.add_qual_score(read_pair[read1_idx])
        self.fastq_out.output_reads(read_pair, read1_idx, pattern_idx)


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

        # introduce seq errors
        self.add_seq_err(read_pair[0])
        self.add_seq_err(read_pair[1])
        
        # introduce quality scores
        self.add_qual_score(read_pair[0])
        self.add_qual_score(read_pair[1])

        # output
        self.fastq_out.output_reads(read_pair)


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


    def create_outdir(self):
        '''create output directory'''
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir, exist_ok=False)
        if not os.path.isdir(self.pkl_dir):
            os.makedirs(self.pkl_dir, exist_ok=False)
        if not os.path.isdir(self.tmp_dir):
            os.makedirs(self.tmp_dir, exist_ok=False)


    def test(self, read_pair):
        return read_pair


    # @property
    # def progress_bar(self):
    #     '''show the progress of read simulaiton'''
    #     with self.countLock:
    #         self.tqdm_count[0] += 2 if self.pair_end else 1

    #     incre_amount = self.tqdm_count[0] - self.tqdm_count[1]
    #     if incre_amount > self.tqdm_count[2]:
    #         self.tqdm_pbar.update(incre_amount)
    #         self.tqdm_count[1] = self.tqdm_count[0]

