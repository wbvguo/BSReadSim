import os
import random
import numpy as np
import subprocess
import multiprocessing
import threading

from Bio import SeqIO
from tqdm import tqdm
from scipy.stats import bernoulli
from typing import Dict, Union, Tuple
from threading import Lock
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from multiprocessing import Queue 

from DataProcessor import DataProcessor
from LockedIterator import LockedIterator
from SetMethylation import SetMethylation
from StreamReads import StreamReads
from StreamHTSIM import StreamHTSIM
from SetExperiment import SetExperiment
from ReadProcessor import ReadProcessor

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
    :param bool  overwrite_bed  : whether to overwrite the bed file for RRBS if it already exists [True]
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
    :param bool  random_error   : whether to generate sequencing error randomly [True]
    :param float error_rate     : sequencing error rate for random error generation[0.005]
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
                 meth_db_path: str = None, asm_file: str = None,                                    # ASM simulation
                 cgmap_file: str = None, cgmap_pool: bool = False, meth_seed: int = None,           # methylation ref
                 beta_params: Dict = {"CG":(0.5, 0.5), "CHG":(0.01, 0.05), "CHH":(0.01, 0.05)},     # beta distribution
                 collect_ch: bool = True, update_boundary: bool = False, overwrite_db: bool = False,# meth_db output
                 vcf_file: str = None,                                                              # genotype reference
                 mut_rate: float = 0.0010, haplo_mode: bool = False, seed: int = -1,                # mutation parameter
                 mut_indel_frac: float = 0.15, indel_ext_prob: float = 0.15,                        # indel parameter
                 insert_mean: int = 400, insert_std: int = 25,                                      # fragment setting
                 insert_min: int = 100, insert_max: int = 1000,
                 read_len: int = 100, depth: int = 20, num_reads: int = None,                       # reads setting
                 propN_cutoff: float = 0.05, conversion_rate: float = 0.998,                 
                 overwrite_bed:bool = True, middle_cuts: bool =True, cover_uniform:bool=True, rrbs_uniform:bool=True,
                 undirectional: bool = False, pair_end: bool = True,                                # sequencing protocol
                 random_error= True, error_rate: float = 0.005, error_model: str = None,                # sequencing error TODO:
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
        
        self.ref_fasta  = ref_fasta
        self.outdir     = outdir
        self.prefix     = prefix
        self.pair_end   = pair_end
        self.n_threads  = n_threads
        self.verbose    = verbose
        
        self.ref_dict   = SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))                            # can be saved
        
        # set experiment
        print('Initiating experiment...')
        self.experiment = SetExperiment(ref_dict=self.ref_dict, outdir=self.outdir, num_reads=num_reads,
                                        collect_ch = collect_ch, overwrite_bed= overwrite_bed, asm_file=asm_file,
                                        vcf_file = vcf_file, middle_cuts=middle_cuts, mut_rate=mut_rate,
                                        haplo_mode = haplo_mode, seed = seed, 
                                        mut_indel_frac= mut_indel_frac,indel_ext_prob= indel_ext_prob, 
                                        insert_mean=insert_mean,insert_std=insert_std,
                                        propN_cutoff=propN_cutoff,conversion_rate=conversion_rate,
                                        undirectional=undirectional, random_error=random_error,error_rate=error_rate,
                                        error_model=error_model,
                                        qual_uniform=qual_uniform,qual_model=qual_model,
                                        cover_uniform=cover_uniform,gc_bias_file=gc_bias_file,
                                        site_dependency=site_dependency,site_model=site_model,
                                        depth=depth, read_len=read_len, pair_end=pair_end, is_uniform=is_uniform,
                                        probe_bed_file=probe_bed_file, cut_site_str=cut_site_str, rrbs_model=rrbs_model,rrbs_uniform=rrbs_uniform,probe_deviation=probe_deviation,
                                        collect_stats=collect_stats,
                                        insert_min=insert_min, insert_max=insert_max)
        self.htsim_path = self.experiment.htsim_path
        self.htsim_opts = self.experiment.htsim_opts
        self.count_dict = self.experiment.count_dict
        self.num_reads  = self.experiment.num_reads # save for progress bar TODO:
        
        # prepare the methylation reference
        print('Initiating methylation profile...')
        self.meth_set   = SetMethylation(ref_dict=self.ref_dict, outdir=outdir,
                                         meth_db_path=meth_db_path,
                                         cgmap_file=cgmap_file, cgmap_pool=cgmap_pool,
                                         asm_file=asm_file,
                                         beta_params=beta_params,
                                         collect_ch=collect_ch,
                                         update_boundary=update_boundary,
                                         overwrite_db=overwrite_db,
                                         verbose=verbose,
                                         seed=meth_seed)
        self.meth_db    = self.meth_set.meth_db
        #self.meth_db.create_share_arr(self.meth_set.arr_max_size)
        
        # prepare output
        self.fastq_out  = StreamReads(outdir=outdir, prefix=prefix, pair_end=pair_end,
                                      shuffle=shuffle, seed_file=ref_fasta, gzip=gzip)
        
        # to hold intermediate data
        self.curr_contig= None
        self.var_profile= None
        # self.processor  = None
        # self.read_queue = None
        # self.write_queue= None
    
    
    def run(self, n_threads=None):
        """Simulating Bisulfite sequencing reads with parallelization:
        # 1. recieve from htsim
        # 2. process and write to disk (multiple threads)
        """
        self.n_threads = n_threads if n_threads else self.n_threads
        cmd_part = [self.htsim_path, self.ref_fasta] + [str(item) for key_val in self.htsim_opts.items() for item in key_val]
        
        print(f'Simulating methylated Reads with {self.n_threads} threads...')
        if self.verbose:
            print(f'[CMD]: {" ".join(cmd_part)}')
            print('[INFO]: #reads/#read pairs for each contig', end=" ")
            print(self.count_dict)
        
        # self.read_queue = Queue(maxsize=10**3)
        # self.write_queue= Queue(maxsize=10**3)
        
        # writer_thread = threading.Thread(target=self.write_fastq)
        # writer_thread.start()
        # with ThreadPoolExecutor(max_workers=self.n_threads) as executor:
        for contig_id in self.count_dict.keys():
            if self.count_dict[contig_id] == 0: # need to test if reads < self.n_threads TODO:
                continue
            sim_cmd  = cmd_part + ['-c', contig_id] + ['-n', str(self.count_dict[contig_id])]
            
            
            read_gen = LockedIterator(StreamHTSIM(sim_cmd=sim_cmd, pair_end=self.pair_end)) # only output 1 header for -c TODO:
            var_contig, sim_data= next(read_gen)                                    # the first element of generator is variants
            self.curr_contig= var_contig                                            # update the profiles
            self.var_profile= self.meth_set.set_var_meth(var_contig, sim_data)      # a dict, can be empty
            self.meth_db.load_contig_share(var_contig)                              # [pos_map, meth_arr, status]
            self.processor  = ReadProcessor(meth_arr= self.meth_db.shared_meth_arr,
                                            pos_map = self.meth_db.shared_pos_map,
                                            var_profile = self.var_profile,
                                            experiment  = self.experiment)
            
            self.data_processor = DataProcessor(processor=self.processor, fastq_out=self.fastq_out)
            self.data_processor.start()
            self.data_processor.spawn_workers()
            for _, read_pair in read_gen:
                self.data_processor.feed_data(read_pair)
            self.data_processor.stop()
            
            # with multiprocessing.Pool(processes=self.n_threads) as pool:
            #     pool.apply_async(self.process_read)
                
            #     for _, read_pair in read_gen:
            #         self.read_queue.put(read_pair)
            #         print(self.read_queue.qsize())
            #         print(self.write_queue.qsize())
            #     self.read_queue.put(None)
        
        # self.write_queue.put(None)
        # writer_thread.join() # wait for the writer thread to complete
        self.fastq_out.close()
        print('Simulation Finished!\n')
    
