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

from SetCyotsineMethylation import SetCytosineMethylation
from StreamWGSIM import StreamWGSIM
from UtilityFunctions import get_wgsim_path, reverse_complement


class SimulateMethylatedReads:
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
    simulation = SimulateMethylatedReads(ref_fasta, outdir, **kwargs)
    simulation.run()
    ```

    """

    def __init__(self, ref_fasta: str = None, outdir: str = None, prefix: str = "sim",          # required arguments
                 meth_db_path: str = None, cgmap_file: str = None, cgmap_pool: bool = False,    # methylation ref
                 asm_sim: bool = False, asm_file: str = None,                                   # ASM simulation
                 beta_params: Dict = {"CG":(0.5, 0.5), "CHG":(0.01, 0.05), "CHH":(0.01, 0.05)}, # beta distribution
                 collect_ch: bool = True, overwrite_db: bool = False,                           # methylation output
                 vcf_file: str = None,                                                          # genetic variant reference
                 mut_rate: float = 0.0010, haplo_mode: bool = False, seed: int = -1,            # mutation parameter
                 update_boundary: bool = False,
                 mut_indel_frac: float = 0.15, indel_ext_prob: float = 0.15,                    # indel parameter
                 insert_mean: int = 400, insert_std: int = 25,                                  # fragment setting
                 read_len: int = 100, depth: int = 20, num_reads: int = None,                   # reads setting
                 seq_err: float = 0.005, conversion_rate: float = 0.998,                        # rates setting
                 undirectional: bool = False, pair_end: bool = True, propN_cutoff: float = 0.05,# sequencing protocol
                 n_threads: int = 4, gzip: bool = True, verbose: bool = True, shuffle: bool = True):
        # not sure collect_sim_stats's role in previous bsbolt simulator
        # check required auguments are provided
        if not os.path.exists(ref_fasta):
            raise ValueError('Cannot find the reference file, please check!')
        if not outdir:
            raise ValueError('Please specify the output directory!')

        self.ref_fasta= ref_fasta
        self.outdir   = outdir
        self.prefix   = prefix
        self.gzip     = gzip
        self.seed     = seed
        self.fastq1   = None
        self.fastq2   = None

        # prepare the methylation reference
        print('Initiating methylation profile:\n')
        self.meth_db  = SetCytosineMethylation(ref_fasta=ref_fasta, outdir=outdir,
                                               meth_db_path=meth_db_path,
                                               cgmap_file=cgmap_file, cgmap_pool=cgmap_pool,
                                               asm_sim=asm_sim, asm_file=asm_file,
                                               beta_params=beta_params,
                                               collect_ch=collect_ch,
                                               overwrite_db=overwrite_db,verbose=verbose,
                                               seed = None if seed < 0 else seed)

        # prepare wgsim command
        wgsim_args    = [get_wgsim_path(), ref_fasta]
        genome_len    = self.meth_db.genome_len
        num_reads     = int(genome_len*depth/read_len/(1+int(pair_end))) if not num_reads else num_reads
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
        self.seq_err        = seq_err
        self.collect_ch     = collect_ch
        self.asm_sim        = asm_sim

        # to hold intermediate data
        self.current_contig = None
        self.pos_map        = None
        self.meth_arr       = None
        self.variant_profile= None
        self.variant_data   = {}
        self.read_count     = 0
        self.read_count_old = 0
        self.tqdm_step_size = num_reads * 0.01
        self.n_threads      = n_threads
        self.lock           = Lock()

        # prepare output
        self.output_obj     = self.get_output_obj
        self.shuffle        = shuffle
        self.verbose        = verbose
        if verbose:
            self.pbar= tqdm(range(self.num_reads), bar_format="{desc:<5.5}{percentage:3.0f}%|{bar:20}{r_bar}")


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
                job_arr = [executor.submit(self.process_read_group, read_pair) for read_pair in read_gen]
                for job in job_arr:
                    job.add_done_callback(self.progress_bar)

        for output in self.output_obj:      # close the fastq object
            output.close()

        if self.verbose:                    # close the progress bar
            self.pbar.close()

        if self.shuffle or self.gzip:       # shuffle or gzip the reads if needed
            print("Shuffle/Gzip the reads...\n")
        for fastq in self.fastq_list:
            if self.shuffle:
                self.shuffle_read(fastq, self.ref_fasta)
            if self.gzip:
                self.gzip_read(fastq)

        print('Simulation Finished!\n')


    def process_read_group(self, read_pair):
        """
        This processing step works as follows:
        1. randomly assign reads to Watson or Crick strand, with corresponding base change pattern
        2. set methylation status according to methylation profile
        3. bisulfite converted and introduce sequencing error
        4. output reads
        """
        pattern_list= [('C', 'T'), ('G', 'A')]
        ref_strand  = random.choice([0, 1]) # randomly select ref strand ['W', 'C']
        sub_pattern = pattern_list[ref_strand]

        if self.pair_end:
            # retrive the methy profile
            if ref_strand: # G to A
                tmp_meth = self.retrive_meth_db(read_rec, sub_pattern)
                tmp_meth = self.set_methylation(read_rec, sub_pattern)
            else: #C to T
                pass
            
        # set methylation
        
        # bisulfite converted
        
        # introduce seq errors
        
        # introduce quality scores
        
        # output
        if self.pair_end:
            if ref_strand:
                sim_data[0], sim_data[1] = sim_data[1], sim_data[0]        #??? need to check
            sim_data[0] = self.treat_bisulfite(sim_data[0], sub_pattern)
            sim_data[1] = self.treat_bisulfite(sim_data[1], sub_pattern)
        else:
            sim_data[0] = self.treat_bisulfite(sim_data[0], sub_pattern)

        # switch subpattern randomly for output if undirectional           #??? need to check
        if self.undirectional and random.choice([0, 1]):
            sub_pattern = ('C', 'T') if sub_pattern[0] == 'G' else ('G', 'A')
            if self.pair_end:
                sim_data[0], sim_data[1] = sim_data[1], sim_data[0]
        self.output_sim_reads(sim_data, sub_pattern[0], ref_strand)    


    def retrive_meth_db(self, read_rec, sub_pattern):
        '''input the cigar string ofs etc, output the values per site'''
        read_ctx = read_rec['ctx']
        if sub_pattern[0] == "C":   # 1, 3, 7
            site_idx = (read_ctx<8) if self.collect_ch else (read_ctx==1)
            read_rec['ctx2'] = np.ma.masked_greater(read_ctx, 8)
        else:                       # 9, 11, 15
            site_idx = (read_ctx>8) if self.collect_ch else (read_ctx==9)
            read_ctx[read_ctx < 8] = 0
            read_rec['ctx2'] = np.ma.masked_less(read_ctx, 8)

        read_pos = read_rec['start'] + np.arange(self.read_len)
        read_meth= np.zeros(self.read_len)

        arr_idx = 2
        if not read_rec['flag_pos']:                # SNP free region match
            match_site = site_idx
            read_meth[match_site]= self.fetch_meth_val(read_pos[match_site], arr_idx, 0)
        else:                                       # contain mutation pos in DNA fragment
            if self.asm_sim:
                arr_idx = 4 if read_rec['flag_mut'] else 3

            match_idx = read_rec['cgr'] == 0        # match
            match_site= match_idx & site_idx
            read_meth[match_site]   = self.fetch_meth_val(read_pos[match_site],arr_idx, 0)
            if read_rec['n_sub']:                   # snp
                snp_idx   = read_rec['cgr'] == 1
                snp_site  = snp_idx & site_idx
                read_meth[snp_site] = self.fetch_meth_val(read_pos[snp_site],  arr_idx, 1)
            if read_rec['n_indel']:                 # indel
                read_pos += read_rec['ofs']
                indel_idx = read_rec['cgr'] == 3
                read_meth[indel_idx]= self.fetch_meth_val(read_pos[indel_idx], arr_idx, 3)

        return read_meth, site_idx


    def fetch_meth_val(self, pos_arr, arr_idx, var_type):
        '''fetch the methylation value from meth_db'''
        if var_type == 1:           # SNP
            return [self.variant_profile[pos][0] for pos in pos_arr]
        elif var_type == 3:         # insertion sites will have the same coordinate
            insert_pos= [x for x in pos_arr if pos_arr.count(x) > 1]
            meth_val = []
            for pos in insert_pos:
                meth_val  += list(self.variant_profile[pos][0])
        else:
            return list(self.meth_arr[self.pos_map[pos_arr], arr_idx])


    def set_meth_state(self, read_rec,  read_meth, site_dependency = False):
        '''set methylation status based on the meth_arr'''
        if site_dependency:
            pass
            # generate methylation pattern according to read_meth and distance
        else:
            meth_states = bernoulli.rvs(read_meth, size = len(read_meth))
            meth_change = np.where(meth_states!=0)
            read_rec['ctx2'][meth_change]+= 1 # methybases increased by 1


    def treat_bisulfite(self, read_rec):
        """ unmethylated C conversion """
        conv_pos = np.where(np.bitwise_and(read_rec['ctx2'], 0x1))
        conv_change = [i for i in bernoulli.rvs(self.conversion_rate, size=len(conv_pos))]
        read_rec['seq'][conv_pos] = np.bitwise_and(read_rec['seq'][conv_change] + 2, 0x3) # C2T, G2A


    def add_seq_err(self, read_rec, random_err = True):
        ''' introduce sequencing error and quality scores'''
        if random_err:
            err_arr = bernoulli.rvs(self.seq_err, size = self.read_len)
            err_idx = np.where(err_arr == 1)
            if err_idx:
                base_set = {0,1,2,3}
                for idx in err_idx:
                    base = read_rec[idx]
                    read_rec['seq'][idx] = np.ramdon.choice(base_set.difference(base), size = 1)
                    read_rec['ctx2'] = 
                
                # change
                # update cigar string
        else:
            pass

    def add_qual_score(self, read_rec, qual_uniform = True):
        '''add quality scores for the read'''
        if qual_uniform:
            qual_arr = np.full(self.read_len, read_rec['qual'])
            read_rec['qual'] = qual_arr
        else:
            pass


    def output_reads(self, read_rec):
        ''' output reads '''
        pass


    def output_sim_reads(self, sim_data, sub_base, ref_strand):
        """Write simulated bisulfite reads"""
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
    def get_output_obj(self):
        """Return io object for fastq writing"""
        self.fastq1 = f'{self.outdir}/{self.prefix}_1.fastq'
        if self.pair_end:
            self.fastq2 = f'{self.outdir}/{self.prefix}_2.fastq'


    @staticmethod
    def shuffle_read(fastq_file: str = None, random_source: str = None):
        '''
        shuffle the reads randomly, if false, the reads will be segemented by contigs
        the number of simulated reads should not exceed the number of bits in the reference genome.
        for HG, should be less than 8*3G < 2.4*10^10 (average DEPTH should be less than read_len*8)
        from link: https://www.biostars.org/p/9764/
        '''
        prefix_split    = os.path.splitext(fastq_file)[0].split("_")
        fastq_shuffle   = "_".join(prefix_split[:-1]) + "_shuffle_" + prefix_split[-1] + ".fastq"
        shuf_cmd_list   = ["awk", "'{OFS=\"\t\"; getline seq; getline sep; getline qual; print $0,seq,sep,qual}'",
                         fastq_file, "|", "shuf --random-source", random_source, "|",
                         "awk", "'{OFS=\"\n\"; print $1,$2,$3,$4}'", ">", fastq_shuffle]
        shuffle_run     = subprocess.Popen(shuf_cmd_list,  stdout=subprocess.PIPE, universal_newlines=True)
        stdout, stderr  = shuffle_run.communicate()
        rename_cmd_list = ["mv", fastq_shuffle, fastq_file]
        rename_run      = subprocess.Popen(rename_cmd_list,stdout=subprocess.PIPE, universal_newlines=True)
        stdout, stderr  = rename_run.communicate()

    @staticmethod
    def gzip_read(fastq_file):
        '''gzip the output reads to fastq.gz format'''
        gzip_cmd_list   = ["gzip", fastq_file]
        gzip_run        = subprocess.Popen(gzip_cmd_list,  stdout=subprocess.PIPE, universal_newlines=True)
        stdout, stderr  = gzip_run.communicate()

    @property
    def progress_bar(self):
        '''show the progress of read simulaiton'''
        if self.verbose:
            with self.lock:
                self.read_count += 2 if self.pair_end else 1

            incre_amount = self.read_count - self.read_count_old
            if incre_amount > self.tqdm_step_size:
                self.pbar.update(incre_amount)
                self.read_count_old = self.read_count
