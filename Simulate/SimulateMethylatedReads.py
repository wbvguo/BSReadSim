import os
import random
import numpy as np

from Bio import SeqIO
from tqdm import tqdm
from scipy.stats import bernoulli
from typing import Dict, Union, Tuple

from Simulate.SetCyotsineMethylation import SetCytosineMethylation
from Simulate.StreamSim import StreamWGSIM
from Utils.UtilityFunctions import get_external_paths, reverse_complement


class SimulateMethylatedReads:
    """Bisulifite read simulation class. The class works as follows:
    
    1. WGSIM (forked and modified version) is called to simulate paired end Illumina reads
        - If run in single end mode, the number of read simulated is double to get desired coverage and the second read
          isn't process as a bisulfite read or output.
    2. Methylation values are set for all methylatable bases (Cytosine and Guanine relative to the reference)
        - Values can be set randomly or taken from a reference file (bsbolt simulation database or CGmap file)
    3. Reads are bisulfite converted (with a conversion rate) and output

    Params:
    
       * *reference_file (str)*: path to reference genome fasta file.
       * *sim_output (str)*: output path/directory
       * *sequencing_error (float)*: simulated sequencing error rate, [0.005]
       * *mutation_rate (float)*: simulated mutation error rate, [0.0010]
       * *mutation_indel_fraction (float)*: fraction of mutations that are INDELs, [0.15]
       * *indel_extension_probability (float)*: probability INDEL length will be extended, [0.15]
       * *random_seed (int)*: random seed for mutation and sequencing error generation, [-1]
       * *paired_end (bool)*: simulate paired end bisulfite sequencing data. [False]
       * *read_length (int)*: length of simulated reads, [100]
       * *read_depth (int)*: average read depth over simulated contigs, [20]
       * *undirectional (bool)*: simulate undirectional (PCR product of Watson and Crick strands), [False]
       * *methylation_reference (str)*: path to previously generated bsbolt reference directory
       * *cgmap (str)*: path to CGmap file to use as methylation reference
       * *ambiguous_base_cutoff (float)*: reference segments where the proportion of ambiguous bases, - or N, greater
                                          than threshold will be skipped, [0.05]
       * *haplotype_mode (bool)*: simulate only homozygous variants, [False]
       * *pe_fragment_size (int)*: maximum fragment size, [400]
       * *insert_deviation (int)*: standard deviation of simulated insert sizes, [25]
       * *mean_inner_dist (int)*: mean insert size, [100]
       * *collect_ch_sites (bool)*: simulated and collect CH methylation sites, [True]
       * *collect_sim_stats (bool)*: output simulated bases for all collected methylatable bases, [False]
       * *verbose (bool)*: verbose output, [True]
       * *overwrite_db (bool)*: overwrite previously generated bsbolt simulated database, [False]

    Usage:
    ```python
    simulation = SimulateMethylatedReads(reference_file, sim_output, **kwargs)
    simulation.run()
    
    simulation = getSimulationObject(reference_file, sim_out, **kwargs)
    sim_meth_db= setMethylationProfile(**kwargs)
    simylation.update(sim_meth_db) # if not set, the output is WGS, not WGBS
    simulation.run()
    ```

    """

    def __init__(self, reference_file: str = None, sim_output: str = None,                               # required arguments
                 cgmap_file: str = None, meth_db_file: str = None,                                       # methylation reference (CGmap.gz)
                 meth_beta_param: Dict = {"CG": (0.5, 0.5), "CHG": (0.01, 0.05), "CHH":(0.01, 0.05)},    # methylation parameter for beta distribution
                 vcf_file: str = None,                                                                   # genetic variant reference (vcf.gz format)
                 mutation_rate: float = 0.0010, haplotype_mode: bool = False, random_seed: int = -1,     # genetic variant parameter
                 mutation_indel_fraction: float = 0.15, indel_extension_probability: float = 0.15,
                 pe_fragment_size: int = 400, insert_deviation: int = 25, mean_inner_dist: int = 100,    # fragment setting
                 read_length: int = 100, read_depth: int = 20, sequencing_error: float = 0.005,          # reads setting
                 undirectional: bool = False, paired_end: bool = True, conversion_rate: float = 0.998,   # sequencing protocol
                 ambiguous_base_cutoff: float = 0.05, 
                 verbose: bool = True, collect_ch_sites: bool = True, collect_sim_stats: bool = False, overwrite_db: bool = False):
        
        # check required auguments are provided
        if not os.path.exists(reference_file):
            raise ValueError('Cannot find the reference file, please check!')
        if not sim_output:
            raise ValueError('Please specify the output directory!')
        
        # parse reference
        self.reference_file = reference_file
        self.reference_dict = SeqIO.to_dict(SeqIO.parse(reference_file, "fasta"))             #??? can it process fa.gz
        
        # prepare wgsim command
        wgsim_args     = [get_external_paths()[1], reference_file] # second element is wgsim path
        genome_length  = sum([len(seq) for key, seq in self.reference_dict.items()])
        number_reads   = int(genome_length * read_depth / read_length / (1 + int(paired_end)))
        wgsim_options  = {'-1': read_length, '-2': read_length, '-N': number_reads,
                          '-r': mutation_rate, '-h': int(haplotype_mode), '-S': random_seed, 
                          '-R': mutation_indel_fraction,'-X': indel_extension_probability,
                          '-d': pe_fragment_size, '-s': insert_deviation, '-I': mean_inner_dist,
                          '-e': 0, '-A': ambiguous_base_cutoff, '-g': vcf_file} #set -e to be 0
        self.wgsim_args     = wgsim_args
        self.wgsim_options  = wgsim_options
        
        # prepare methylation profile
        self.meth_options   = {'cgmap': cgmap_file, 'meth_ref': meth_db_file, 'beta_param': meth_beta_param}
        
        # sequencing settings
        self.paired_end     = paired_end
        self.undirectional  = undirectional
        self.conversion_rate= conversion_rate
        self.seq_error      = sequencing_error
        
        # simulation process setting
        self.tqdm_disable   = not verbose
        self.overwrite_db   = overwrite_db
        self.collect_sim_stats = collect_sim_stats
        self.collect_ch_sites  = collect_ch_sites
        
        # to hold intermediate data
        self.current_contig = None
        self.contig_variant = None
        self.contig_profile = None
        self.contig_values  = {}
        self.variant_data   = {}
        
        
        # prepare output
        self.sim_output = sim_output
        self.output_objects = self.get_output_objects
        
        
    def run(self):
        """Simulated bisulfite sequencing reads"""
        print('Generating methylation profile:\n')
        self.sim_meth_db    = SetCytosineMethylation(reference_file=self.reference_file,
                                                     sim_output=self.sim_output,
                                                     meth_ref=self.meth_options['meth_ref'],
                                                     cgmap=self.meth_options['cgmap'],
                                                     beta_param = self.meth_options['beta_param'],
                                                     collect_ch_sites=self.collect_ch_sites,
                                                     overwrite_db=self.overwrite_db)
        
        # simulation
        print('Simulating methylated Reads:\n')
        self.sim_command    = self.wgsim_args + [str(item) for key_val in self.wgsim_options.items() for item in key_val]
        print(f'[CMD]: {" ".join(self.sim_command)}\n')
        self.simulate_methylated_reads()
        print('Simulation Finished!')
        
        
    def simulate_methylated_reads(self):
        """
        Read processing loop
        """
        for variant_contig, sim_data in tqdm(StreamWGSIM(sim_command=self.sim_command),
                                             desc='Simulating Bisulfite Converted Read Pairs',
                                             position=0, leave=True,
                                             disable=self.tqdm_disable):
            if variant_contig:
                self.current_contig = variant_contig
                self.contig_variant = sim_data
                self.contig_profile = self.get_methylation_reference(variant_contig, sim_data) #sim_data can be none
                #TODO: self.save_contig_data()
            else:
                assert sim_data[0][0]['chrom'] == self.current_contig
                self.process_read_group(sim_data)
        
        # close the fastq object
        for output in self.output_objects:
            output.close()
        
        
    def process_read_group(self, sim_data, stranded_capture=False):
        """Set read methylation values, randomly assign reads to Watson or Crick strand, bisulfite conversion"""
        ref_strand  = 'W' # watson
        sub_pattern = ('C', 'T')
        strand_swtich = (not stranded_capture) and bernoulli.rvs(0.5) # randomly select reference strand
        
        if strand_swtich:
            ref_strand  = 'C' # crick
            sub_pattern = ('G', 'A')
            sim_data[0], sim_data[1] = sim_data[1], sim_data[0]                         #??? need to check
        
        # set read methylation, bisulfite conversion, add sequencing error
        sim_data[0] = self.set_read_methylation_bisulfite(sim_data[0], sub_pattern)
        if self.paired_end:
            sim_data[1] = self.set_read_methylation_bisulfite(sim_data[1], sub_pattern)
        
        # switch subpattern randomly for output if undirectional                        #??? need to check
        if self.undirectional and bernoulli.rvs(0.5):
            sub_pattern = ('C', 'T') if sub_pattern[0] == 'G' else ('G', 'A')
            if self.paired_end:
                sim_data[0], sim_data[1] = sim_data[1], sim_data[0]
        self.output_sim_reads(sim_data, sub_pattern[0], ref_strand)
        
        
    def set_read_methylation_bisulfite(self, read, sub_pattern):
        """
        1. Set methylation according to sim value, variants can be methylated
        2. bisulfite conversion
        3. add sequencing error (independent from bisulfite conversion)
        """
        read_seq = read[1]
        read_info= read[0]

        self.seq_array  = np.array(list(read_seq))
        self.cigar_array= np.array(list(read_info['cigar']))
        self.start  = read_info['start'] 

        sub_base = sub_pattern[0]
        if sub_base == "C":
            meth_base_info = read_info['c_base_info']
        else:
            meth_base_info = read_info['g_base_info']
        

        ref_pos_dict= dict()
        # if there is methylable base
        if len(meth_base_info):
            seq_pos_offset_array = np.array([site.split("_") for site in meth_base_info.split(",")[:-1]]).transpose().astype(int)
            seq_pos_array = seq_pos_offset_array[0]                    # methylable base position on the read
            offset_array  = seq_pos_offset_array[1]
            ref_pos_array = self.start + seq_pos_array + offset_array # methylable base position w.r.t the chromosome
            ref_pos_dict  = dict(zip(seq_pos_array, ref_pos_array))

        
        # handle match
        cigar_M_idx = np.where(self.cigar_array == "M")[0] # match
        cigar_X_idx = np.where(self.cigar_array == "X")[0] # substitution
        cigar_I_idx = np.where(self.cigar_array == "I")[0] # insertion
        seq_pos_array  = np.array(list(ref_pos_dict.keys()))
        
        cigar_M_cg_idx = np.intersect1d(cigar_M_idx, seq_pos_array)
        cigar_X_cg_idx = np.intersect1d(cigar_X_idx, seq_pos_array)
        cigar_I_cg_idx = np.intersect1d(cigar_I_idx, seq_pos_array)
        
        for idx in cigar_M_cg_idx:
            meth_pos = ref_pos_dict[idx]
            meth_status, meth_context = self.set_base_methylation(meth_pos)
            if meth_status:
                cigar_status = 'C' if meth_context == "CG" else "Y"
            else:
                cigar_status = 'c' if meth_context == "CG" else "y"
            self.cigar_array[idx] = cigar_status
        
        # for idx in cigar_X_idx:
        #     if self.seq_array[idx] == sub_base:
        #         ref_pos  = ref_pos_dict[idx]
        #         ref_base = self.reference_dict[self.current_contig][ref_pos]
        #         seq_base = self.seq_array[idx]
        #         meth_pos = f'{ref_pos}_{ref_base}_{seq_base}'
        #         cigar_status = self.set_variant_methylation(meth_pos)
        #         self.cigar_array[idx] = cigar_status
        #     elif self.seq_array[idx] == "T" and sub_base == "C":
        #         self.cigar_array[idx] = 'V'
        #     elif self.seq_array[idx] == "A" and sub_base == "G":
        #         self.cigar_array[idx] = 'V'
        #     else:
        #         pass
        
        # for idx in cigar_I_cg_idx:
        #     meth_pos = f'{ref_pos}_{self.cigar_array[idx]}'
        #     cigar_status = self.set_variant_methylation(meth_pos, variant_type='I')
        #     self.cigar_array[idx] = cigar_status
            
        
        # bisulfite conversion 
        self.read_bisulfite_convert(sub_pattern)
        # add sequencing errors
        self.add_sequencing_error()
        
        # convert to string for output
        read[1] = ''.join(list(self.seq_array))
        read[0]['cigar'] = ''.join(list(self.cigar_array))
        read[0]['sub_base'] = sub_base
        return read
        
        
    def set_base_methylation(self, meth_pos) -> Tuple[str, str]:
        """set methylation satus for the matched methyable base"""
        # nucleotide, methylation_level, context, 0, 0, 0
        try:
            meth_info = self.contig_profile[meth_pos]
        except KeyError:
            return False, 0
        meth_status = bernoulli.rvs(meth_info[2])
        meth_context= meth_info[1]
        return meth_status, meth_context
        
        
    def set_variant_methylation(self, meth_pos, variant_type='X') -> Tuple[str, str]:
        cigar_status = 'Z' if variant_type == 'X' else 'R'
        meth_status, meth_context = self.set_base_methylation(meth_pos)
        if meth_status:
            return cigar_status
        return cigar_status.lower()
        
        
    def read_bisulfite_convert(self, sub_pattern):
        '''substitute the pattern by a probability'''
        candidate_pos_cg = np.where(self.cigar_array == "c")[0]
        candidate_pos_ch = np.where(self.cigar_array == "y")[0]
        candidate_pos = np.append(candidate_pos_cg, candidate_pos_ch)

        self.cigar_array[candidate_pos_cg] = "b"
        self.cigar_array[candidate_pos_ch] = "d"

        converted_res = [sub_pattern[i] for i in bernoulli.rvs(self.conversion_rate, size=len(candidate_pos))]
        self.seq_array[candidate_pos] = converted_res
    
    
    def add_sequencing_error(self):
        error_flag_list= np.array(bernoulli.rvs(self.seq_error, size = len(self.seq_array)))
        error_idx_list = np.where(error_flag_list == 1)[0]
        base_set = set(['A', 'G', 'C', 'T'])
        for idx in error_idx_list:
            candidate_base_set = base_set - set(self.seq_array[idx])
            self.seq_array[idx] = random.sample(candidate_base_set, k=1)[0]
            self.cigar_array[idx] = "E"
        
        
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
        self.output_objects[0].write(read)
        if self.paired_end:
            read_label = f'@{sim_data[1][0]["read_id"]}_{sim_data[1][0]["chrom"]}/2'
            read_comment = f'+{sim_data[1][0]["chrom"]}:{sim_data[1][0]["start"]+1}:' \
                           f'{sim_data[1][0]["end"]}:{sim_data[1][0]["cigar"]}:{ref_strand}{conversion_2}'
            read = f'{read_label}\n{sim_data[1][1].upper()}\n{read_comment}\n{sim_data[1][3]}\n'
            self.output_objects[1].write(read)
    
    
    def get_methylation_reference(self, contig: str, variant_data:  Dict= False) -> Dict:
        """ Set variant methylation if variant data provided and return methylation profile else
        return methylation profile

        Params:

        * *contig (str)*: contig id
        * *variant_data (dict)*: simulated variant information

        Returns:

        * *contig_profile (Dict[str, float])*: methylation reference values"""
        #if variant_data:
        self.contig_profile = self.sim_meth_db.get_contig_methylation(contig)
            # self.sim_meth_db.set_variant_methylation(variant_data, self.contig_profile, self.current_contig)
        return self.contig_profile
        #return self.sim_meth_db.get_contig_methylation(contig)
    
    
    @property
    def get_output_objects(self):
        """Return io object for fastq writing"""
        output_list = [open(f'{self.sim_output}_1.fq', 'w')]
        if self.paired_end:
            output_list.append(open(f'{self.sim_output}_2.fq', 'w'))
        return output_list

