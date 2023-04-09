#!/usr/bin/env python
# coding: utf-8

# In[1]:


import gzip
import io
import pysam
import numpy as np
import argparse
from typing import List, Union, Dict, NamedTuple, Tuple
from tqdm import tqdm
from collections import defaultdict

class OpenFastq:
    """ Simple class to simplify iterating through fastq files. The script yields a tuple for every four lines
    in a fastq file
    ------------------------------------------------------------------------------------
    KeywordArguments:
        fastq: str = path to fastq file
    """

    def __init__(self, fastq: str = None):
        if fastq.endswith(".gz"):
            self.f = io.BufferedReader(gzip.open(fastq, 'rb'))
        else:
            self.f = open(fastq, 'r')

    def __iter__(self) -> List[str]:
        with self.f as fastq:
            while True:
                line1 = fastq.readline()
                if not line1.strip():
                    break
                line1 = self.process_line(line1)
                line2 = self.process_line(fastq.readline())
                line3 = self.process_line(fastq.readline())
                line4 = self.process_line(fastq.readline())
                yield [line1, line2, line3, line4]

    @staticmethod
    def process_line(line: Union[bytes, str]):
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '')
        else:
            return line.replace('\n', '')


# In[3]:


def get_correct_bsseeker_name(read, undirectional=True):
    name = read.qname
    bs_id, conversion = None, None
    strand = 'C' if read.is_reverse else 'W'
    for tag in read.tags:
        if tag[0] == 'XO':
            bs_id = tag[1]
            break
    if bs_id == '+FR' or bs_id == '-FR':
        name = f'{name.split("/")[0]}/1' if read.is_read1 else f'{name.split("/")[0]}/2'
    return name, f'{strand}'


# In[67]:


def get_read_reference_info(fastq_files: List[str] = None):
    reference_info = {}
    fastq_iterators = (OpenFastq(file) for file in fastq_files)
    for fastq in fastq_iterators:
        for line in fastq:
            read_info = line[0].replace('@', '').strip()
            read_name = read_info
            read_string = line[1].strip()
            reference_info[read_name] = read_info + f'string{read_string}'
    return reference_info

def parse_alignment_comment(read_comment):
    read_comment, read_string = read_comment.split("string")
    read_comment = read_comment.rsplit(":", 1)
    chrom = read_comment[0].split("_", maxsplit = 1)[1]
    start, end = read_comment[1].split("_")[0].rsplit("-", 1)
    ref_strand = "C"
    return dict(chrom=chrom, start=int(start),
                end=int(end), ref_strand=ref_strand, read_len = len(read_string))


class AlignmentEvaluator:
    """Evaluate alignment against simulated bisulfite sequencing data.

    Params:

    * *duplicated_regions (dict)*: regions duplicated in the simulation reference, [None]
    * *matching_target_prop (float)*: proportion of alignment that most overlap with target region
                                     for a valid alignment to be called, [0.95]
    """

    def __init__(self, duplicated_regions: Dict[str, Tuple[int, int]] = {}, matching_target_prop: float = 0.95,
                 verbose: bool = False):
        self.duplicated_regions = duplicated_regions
        self.matching_target_prop = matching_target_prop
        self.verbose = verbose

    def evaluate_alignment(self, alignment_file: str, 
                           fastq_files: List[str] = None,
                           tool: str = 'BSBolt',
                           MAPQ_threshold = 30) -> Dict[str, int]:
        """

        Params:

        * *alignment_file (str)*: path to alignment file
        * *fastq_files (list)*: list of paths to fastq files

        Returns:

        * *alignment_evaluations (dict)*: target alignment stats
        """
        alignment_evaluations = defaultdict(int)
        off_target_reads = {}
        reads_observed = []
        reference_info = get_read_reference_info(fastq_files=fastq_files)
        print("processing aligned reads")
        for alignment in tqdm(pysam.AlignmentFile(alignment_file, "rb"), disable=True if not self.verbose else False):
            strand = None
            read_name = alignment.qname
            if tool == 'bsseeker':
                read_name, strand = get_correct_bsseeker_name(alignment)
            reads_observed.append(read_name)
            alignment_evaluations['Observed'] += 1
            alignment_info = parse_alignment_comment(reference_info[read_name])
            dup_region = False
            if alignment_info['chrom'] in self.duplicated_regions:
                duplicate_range = self.duplicated_regions[alignment_info['chrom']]
                if duplicate_range[0] <= alignment_info['start'] <= duplicate_range[1]:
                    dup_region = True
            dup = 'Dup' if dup_region else 'NoDup'
            if alignment.is_qcfail:
                alignment_evaluations['QC Fail'] += 1
            elif not alignment.is_unmapped:
                if not strand:
                    strand = 'C' if alignment.is_reverse else 'W'
                chrom_match, matching_prop = self.assess_alignment(alignment, alignment_info)
                on_target = 'On' if chrom_match and matching_prop >= self.matching_target_prop else 'Off'
                secondary = 'Sec' if any([alignment.is_secondary, alignment.is_supplementary]) else 'Prim'
                pair = 'PropPair' if alignment.is_proper_pair else 'Discord'
                #strand_match = 'StrandMatch' if strand == alignment_info['ref_strand'] else 'StrandMismatch'
                strand_match = 'StrandMatch'
                unique = 'Unique' if alignment.mapping_quality >= MAPQ_threshold else 'Multiple'
                if on_target == 'On':
                    on_target = 'On' if strand_match == 'StrandMatch' else 'Off'
                alignment_evaluations[f'{on_target}_{secondary}_{unique}_{pair}_{dup}_{strand_match}'] += 1
            else:
                off_target_reads[read_name] = f'unaligned_{dup}'
                alignment_evaluations['Unaligned'] += 1
        print("summarizing results")
        alignment_evaluations['ReadObserved'] = len(reads_observed)
        alignment_evaluations['SimulatedReads'] = len(reference_info.keys())
        alignment_evaluations['Unobserved'] = alignment_evaluations['SimulatedReads'] - alignment_evaluations['ReadObserved']
        alignment_evaluations['PercentageUniquelyAlignedPrim'] = 100*(np.round(alignment_evaluations['On_Prim_Unique_PropPair_NoDup_StrandMatch']/alignment_evaluations["Observed"], 2))
        return alignment_evaluations, off_target_reads

    @staticmethod
    def assess_alignment(alignment: pysam.AlignedSegment, alignment_info: Dict):
        """ Compare alignment against reference alignment"""
        chrom_match = alignment.reference_name == alignment_info['chrom']
        # assess reference bases that match between the two reads
        if not chrom_match:
            return chrom_match, 0.0
        else:
            matching_prop = alignment.get_overlap(alignment_info['start'], alignment_info['end'])/alignment_info['read_len']
            return chrom_match, matching_prop


# In[5]:


parser = argparse.ArgumentParser(description='Evaluate read correct mapping statistics')


# In[80]:


parser.add_argument('-d', '--duplicated_regions', default = {}, metavar = '', help = "Dict[str, Tuple[int, int]], [{}]")
parser.add_argument('-m', '--matching_target_prop', type = float, default = 0.95, 
                    help = "proportion of alignment that must overlap with target region for a valid alignment to be called, [0.95]", 
                    metavar = '')
parser.add_argument('-v', '--verbose', default = False, action = "store_true", help = "verbose evaluation")
parser.add_argument('-f', '--fastq_files', type = lambda file: [file_path for file_path in file.split(',')], required = True,
                   help = "comma separated list of fastq files", metavar = '')
parser.add_argument('-a', '--alignment_file', type = str, required = True, help = 'path for alignment .bam file'
                   , metavar = '')
parser.add_argument('-t', '--tool', type = str, help = 'aligner name, one of bwa_meth, BSBolt, biscuit, bsseeker, bismark, [BSBolt]', 
                    default = "BSBolt", metavar = '')
parser.add_argument('-q','--mapq_threshold', type = float, default = 30, 
                    help = 'MAPQ threshold for an alignment to be considered unique, [30]', metavar = '')


# In[85]:


args = parser.parse_args()
duplicated_regions = args.duplicated_regions
matching_target_prop = args.matching_target_prop
verbose = args.verbose
fastq_files = args.fastq_files
alignment_file = args.alignment_file
tool = args.tool
MAPQ_threshold = args.mapq_threshold


# In[57]:


evaluator = AlignmentEvaluator(duplicated_regions=duplicated_regions,
                              matching_target_prop=matching_target_prop,
                              verbose=verbose)
alignment_evaluations, off_target_reads = evaluator.evaluate_alignment(alignment_file = alignment_file,
                                                                      fastq_files = fastq_files,
                                                                       tool = tool,
                                                                       MAPQ_threshold = MAPQ_threshold)


# In[86]:


file = open(f"{tool}_Alignment_statistics.txt", "w")
for key, value in alignment_evaluations.items():
    file.write('%s:%s\n' % (key, value))
file.close()
file = open(f"{tool}_Off_target_reads.txt", "w")
for key, value in off_target_reads.items():
    file.write('%s:%s\n' % (key, value))
file.close()


# In[68]:


# evaluator = AlignmentEvaluator(duplicated_regions={},
#                               matching_target_prop=0.9,
#                               verbose=True)


# In[69]:


# alignment_evaluations, off_target_reads = evaluator.evaluate_alignment(alignment_file = "/u/scratch/h/hz991224/sherman_reads/simulated/simulated_1_bismark_bt2_pe.bam",
#                                                                       fastq_files = ["/u/scratch/h/hz991224/sherman_reads/simulated_1.fastq", "/u/scratch/h/hz991224/sherman_reads/simulated_2.fastq"],
#                                                                        tool = "bismark",
#                                                                        MAPQ_threshold = 30)
