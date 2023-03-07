import pickle
import re
import pysam
import numpy as np
from typing import Any, Dict, List, Tuple, Union

class CallMethVector:
    """
     return vector of methylated sites
    """

    def __init__(self, bam_file: str = None, genome_db: str = None, 
                 contig: str = None, start=None, end=None, cg_only: bool = False, 
                 min_base_quality: int = 20, min_mapping_quality: int = 30, 
                 return_queue=None, filter_duplicates=True, fully_contain=True):
        
        self.bam_file  = bam_file
        self.input_bam = pysam.AlignmentFile(self.bam_file, 'rb', require_index=True)
        self.genome_db = f'{genome_db}/'
        self.mate_flag_dict = {67: 131, 323: 387, 115: 179, 371: 435, 
                               131: 67, 387: 323, 179: 115, 435: 371}
        self.meth_state_dict= {('C','C'): [True, 1], ('C','T'): [True, 0], 
                               ('G','G'): [True, 1], ('G','A'): [True, 0]}
        
        self.contig = contig
        self.start  = start
        self.end    = end
        self.cg_only= cg_only
        self.fully_contain = fully_contain
        
        self.filter_duplicates  = filter_duplicates
        self.min_base_quality   = min_base_quality
        self.min_mapping_quality= min_mapping_quality
        
        self.chunk_size = 10000
        self.return_queue   = return_queue


    def call_meth(self):
        try:
            chrom_seq = self.get_ref_sequence(f'{self.genome_db}{self.contig}.pkl')
        except FileNotFoundError:
            self.return_queue.put([])
            print(f'{self.contig} not found in BSBolt DB, Methylation Calls for {self.contig} skipped. Methylation '
                  f'values should be called using the same DB used for alignment.')
            self.return_queue.put([])
        else:
            self.call_contig(chrom_seq)


    def align_filter(self, align):
        is_mapped  = (not align.is_unmapped)
        is_primary  = not (align.is_secondary or align.is_supplementary)
        is_qualified = (not align.is_qcfail) and align.mapping_quality >= self.min_mapping_quality
        not_duplicate = not (align.is_duplicate and self.filter_duplicates)
        is_fully_align = all(op in (0, 7, 8) for op, length in align.cigar)
        return is_mapped and is_primary and is_qualified and is_fully_align and not_duplicate


    def call_contig(self, chrom_seq: str):
        """Iterates through bam reads call methylation along vectors. When a read overlaps a site where methylation is
         called the site with a higher quality is taken. If overlapping sites with the same quality are observed
         the first observed site is reported. If an overlapping site is reported as a mismatch only the site with
         a methylation call is reported (this should be extremely rare but is observed in the test cases)
        """
        # set search pattern for only CG sites or all Cs
        pattern_list = ('CG', 'CG') if self.cg_only else ('C', 'G')
        # iterate through pileup
        contig_chunk = []
        meth_vectors = {}
        
        for aligned_read in self.input_bam.fetch(contig=self.contig, start=self.start, end=self.end, multiple_iterators=True):
            if not self.align_filter(aligned_read):
                continue
            if self.fully_contain and (aligned_read.reference_start < self.start or aligned_read.reference_end > self.end):
                continue # exclude the partially overlapped reads
            
            # get sequence around pileup site
            ref_seq = chrom_seq[aligned_read.reference_start - 1: aligned_read.reference_end + 2].upper()
            if aligned_read.get_tag('YS') == 'W_C2T': # only for directional library
                search_pattern, ref_base, strand = pattern_list[0], 'C', 'Watson'
                offset = 0
            else:
                search_pattern, ref_base, strand = pattern_list[1], 'G', 'Crick'
                offset = 1 if self.cg_only else 0
            
            search_res = [match.start() + offset + aligned_read.reference_start - 1 for match in 
                          re.finditer(search_pattern, ref_seq)]
            
            # if reference sequence not found proceed to next sequence
            if not search_res:
                continue
            
            meth_calls = self.call_vector(aligned_read, set(search_res), ref_base)
            
            if not meth_calls[0]:
                continue
            
            processed_vector = self.process_meth_vector(aligned_read, meth_calls, strand, meth_vectors)
            
            if processed_vector:
                contig_chunk.append(processed_vector)
            
            if len(contig_chunk) == self.chunk_size:
                self.return_queue.put((self.contig, contig_chunk))
                contig_chunk = []

        # process reads that didn't have a pair with a observed methylation site
        for call in meth_vectors.values():
            contig_chunk.append((call['read_name'], call['calls'][1][0], call['calls'][1][-1],
                                 np.asarray(call['calls'][0]), np.asarray(call['calls'][1]), call['flag'],
                                 self.mate_flag_dict[call['flag']], call['strand']))
        if contig_chunk:
            self.return_queue.put((self.contig, contig_chunk))


    def call_vector(self, aligned_read: pysam.AlignedRead, position_set: set, ref_base: str) -> List:
        meth_calls   = [[], [], []]
        reference_consumers = {0, 2, 3, 7, 8}
        query_consumers     = {0, 1, 4, 7, 8}
        # set relative to genomic position so add reference start and one since capturing the first base
        ref_pos     = aligned_read.reference_start
        query_seq   = aligned_read.query_sequence
        query_qual  = aligned_read.query_qualities
        query_pos   = 0
        
        for cigar_type, cigar_count in aligned_read.cigartuples:
            if cigar_type in reference_consumers and cigar_type in query_consumers:
                for _ in range(cigar_count):
                    if ref_pos in position_set:
                        pos_qual = query_qual[query_pos]
                        if pos_qual > self.min_base_quality:
                            query_base = query_seq[query_pos]
                            call_made, meth_state = self.get_meth_call(ref_base, query_base)
                            if call_made:
                                meth_calls[0].append(meth_state)
                                meth_calls[1].append(ref_pos)
                                meth_calls[2].append(pos_qual)
                    ref_pos += 1
                    query_pos += 1
            elif cigar_type in reference_consumers and cigar_type not in query_consumers:
                ref_pos += cigar_count
            elif cigar_type in query_consumers and cigar_type not in reference_consumers:
                query_pos += cigar_count
        return meth_calls

    def process_meth_vector(self, aligned_read: pysam.AlignedRead, meth_calls: List,
                            strand: str, meth_vectors: Dict[str, Any]) -> Union[None, Tuple]:
        if aligned_read.is_proper_pair:
            mate_flag = self.mate_flag_dict[aligned_read.flag]
            mate_pair_label = f'{aligned_read.query_name}_{mate_flag}_{aligned_read.next_reference_start}'
            if mate_pair_label in meth_vectors:
                paired_calls = meth_vectors.pop(mate_pair_label)['calls']
                paired_calls[0].extend(meth_calls[0])
                paired_calls[1].extend(meth_calls[1])
                paired_calls[2].extend(meth_calls[2])
                paired_calls = self.clean_overlap(paired_calls)
                return (aligned_read.query_name, paired_calls[1][0], paired_calls[1][-1],
                        np.array(paired_calls[0]), np.array(paired_calls[1]),
                        aligned_read.flag, mate_flag, strand)
            else:
                vector_label = f'{aligned_read.query_name}_{aligned_read.flag}_{aligned_read.reference_start}'
                meth_vectors[vector_label] = {'calls': meth_calls,
                                                     'read_name': aligned_read.query_name,
                                                     'flag': aligned_read.flag, 'strand': strand}
                return None
        else:
            return (aligned_read.query_name, meth_calls[1][0], meth_calls[1][-1],
                    np.array(meth_calls[0]), np.array(meth_calls[1]), aligned_read.flag, None, strand)


    def get_meth_call(self, ref_base: str, base_call: str) -> Tuple[bool, Union[None, int]]:
        """
        Methylation for a C relative to the sense strand of the reference can only be called using watson reads,
        and G with crick reads
        Arguments
            nucleotide (str): reference nucleotide
            base_call (collections.Counter): watson nucleotides are Uppercase and crick nucleotides lowercase
        Returns:
            methylation call dictionary
        """
        # call cytosine with watson
        try:
            tmp_call = self.meth_state_dict[(ref_base, base_call)]
        except:
            return False, None
        else:
            return tmp_call


    @staticmethod
    def clean_overlap(meth_calls: List) -> List:
        cleaned_calls = {}
        for meth_call, pos, qual in zip(meth_calls[0], meth_calls[1], meth_calls[2]):
            if pos not in cleaned_calls:
                cleaned_calls[pos] = (meth_call, qual)
            else:
                if qual > cleaned_calls[pos][1]:
                    cleaned_calls[pos] = (meth_call, qual)
        return [[call[0] for call in cleaned_calls.values()], list(cleaned_calls.keys())]


    @staticmethod
    def get_ref_sequence(path: str) -> str:
        """load serialized reference file from path
        """
        with open(path, 'rb') as genome_file:
            return pickle.load(genome_file)
