import pysam
import threading
from collections import defaultdict


class ReadPair:
    """
    return the read pairs 
    """
    
    def __init__(self, bam_file: str=None, contig: str= None, start:int = None, end:int =None):
        self.bam_file= bam_file
        self.contig = contig
        self.start = start
        self.end  = end
        self.it = iter(self.pair_generator())
        self.lock = threading.Lock()

    def __iter__(self):
        return self

    def next(self):
        with self.lock:
            return next(self.it)
    
    def pair_generator(self):
        """
        Generate read pairs in a BAM file or within a region.
        Reads are added to read_dict until a pair is found.
        """
        read_dict = defaultdict(lambda: [None, None])
        bam = pysam.AlignmentFile(self.bam_file, "rb")
        #bam.reset()
        
        for align in bam.fetch(until_eof=True, contig=self.contig, start=self.start, stop=self.end):
            if not align.is_proper_pair or align.is_secondary or align.is_supplementary:
                continue
            qname = align.query_name
            if qname not in read_dict.keys():
                read_dict[qname][int(align.is_read2)] = align
            else:
                read_dict[qname][int(align.is_read2)] = align
                yield read_dict[qname]
                del read_dict[qname]

