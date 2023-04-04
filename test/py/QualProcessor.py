import numpy as np
from collections import defaultdict, Counter
from ReadPair import ReadPair
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import multiprocessing

class QualProcessor:
    
    def __init__(self, bam_file: str = None, read_len: int = 151, n_workers: int = 8,
                 contig: str = None, start:int = None, end:int=None):
        self.bam_file = bam_file
        self.contig = contig
        self.start = start
        self.end = end
        
        self.n_workers  = n_workers
        self.read_len   = read_len
        self.qual_pos   = [[Counter({})]*self.read_len]*2
        self.qual_trans = [Counter({})]*2


    def collect_qual_trans_read(self, read1, read2):
        # record quality along the read, and quality transition matrix
        read1_5to3 = 1 if read1.get_tag('YS') == "W_C2T" else 0 # note, only work for directional library
        read1_qual = read1.get_forward_qualities() if read1_5to3 else np.flip(read1.get_forward_qualities())
        read2_qual = np.flip(read2.get_forward_qualities()) if read1_5to3 else read2.get_forward_qualities()
        
        qual_pos   = [(0, ix, qual) for ix, qual in enumerate(read1_qual)] + [(1, ix, qual) for ix, qual in enumerate(read2_qual)]
        qual_trans = [(0, read1_qual[ix-1], read1_qual[ix]) for ix in range(1, len(read1_qual))] + [(1, read2_qual[ix-1], read2_qual[ix]) for ix in range(1, len(read2_qual))]
        return [qual_pos, qual_trans]

    def update_counter(self, results, counter1, counter2, lock):
        with lock:
            for res in results:
                for res_i in res[0]:
                    if res_i in counter1:
                        counter1[res_i] += 1
                    else:
                        counter1[res_i] = 1
                for res_i in res[1]:
                    if res_i in counter1:
                        counter2[res_i] += 1
                    else:
                        counter2[res_i] = 1
        

    def collect_qual_trans_parallel(self):
        read_gen = ReadPair(self.bam_file, contig=self.contig, start=self.start, end=self.end)
        manager  = multiprocessing.Manager()
        counter1 = manager.dict()  # Shared counter
        counter2 = manager.dict()  # Shared counter
        lock     = manager.Lock()  # Shared lock
        with multiprocessing.Pool(processes=self.n_workers) as pool:
            results = [pool.apply_async(self.collect_qual_trans_read, args=(counter, items)) for items in items_list]
            
            for res in results:
                self.update_counter(res, counter1, counter2, lock)
        return counter1, counter2


    def process_contigs(self):
        """Launches a processing pool to call methylation values across the input file contigs
        """
        # initialize manager
        manager = multiprocessing.Manager()
        # get return dictionary
        self.return_queue = manager.Queue(maxsize=20)
        self.completed_contigs = manager.list()
        # threads for methylation calling, if one thread use thread for calling and watching
        pool_threads = self.threads - 1 if self.threads != 1 else 1
        # start pool
        self.pool = multiprocessing.Pool(processes=pool_threads)
        # for contig call methylation and return values to dict
        for contig in self.contigs:
            contig_kwargs = dict(self.call_methylation_kwargs)
            contig_kwargs.update(dict(contig=contig, return_queue=self.return_queue))
            self.pool.apply_async(call_contig_methylation,
                                  args=[self.completed_contigs, contig_kwargs],
                                  error_callback=self.methylation_process_error)
        self.pool.close()