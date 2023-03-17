import numpy as np
from collections import defaultdict, Counter
from ReadPair import ReadPair

class QualProcessor:
    
    def __init__(self, bam_file: str = None, contig: str = None, start:int = None, end:int=None):
        self.bam_file = bam_file
        self.contig = contig
        self.start = start
        self.end = end
        
        self.qual_pos_1  = [Counter({})]*151
        self.qual_pos_2  = [Counter({})]*151

        self.qual_trans_1= Counter({})
        self.qual_trans_2= Counter({})
    
    
    def collect_qual_trans(self):
        # assumed used directional library
        for read1, read2 in ReadPair(self.bam_file, contig=self.contig, start=self.start, stop=self.end):
            # record quality along the read, and quality transition matrix
            read1_5to3 = 1 if read1.get_tag('YS') == "W_C2T" else 0
            read1_qual = read1.get_forward_qualities() if read1_5to3 else np.flip(read1.get_forward_qualities())
            read2_qual = np.flip(read2.get_forward_qualities()) if read1_5to3 else read2.get_forward_qualities()

            for ix, base_q in enumerate(read1_qual):
                self.qual_pos_1[ix].update([base_q])
                if ix:
                    self.qual_trans_1.update([(read1_qual[ix-1], base_q)])

            for ix, base_q in enumerate(read2_qual):
                self.qual_pos_2[ix].update([base_q])
                if ix:
                    self.qual_trans_2.update([(read2_qual[ix-1], base_q)])
                    
    def counter_update:
        # dedicated file writing task


def file_writer(filepath, queue):
    # run forever
    while True:
        # get a line of text from the queue
        line = queue.get()
        # write it to file
        
        
        # flush the buffer
        file.flush()
        # mark the unit of work complete
        queue.task_done()
        



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