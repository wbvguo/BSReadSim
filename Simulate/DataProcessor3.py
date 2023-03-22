import multiprocessing
import queue
import threading

from multiprocessing.pool import ThreadPool
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from ReadProcessor import ReadProcessor
from StreamReads import StreamReads

class DataProcessor:
    def __init__(self, read_gen, n_workers=4, 
                 processor: ReadProcessor = None, fastq_out: StreamReads = None):
        self.read_gen   = read_gen
        self.n_workers  = n_workers
        self.processor  = processor
        self.fastq_out  = fastq_out

    def start_processing_2q(self):
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            for _, read_pair in self.read_gen:
                future = executor.submit(self.worker, read_pair)
                #future.add_done_callback(lambda f: self.fastq_out.output_reads_lock(f.result()))
            executor.shutdown(wait=True)


    def worker(self, read_pair):
        result = self.processor.process_read_pair2(read_pair)
        self.fastq_out.output_reads_lock(result)

