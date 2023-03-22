import queue
import threading

from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from ReadProcessor import ReadProcessor
from StreamReads import StreamReads
from multiprocessing.pool import ThreadPool

class DataProcessor:
    def __init__(self, read_gen, n_workers=4, 
                 processor: ReadProcessor = None, fastq_out: StreamReads = None):
        self.read_gen   = read_gen
        self.n_workers  = n_workers
        self.processor  = processor
        self.fastq_out  = fastq_out


    def start_processing_mt(self):
        self.out_queue  = queue.Queue()
        self.writer_thread = threading.Thread(target=self.write_fastq)
        self.writer_thread.start()

        self.pool = ThreadPool(self.n_workers)
        for result in self.pool.imap_unordered(self.processor.process_read_pair3, self.read_gen):
            # put the result into the queue for the writer thread to process
            self.out_queue.put(result)
        
        # self.pool = ThreadPoolExecutor(max_workers=self.n_workers)
        # # Generate and process data in parallel using the worker threads
        # for _, read_pair in self.read_gen:
        #     future = self.pool.submit(self.processor.process_read_pair2, read_pair)
        #     future.add_done_callback(lambda f: self.out_queue.put(f.result()))
        
    def write_fastq(self):
        while True:
            read_pair = self.out_queue.get()
            if read_pair is None:
                break
            self.fastq_out.output_reads(read_pair)


    def stop(self):
        # Wait for all worker threads to finish before terminating the write thread
        #self.pool.shutdown(wait=True)
        self.pool.close()
        self.pool.join()
        self.out_queue.put(None)
        self.writer_thread.join()

