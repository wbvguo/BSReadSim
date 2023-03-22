import multiprocessing
import queue
import threading

from multiprocessing.pool import ThreadPool
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from ReadProcessor import ReadProcessor
from StreamReads import StreamReads

class DataProcessor:
    def __init__(self, read_gen, n_workers=4, 
                 processor: ReadProcessor = None, fastq_out: StreamReads = None,
                 max_qsize: int = 1000):
        self.read_gen   = read_gen
        self.n_workers  = n_workers
        self.processor  = processor
        self.fastq_out  = fastq_out
        self.max_qsize  = max_qsize
        self.end_signal = False

    def start_processing_2q(self):
        self.out_queue  = queue.Queue()
        self.writer_thread = threading.Thread(target=self.write_fastq)
        self.writer_thread.start()
        
        #self.pool = ThreadPool(self.n_workers)
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            for _, read_pair in self.read_gen:
                future = executor.submit(self.worker, read_pair)
                #future.add_done_callback(lambda f: self.output_queue.put(f.result()))

            executor.shutdown(wait=True)
        self.end_signal = True


    def write_fastq(self):
        while True:
            if self.out_queue.qsize() > self.max_qsize or self.end_signal:
                chunk_size = self.out_queue.qsize() if self.end_signal else self.max_qsize
                self.fastq_out.output_reads_chunk([self.out_queue.get() for _ in range(chunk_size)])
                if self.end_signal:
                    break


    def stop(self):
        # self.pool.close()
        # self.pool.join()
        self.out_queue.put(None)
        self.writer_thread.join()


    def worker(self, read_pair):
        result = self.processor.process_read_pair2(read_pair)
        self.out_queue.put(result)

    # def spawn_workers(self):
    #     for i in range(self.n_workers):
    #         self.pool.apply_async(self.worker)

