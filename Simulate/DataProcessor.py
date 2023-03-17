import multiprocessing
import queue
import threading

from ReadProcessor import ReadProcessor
from StreamReads import StreamReads

class DataProcessor:
    def __init__(self, num_processes=4, processor: ReadProcessor = None, fastq_out: StreamReads = None):
        self.num_processes = num_processes
        self.in_queue   = multiprocessing.Queue()
        self.out_queue  = multiprocessing.Queue()
        self.pool = multiprocessing.Pool(self.num_processes)
        self.writer_thread = threading.Thread(target=self.write_fastq, args=())
        self.processor  = processor
        self.fastq_out  = fastq_out


    def start(self):
        self.writer_thread.start()


    def write_fastq(self):
        while True:
            try:
                read_pair = self.out_queue.get(timeout=1)
            except queue.Empty:
                if not any(p.is_alive() for p in self.pool._pool):
                    break
                else:
                    continue
            if read_pair is None:
                break
            self.fastq_out.output_reads(read_pair)


    def feed_data(self, data):
        self.in_queue.put(data)


    def stop(self):
        self.pool.close()
        self.pool.join()
        self.out_queue.put(None)
        self.writer_thread.join()


    def worker(self):
        while True:
            data = self.in_queue.get()
            if data is None:
                break
            result = self.processor.process_read_pair2(data)
            self.out_queue.put(result)


    def spawn_workers(self):
        for i in range(self.num_processes):
            self.pool.apply_async(self.worker)

