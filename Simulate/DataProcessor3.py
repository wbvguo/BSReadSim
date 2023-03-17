import multiprocessing
import queue
import threading

from ReadProcessor import ReadProcessor
from StreamReads import StreamReads

class DataProcessor:
    def __init__(self, read_gen, n_process=4, processor: ReadProcessor = None, fastq_out: StreamReads = None):
        self.read_gen   = read_gen
        self.n_process  = n_process
        self.out_queue  = multiprocessing.Queue()
        self.pool = multiprocessing.Pool(self.n_process)
        self.writer_thread = threading.Thread(target=self.write_fastq, args=())
        self.processor  = processor
        self.fastq_out  = fastq_out

    def start_processing_mp(self):
        # start the thread for writing to the file
        writer_thread = threading.Thread(target=self.write_fastq)
        writer_thread.start()

        # use imap_unordered to apply the processing function to each piece of data
        for result in self.pool.imap_unordered(self.processor.process_read_pair3, self.read_gen):
            # put the result into the queue for the writer thread to process
            self.out_queue.put(result)

        # add None to the queue to signal the writer thread to exit
        # self.out_queue.put(None)
        # writer_thread.join()


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


    def stop(self):
        self.pool.close()
        self.pool.join()
        self.out_queue.put(None)
        self.writer_thread.join()
