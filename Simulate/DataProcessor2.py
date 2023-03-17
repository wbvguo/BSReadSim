import multiprocessing
import queue
import threading
import numpy as np
import pickle

from ReadProcessor import ReadProcessor
from StreamReads import StreamReads
from StreamMethDB import StreamMethDB

class DataProcessor:
    def __init__(self, contig_id: str = None, read_gen, n_process=4, 
                 processor: ReadProcessor = None, arr_max_size: List = None,
                 fastq_out: StreamReads = None, meth_db: StreamMethDB = None):
        self.read_gen   = read_gen
        self.n_process  = n_process
        self.out_queue  = multiprocessing.Queue()
        self.pool = multiprocessing.Pool(self.n_process)
        self.writer_thread = threading.Thread(target=self.write_fastq, args=())
        self.processor  = processor
        self.fastq_out  = fastq_out
        self.meth_db    = meth_db
        self.shared_pos_map = None
        self.shared_meth_arr= None
        self.contig_id  = contig_id
        self.arr_max_size= arr_max_size

    def start_processing_mp(self, ):
        # start the thread for writing to the file
        self.create_share_arr(self.arr_max_size)
        self.load_contig_share(self.contig_id)
        
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


    def create_share_arr(self, max_arr_size):
        shared_data1 = multiprocessing.RawArray('I', max_arr_size[0])
        shared_data2 = multiprocessing.RawArray('f', max_arr_size[1]*5)
        self.shared_pos_map = np.frombuffer(shared_data1, dtype=np.uint32)
        self.shared_meth_arr= np.frombuffer(shared_data2, dtype=np.float32).reshape((max_arr_size[1], 5))


    def load_contig_share(self, contig_id):
        '''load the contig profiles'''
        input_file = f'{self.meth_db.tmp_dir}/{contig_id}_values.pkl'
        
        try:
            with open(input_file, 'rb') as FILE:
                contig_profile = pickle.load(FILE)
        except FileNotFoundError:
            print(f'{contig_id}: methylation profile not found in {self.meth_db.tmp_dir}')
            return None
        else:
            # copy the data from the original numpy array to the shared numpy array
            self.shared_pos_map[:]   = np.zeros_like(self.shared_pos_map)
            self.shared_meth_arr[:,:]= np.zeros_like(self.shared_meth_arr)
            self.shared_pos_map[:contig_profile[0].shape[0]]   = contig_profile[0][:]
            self.shared_meth_arr[:contig_profile[1].shape[0],:]= contig_profile[1][:,:]

