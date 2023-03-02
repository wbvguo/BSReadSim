import subprocess
from threading import Lock
import os
import numpy as np
import warnings

class StreamReads:
    '''stream reads and write into fastq files'''
    def __init__(self, outdir: str = None, prefix: str = None, pair_end: bool = True,
                 write_lock: Lock = None, shuffle: bool = False, ref_fasta: str = None, gzip: bool = False):
        self.outdir   = outdir
        self.prefix   = prefix
        self.pair_end = pair_end
        self.shuffle  = shuffle
        self.ref_fasta= ref_fasta
        self.gzip     = gzip
        self.writeLock= write_lock
        self.create_fastq()
        
        # column context, row cigar
        self.cigar_table = np.array([['M', 'c', 'C', 'b',  'B', '-', '-', 'a',  'A', 'c', 'C', 'b',  'B', '-', '-', 'a',  'A'], # match
                                     ['x', 'x', 'X', 'x',  'X', '-', '-', 'x',  'X', 'x', 'X', 'x',  'X', '-', '-', 'x',  'X'], # snp
                                     ['-', '-', '-', '-',  '-', 'e', 'E', '-',  '-', '-', '-', '-',  '-', 'e', 'E', '-',  '-'], # seq err
                                     ['i', 'i', 'I', 'i',  'I', '-', '-', 'i',  'I', 'i', 'I', 'i',  'I', '-', '-', 'i',  'I'], # insert
                                     ['-', '#', '-', '#',  '-', '-', '-', '#',  '-', '#', '-', '#',  '-', '-', '-', '#',  '-']])# convert failed


    def create_fastq(self):
        '''Return io object for fastq writing'''
        fastq1 = f'{self.outdir}/{self.prefix}_1.fastq'
        fastq2 = f'{self.outdir}/{self.prefix}_2.fastq'
        self.fastq_list = [fastq1, fastq2] if self.pair_end else [fastq1]
        
        # check the existence of fastq and create object
        self.output = []
        for fastq_file in self.fastq_list:
            if os.path.exists(fastq_file):
                warnings.warn(f'Fastq file exists, will overwrite... {fastq_file}')
                os.remove(fastq_file)           
            self.output.append(open(fastq_file, 'a'))


    def output_reads(self, read_pair, read1_idx, pattern_idx):
        # take read_pairs, read1_idx:{0,1}; pattern_idx:{0,1}
        read_lines = []

        read_lines.append(self.gen_fastq_lines(read_pair[read1_idx], 1, read1_idx, pattern_idx))
        if self.pair_end:
            read_lines.append(self.gen_fastq_lines(read_pair[1-read1_idx], 2, read1_idx, 1- pattern_idx))
        
        self.write_file(read_lines)


    def gen_fastq_lines(self, read_rec, pair_idx, read1_idx, pattern_idx):
        id  = f'{read_rec["read_id"]}/{pair_idx}'
        seq = ''.join("ACGTN"[i] for i in read_rec['seq'])
        
        strand  = ["W","C"][read1_idx]
        pattern = ["C2T", "G2A"][pattern_idx]
        
        cmt = "+"
        for ix, ctx in enumerate(read_rec["ctx"].filled()):
            cmt += self.cigar_table[read_rec["cgr"][ix], ctx]
        
        cmt += f':{strand}_{pattern}'
        qual= ''.join([chr(i) for i in read_rec['qual']])
        return f'{id}\n{seq}\n{cmt}\n{qual}\n'


    def write_file(self, read_lines):
        """write bisulfite reads to disk"""
        with self.writeLock:
            for idx, line in enumerate(read_lines):
                self.output[idx].write(line)


    def close(self):
        '''close the fastq object and shuffle or gzip reads'''
        for obj in self.output:
            obj.close()

        if self.gzip or self.shuffle:
            print("Shuffling/Compressing reads\n")
        for fastq in self.fastq_list:
            if self.shuffle:
                self.shuffle_read(fastq, self.ref_fasta)
            if self.gzip:
                self.gzip_read(fastq)


    @staticmethod
    def shuffle_read(fastq_file: str = None, random_source: str = None):
        '''
        shuffle the reads randomly, if false, the reads will be segemented by contigs
        the number of simulated reads should not exceed the number of bits in the reference genome.
        for HG, should be less than 8*3G < 2.4*10^10 (average DEPTH should be less than read_len*8)
        from link: https://www.biostars.org/p/9764/
        '''
        print(f"shuffle reads {fastq_file}")
        
        prefix_split    = os.path.splitext(fastq_file)[0].split("_")
        fastq_shuffle   = "_".join(prefix_split[:-1]) + "_shuffle_" + prefix_split[-1] + ".fastq"
        shuf_cmd_list   = ["awk", "'{OFS=\"\t\"; getline seq; getline sep; getline qual; print $0,seq,sep,qual}'",
                           fastq_file, "|", "shuf --random-source", random_source, "|",
                           "awk", "'{OFS=\"\n\"; print $1,$2,$3,$4}'", ">", fastq_shuffle]
        shuffle_run     = subprocess.Popen(shuf_cmd_list,  stdout=subprocess.PIPE, universal_newlines=True)
        stdout, stderr  = shuffle_run.communicate()
        
        rename_cmd_list = ["mv", fastq_shuffle, fastq_file]
        rename_run      = subprocess.Popen(rename_cmd_list,stdout=subprocess.PIPE, universal_newlines=True)
        stdout, stderr  = rename_run.communicate()


    @staticmethod
    def gzip_read(fastq_file):
        '''gzip the output reads to fastq.gz format'''
        print(f"gzip reads {fastq_file}")
        
        gzip_cmd_list   = ["gzip", fastq_file]
        gzip_run        = subprocess.Popen(gzip_cmd_list,  stdout=subprocess.PIPE, universal_newlines=True)
        stdout, stderr  = gzip_run.communicate()

