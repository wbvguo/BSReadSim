import subprocess
import threading
import os
import numpy as np


class StreamReads:
    '''stream reads and write into fastq files'''
    def __init__(self, outdir: str = None, prefix: str = None, pair_end: bool = True,
                 shuffle: bool = True, ref_fasta: str = None, gzip: bool = True):
        self.outdir   = outdir
        self.prefix   = prefix
        self.pair_end = pair_end
        self.shuffle  = shuffle
        self.ref_fasta= ref_fasta
        self.gzip     = gzip
        self.writeLock= threading.Lock()

        fastq1 = f'{self.outdir}/{self.prefix}_1.fastq'
        fastq2 = f'{self.outdir}/{self.prefix}_2.fastq'
        if self.pair_end:
            self.fastq_list = [fastq1, fastq2]
        else:
            self.fastq_list = [fastq1]
        self.create_fastq()
        
        # column context, row cigar
        self.cigar_table = np.array([['M', 'c', 'C', 'b',  'B', '-', '-', 'a',  'A', 'c', 'C', 'b',  'B', '-', '-', 'a',  'A'], # match
                                     ['x', 'x', 'X', 'x',  'X', '-', '-', 'x',  'X', 'x', 'X', 'x',  'X', '-', '-', 'x',  'X'], # snp
                                     ['-', '-', '-', '-',  '-', '-', '-', '-',  '-', '-', '-', '-',  '-', '-', '-', '-',  '-'], # empty
                                     ['M', 'i', 'I', 'i',  'I', '-', '-', 'i',  'I', 'i', 'I', 'i',  'I', '-', '-', 'i',  'I'], # insert
                                     ['-', '#', '-', '#',  '-', '-', '-', '#',  '-', '#', '-', '#',  '-', '-', '-', '#',  '-']])# convert failed
        self.qual_str    = '''!"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'''


    def create_fastq(self):
        '''Return io object for fastq writing'''
        self.output_obj = []
        for fastq_file in self.fastq_list:
            # check if the file exists
            self.output_obj.append(open(fastq_file, 'a'))


    def output_reads(self, read_pair, read_flip, read_sub): 
        # take into the read_pairs, read_flip:[0,1]; read_sub:[0,1]
        read_lines = []
        for idx in range(1+int(self.pair_end)):
            id  = f'{read_pair[idx]["read_id"]}/{idx+1}'
            seq = ''.join("ACGT"[i] for i in read_pair[idx]['seq'])
            strd= read_flip[idx]
            sub = read_sub[idx]
            if strd:
                pass
            
            cmt = '+'
            for ix, ctx in enumerate(read_pair[idx]["ctx"]):
                cmt += self.cigar_table[read_pair[idx]["cgr"][ix], ctx]
            cmt += f':{["W","C"][strd]}{["C2T", "G2A"][sub]}'
            qual = '\n'
            #qual= ''.join([self.qual_str[i] for i in read_pair[idx]['qual']])

            read_lines.append(f'{id}\n{seq}\n{cmt}\n{qual}\n')
        self.write_file(read_lines)


    def write_file(self, read_lines):
        """write bisulfite reads to disk"""
        with self.writeLock:
            for idx, line in enumerate(read_lines):
                self.output_obj[idx].write(line)


    def close(self):
        '''close the fastq object and shuffle or gzip reads'''
        for obj in self.output_obj:
            obj.close()

        if self.gzip or self.shuffle:
            print("Shuffling/Compressing reads\n")
        for fastq in self.fastq_list:
            if self.shuffle and self.ref_fasta:
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

