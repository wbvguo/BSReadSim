import subprocess
import threading
import os


class StreamReads:
    '''stream reads and write into fastq files'''
    def __init__(self, outdir: str = None, prefix: str = None, 
                 ref_fasta: str = None, pair_end: bool = True,
                 shuffle: bool = True, gzip: bool = True):
        self.outdir   = outdir
        self.prefix   = prefix
        self.ref_fasta= ref_fasta
        self.pair_end = pair_end
        self.shuffle  = shuffle
        self.gzip     = gzip
        self.output_fastq = []
        self.writeLock= threading.Lock()

        fastq1 = f'{self.outdir}/{self.prefix}_1.fastq'
        fastq2 = f'{self.outdir}/{self.prefix}_2.fastq'
        if self.pair_end:
            self.fastq_list = [fastq1, fastq2]
        else:
            self.fastq_list = [fastq1]


    def creata_fastq(self):
        '''Return io object for fastq writing'''
        for fastq_file in self.fastq_list:
            # check if the file exists
            self.output_fastq.append(open(fastq_file, 'a'))


    def output_reads(self, read_pair, read_flip, read1_sub): # take into the read_pairs
        read_lines = []
        for idx in range(1+int(self.pair_end)):
            read_id  = f'{read_pair[idx]["read_id"]}/{idx+1}'
            read_seq = ''.join("ACGT"[i] for i in read_pair[idx]['seq'])
            read_sub = ('C2T', 'G2A')[read_pair[idx]['sub']]
            read_strd= ('W','C')[read_pair[idx]['flip']]
            read_cmt = f'+{read_pair[idx]["cgr_str"]}:{read_strd}_{read_sub}'
            read_qual= ''.join(read_pair[idx]['qual'])

            read_lines.append(f'{read_id}\n{read_seq}\n{read_cmt}\n{read_qual}\n')
        self.write_file(read_lines)


    def write_file(self, read_lines):
        """write bisulfite reads to disk"""
        with self.writeLock:
            for idx, line in enumerate(read_lines):
                self.output_fastq[idx].write(read_lines[idx])


    def close(self):
        '''close the fastq object and shuffle or gzip reads'''
        for output in self.output_fastq:
            output.close()

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




# column context, row cigar
cigar_table = np.array([['M', 'c', 'C', 'b',  'B', '-', '-', 'a',  'A', 'c', 'C', 'b',  'B', '-', '-', 'a',  'A'], # match
                        ['x', 'x', 'X', 'x',  'X', '-', '-', 'x',  'X', 'x', 'X', 'x',  'X', '-', '-', 'x',  'X'], # snp
                        ['-', '-', '-', '-',  '-', '-', '-', '-',  '-', '-', '-', '-',  '-', '-', '-', '-',  '-'], # empty
                        ['M', 'i', 'I', 'i',  'I', '-', '-', 'i',  'I', 'i', 'I', 'i',  'I', '-', '-', 'i',  'I'], # insert
                        ['-', '#', '-', '#',  '-', '-', '-', '#',  '-', '#', '-', '#',  '-', '-', '-', '#',  '-']])# convert failed
