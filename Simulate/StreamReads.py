import subprocess


class StreamReads:
    '''
    :param bool shuffle: whether to shuffle the reads (Default: reads are segemented by contig_id)
    '''
    
    def __init__(self, shuffle:bool = True, pair_end:bool = True, fasta:str = None,
                 outdir: str = None, prefix: str = "sim") -> None:
        self.shuffle    = shuffle
        self.fasta      = fasta
        self.pair_end   = pair_end
        self.outdir     = outdir
        self.prefix     = prefix
        self.fastq_list = None
        self.fastq_obj  = None
    
    @property
    def creata_fastq(self):
        """Return io object for fastq writing"""
        fastq1 = f'{self.outdir}/{self.prefix}_1.fastq'
        fastq2 = f'{self.outdir}/{self.prefix}_2.fastq'
        if self.pair_end:
            self.fastq_list = [fastq1, fastq2]
        else:
            self.fastq_list = [fastq1]
        
        for fastq in self.fastq_list:
            # check if the file exists
            self.fastq_obj.append(open(fastq, 'a'))
        

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