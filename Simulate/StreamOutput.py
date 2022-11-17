import os
import pickle
import subprocess


class StreamOutput:
    '''
    write simulation values/variants/reads to disk
    :param str  outdir: path to the simulation folder
    :param bool shuffle: whether to shuffle the reads (Default: reads are segemented by contig_id)
    :rtype None
    '''

    def __init__(self, outdir: str = None, overwrite_db: bool = False):
        self.outdir = outdir
        self.pkl_dir= f'{self.outdir}/pkl/'
        self.overwrite_db = overwrite_db

    def create_outdir(self):
        '''create output directory'''
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir, exist_ok=False)
            os.makedirs(self.pkl_dir, exist_ok=False)

    def check_outdir(self):
        '''check if we have existence and permission'''
        if not os.path.isdir(self.outdir):
            print(f"No such folder: {self.outdir}")
        if not os.path.isdir(self.pkl_dir):
            print(f"No such folder: {self.pkl_dir}")


    def output_contig(self, contig_id, contig_profile, is_variant=False):
        '''output methylation or variants'''
        if is_variant:
            contig_label = f'{contig_id}_variants'
        else:
            contig_label = f'{contig_id}_values'

        output_file = f'{self.pkl_dir}/{contig_label}.pkl'

        if not self.overwrite_db and os.path.exists(output_file):
            raise ValueError("Output file exists but overwrite_db is false, please check")

        with open(output_file, 'wb') as file:
            pickle.dump(contig_profile, file)


    def load_contig(self, contig_id, is_variant=False):
        '''load the contig profiles'''
        if is_variant:
            contig_label = f'{contig_id}_variants'
        else:
            contig_label = f'{contig_id}_values'

        try:
            with open(f'{self.pkl_dir}/{contig_label}.pkl', 'rb') as file:
                contig_profile = pickle.load(file)
        except FileNotFoundError:
            profile_type = 'variant' if is_variant else 'methylation values'
            print(f'{contig_id}: {profile_type} profile not found in {self.pkl_dir}')
            return None
        else:
            return contig_profile


    def output_reads(self, reads):
        '''write reads to disk'''
        pass


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