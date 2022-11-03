import os
import pickle


class StreamOutput:
    '''
    write simulation values/variants/reads to disk
    :param str  outdir: path to the simulation folder
    :param bool shuffle: whether to shuffle the reads (Default: reads are segemented by contig_id)
    :rtype None
    '''

    def __init__(self, outdir: str = None):
        self.outdir = outdir
        self.pkl_dir= f'{self.outdir}/pkl/'

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

        with open(f'{self.pkl_dir}/{contig_label}.pkl', 'wb') as file:
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
