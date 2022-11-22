import os
import pickle


class StreamMethDB:
    '''
    write methylation values/variants to disk
    :param str  outdir      : path to the simulation folder
    :param str  pkl_dir     : subfolder of the meth_db pkl file
    :param bool overwrite_db: if we should overwrite exisiting files
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


    def save_ref(self, ref_dict):
        if not self.overwrite_db:
            raise ValueError("cannot save ref_dict, overwrite is set to be false\n")
        with open(self.pkl_dir + '/ref_dict.pkl', 'wb') as file:
            pickle.dump(ref_dict, file)


    def load_ref(self):
        with open(self.pkl_dir + '/ref_dict.pkl', 'rb') as file:
            return pickle.load(file)


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


