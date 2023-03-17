import os
import pickle
from typing import Dict

class StreamMethDB:
    '''
    write/load methylation values or variants to/from disk
    :param str  outdir      : path to the simulation folder
    :param bool overwrite_db: if we should overwrite exisiting files
    :rtype None
    '''

    def __init__(self, outdir: str = None, overwrite_db: bool = False):
        self.outdir = outdir
        self.pkl_dir= f'{self.outdir}/pkl/'
        self.tmp_dir= f'{self.outdir}/tmp/'
        self.ref_pkl= f'{self.pkl_dir}/ref_dict.pkl'
        self.overwrite_db = overwrite_db
        
        self.create_outdir()


    def create_outdir(self):
        '''create output directory'''
        if not os.path.isdir(self.outdir):
            os.makedirs(self.outdir, exist_ok=False)
        if not os.path.isdir(self.pkl_dir):
            os.makedirs(self.pkl_dir, exist_ok=False)
        if not os.path.isdir(self.tmp_dir):
            os.makedirs(self.tmp_dir, exist_ok=False)


    def check_outdir(self):
        '''check if we have existence and permission'''
        if not os.path.isdir(self.outdir):
            print(f"No such folder: {self.outdir}")
        if not os.path.isdir(self.pkl_dir):
            print(f"No such folder: {self.pkl_dir}")
        if not os.path.isdir(self.tmp_dir):
            print(f"No such folder: {self.tmp_dir}")


    def save_ref(self, ref_dict: Dict = None, overwrite_db: bool = False):
        if not self.overwrite_db and not overwrite_db and os.path.exists(self.ref_pkl):
            raise ValueError("cannot save reference dict, overwrite_db is set to be false\n")
        with open(self.ref_pkl, 'wb') as FILE:
            pickle.dump(ref_dict, FILE)


    def load_ref(self):
        with open(self.ref_pkl, 'rb') as FILE:
            return pickle.load(FILE)


    def output_contig(self, contig_id, contig_profile, is_variant=False, overwrite_db: bool = False):
        '''output methylation or variants'''
        type_label = 'variants' if is_variant else 'values'
        output_file= f'{self.pkl_dir}/{contig_id}_{type_label}.pkl'
        
        if not self.overwrite_db and not overwrite_db and  os.path.exists(output_file):
            raise ValueError("Output file exists but overwrite_db is false, please check")
        with open(output_file, 'wb') as FILE:
            pickle.dump(contig_profile, FILE)


    def load_contig(self, contig_id, is_variant=False):
        '''
        load the contig_profile: pos_map, meth_arr, flag (0 for not updated, 1 for updated by variants)
        '''
        type_label = 'variants' if is_variant else 'values'
        input_file = f'{self.pkl_dir}/{contig_id}_{type_label}.pkl'
        
        try:
            with open(input_file, 'rb') as FILE:
                contig_profile = pickle.load(FILE)
        except FileNotFoundError:
            profile_type = 'variant' if is_variant else 'methylation'
            print(f'{contig_id}: {profile_type} profile not found in {self.pkl_dir}')
            return None
        else:
            return contig_profile

