import os
import pickle
import numpy as np

from typing import Dict, List


class StreamMethDB:
    '''
    write methylation values/variants to disk
    :param str  outdir      : path to the simulation folder
    :param str  pkl_dir     : subfolder of the meth_db pkl file
    :param bool overwrite_db: if we should overwrite exisiting files
    :rtype None
    '''

    def __init__(self, outdir: str = None, meth_db_path: str = None,
                 overwrite_db: bool = False, ref_dict: Dict = None):
        self.outdir = outdir
        self.pkl_dir= f'{self.outdir}/pkl/'
        self.overwrite_db = overwrite_db
        if ref_dict:
            self.save_ref(ref_dict)
        else:
            self.load_ref()
        self.genome_len = sum([len(seq) for _, seq in self.ref_dict.items()])

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
        pass

    def load_ref(self):
        pass

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


    def set_var_meth(self, contig_id, sim_data, update_boundary=True) -> Dict[str, List]:
        '''set random methylation due to variants are random, update the meth_arr on the boundary'''
        if not sim_data: # can have no variant
            return None
        self.update_boundary = update_boundary
        if update_boundary:
            self.pos_map, self.meth_arr, _ = self.load_contig(contig_id)
        var_meth_dict = {}

        seq = self.ref_dict[contig_id].seq.upper()
        seq_len = len(seq)
        for pos, variant_info in sim_data.items():
            if pos<2 or pos>(seq_len-2):
                continue
            if variant_info['indel'] == -1:  # deletion starts at pos
                offset = variant_info['offset']
                local_seq = f'{seq[(pos-2):(pos)]}{seq[(pos+offset):(pos+offset+2)]}'
                pos_list  = [pos-2, pos-1, pos+offset, pos+offset+1]
                self.handle_boundary(pos_list, local_seq)
                continue
            elif variant_info['indel'] == 1:  # insertion starts at pos
                offset = variant_info['offset']
                local_seq = f'{seq[(pos-2):(pos)]}{variant_info["alt"]}{seq[(pos):(pos+2)]}'
                pos_list  = [pos-2, pos-1, pos, pos+1]
                self.handle_boundary(pos_list, local_seq)

                ins_meth_arr = np.zeros(offset)
                ins_ctx_arr  = np.zeros(offset)
                for ins_idx, base in enumerate(variant_info['alt']):
                    if base not in {'C', 'G'}:
                        continue
                    updown  = 1 if base == "C" else -1  # whether to go upstream or downstream
                    base_d1 = local_seq[1*updown+ins_idx+2]
                    base_d2 = local_seq[2*updown+ins_idx+2]
                    context = self.get_cg_context(base, base_d1, base_d2)
                    ins_meth_arr[ins_idx]= self.simu_beta_dist(context=context)[0] # base,3,context
                    ins_ctx_arr[ins_idx] = context

                if np.any(ins_ctx_arr):
                    variant_info['meth'] = ins_meth_arr
                    variant_info['ctx']  = ins_ctx_arr
                    var_meth_dict[pos]   = (ins_meth_arr, ins_ctx_arr)
            else: # substitution
                base = variant_info['alt']
                local_seq = f'{seq[(pos-2):pos]}{variant_info["alt"]}{seq[(pos+1):(pos+3)]}'
                pos_list  = [pos-2, pos-1, pos+1, pos+2]
                self.handle_boundary(pos_list, local_seq)

                if base not in {'C', 'G'}:
                    continue
                updown  = 1 if base == "C" else -1
                base_d1 = local_seq[1*updown+2]
                base_d2 = local_seq[2*updown+2]
                context = self.get_cg_context(base, base_d1, base_d2)
                snp_meth= self.simu_beta_dist(context=context)[0]
                variant_info['meth'] = snp_meth
                variant_info['ctx']  = context
                var_meth_dict[pos]   = (snp_meth, context)
            sim_data[pos] = variant_info
        self.output_contig(contig_id, sim_data, is_variant=True)
        if self.update_boundary:
            self.output_contig(contig_id, [self.pos_map, self.meth_arr, 1], is_variant=False)
        return var_meth_dict


    def handle_boundary(self, pos_list, local_seq):
        '''accomandate the boundary of the mutations'''
        if self.update_boundary:
            ptr_list = [0, 1, -2, -1]
            assert len(pos_list) == 4 and len(local_seq) >= 4
            for idx, ptr in enumerate(ptr_list):
                base = local_seq[ptr]
                if ptr >= 0 and base == "C":
                    base_d1 = local_seq[ptr+1]
                    base_d2 = local_seq[ptr+2]
                elif ptr <0 and base == "G":
                    base_d1 = local_seq[ptr-1]
                    base_d2 = local_seq[ptr-2]
                else:
                    continue
                context     = self.get_cg_context(base, base_d1, base_d2)
                change_pos  = pos_list[idx]
                self.meth_arr[self.pos_map[change_pos], 4] = self.simu_beta_dist(context=context)[0]


    @classmethod
    def simu_beta_dist(self, context = "CG", size = 1):
        '''output the values accordig to the context using beta distribution'''
        if isinstance(context, int):
            context = self.context_dict[context]
        return beta.rvs(a=self.beta_params[context][0],
                        b=self.beta_params[context][1],
                        size=size, random_state=self.seed).astype(np.float16)


    @classmethod
    def get_cg_context(self, base, base_d1, base_d2):
        '''input the base and surrounding, output context'''
        if base == "C":
            flag_d1 = int(base_d1 == "G")
            flag_d2 = int(base_d2 == "G")
        elif base == "G":
            flag_d1 = int(base_d1 == "C")
            flag_d2 = int(base_d2 == "C")
        else:
            return None
        return self.base_context_table[base][flag_d1, flag_d2]
