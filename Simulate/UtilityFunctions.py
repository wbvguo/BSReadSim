from distutils.version import LooseVersion
import os
import io
import sys
import gzip
import pysam
from typing import List, Union


class parseCGmap:
    '''turn a CGmap file into a generator'''
    def __init__(self, cgmap_file: str = None, contig_id: str = None, collect_ch: bool = False):
        self.cgmap_file = cgmap_file
        self.contig_id  = contig_id
        self.collect_ch = collect_ch

        if cgmap_file.endswith('.gz'):
            self.file_obj = io.BufferedReader(gzip.open(cgmap_file, 'rb'))
        else:
            self.file_obj = open(cgmap_file, 'r')


    def __iter__(self):
        with self.file_obj as cg:
            while True:
                line = cg.readline()
                if not line:
                    break
                chr_id, base, pos, context, diN, meth_level, meth_count, tot_count = self.process_line(line)
                if self.contig_id:
                    if chr_id != self.contig_id:
                        continue
                    if not self.collect_ch and context != 'CG':
                        continue
                yield [chr_id, base, pos, context, meth_level]

    @staticmethod
    def process_line(line: Union[str, bytes]) -> List[str]:
        '''parse the CGmap lines'''
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '').split('\t')
        else:
            return line.replace('\n', '').split('\t')


class parseASM:
    '''parse the ASM file'''
    def __init__(self, asm_file: str = None, contig_id: str = None, collect_ch: bool = False):
        self.asm_file   = asm_file
        self.contig_id  = contig_id
        self.collect_ch = collect_ch

        if asm_file.endswith('.gz'):
            self.file_obj = io.BufferedReader(gzip.open(asm_file, 'rb'))
        else:
            self.file_obj = open(asm_file, 'r')


    def __iter__(self):
        with self.file_obj as asm:
            while True:
                line = asm.readline()
                if not line:
                    break
                # chr, base, pos, context, dinucleotide, meth, snp_pos, REF, ALT, ref_meth, alt_meth, fold_change, pval, comment
                # comment can be meth_count/tot_count;ref_meth_count/ref_tot_count;alt_meth_count/alt_tot_count
                chr_id, base, pos, context, diN, tot_meth, snp_pos, ref, alt, ref_meth, alt_meth, fold_change, p_val, comment = self.process_line(line)
                if chr_id != self.contig_id:
                    continue
                if not self.collect_ch and context != 'CG':
                    continue
                yield [chr_id, base, pos, context, tot_meth, ref_meth, alt_meth]


    @staticmethod
    def process_line(line: Union[str, bytes]) -> List[str]:
        '''parse the ASM lines'''
        if isinstance(line, bytes):
            return line.decode('utf-8').replace('\n', '').split('\t')
        else:
            return line.replace('\n', '').split('\t')


def get_htsim_path():
    """Get paths of dependencies. Print warning if setup.py not run and dependencies not compiled.
    rtype: str htsim_path: path to htsim executable htsim
    """
    utility_dir = os.path.dirname(os.path.realpath(__file__))
    base_dir    = os.path.dirname(utility_dir)
    htsim_path  = f'{base_dir}/HTSIM/htsim'
    if not os.path.exists(htsim_path):
        raise ValueError("[ERROR] Executable htsim not found, please check!")
    return htsim_path


def complement(sequence):
    """
    Params:

    * *sequence(str) *: DNA sequence, non ATGC nucleotide will be returned unaltered

    Returns:

    * *sequence.translate(_rc_trans)(str) *: complement of input sequence
    """
    _rc_trans = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
    return sequence.translate(_rc_trans)


def reverse_complement(sequence):
    """
    Params:

    * *sequence (str)*: DNA sequence, non ATGC nucleotide will be returned unaltered

    Returns:

    * *reversed_string.translate(_rc_trans) (str)*: reverse complement of input sequence
    """
    # reverse string
    reversed_string = sequence[::-1]
    return complement(reversed_string)


def retrieve_iupac(nucleotide):
    """
    Params:

    * *nucleotide (str)*: single character

    Returns:

    * *iupac_tuple (tuple)*: tuple of strings with possible bases
    """
    iupac_key = {'R': ('A', 'G'), 'Y': ('C', 'T'), 'S': ('G', 'C'), 'W': ('A', 'T'), 'K': ('G', 'T'), 'M': ('A', 'C'),
                 'B': ('C', 'G', 'T'), 'D': ('A', 'G', 'T'), 'H': ('A', 'C', 'T'), 'V': ('A', 'C', 'G'),
                 'N': ('A', 'C', 'G', 'T')}
    try:
        iupac_tuple = iupac_key[nucleotide]
    except KeyError:
        iupac_tuple = tuple(nucleotide)
    return iupac_tuple


def import_package_check(package_name):
    try:
        package = __import__(package_name)
    except ModuleNotFoundError:
        print(f'{package_name} not found, please install {package_name}')
        print(f'pip3 install {package_name} --user')
        sys.exit()
    else:
        return package.__version__


def check_package_version():
    """Check third party packages to make sure version is greater than requirement"""
    packages = dict(pysam='0.9.1', tqdm='4.31.1', numpy='1.16.0')
    all_packages = True
    for package, required_version in packages.items():
        installed_version = import_package_check(package)
        if LooseVersion(installed_version) < LooseVersion(required_version):
            print(f'{package} {installed_version} < {package} {required_version}, please update {package}')
            print(f'pip3 install {package} --upgrade')
            all_packages = False
    return all_packages


def check_python_version():
    if sys.version_info < (3, 6):
        print('Python must be >= 3.6.0')
        raise OSError





def propagate_error(error):
    raise error


def sort_bam(bam_output=None, bam_input=None):
    """ Sort bam file

    Params:
    * *bam_output (str)*: output path for sorted bam file
    * *bam_input (str)*: input bam file"""
    pysam.sort('-o', bam_output, bam_input)


def index_bam(bam_input=None):
    """
    Index bam file

    Params:
    * *bam_input (str)*: input bam files
    """
    pysam.index(bam_input)
