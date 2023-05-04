conda create -n bsb python=3.8
conda activate bsb
conda install -c anaconda entrypoints
conda install requests ipykernel nb_conda nb_conda_kernels 
conda install jupyter_nbextensions_configurator
conda install numpy==1.23.5 pandas scipy matplotlib tqdm
conda install -c bioconda pysam pybedtools
conda install -c conda-forge biopython
conda install -c pytorch pytorch torchvision cudatoolkit=10.0
conda install -c conda-forge boto3 filelock tokenizers regex sentencepiece sacremoses 
conda install -c conda-forge scikit-learn tensorflow
conda install -c conda-forge pympler
conda install -c anaconda statsmodels
