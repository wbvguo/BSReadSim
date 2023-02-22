#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=2G,h_rt=1:00:00
#$ -pe shared 8
#$ -r n # job is NOT rerunable
#$ -m a # Email on abort
#$ -o joblog.$JOB_ID

source ~/.bash_profile
source /u/local/Modules/default/init/modules.sh
source /u/home/w/wbguo/.bash_profile

module load gcc/10.2.0
module load anaconda3
conda activate bsb


working_path=/u/home/w/wbguo/project-pellegrini/test/

bsbolt Index -G $working_path/chr21.fa -DB $working_path/idx/bsbolt
