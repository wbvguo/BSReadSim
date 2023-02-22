#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -j y # Error stream is merged with the standard output
#$ -l h_data=2G,h_rt=4:00:00
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


# set path
working_path=/u/home/w/wbguo/project-pellegrini/test/
DB_path=$working_path/idx/bsbolt

# call methylation
bsbolt CallMethylation -I $working_path/ERR2359938.chr21.10.sorted.mdup.fix.bam -O $working_path/ERR2359938 -DB $DB_path -t 8 -min 10

