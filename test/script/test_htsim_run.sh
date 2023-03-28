# cd in HTSIM folder

# input fasta alone
./htsim ~/iproject/BSReadSim/test/data/ref/BSB_test.fa | more

./htsim ~/iproject/BSReadSim/test/data/ref/BSB_test.fa -N 100 -R 0 -E 0 1> data 2>log
#################### number
####-d 
####sort by coordinate first
# -N segmentation fault when N 100,000, it also not work as expected
# -N -c
# -n -c 
# -n 
# -c

#test fa.gz 

# input VCF file

#################### SNP

#################### methylation





