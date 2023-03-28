## check the output
./wgsim -g ~/iproject/BSReadSim/test/vcf/SNP220624.vcf ~/iproject/BSReadSim/test/ref/BSB_test.fa | more



/home/wbguo/iproject/BSBolt/bsbolt/External/WGSIM/wgsim -1 100 -2 100 -e 0.005 -d 400 -s 25 -r 0.001 -R 0.15 -X 0.15 -S -1 -A 0.05 -I 100 -h 0 -N 392320 /home/wbguo/iproject/BSBolt/tests/TestData/BSB_test.fa 2>/dev/null | more

To separate the output file into variant file and the merged fastq file
sed '/variant start/,/variant end/d' test > fastq 
sed '/variant start/,/variant end/d' test | awk 'NR%2==0' > fastq.1
sed '/variant start/,/variant end/d' test | awk 'NR%2==1' > fastq.2

## extract WGSIM
sed -n '/variant start/, /variant end/{ /variant start/! { /variant end/! p } }' test > variants
awk '/variant start/{ f = 1; next } /variant end/{ f = 0 } f' test > variants


#### randomly shuffle the fq1 and fq2 file
# https://www.biostars.org/p/9764/
# the number of simulated reads should not exceed the number of bits in the reference genome for HG, should be less than 8*3G< 2.4 * 10^10 (average DEPTH should be less than read_len * 8)
time awk '{OFS="\t"; getline seq; getline sep; getline qual; print $0,seq,sep,qual}' reads_1.fq | shuf --random-source ../ref/BSB_test.fa | awk '{OFS="\n"; print $1,$2,$3,$4}' > reads1_shuffle.fq
time awk '{OFS="\t"; getline seq; getline sep; getline qual; print $0,seq,sep,qual}' reads_2.fq | shuf --random-source ../ref/BSB_test.fa | awk '{OFS="\n"; print $1,$2,$3,$4}' > reads2_shuffle.fq

