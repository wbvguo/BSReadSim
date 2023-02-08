##### find all the eligible regions
cd PGP-UK
samtools view -bq 30 -h -F 0x4 -F 0x100 -F 0x200 -F 0x400 -F 0x800 bam/ERR2359938.sorted.bam chr21 | bedtools genomecov -ibam - -bg > ERR2359938.chr21.bed
awk '{ if ($4 > 20) { print } }' ERR2359938.chr21.bed |  bedtools merge -i -  | awk '{ if ($3 - $2 > 500) { print } }' > ERR2359938.chr21.d20.l500.merged.bed

# filter based on the mapping quality 20
# exclude ummaped
# exclude not primary alignment
# exclude failed qualify check
# exclude duplicates
# exclude supplementary

