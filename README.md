# 1. Introduction
**BSReadSim** is a efficient and flexible <ins>**B**</ins>isulfite <ins>**S**</ins>equencing <ins>**Read**</ins> <ins>**Sim**</ins>ulator that allows for profile-based bisulfite sequencing reads simulation. To our best knowledge, BSReadSim is the only simulator that allows user to incooperate both genetic variant and methylation profile with high fidelity. As a result, it can generate more realistic sequencing reads compared to other bisulfite read simulators. The following table summarized the difference between BSReadSim and other tools:

Features | Sherman | BSBolt | BSSim |  BSReadSim
----| ---- | ---- | ---- | ---- | 
sequencing technology| WGBS | WGBS | WGBS |  WGBS/RRBS/TBS
adjustable bisulfite conversion rate | yes | no | yes | yes
support multi-thread | no | no | yes | yes
allow genetic variant input | no   | no | $\textsf{yes}^*$ | yes
allow methylation profile input | no |   $\textsf{yes}^*$ | no | yes
allow alleclic-specific methylation | no | no | no | yes
allow GC bias/nonuniform coverage| no | no | no | yes
haplotype-aware | no   | no | no | yes
site-site dependency | no | no | no| [todo]

$\textsf{yes}^*$: Limited support
- BSSim only accept snps input using a frequency table, it cannot faithfully simulate given genotypes, preserve haplotype information, nor handle indel variants
- BSBolt allows users to input the methylation reference. But during simulation, it **randomly** picks a value from this methylation reference pool. As a result, for a paticular CG site, the simulated data and the reference profile will likely have different methylation level


# 2. Installation
## Dependency
- Python 3.8 or later
- Biopython
- numpy, scipy
- tqdm
- GCC >= 9.4.0

## Install from conda
```
Toadd
```
## Install from github
```
git clone https://github.com/wbvguo/BSReadSim.git
cd BSReadSim
pip install -e .
```


# 3. Usage
## Command
```
bsreadsim ... [to add]
```

## Input file format  
click here to see the [example](./data/example/) data
- reference genome: [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format ( `.fasta`, `.fa` or `.fa.gz` ) 
- methylation profile: [CGmap](https://bsbolt.readthedocs.io/en/latest/methylation_calling/) format ( `.CGmap` or `.CGmap.gz` ) 
- genetic variants: [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format, v4.0 or later, **sorted** by chromosome id and position ( `.vcf` or `.vcf.gz` )
- probes for targeted sequencing: [BED](https://genome.ucsc.edu/FAQ/FAQformat.html) format, contains at least 6 columns, **sorted** by chromosome id and postion ( `bed` or `.bed.gz` )


## 3.1 Simulate WGBS data
1. with randomly generated genetic variants and methylation levels
```
```

2. with specified genetic variant profiles
```
```

3. with specified methylation profiles
```
```

4. simulate site-site dependency
```
```

## 3.2 Simulate RRBS data


## 3.3 Simulate TBS data


# 4. Contact
Please raise up issues through the github [issue](https://github.com/wbvguo/BSReadSim/issues) page

# 5. Reference
[1] Sherman: https://www.bioinformatics.babraham.ac.uk/projects/sherman/

[2] Farrell, C., Thompson, M., Tosevska, A., Oyetunde, A., & Pellegrini, M. (2021). BiSulfite Bolt: A bisulfite sequencing analysis platform. GigaScience, 10(5), giab033.

[3] Xie, Q., Liu, Q., Mao, F., Cai, W., Wu, H., You, M., ... & Wu, J. (2014). A Bayesian framework to identify methylcytosines from high-throughput bisulfite sequencing data. PLoS computational biology, 10(9), e1003853.