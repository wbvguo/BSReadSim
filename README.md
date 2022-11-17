# Introduction
BSReadSim (<u>B</u>isulfite <u>S</u>equencing <u>Read</u> <u>Sim</u>ulator) is a bisulfite sequencing reads simulator that allows for profile-based read simulation. To our best knowledge, BSReadSim is the only bisulfite sequencing read simulator that allows user to incooperate the genetic variant and methylation profile with high fidelity, thus it can generate more realistic bisulfite sequencing reads compared to other bisulfite read simulators. The following table summarized the difference between BSReadSim and other tools:

Features | Sherman | BSBolt | BSReadSim
----| ---- | ---- | ---- |
bisulfite conversion rate | yes | no | yes
methylation profile input | no | $\textsf{yes}^*$ | yes
genetic variant input | no   | no | yes
alleclic-specific methylation simulation | no | no | [todo]
site-site dependency | no | no | [todo]
multi-thread support | no | no | [todo]

( $\textsf{yes}^*$ : BSBolt allows users to input the methylation profile, but used the profile as a reference pool. During simulation, it **randomly** picks a methylation value from this pool. As a result, for a paticular CG site, the simulated data and the reference profile will likely not have the same methylation level)

# Installation
## Dependency
- Python 3.8 or later
- Biopython
- tqdm
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


# Usage
## Command
```
bsreadsim ... [to add]
```

## Input file format
1. reference genome: [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format (`.fasta`, `.fa` or `.fa.gz`)
2. methylation profile: [CGmap](https://bsbolt.readthedocs.io/en/latest/methylation_calling/) format (`.CGmap` or `.CGmap.gz`)
3. genetic variants: [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format, v4.2 or later (`.vcf` or `.vcf.gz`)


# Reference
[1] Sherman: https://www.bioinformatics.babraham.ac.uk/projects/sherman/

[2] Farrell, C., Thompson, M., Tosevska, A., Oyetunde, A., & Pellegrini, M. (2021). BiSulfite Bolt: A bisulfite sequencing analysis platform. GigaScience, 10(5), giab033.