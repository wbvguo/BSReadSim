# 1. Introduction
BSReadSim is a fast and flexible <u>B</u>isulfite <u>S</u>equencing <u>Read</u> <u>Sim</u>ulator that allows for profile-based read simulation. To our best knowledge, BSReadSim is the only bisulfite sequencing read simulator that allows user to incooperate the genetic variant and methylation profile with high fidelity. As a result, it can generate more realistic bisulfite sequencing reads compared to other bisulfite read simulators. The following table summarized the difference between BSReadSim and other tools:

Features | Sherman | BSBolt | BSReadSim
----| ---- | ---- | ---- |
variable bisulfite conversion rate | yes | no | yes
methylation profile input | no | $\textsf{yes}^*$ | yes
genetic variant input | no   | no | yes
alleclic-specific methylation simulation | no | no | yes
site-site dependency | no | no | [todo]
multi-thread support | no | no | yes

( $\textsf{yes}^*$ : BSBolt allows users to input the methylation reference. But during simulation, it **randomly** picks a methylation value from this reference pool. As a result, for a paticular CG site, the simulated data and the reference profile will likely have different methylation level)

# 2. Installation
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


# 3. Usage
## Command
```
bsreadsim ... [to add]
```

## Input file format
- reference genome: [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format (`.fasta`, `.fa` or `.fa.gz`)
- methylation profile: [CGmap](https://bsbolt.readthedocs.io/en/latest/methylation_calling/) format (`.CGmap` or `.CGmap.gz`)
- genetic variants: [VCF](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format, v4.0 or later(`.vcf` or `.vcf.gz`)

# 4. Contact
Please raise up issues through the github [issue](https://github.com/wbvguo/BSReadSim/issues) page

# 5. Reference
[1] Sherman: https://www.bioinformatics.babraham.ac.uk/projects/sherman/

[2] Farrell, C., Thompson, M., Tosevska, A., Oyetunde, A., & Pellegrini, M. (2021). BiSulfite Bolt: A bisulfite sequencing analysis platform. GigaScience, 10(5), giab033.