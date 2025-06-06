# Installation

## Docker
We provide a docker image to make the installation steps more convenient:
```bash
docker pull osvolo/anasfv:latest
docker container run -it osvolo/anasfv /bin/bash
```	 
## Install requirements in conda environment and install ANASFV via PyPI
Requirements:

1. python: 3.11 (tested). Most Python 3 versions should work.
   
2. Software versions tested:
 	  \- Samtools: 1.17
  	 \- BEDTools: 2.26.0
  	 \- Minimap2: 2.17
  	 \- Prodigal: 2.6.3
  	 \- Exonerate: 2.4.0
  	 \- BLAST: 2.12.0
  	 \- MUSCLE: 5.1
    \- Medaka: 1.11.3
    \- Homopolish: 0.4.1	
  	 \- uDance: 1.6.4
```bash
conda create -n anasfv -c conda-forge python=3.11 -y
conda activate anasfv
conda install -c bioconda samtools bedtools minimap2 prodigal exonerate blast muscle -y
pip install anasfv
```

If you need to use medaka and homopolish for polish, you need to create their corresponding conda environments and install them, because there will be some conflicts if you install them directly in the ANASFV runtime environment.
```bash
conda create -n medaka -c bioconda -c conda-forge medaka=1.11.3 -y
conda config --set channel_priority flexible
conda create -n homopolish -c conda-forge -c bioconda -c defaults more-itertools=8.4.0 homopolish=0.4.1 -y
```

The tree building process uses uDance. For uDance installation refer to [uDance](https://github.com/balabanmetin/uDance)
