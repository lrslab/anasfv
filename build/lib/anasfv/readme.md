## ANASFV
The **ANASFV** project focuses on analyzing nanopore-sequenced data of PCR-amplified African Swine Fever Virus (ASFV). It consists of 4 tasks:

Task 1: Using ONT reads of PCR-amplified ASFV to assemble a genome.

Task 2: Analyzing the completeness of the assembled ASFV genome.

Task 3: Checking for any evidence of recombination between genotypes I and II.

Task 4: Constructing a phylogenetic tree.

Full documentation is available at [read the docs](https://anasfv.readthedocs.io/en/latest/).


<a href="https://pypi.python.org/pypi/anasfv" rel="pypi">![PyPI](https://img.shields.io/pypi/v/anasfv?color=green) </a>


## Docker:
We provide a docker image. In case you find the installation too troublesome:
```
docker pull osvolo/anasfv:latest
docker container run -it osvolo/anasfv /bin/bash
```
## Installation:
Requirements：
1. python: 3.11 (tested). Most Python 3 versions should work.
2. Software versions tested:
 	 \- Samtools: 1.17
  	 \- BEDTools: 2.26.0
  	 \- Minimap2: 2.17-r941
  	 \- Prodigal: 2.6.3
  	 \- Exonerate: 2.4.0
  	 \- blast: 2.12.0
  	 \- MUSCLE: 5.1
   	 \- Medaka: 1.11.3
  	 \- Homopolish: 0.4.1
  	 \- uDance: 1.6.5
 
Install requirements in conda environment and install ANASFV via PyPI:
```
conda create -n anasfv -c conda-forge python=3.11 -y
conda activate anasfv
conda install -c bioconda samtools bedtools minimap2 prodigal exonerate blast muscle -y
pip install anasfv
```

If you need to use medaka and homopolish for polish, you need to create their corresponding conda environments and install them, because there will be some conflicts if you install them directly in the ANASFV runtime environment.
```
conda create -n medaka -c bioconda -c conda-forge medaka -y
conda config --set channel_priority flexible
conda create -n homopolish -c conda-forge -c bioconda -c defaults more-itertools=8.4.0 homopolish=0.4.1 -y
```

The tree building process uses uDance. For uDance installation refer to [uDance](https://github.com/balabanmetin/uDance)

## A Quick Example:
### Prepare data：
1. Test data: Downloads test_data.fasta to the working directory
```
wget https://github.com/lrslab/anasfv/releases/download/test_data.fasta/test_data.fasta
```
2. Other ASFV genomes. These genomes are used for mapping assembly and tree building. You can directly use the [single_fasta](https://github.com/nimua/single_fasta.git), which contains 312 downloaded ASFV genomes, or you can run download_asfv_genome.py, which will create a directory name "single_fasta" and download all the latest ASFV genomes on NCBI to the directory.
```
download_asfv_genome.py
```
### Task 1 (Assembling a genome):
Finding nearest genome from "./single_fasta" as ref to perform mapping assebly. Then two rounds of polish.
```
mapping_assembly.py -p 4 -r single_fasta -i test_data.fasta -o genome.fasta --medaka r941_min_high_g303 --homopolish R9.4.pkl 
```
### Task 2 (Genome completeness evaluation):
We only established consensus gene sets for genotype I and genotype II. Using -c to assign consensus gene sets.
Using OQ504956.1.fasta as example：
```
completeness.py single_fasta/OQ504956.1.fasta -c II > OQ504956.1_completeness.tsv
```
### Task 3 (Recombination test):
Using OQ504956.1 as example：
```
recombination_test.py single_fasta/OQ504956.1.fasta > OQ504956.1_recombination_test.tsv
```
### Task 4 (Constructing a tree):
Building the tree with the following command. It will use all genome files from "./single_fasta" and get a tree in Newick format.
```
make_tree.py -f single_fasta -o tree
```


