## ANASFV
The **ANASFV** project focuses on analyzing nanopore-sequenced data of PCR-amplified African Swine Fever Virus (ASFV). It consists of 4 parts:

Part 1: Using ONT reads of PCR-amplified ASFV to assemble a genome.

Part 2: Analyzing the completeness of the assembled ASFV genome.

Part 3: Checking for any evidence of recombination between genotypes I and II.

Part 4: Constructing a phylogenetic tree.

Full documentation is available at [read the docs](https://anasfv.readthedocs.io/en/latest/).


<a href="https://pypi.python.org/pypi/anasfv" rel="pypi">![PyPI](https://img.shields.io/pypi/v/anasfv?color=green) </a>


## Installation:

Requirements：
1. python: 3.11 (tested). Most Python 3 versions should work.
2. Software versions tested:
	 \- Medaka: 1.11.3
 	 \- Samtools: 1.17
  	 \- BEDTools: 2.26.0
  	 \- Minimap2: 2.17-r941
  	 \- NanoFilt: 2.8.0
  	 \- Homopolish: 0.4.1
  	 \- Prodigal: 2.6.3
  	 \- Exonerate: 2.4.0
  	 \- blast: 2.12.0
  	 \- MUSCLE: 5.1

  	 
 
Install requirements:
```
conda install -c bioconda medaka samtools bedtools minimap2 nanofilt prodigal exonerate blast muscle
conda install -c conda-forge -c bioconda homopolish=0.4.1=pyhdfd78af_1
```
Install ANASFV by PyPI:
```
pip install anasfv
```

It may not be easy to successfully install all the dependent software. So we provide a docker installation:
```
docker.....
```

## A Quick Example:
### Prepare data：
1. Test data: One line of command downloads test.fq to the working directory
```
download xxxxxxx
```
2. Other ASFV genomes. These genomes are used for mapping assembly and tree building. You can directly use the single_fasta directory in this project, which contains 312 downloaded ASFV genomes, or you can run download_asfv_genome.py, which will create a single_fasta directory in the working directory and download all the latest ASFV genomes on NCBI to the single_fasta directory.
```
download_asfv_genome.py
```
### Part 1 (Assembling a genome):
Finding nearest genome from "./single_fasta" as ref to perform mapping assebly. Then two rounds of polish.
```
mapping_assembly.py -p 4 -r single_fasta -i test.fq -o genome.fasta --medaka r941_min_high_g303 --homopolish R9.4.pkl 
```
### Part 2 (Genome completeness evaluation):
We only established consensus gene sets for genotype I and genotype II. Using -c to assign consensus gene sets.
Using OQ504956.1.fasta as example：
```
completeness.py single_fasta/OQ504956.1.fasta -c II > OQ504956.1_completeness.tsv
```
### Part 3 (Recombination test):
Using OQ504956.1 as example：
```
recombination_test.py single_fasta/OQ504956.1.fasta > OQ504956.1_recombination_test.tsv
```
### Part 4 (Constructing a tree):
We installed uDance in the docker environment. If you are using docker, you can directly build the tree with the following command. It will use all genome files from "./single_fasta" and get a tree in Newick format
```
make_tree.py -f single_fasta -o tree
```

If you do not use docker. You can use the following command to find CDS in all genome files from "./single_fasta" and get an "./aligenments" directory as input for uDance. Then perform a tree construction in de-novo mode and an iterative tree construction in tree mode. Refer to [uDance](https://github.com/balabanmetin/uDance)
```
get_cds_alignments.py -f single_fasta
```






