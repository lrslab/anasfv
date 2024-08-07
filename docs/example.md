# Example Projects
## Docker
Using Docker is more convenient, eliminating the troublesome of installation steps. So we provide a docker image:
```bash
docker pull osvolo/anasfv:latest
docker container run -it osvolo/anasfv /bin/bash
```	 
## Prepare your data
The test files are already included in the docker container. If you don't use Docker, you can get it by the following method.

① ONT reads (fasta or fastq file). You can download our test file using the following command.
```bash
wget https://github.com/lrslab/anasfv/releases/download/test_data.fasta/test_data.fasta
```

② Other ASFV genomes. These genomes are used for mapping assembly and tree building. You can directly use the single_fasta directory in this project, which contains 312 downloaded ASFV genomes.
```bash
git clone https://github.com/nimua/single_fasta.git
```
Or you can run download_asfv_genome.py, which will create a single_fasta directory in the working directory and download all the latest ASFV genomes on NCBI to the single_fasta directory.
```bash
download_asfv_genome.py
```

### Task 1 (Assembling a genome):
Finding nearest genome from "./single_fasta" as ref to perform mapping assebly. Then two rounds of polish with medaka and homopolish.
```bash
mapping_assembly.py -p 4 -r single_fasta -i test_data.fasta -o genome.fasta --medaka r941_min_high_g303 --homopolish R9.4.pkl 
```
### Task 2 (Genome completeness evaluation):
We only established consensus gene sets for genotype I and genotype II. Using -c to assign consensus gene sets.
Using OQ504956.1.fasta as example：
```bash
completeness.py single_fasta/OQ504956.1.fasta -c II > OQ504956.1_completeness.tsv
```
### Task 3 (Recombination test):
Using OQ504956.1 as example：
```bash
recombination_test.py single_fasta/OQ504956.1.fasta > OQ504956.1_recombination_test.tsv
```
### Task 4 (Constructing a tree):
Using all genome files from "./single_fasta" and get a final file "tree.nwk" in Newick format
```bash
make_tree.py -f single_fasta -o tree
```
