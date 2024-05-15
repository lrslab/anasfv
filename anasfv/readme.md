
The **anasfv** project focuses on analyzing nanopore-sequenced data of PCR-amplified African Swine Fever Virus (ASFV). 

The workflow consists of 4 parts:

Part 1: Using nanopore-sequenced data of PCR-amplified ASFV to assemble a genome.

Part 2: Analyzing the completeness of the assembled ASFV genome.

Part 3: Checking for any evidence of recombination between genotypes I and II.

Part 4: Constructing a phylogenetic tree.

## Requirements:

1. python: 3.11 (tested). Most Python 3 versions should work.
2. Software versions tested:
	 \- Biopython: 1.8.1
	 \- Pandas: 2.0.3
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
pip install anasfv
```
## Workflow Example:
Get the number of available processors.
```
NPROC=$(nproc)
```
### Part 1 (Assembling a genome):
1. Download all ASFV genomes from NCBI to the "./single_fasta" directory.
```
download_asfv_genome.py
```
2. Trimming.
```
NanoFilt your_asfv_reads.fastq -q 10 -l 1000 --maxlength 200000 --headcrop 50 > all_trimmed.fq
```
3. Find nearest genome as ref.
```
find_near_ref.py -r single_fasta -f all_trimmed.fq > near.fasta
```
4. Use near.fasta as ref to generate sam file.
```
minimap2 -a near.fasta ./all_trimmed.fq > all-alignment.sam
```
5. Generate consensus file.
```
samtools view -b -F 4 all-alignment.sam > all-alignment.bam
samtools sort -@ ${NPROC} -o all-sorted_alignment.bam all-alignment.bam
samtools consensus -f fasta all-sorted_alignment.bam -o all-assembled.fa
```
6. Polish with medaka (For model selection, please refer to [medaka](https://github.com/nanoporetech/medaka#Models) ).
```
medaka_consensus -i all_trimmed.fq -d all-assembled.fa -o all-assembly_medaka_result -m <suitable_model> -t ${NPROC} > medaka.log
```
7. Polish with homopolish (model selection: R9.4.pkl/R10.3.pkl).
```
homopolish polish -a ./all-assembly_medaka_result/consensus.fasta -l ./near.fasta -m <suitable_model> -o homopolish-output
```
8. Rename final genome file and move it to the "./single_fasta" directory for Part 2 analysis.
```
#Rename final genome file and move it to the "./single_fasta" directory
cp ./homopolish_output/consensus_homopolished.fasta ./single_fasta/strain_name.fasta
#Change sequence ID
sed -i '1s/.*/>strain_name/' ./single_fasta/strain_name.fasta
```
### Part 2 (Genome completeness evaluation):
We only established consensus gene sets for genotype I and genotype II. Using -c to assign consensus gene sets.
Using OQ504956.1 as example£º
```
completeness.py ./single_fasta/OQ504956.1.fasta -c II > OQ504956.1_completeness.tsv
```
### Part 3 (Recombination test):
Using OQ504956.1 as example£º
```
recombination_test.py ./single_fasta/OQ504956.1.fasta > OQ504956.1_recombination_test.tsv
```
### Part 4 (Constructing a tree):
1. Download all ASFV genomes from NCBI to the "./single_fasta" directory (If it has already been downloaded in Part 1, please ignore this step).
```
download_asfv_genome.py
```
2. Get aligenments for uDance.
( find cds in all genome files from "./single_fasta" and get a "./aligenments" directory as input for uDance )
```
get_cds_alignments.py -f single_fasta
```
3. Build a tree using uDance.

Perform a tree construction in de-novo mode and an iterative tree construction in tree mode.
Refer to [uDance](https://github.com/balabanmetin/uDance)
