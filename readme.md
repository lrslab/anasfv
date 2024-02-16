
The **anasfv** project focuses on analyzing nanopore-sequenced data of PCR-amplified African Swine Fever Virus (ASFV). The workflow accepts fastq files as input and produces three main outputs: an assembled genome, a genome completeness evaluation similar to BUSCO, and a phylogenetic tree.

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
pip install biopython
pip install pandas
conda install -c bioconda medaka
conda install -c bioconda samtools
conda install -c bioconda bedtools
conda install -c bioconda minimap2
conda install -c bioconda nanofilt 
conda install -c conda-forge -c bioconda homopolish=0.4.1=pyhdfd78af_1
conda install -c bioconda prodigal
conda install -c bioconda exonerate
conda install -c bioconda blast
conda install -c bioconda muscle
```
## Workflow Example:

1. Download all ASFV genomes from NCBI to the "./single_fasta" directory
```
./download_asfv_genome.py
```
2. Merge asfv genome files to one file
( If you already have assembled genome, you can skip steps 2-8 )
```
cat single_fasta/*.fasta > allde.fasta
```
3. Trimming
```
NanoFilt all.fastq -q 10 -l 1000 --maxlength 200000 --headcrop 50 > all_trimmed.fq
```
4. Find nearest genome as ref
```
./find_near_ref.py -r allde.fasta -f all_trimmed.fq > near.fasta
```
5. Use near.fasta as ref to generate sam file
```
minimap2 -a near.fasta ./all_trimmed.fq > all-alignment.sam
```
6. Generate consensus file
```
samtools view -b -F 4 all-alignment.sam > all-alignment.bam
samtools sort -@ 8 -o all-sorted_alignment.bam all-alignment.bam
samtools consensus -f fasta all-sorted_alignment.bam -o all-assembled.fa
```
7. Polish with medaka
```
medaka_consensus -i all_trimmed.fq -d all-assembled.fa -o all-assembly_medaka_result -m r941_min_fast_g303 -t 2 > medaka.log
```
8. Polish with homopolish
( get a final genome)
```
homopolish polish -a ./all-assembly_medaka_result/consensus.fasta -l ./near.fasta -m R9.4.pkl -o homopolish-output
```
9. Rename final_genome file and move it to the "./single_fasta" directory.
```
cp ./homopolish_output/consensus_homopolished.fasta ./single_fasta/strain_name.fasta
sed -i '1s/.*/>strain_name/' ./single_fasta/strain_name.fasta
``` 
10. Genome completeness evaluation
```
./completeness.py II > completeness.tsv
```
11. Get aligenments for uDance
( find cds in all genome files from "./single_fasta" and get a "./aligenments" directory as input for uDance )
```
get_cds_alignments.py
```
12. Build a tree using uDance
Perform a tree construction in de-novo mode and an iterative tree construction in tree mode. For the configuration used, see "./uDance_config_example"
Refer to [uDance](https://github.com/balabanmetin/uDance)
