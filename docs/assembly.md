#  Mapping Assembly
We provide the [mapping_assembly.py](#mapping_assemblypy) script to perform the entire process from reads to assembled genomes. In addition, we provide [polish_asfv.py](#polish_asfvpy) to polish the homopolymers, [download_asfv_genome.py](#download_asfv_genomepy) to download all latest ASFV genomes from NCBI, and [find_near_ref.py](#find_near_refpy) to select the ASFV genome with the most mapped reads as the reference genome.
## mapping_assembly.py
### Description
The task requires the ONT reads and a reference sequence. If a reference sequence is specified directly through the -r parameter, it will be used directly. Alternatively, a directory containing multiple reference sequences can also be specified through the -r parameter. We have provided a directory "single_fasta" containing 406 genomes on [github](https://github.com/nimua/single_fasta.git). You can also download the latest ASFV genome by executing download_asfv_genome.py. If a directory is specified through the -r parameter, the genome exhibiting the highest read mapping coverage among all available ASFV genomes will be selected as the reference genome. This reference genome will then be utilized for mapping assembly of the sequenced data. The alignment is performed using the minimap2 aligner with the -a option, and ONT reads as input. A consensus sequence is generated using samtools. The consensus sequence is further polished with medaka.
### Arguments
| Argument name	  | Required | Description |
| --------------  | ----- | -------- |
| -p, --processes |  No  | number of processes (default = 4)   |
| -i,	--input   |  Yes  | input ASFV reads (fasta or fastq file) |
| -r, --ref     |  Yes  | a fasta file of ASFV genome or a folder containing multiple ASFV genomes (if the it is a folder, the program will automatically select the nearest one as the reference)|
| -o, --output   |  Yes  | file name of the output of assembled ASFV genome  |
| --medaka      |  No  | medaka model  |

### Example
```bash
mapping_assembly.py -p 4 -r single_fasta -i test_data.fasta -o asfv_genome.fasta --medaka r941_min_high_g303
```
### Output
A fasta file of the assembled ASFV genome.

------------------------------------
## polish_asfv.py
### Description
Polish the homopolymers by homopolish. It is recommended to choose the closest non-ONT sequenced ASFV genome as the reference genome in NCBI by blastn.
### Arguments
| Argument name	  | Required | Description |
| --------------  | ----- | -------- |
| -i,	--input   |  Yes  | fasta file of input ASFV genome |
| -r, --ref     |  Yes  | fasta file of reference ASFV genome |
| -m, --model   |  Yes  | model used in homopolish |

### Example
```bash
polish_asfv.py -i single_fasta/MN194591.1.fasta -r single_fasta/OR180113.1.fasta -m R9.4.pkl
```
### Output
Polished ASFV genome in the current working directory.

------------------------------------
## download_asfv_genome.py
### Description
Download all ASFV sequences with a length ranging from 160,000 to 250,000 from NCBI.
### Example
```bash
download_asfv_genome.py
```
### Output
A directory "single_fasta" containing all ASFV genomes from NCBI in the current working directory.

------------------------------------
## find_near_ref.py
### Description
From multiple references, get the nearest reference (The genome exhibiting the highest read mapping coverage among all available ASFV genomes). This task has been integrated into mapping_assembly.py. If you want to use this task separately, you can use this script.
### Arguments
| Argument name	  | Required | Description |
| --------------  | ----- | -------- |
| -f, --file |  Yes  |  input ASFV reads (fasta or fastq file)  |
| -r, --reference |  Yes  |  dir of the reference files  |
| -c, --core |  No  |  number of processes (default = 32)   |
| -q, --qscore |  No  |  the qscore used to filter the bam file (default= 0)  |
| -m, --mapper |  No  |  the mapper used to map, can be minimap2 or bwa (default="minimap2") |

### Example
```bash
find_near_ref.py -r ./single_fasta -f read.fastq > near.fasta 
```
### Output
A genome file exhibiting the highest read mapping coverage among all available ASFV genomes.
