#  Mapping Assembly
We provide the [mapping_assembly.py](#mapping_assemblypy) script to perform the entire process from reads to assembled genomes. In addition, we provide [download_asfv_genome.py](#download_asfv_genomepy) to download all latest ASFV genomes from NCBI, and [find_near_ref.py](#find_near_refpy) to select the ASFV genome with the most mapped reads as the reference genome.
## mapping_assembly.py
### Description
The task requires the ONT reads and a reference sequence. If a reference sequence is specified directly through the -r parameter, it will be used directly. Alternatively, a directory containing multiple reference sequences can also be specified through the -r parameter. We have provided a directory "single_fasta" containing 312 genomes on [github](https://github.com/lrslab/anasfv). You can also download the latest ASFV genome by executing download_asfv_genome.py. If a directory is specified through the -r parameter, the genome exhibiting the highest read mapping coverage among all available ASFV genomes will be selected as the reference genome. This reference genome was then utilized for mapping assembly of the sequenced data. The alignment was performed using the minimap2 aligner with the -a option, and ONT reads as input. A consensus sequence was generated using samtools. The consensus sequence was further polished using two rounds of error correction. First, medaka was used. The ONT reads and the initial consensus sequence were provided as input to the medaka_consensus program. Next, homopolish was utilized to perform additional error correction.
### Arguments
| Argument name	  | Required | Description |
| --------------  | ----- | -------- |
| -p, --processes |  No  | number of processes (default = 4)   |
| -i,	--input   |  Yes  | input ASFV reads (fasta or fastq file) |
| -r, --ref     |  Yes  | a fasta file of ASFV genome or a folder containing multiple ASFV genomes (if the it is a folder, the program will automatically select the nearest one as the reference)|
| -o, --output   |  Yes  | file name of the output of assembled ASFV genome  |
| --medaka      |  No  | medaka model  |
| --homopolish   |  No  | homopolish model  |

### Example
```
mapping_assembly.py -p 4 -r single_fasta -i test_data.fasta -o asfv_genome.fasta --medaka r941_min_high_g303 --homopolish R9.4.pkl
```
### Output
A fasta file of the assembled ASFV genome.

## download_asfv_genome.py
### Description
Download all ASFV sequences with a length ranging from 160,000 to 250,000 from NCBI.
### Example
```
download_asfv_genome.py
```
### Output
A directory "single_fasta" containing all ASFV genomes from NCBI in the current working directory.

## find_near_ref.py
### Description
From multiple references, get the nearest reference. This task has been integrated into mapping_assembly.py. If you want to use this task separately, you can use this script.
### Arguments
| Argument name	  | Required | Description |
| --------------  | ----- | -------- |
| -f, --file |  Yes  |  input ASFV reads (fasta or fastq file)  |
| -r, --reference |  Yes  |  dir of the reference files  |
| -c, --core |  No  |  number of processes (default = 32)   |
| -q, --qscore |  No  |  the qscore used to filter the bam file (default= 0)  |
| -m, --mapper |  No  |  the mapper used to map, can be minimap2 or bwa (default="minimap2") |

### Example
```
find_near_ref.py -r ./single_fasta -f read.fastq > near.fasta 
```
### Output
A genome file exhibiting the highest read mapping coverage among all available ASFV genomes.
