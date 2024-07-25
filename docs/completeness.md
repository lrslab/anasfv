# Genome Completeness Evaluation
BUSCO is a computational tool used to assess the completeness and quality of genome assemblies(23). ASFV is not included in the BUSCO’s database (OrthoDB). Therefore, we developed a completeness evaluation system to generate a BUSCO-like notation. 
## completeness.py
### Description
We generated consensus sequences from genotype I and II isolates. The CDS predicted from input ASFV genome using Prodigal were compared to consensus sequences using BLASTN (e-value 1e-5). This comparison yielded a BUSCO like genome completeness evaluation.
### Input
An ASFV genome file.
### Arguments
| Argument name	  | Required | Description |
| --------------  | ----- | -------- |
| -c |  No  | consensus sequences to be use: I or II (default='II')  |

### Example
```
completeness.py single_fasta/OM966717.1.fasta -c II > OM966717.1_completeness.tsv
```
### Output
A BUSCO like genome completeness evaluation, with C:complete [D:duplicated], F:fragmented, M:missing, n:number of genes used. If a consensus CDS term can find a mapping from the predicted gene sequences with an identity larger than 90%, together with a unique mapping length longer than 90%, it can be considered "complete." A consensus CDS term that cannot find a valid mapping (with an identity greater than 30% and a mapping length greater than 30%) in the predicted genes is considered "missing". The other consensus CDS with partial hit is termed "fragmented". The "duplicated" term means that there is more than one 'complete' hit in the predicted gene sequences.
