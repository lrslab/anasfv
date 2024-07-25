# Recombination Test
With the emergence of recombinant ASFV of genotype I and genotype II, we use recombination_test.py to quickly identify whether an ASFV genome is recombinant and what its recombination pattern is.
## recombination_test.py
### Description
We generated consensus CDS sequences from genotype I and II isolates. The CDS predicted from input ASFV genome using Prodigal were compared to consensus sequences using BLASTN (e-value 1e-5). Output a table indicating whether each CDS is closer to genotype I or genotype II.
### Input
An ASFV genome file to be tested.
### Example
```
recombination_test.py single_fasta/OQ504956.1.fasta > OQ504956.1_recombination_test.tsv
```
### Output
A 6-column tsv table, the 6 fields are **CDS name**, **Start**, **End**, **Similarity with GenotypeI**, **Similarity with GenotypeII** and **Conclusion**.
The following table is a partial example of a recombinant

| CDS | Start | End | GenotypeI | GenotypeII | Conclusion |
| --- | ----- | --- | --------- | ---------- | ---------- |
|F165R |	54963	| 55373 |	100 |	97.567 |	I|
|F1055L |	55360 |	58527 |	99.937 |	97.277 |	I|
|K205R |	58696 |	59313 |	100 |	97.896 |	I|
|K78R |	59402 |	59638 |	99.578 |	100 |	II|
|K196R |	59635 |	60225 |	98.816 |	100 |	II|
|K145R |	60241 |	60678 |	98.63 |	100 |	II|
|K421R |	60712 |	61980 |	98.503 |	100 |	II|
|EP1242L |	62021 |	65749 |	99.142 |	98.713 |	I|
|EP84R |	65828 |	66082 |	99.608 |	95.686 |	I|
|EP424R | 66118 |	67377 |	97.536 |	98.73 |	II|
|EP152R	| 67872 |	68330 |	97.386 |	100 |	II|
