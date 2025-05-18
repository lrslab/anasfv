from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

'''
This script generates a random simulated recombinant ASFV genome for Supplymentary Figure S8.
i_ii.fasta contains genotype I genome (NC_044957.1) and genotype II genome (NC_044959.2), and the aligned fasta file is obtained through mafft:
mafft --auto i_ii.fasta > i_ii.aln
'''

fasta_file = "i_ii.aln"
recombinant_events=15
output_file = "sim.fasta"

# 1. Reading aligned file
with open(fasta_file, "r") as handle:
    l=[]
    for record in SeqIO.parse(handle, "fasta"):
        sequence_length = len(record.seq)
        l.append(str(record.seq))
gi,gii=l

# 2. Generate random numbers
sequence_length=len(gi)
random_numbers = random.sample(range(sequence_length), recombinant_events)
random_numbers=sorted(random_numbers)

# 3. Generating recombinant sequences
if recombinant_events%2==0:   
    recomb=gi[0:random_numbers[0]]
    print('i',0,'to',random_numbers[0])
else:
    recomb=gii[0:random_numbers[0]]
    print('ii',0,'to',random_numbers[0])

for i in range(len(random_numbers)-1):
    if i%2==0:
        recomb+=gii[random_numbers[i]:random_numbers[i+1]]
        print('ii',random_numbers[i],'to',random_numbers[i+1])
    else:
        recomb+=gi[random_numbers[i]:random_numbers[i+1]]
        print('i',random_numbers[i],'to',random_numbers[i+1])

if recombinant_events%2==0:    
    recomb+=gi[random_numbers[-1]:]
    print('i',random_numbers[-1],'to','end')
else:
    recomb+=gii[random_numbers[-1]:]
    print('ii',random_numbers[-1],'to','end')
    
recomb=recomb.replace('-','').upper()

# 4. Save
seq_record = SeqRecord(Seq(recomb), id="sim", description="Simulated Recombinant")
SeqIO.write(seq_record, output_file, "fasta")
