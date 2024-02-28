#!/usr/bin/env python
# coding: utf-8
from Bio import SeqIO
import os
import subprocess
import pandas as pd
from completeness import run_blast


def like_I_or_II(single_fasta,prodigal_file):

    run_blast(single_fasta,'Imore')
    run_blast(single_fasta,'IImore')
    
    print(f"CDS\tStart\tEnd\tGenotypeI\tGenotypeII\tConclusion\n")
    for record in prodigal_file:
        
        SeqIO.write(record, f"./prodigal_result/{record.id}.fasta", "fasta")

        cds=None
        start=record.description.split(' # ')[1]
        end=record.description.split(' # ')[2]
        

        
        try:
            df_query=pd.read_csv(f'typeImore_blast/{record.id}.blast',sep='\t', header=None)
        except:
            typei=None
        else:
            typei=df_query[2][0]
            cds=df_query[1][0]
            
        try:
            df_query=pd.read_csv(f'typeIImore_blast/{record.id}.blast',sep='\t', header=None)
        except:
            typeii=None
        else:
            typeii=df_query[2][0]
            if not cds:
                cds=df_query[1][0]
        
        conclusion=None
        if typei and typeii:
            if typei>typeii:
                conclusion='I'
            elif typei<typeii:
                conclusion='II'             
            
        if not cds:
                cds=record.id
                
        print(f"{cds}\t{start}\t{end}\t{typei}\t{typeii}\n")
        os.remove(f'./prodigal_result/{record.id}.fasta')



if __name__ =='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help='path of an ASFV genome file')
    args = parser.parse_args()
    input_path = args.input
    single_fasta=os.path.basename(input_path)
    
    os.system('rm -rf prodigal_result')
    os.makedirs('prodigal_result')

    prodigal_cmd = f'prodigal -i {input_path} -o ./prodigal_result/{single_fasta}.gff -f gff -a ./prodigal_result/{single_fasta}.faa -d ./prodigal_result/{single_fasta}.fna'
    os.system(prodigal_cmd)

    prodigal_file=list(SeqIO.parse(f'./prodigal_result/{single_fasta}.fna', "fasta"))
    like_I_or_II(single_fasta,prodigal_file)