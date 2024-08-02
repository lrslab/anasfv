#!/usr/bin/env python
# coding: utf-8
from Bio import SeqIO
import os
import subprocess
import pandas as pd
from completeness import run_blast
import pkgutil



def like_I_or_II(single_fasta,prodigal_file):
    print("CDS\tStart\tEnd\tGenotypeI\tGenotypeII\tConclusion")
    run_blast(single_fasta,'Imore')
    run_blast(single_fasta,'IImore')

    d_i_seq_len = {}
    for record in SeqIO.parse('./consensus_ffn/typeImore_consensus.ffn', "fasta"):
        d_i_seq_len[record.id]=len(record.seq)
    
    d_ii_seq_len = {}
    for record in SeqIO.parse('./consensus_ffn/typeIImore_consensus.ffn', "fasta"):
        d_ii_seq_len[record.id]=len(record.seq)
        
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
            if df_query[3][0]/d_i_seq_len[df_query[1][0]]>0.8:
                typei=df_query[2][0]
                cds=df_query[1][0]
            else:
                typei = None
            
        try:
            df_query=pd.read_csv(f'typeIImore_blast/{record.id}.blast',sep='\t', header=None)
        except:
            typeii=None
        else:
            if df_query[3][0]/d_ii_seq_len[df_query[1][0]]>0.8:
                typeii = df_query[2][0]
                cds = df_query[1][0]
            else:
                typeii = None
        
        conclusion=None
        if typei and typeii:
            if typei>typeii:
                conclusion='I'
            elif typei<typeii:
                conclusion='II'             
            
        if not cds:
                cds=record.id
        
        conclusion=''
        if typei and typeii:
            if typei > typeii:
                conclusion='I'
            elif typei < typeii:
                conclusion='II'
        elif typei and typei==100:
            conclusion='I'
        elif typeii and typeii==100:
            conclusion='II'
                
        print(f"{cds}\t{start}\t{end}\t{typei}\t{typeii}\t{conclusion}")
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
    
    
    if not os.path.exists('consensus_ffn'):
        os.mkdir('consensus_ffn')
    for suffix in ['ffn','ndb','nhr','nin','not','nsq','ntf','nto']:
        for consensus_use in ['Imore','IImore']:
            data = pkgutil.get_data('anasfv', f'consensus_ffn/type{consensus_use}_consensus.{suffix}')
            with open(f'consensus_ffn/type{consensus_use}_consensus.{suffix}', 'wb') as f:
                f.write(data)
    
    like_I_or_II(single_fasta,prodigal_file)
