#!/usr/bin/env python
# coding: utf-8
from Bio import SeqIO
import os
import subprocess
import pandas as pd
import pkgutil
from collections import Counter

def cal_BUSCO(consensus_use, prodigal_file, with_MGF=True):       
    consensusfile = f'./consensus_ffn/type{consensus_use}_consensus'
    df_querys = pd.DataFrame(
        columns=['name', 'Subject_id', 'identity', 'alignment_length', 'mismatches', 'gap_openings', 'q.start', 'q.end',
                 's.start', 's.end', 'e-value', 'bit_score'])
    total_genes = [i.id for i in list(SeqIO.parse(consensusfile + '.ffn', "fasta"))]
    if not with_MGF:
        total_genes = [i for i in total_genes if 'MGF' not in i]

    total = len(total_genes)
    subject_cds=list(SeqIO.parse(f'{consensusfile}.ffn', "fasta"))
    d_subject_cds={}
    for cds in subject_cds:
        d_subject_cds[cds.id]=len(cds.seq)
        
    for record in prodigal_file:
        try:
            df_query = pd.read_csv(f'./type{consensus_use}_blast/{record.id}.blast', sep='\t', header=None)
        except:
            pass
        else:
            df_query.rename(columns={0: 'name', 1: 'Subject_id', 2: 'identity', 3: 'alignment_length', 4: 'mismatches',
                                     5: 'gap_openings', 6: 'q.start', 7: 'q.end', 8: 's.start', 9: 's.end',
                                     10: 'e-value', 11: 'bit_score'}, inplace=True)
            df_querys = pd.concat([df_querys, df_query.loc[:0]])

    l_Complete = []
    l_Fragmented = []

    df_querys = df_querys.reset_index()
    
    for i in range(df_querys.shape[0]):
        if not with_MGF and 'MGF' in df_querys['Subject_id'][i]:
            pass
        elif df_querys['identity'][i]>=90 and df_querys['alignment_length'][i]/d_subject_cds[df_querys['Subject_id'][i]] >0.9:
            l_Complete.append(df_querys['Subject_id'][i])
        elif df_querys['identity'][i]>=30 and df_querys['alignment_length'][i]/d_subject_cds[df_querys['Subject_id'][i]] >0.3:
            l_Fragmented.append(df_querys['Subject_id'][i])

    l_Fragmented = list(set(l_Fragmented) - set(l_Complete))
    l_Duplicate= [element for element, count in Counter(l_Complete).items() if count > 1]
    l_Missing=[i for i in total_genes if not i in l_Fragmented+l_Complete]
    nComplete = len(set(l_Complete))
    nDuplicate = len(l_Complete) - nComplete
    nFragmented = len(set(l_Fragmented))
    nMissing = total - nComplete - nFragmented
    complete_percentage = round((nComplete / total) * 100, 2)
    duplicate_percentage = round((nDuplicate / total) * 100, 2)
    fragmented_percentage = round((nFragmented / total) * 100, 2)
    missing_percentage = round((nMissing / total) * 100, 2)
    
    busco = f'C:{complete_percentage}%[D:{duplicate_percentage}%],F:{fragmented_percentage}%,M:{missing_percentage}%,n:{total}'
    return busco,', '.join(l_Fragmented),', '.join(l_Duplicate),', '.join(l_Missing)

def run_blast(strain, consensus_use):

    file = f'./prodigal_result/{strain}.fna'
    prodigal_file = list(SeqIO.parse(file, "fasta"))
    consensusfile = f'./consensus_ffn/type{consensus_use}_consensus'
    
    file = f'./prodigal_result/{strain}.fna'
    prodigal_file = list(SeqIO.parse(file, "fasta"))
    
    try:
        os.makedirs(f'type{consensus_use}_blast')
    except FileExistsError:
        pass
    for record in prodigal_file:
        SeqIO.write(record, f"./type{consensus_use}_blast/{record.id}.fasta", "fasta")
        os.system(f"blastn -query ./type{consensus_use}_blast/{record.id}.fasta -out ./type{consensus_use}_blast/{record.id}.blast -db {consensusfile} -outfmt 6 -evalue 1e-5 -num_threads 12")


if __name__ =='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='The consensus cds of genotype II is used by default.')
    parser.add_argument('input', type=str, help='input an ASFV genome file')
    parser.add_argument('-c', type=str, default='II', help='I or II')
    
    args = parser.parse_args()
    input_path = args.input
    consensus_use = args.c
    single_fasta=os.path.basename(input_path)

    os.system('rm -rf prodigal_result')
    os.makedirs('prodigal_result')
    
    
    prodigal_cmd = f'prodigal -i {input_path} -o ./prodigal_result/{single_fasta}.gff -f gff -a ./prodigal_result/{single_fasta}.faa -d ./prodigal_result/{single_fasta}.fna'
    os.system(prodigal_cmd)

    if not os.path.exists('consensus_ffn'):
        os.mkdir('consensus_ffn')
    for suffix in ['ffn','ndb','nhr','nin','not','nsq','ntf','nto']:
        data = pkgutil.get_data('anasfv', f'consensus_ffn/type{consensus_use}_consensus.{suffix}')
        with open(f'consensus_ffn/type{consensus_use}_consensus.{suffix}', 'wb') as f:
            f.write(data)
    
    run_blast(single_fasta,consensus_use)

    print('file_name\tsize\tprodigal_gene_num\twith_MGF\twithout_MGF\tduplicate_genes\tfragmented_genes\tmissing_genes')

    file = f'./prodigal_result/{single_fasta}.fna'
    prodigal_file = list(SeqIO.parse(file, "fasta"))
    size = len(SeqIO.read(input_path, "fasta"))
    busco,fragmented_genes,duplicate_genes,missing_genes = cal_BUSCO(consensus_use, prodigal_file)
    busco_without_MGF = cal_BUSCO(consensus_use, prodigal_file, with_MGF=False)[0]
    print(f'{single_fasta}\t{size}\t{len(prodigal_file)}\t{busco}\t{busco_without_MGF}\t{duplicate_genes}\t{fragmented_genes}\t{missing_genes}')
