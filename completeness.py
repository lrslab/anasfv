#!/usr/bin/env python
# coding: utf-8
from Bio import SeqIO
import os
import subprocess
import pandas as pd
import multiprocessing


def cal_BUSCO(consensus_use, prodigal_file, with_MGF=True):
    consensusfile = f'./consensus_ffn/type{consensus_use}_consensus'
    df_querys = pd.DataFrame(
        columns=['name', 'Subject_id', 'identity', 'alignment_length', 'mismatches', 'gap_openings', 'q.start', 'q.end',
                 's.start', 's.end', 'e-value', 'bit_score'])
    total_genes = [i.id for i in list(SeqIO.parse(consensusfile + '.ffn', "fasta"))]
    if not with_MGF:
        total_genes = [i for i in total_genes if 'MGF' not in i]

    total = len(total_genes)
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
        nMissing = 0

    df_querys = df_querys.reset_index()
    for i in range(df_querys.shape[0]):
        if not with_MGF and 'MGF' in df_querys['Subject_id'][i]:
            pass
        elif df_querys['identity'][i] >= 90 and df_querys['alignment_length'][i] / (
                df_querys['s.end'][i] - df_querys['s.start'][i] + 1) > 0.9:
            l_Complete.append(df_querys['Subject_id'][i])
        elif df_querys['identity'][i] >= 30:
            l_Fragmented.append(df_querys['Subject_id'][i])

    l_Fragmented = list(set(l_Fragmented) - set(l_Complete))
    nComplete = len(set(l_Complete))
    nDuplicate = len(l_Complete) - nComplete
    nFragmented = len(set(l_Fragmented))
    nMissing = total - nComplete - nFragmented
    complete_percentage = round((nComplete / total) * 100, 2)
    duplicate_percentage = round((nDuplicate / total) * 100, 2)
    fragmented_percentage = round((nFragmented / total) * 100, 2)
    missing_percentage = round((nMissing / total) * 100, 2)

    busco = f'C:{complete_percentage}%[D:{duplicate_percentage}%],F:{fragmented_percentage}%,M:{missing_percentage}%,n:{total}'
    return busco

def run_blast(strain, consensus_use):
    file = f'./prodigal_result/{strain}.fna'
    prodigal_file = list(SeqIO.parse(file, "fasta"))
    consensusfile = f'./consensus_ffn/type{consensus_use}_consensus'
    for record in prodigal_file:
        SeqIO.write(record, f"./type{consensus_use}_blast/{record.id}.fasta", "fasta")
        os.system(f"blastn -query ./type{consensus_use}_blast/{record.id}.fasta -out ./type{consensus_use}_blast/{record.id}.blast -db {consensusfile} -outfmt 6 -evalue 1e-5 -num_threads 12")


if __name__ =='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='The consensus cds of genotype II is used by default.')
    parser.add_argument('consensus', type=str, default='II', help='I or II')
    args = parser.parse_args()
    consensus_use = args.consensus

    single_fasta_list = [i[:-6] for i in os.listdir("./single_fasta") if i.endswith('.fasta')]
    os.system('rm -rf prodigal_result')
    os.makedirs('prodigal_result')

    for single_fasta in single_fasta_list:
        prodigal_cmd = f'prodigal -i ./single_fasta/{single_fasta}.fasta -o ./prodigal_result/{single_fasta}.gff -f gff -a ./prodigal_result/{single_fasta}.faa -d ./prodigal_result/{single_fasta}.fna'
        os.system(prodigal_cmd)

    try:
        os.makedirs(f'type{consensus_use}_blast')
    except FileExistsError:
        pass

    p = multiprocessing.Pool(60)
    for strain in [i[:-4] for i in os.listdir('./prodigal_result') if i.endswith(".fna")]:
        p.apply_async(run_blast, args=(strain, consensus_use))
    p.close()
    p.join()

    print('accession_ID\tsize\tprodigal_gene_num\twith_MGF\twithout_MGF')
    for strain in [i[:-4] for i in os.listdir('./prodigal_result') if i.endswith(".fna")]:
        file = f'./prodigal_result/{strain}.fna'
        prodigal_file = list(SeqIO.parse(file, "fasta"))
        size = len(SeqIO.read(f'./single_fasta/{strain}.fasta', "fasta"))
        busco = cal_BUSCO(consensus_use, prodigal_file)
        busco_without_MGF = cal_BUSCO(consensus_use, prodigal_file, with_MGF=False)
        print(f'{strain}\t{size}\t{len(prodigal_file)}\t{busco}\t{busco_without_MGF}')
