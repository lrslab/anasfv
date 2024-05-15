#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import Entrez, SeqIO
from collections import defaultdict
import os


def download_genome(search_term,download_path):
    """
    download asfv genomes from NCBI
    """
    Entrez.email = "your_email@xxx.com"
    handle = Entrez.esearch(db="nucleotide", term=search_term, retmax=1000)
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]
    for record_id in id_list:
        handle = Entrez.efetch(db="nucleotide", id=record_id, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        with open(download_path + record_id, "w") as file:
            file.write(fasta_data)
        records = SeqIO.parse(download_path + record_id, "fasta")
        record = next(records)
        accession_id=record.id.split()[0]
        with open(download_path+accession_id+".fasta", "w") as f:
            f.write('>'+accession_id+"\n")
            f.write(f"{str(record.seq)}\n")
        os.remove(download_path+record_id)
        print(accession_id+".fasta")
    print('Download completed')

def find_duplicate_sequences(fasta_files):
    """
    find identical genome
    """

    sequence_dict = defaultdict(list)
    duplicate_ids = []
    for file in fasta_files:
        for record in SeqIO.parse(file, "fasta"):
            sequence = str(record.seq)
            sequence_dict[sequence].append(record.id)
    for sequence, ids in sequence_dict.items():
        if len(ids) > 1:
            duplicate_ids.append(ids)

    return duplicate_ids

if __name__ == "__main__":

    try:
        os.makedirs('single_fasta')
    except FileExistsError:
        print("Note: the single_fasta folder already exists and may already contain previously downloaded content.")

    search_term = '((txid10497[organism:exp] AND biomol_genomic[prop])) AND 160000:250000[Sequence Length]'
    download_path='single_fasta/'
    download_genome(search_term,download_path)


    # # Delete those with poor quality.
    # for file in ['MN194591.1.fasta','MN318203.3.fasta','OF448913.1.fasta']:
    #     os.remove(download_path+file)
    # (These three genome didn't look that good, but in the end I kept them.)


    fasta_files = [download_path+i for i in os.listdir(download_path) if i.endswith('.fasta')]
    duplicate_ids = find_duplicate_sequences(fasta_files)

    # If "NC" exists, keep "NC", otherwise keep the first one
    for ids in duplicate_ids:
        for accession_id in ids:
            if accession_id.startswith('NC'):
                ids.remove(accession_id)
                break
            if accession_id == ids[-1]:
                del ids[0]

    # Remove remaining duplicates
    for ids in duplicate_ids:
        for accession_id in ids:
            os.remove(f'{download_path}{accession_id}.fasta')
            print(f'delete {accession_id}.fasta (identical to other genome)')


    # Even without performing duplicate removal in this step, uDance will only keep one instance of identical sequences when constructing the tree.
    # I want to keep the sequence with an accession ID starting with "NC", so I have removed duplicate sequences here.
