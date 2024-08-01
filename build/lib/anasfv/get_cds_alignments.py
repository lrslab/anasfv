#!/usr/bin/env python
# coding: utf-8



import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO, SeqIO, SeqUtils
try:
    from cStringIO import StringIO
except ImportError:
    from io import StringIO
import subprocess
import sys
import signal
from collections import OrderedDict
import os
import pandas as pd
import Bio.SeqIO
import pkgutil



def fasta2dic(fastafile):
        if ".gz" in fastafile:
            handle=gzip.open(fastafile, "r")
        else:
            handle=open(fastafile, "r")
        record_dict=SeqIO.to_dict(SeqIO.parse(handle,"fasta"))
        handle.close()
        return record_dict

def myexe(cmd, timeout=0):
    """
    a simple wrap of the shell
    mainly used to run the bwa mem mapping and samtool orders
    """
    def setupAlarm():
        signal.signal(signal.SIGALRM, alarmHandler)
        signal.alarm(timeout)

    def alarmHandler(signum, frame):
        sys.exit(1)

    proc=subprocess.Popen(cmd, shell=True, preexec_fn=setupAlarm,
                 stdout=subprocess.PIPE, stderr=subprocess.PIPE,cwd=os.getcwd())
    out, err=proc.communicate()
    print (err, "Run finished with return code:", proc.returncode)
    return out


##################CDS part
def exonerate_wrapper(query, target, outfile=None, geneticcode=5, score=100, bestn=1):
    """
    --geneticcode 5

    return is a outfile name in relative path
    todo: using stringIO to hinder the file IO
    """
    if bestn =="auto":
        bestn=len(fasta2dic(target)) # default, output one region for one query

    exonerate_cmd="""exonerate {query} {target} \
                   --geneticcode {geneticcode} \
                   --score {score} \
                   --bestn {bestn} \
                   --model ungapped
                   """.format(
                        query=query, target=target,
                        geneticcode=geneticcode,
                        score=score,
                        bestn=bestn)
    out=myexe(exonerate_cmd)

    ## trigger to write the outfile to disk
    if outfile:
        outname=query.split("/")[-1].split(".")[0]+".exonerate"
        with open(outname, "w") as fw:
            fw.write(outname)

    return out


def exonerate_parser(exonerate_file):
    """
    parser the exonerate result, and return the position of the feather in 4-col bed format
    4 col bed4: [chro, start,end, name], example ["seq1", 1, 55, "trnP"]
    :param query:
    :param exonerate_file:
    :param prefix:
    :return: list of bed4
    """
    #fw=open(tbl_outname, "w") # change IO to list store
    bed4=[]

    texts=SearchIO.parse(StringIO(exonerate_file), format="exonerate-text")
    try:
        for record in texts:
            for hsp in record:
                for s in hsp:
                    # the biopython.SearchIO interval is 0 based [start, end), so start+1, end+0 to get 1 based coords
                    table_4=[s.fragment.query_id, s.fragment.query_start+1, s.fragment.query_end,s.fragment.hit_id]
                    bed4.append(table_4)
    except ValueError as e:
        pass
    bed4.sort()
    return bed4


def exonerate_parser_write(query, exonerate_file, prefix=None):
    """
    parser the exonerate result, and return the protein and cds file
    modification: add the position to the name of the cds and pro, use space to add interval
    :param query:
    :param exonerate_file:
    :param prefix:
    :return:
    """

    
    ########### functional code
    ref_dict=fasta2dic(query)

    if prefix is None:
        prefix=query.split(".")[0]

    p_outname=(prefix+"_exonerate_p.fa")
    cds_outname=(prefix + "_exonerate_cds.fa")

    fw_p=open(p_outname, "w")
    fw_cds=open(cds_outname, "w")

    texts=SearchIO.parse(StringIO(exonerate_file.decode('ascii')), format="exonerate-text")
    
    # to avoid write duplications
    used_name=set()
    
    for record in texts:
        for hsp in record:
            for s in hsp:
                #print(s.fragment.hit_id)
                name_str=">"+s.fragment.hit_id
                cds_str=str(s.fragment.hit.seq)


                #return a species#gene fasta
                single_name=s.fragment.query_id
                if single_name not in used_name:
#                     fw_p.write(name_str+"#"+name_cds+"\n"+p_str+"\n")
                    fw_cds.write(name_str+"#"+single_name+"\n"+cds_str+"\n")
                    used_name.add(single_name)

    return cds_outname, p_outname


def flow_exon(query, target, outfile=None, geneticcode=5, prefix=None):
    outfile=exonerate_wrapper(query, target, outfile, geneticcode)
    cds_out, p_out=exonerate_parser_write(query, outfile, prefix)
    
    
def get_unaligned_cds(all_cds,output_dir):
    """
    
    Parser the all_cds file. Put the same genes of each different strain in the same file.
    
    """
    records = SeqIO.parse(all_cds, "fasta")
    gene_sequences = {}
    l_species=[]
    for record in records:
        description = record.description
        species, gene = description.split("#")
        l_species.append(species)
        if gene not in gene_sequences:
            gene_sequences[gene] = []

        record.id=record.id.split('#')[0]
        record.description=''
        gene_sequences[gene].append(record)

    for gene, sequences in gene_sequences.items():
        with open(output_dir+gene + ".fasta", "w") as f:
            for sequence in sequences:
                f.write(f">{sequence.id}\n")
                f.write(f"{str(sequence.seq)}\n")



if __name__=="__main__":
    import argparse
    
    if not os.path.exists('ref_cds'):
        os.mkdir('ref_cds')
    data = pkgutil.get_data('anasfv', f'ref_cds/asfv_cds.fa')
    with open(f'ref_cds/asfv_cds.fa', 'wb') as f:
        f.write(data)
    cds_file = 'ref_cds/asfv_cds.fa'
    
    example_text = '''example:
        get_cds_alignments.py -f ./single_fasta
        '''

    parser = argparse.ArgumentParser(prog='get_cds_alignments.py',
                                     description='Run bash cmd lines for files',
                                     epilog=example_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    
    parser.add_argument("-f", "--file", help="dir of the fasta files")
    parser.add_argument("-c", "--core", help="the core", default= 32)
    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    
    single_fasta_list=[i for i in os.listdir(args.file) if i.endswith('.fasta')]
    for single_fasta in single_fasta_list:
        flow_exon(cds_file, f"{args.file}/"+single_fasta, geneticcode=11, prefix=single_fasta[:-6])

    os.system("rm single_process -rf")
    os.makedirs("single_process")
    os.system("mv *exonerate* single_process/")
    os.system("cat single_process/*exonerate_cds.fa > all_cds.fa")
    os.system("rm unaligned_cds alignments -rf")
    os.makedirs("unaligned_cds")
    os.makedirs("alignments")    
    
    all_cds="all_cds.fa"
    output_dir='unaligned_cds/'
    get_unaligned_cds(all_cds,output_dir)
    
    
    for file in os.listdir('./unaligned_cds'):
        os.system(f'muscle -super5 ./unaligned_cds/{file} -output ./alignments/{file} -nt -threads {args.core}')
