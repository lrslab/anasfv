#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Time    : 3/2/2021 5:32 PM
# @Author  : Runsheng
# @File    : get_near_ref.py

"""
from multiple references, get the nearest reference for further polish
mostly used for RNA virus reference choosing
"""

from __future__ import print_function
import logging
import operator
import os
from anasfv.utils import myexe, fasta2dic, chr_select

# logger
from anasfv.logger import init_logger, get_logger
LOGGER=get_logger()


def wrapper_run_get_bedfile(ref, fastq,  core=32, qscore=0, mapper="minimap2"):
    """
    :param ref:
    :param fastq:
    :param core:
    :param qscore:
    :param mapper: minimap2 or bwa, bwa will use bwa mem
    :return:
    """
    if mapper=="minimap2":
        cmd_minimap2="""
        minimap2 -ax map-ont -t {core} {ref} {fastq} > map.sam
        """.format(ref=ref, fastq=fastq, core=core, qscore=qscore)
        cmd_use=cmd_minimap2

    if mapper=="bwa":
        cmd_bwa="""
        bwa mem -t {core} {ref} {fastq} | samtools view -h  -F 260 -q {qscore} -o map.sam
        """.format(ref=ref, fastq=fastq, core=core, qscore=qscore)
        cmd_use=cmd_bwa

    LOGGER.warning(cmd_use)
    myexe(cmd_use)

    prefix=ref.split("/")[-1].split(".")[0]

    cmd_samtools="""
    samtools view -h -F 260 -q {qscore} -b map.sam > map.bam
    samtools sort map.bam > maps.bam
    samtools index maps.bam
    bedtools genomecov -ibam maps.bam -bga > {prefix}.bed 
    rm map.sam
    rm map.bam
    """.format(qscore=qscore, prefix=prefix)
    LOGGER.warning(cmd_samtools)
    myexe(cmd_samtools)

    return prefix+".bed"


def bed_parser_get_higest_coverage(bedfile):
    """
    parser the bed file and get the refname which has higest coverage
    :param bedfile:
    :return:
    """
    cov_sum={}

    f=open(bedfile, "r")
    for line in f.readlines():
        line_l=line.strip().split("\t")
        name, start, end, coverage=line_l
        try:
            cov_sum[name]+=(int(end)-int(start)) * int(coverage)
        except KeyError:
            cov_sum[name] = (int(end) - int(start)) * int(coverage)
    sorted_d = sorted(cov_sum.items(), key=operator.itemgetter(1), reverse=True)
    LOGGER.warning(sorted_d[:10]) # show top 10 in logger
    chrname=str(sorted_d[0][0])
    LOGGER.warning("the highest coverage reference is: "+chrname)

    f.close()
    return chrname


if __name__=="__main__":
    import sys
    import argparse

    example_text = '''example:
        find_near_ref.py -r ./single_fasta -f read.fastq > near.fasta 
        '''

    parser = argparse.ArgumentParser(prog='get_near_ref.py',
                                     description='Run bash cmd lines for files',
                                     epilog=example_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("-f", "--file", help="the fastq file")
    parser.add_argument("-r", "--reference", help="dir of the reference files")
    parser.add_argument("-c", "--core", help="the core", default= 32)
    parser.add_argument("-q", "--qscore", help="the qscore used to filter the bam file", default= 0)
    parser.add_argument("-m", "--mapper", help="the mapper used to map, can be minimap2 or bwa", default="minimap2")

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])

    os.system(f'cat {args.reference}/*.fasta > allde.fasta')
    ref_bed=wrapper_run_get_bedfile('allde.fasta', args.file, args.core, args.qscore, args.mapper)
    chrname=bed_parser_get_higest_coverage(ref_bed)
    fa_dic=fasta2dic('allde.fasta')
    # name, seq=chr_select(fa_dic, chrname)
    seq = str(fa_dic[chrname].seq)
    print(">{name}\n{seq}\n".format(name=chrname, seq=seq))
