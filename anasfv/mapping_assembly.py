#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
mapping assembly and polish
"""

import subprocess
import os

def execute_command_in_conda_env(conda_env, command):
    full_command = f'conda run -n {conda_env} {command}'
    subprocess.run(full_command, shell=True)

if __name__=="__main__":
    import sys
    import argparse

    example_text = '''example:
        mapping_assembly.py -p 4 -r single_fasta -i test_data.fasta -o genome.fasta --medaka r941_min_high_g303
        '''

    parser = argparse.ArgumentParser(prog='mapping_assembly.py',
                                     description='Run bash cmd lines for files',
                                     epilog=example_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    
    parser.add_argument("-p", "--processes", help="number of processes", default= 4)
    parser.add_argument("-i", "--input", help="input fasta or fastq file")
    parser.add_argument("-r", "--ref", help="a folder containing multiple ASFV genomes, or a single reference sequence file")
    parser.add_argument("-o", "--output", help="assembled ASFV genome")
    parser.add_argument("--medaka", help="medaka model",default=None)
    # parser.add_argument("--homopolish", help="homopolish model",default=None)

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    num_processes = args.processes
    input_reads = args.input
    ref = args.ref
    output = args.output
    medaka_model = args.medaka
    # homopolish_model = args.homopolish
    
    if os.path.isdir(ref):
        subprocess.run(f'find_near_ref.py -r {ref} -f {input_reads} -c {num_processes} > near.fasta', shell=True)
        ref = 'near.fasta'
        subprocess.run('rm allde.bed allde.fasta maps.bam maps.bam.bai', shell=True)
    elif os.path.isfile(ref):
        pass
    else:
        raise ValueError(f"{ref} is not a valid directory or file")
    
    command = f'minimap2 -a {ref} {input_reads} | samtools view -b -F 4 | samtools sort -@ {num_processes} -o - | samtools consensus -f fasta - -o {output}'
    subprocess.run(command, shell=True)
    
    strain = output
    if strain.endswith('.fasta'):
        strain = strain[:-6]
    elif strain.endswith('.fa'):
        strain = strain[:-3]
    subprocess.run(f"sed -i '1s/.*/>{strain}/' {output}", shell=True)

    if medaka_model:
        command = f'medaka_consensus -i {input_reads} -d {output} -o medaka_result -m {medaka_model} -t {num_processes} > medaka.log'
        conda_env = 'medaka'
        execute_command_in_conda_env(conda_env, command)
        subprocess.run(f'cp ./medaka_result/consensus.fasta {output}', shell=True)
        subprocess.run(f'rm medaka.log {output}.fai {output}.map-ont.mmi', shell=True)

    # if homopolish_model:
    #     command = f'homopolish polish -a {output} -l {ref} -m {homopolish_model} -o homopolish-output'
    #     conda_env = 'homopolish'
    #     execute_command_in_conda_env(conda_env, command)
    #     subprocess.run(f'cp ./homopolish-output/{strain.split('.')[0]}_homopolished.fasta {output}', shell=True)
    
    
    
    
