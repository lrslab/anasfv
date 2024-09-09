#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
polish 
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
        polish_asfv.py -i single_fasta/MN194591.1.fasta -r single_fasta/OR180113.1.fasta -m R9.4.pkl
        '''

    parser = argparse.ArgumentParser(prog='polish_asfv.py',
                                     description='Run bash cmd lines for files',
                                     epilog=example_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    
    parser.add_argument("-i", "--input", help="fasta file of input ASFV genome")
    parser.add_argument("-r", "--ref", help="fasta file of reference ASFV genome")
    parser.add_argument("-m","--model", help="homopolish model",default='R9.4.pkl')

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    input = args.input
    ref = args.ref
    homopolish_model = args.model
    
    command = f'homopolish polish -a {input} -l {ref} -m {homopolish_model} -o {os.path.basename(input)}-homopolished'
    conda_env = 'homopolish'
    execute_command_in_conda_env(conda_env, command)
    

    
    
    
