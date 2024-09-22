#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
using udance to make a tree 
"""

import subprocess
import os
import yaml


def execute_command_in_conda_env(conda_env, command):
    full_command = f'conda run -n {conda_env} {command}'
    subprocess.run(full_command, shell=True)


def get_abs_path(path, current_path):
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(current_path, path)

if __name__=="__main__":
    import sys
    import argparse

    example_text = '''example:
        make_tree.py -p 4 -f single_fasta -o tree --udance ./uDance
        '''

    parser = argparse.ArgumentParser(prog='make_tree.py',
                                     description='Run bash cmd lines for files',
                                     epilog=example_text,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    
    
    parser.add_argument("-p", "--processes", help="number of processes", type=int, default= 4)
    parser.add_argument("-f", "--file", help="dir of the fasta files")
    parser.add_argument("-o", "--output", help="output dir")
    parser.add_argument("--udance", help="path to udance dir")
    parser.add_argument("--iteration", help="Perform an iterative tree construction", action="store_true")

    args = parser.parse_args(args=None if sys.argv[1:] else ['--help'])
    num_processes = args.processes
    input_folder = args.file
    output_folder = args.output
    udance_folder = args.udance
    

    current_path = os.getcwd()
    output_folder = get_abs_path(output_folder, current_path)
    input_folder = get_abs_path(input_folder, current_path)
    udance_folder = get_abs_path(udance_folder, current_path)

    try:
        os.makedirs(output_folder)
    except FileExistsError:
        pass
    
    os.chdir(output_folder)
    subprocess.run(f'get_cds_alignments.py -f {input_folder}', shell=True)
    
    n_taxa = len([file for file in os.listdir(input_folder) if file.endswith('.fasta')])

    os.chdir(udance_folder)
    with open('config.yaml', 'r') as file:
        data = yaml.safe_load(file)
        data['workdir'] = output_folder
        data['chartype'] = 'nuc'
        data['backbone'] = 'de-novo'
        data['resources']['cores'] = num_processes
        data['mainlines_config']['n'] = n_taxa
        data['mainlines_config']['length'] = 50000
    with open('config1.yaml', 'w') as file:
        yaml.safe_dump(data, file)

    command = f'snakemake -c {num_processes} --configfile config1.yaml --snakefile udance.smk all'
    execute_command_in_conda_env('udance', command)
    
    if args.iteration:
        # round2
        subprocess.run(f'cp {output_folder}/output/udance.maxqs.nwk {output_folder}/backbone.nwk', shell=True)
        subprocess.run(f'rm {output_folder}/udance.log', shell=True)
        subprocess.run(f'rm -rf {output_folder}/output', shell=True)
        
        data['backbone'] = 'tree'
        with open('config2.yaml', 'w') as file:
            yaml.safe_dump(data, file)
        command = f'snakemake -c {num_processes} --configfile config2.yaml --snakefile udance.smk all'
        execute_command_in_conda_env('udance', command)
    
    subprocess.run(f'mv {output_folder}/output {output_folder}/udance_output', shell=True)
    subprocess.run(f'cp {output_folder}/udance_output/udance.maxqs.nwk {output_folder}/tree.nwk', shell=True)
            
    rm_command='rm -rf all_cds.fa ref_cds single_process udance.log'
    subprocess.run(rm_command, shell=True)
