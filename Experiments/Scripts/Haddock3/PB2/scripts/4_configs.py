#!/bin/python3

import configparser
import pandas as pd
import os
import glob

## Define functions for creating ambiguous AIR files
def write_ambig_air_file(active1, passive1, active2, passive2, segid1='A', segid2='B', output_file="ambig.tbl"):
    with open(output_file, "w") as output_file:
        ## Convert residues to integers
        all1 = active1 + passive1
        all2 = active2 + passive2

        ## Write lines from the active1 list
        for resi1 in active1:
            output_file.write('assign (resi {:d} and segid {:s})'.format(resi1, segid1) + '\n')
            output_file.write('(\n')
            c = 0
            for resi2 in all2:
                output_file.write('       (resi {:d} and segid {:s})'.format(resi2, segid2) + '\n')
                c += 1
                if c != len(all2):
                    output_file.write('        or\n')
            output_file.write(') 2.0 2.0 0.0\n\n')

        ## Write lines from the active2 list
        for resi2 in active2:
            output_file.write('assign (resi {:d} and segid {:s})'.format(resi2, segid2) + '\n')
            output_file.write('(\n')
            c = 0
            for resi1 in all1:
                output_file.write('       (resi {:d} and segid {:s})'.format(resi1, segid1) + '\n')
                c += 1
                if c != len(all1):
                    output_file.write('        or\n')
            output_file.write(') 2.0 2.0 0.0\n\n')



def create_config(protein1, protein2, output):
    
    config = configparser.ConfigParser()

    ## Read the configuration file
    config.read('protein_protein_config.cfg')

    ## Update the configuration
    config['main'] = {'run_dir': f'"{output}output"',
                        'mode': '"local"',
                        'ncores': 16,
                        'molecules': [
                            protein1,
                            protein2
                            ], 
                        'concat': 5,
                        'queue_limit': 100}
    
    config['rigidbody'] = {'tolerance': 5,
                            'ambig_fname': f'"{output}ambig.tbl"',
                            'sampling': 100}
    
    config['flexref'] = {'tolerance': 5,
                            'ambig_fname': f'"{output}ambig.tbl"'}
    
    config['emref'] = {'tolerance': 5,
                        'ambig_fname': f'"{output}ambig.tbl"'}


    ## Write the configuration to a file
    output_file = f"{output}config.cfg"
    with open(output_file, 'w') as configfile:
        config.write(configfile)

    ## Replace specific lines in config file (HACKY FIX)
    with open(output_file, 'r') as configfile:
      cfgdata = configfile.read()
    cfgdata = cfgdata.replace('[main]', '## Protein-Glycan Docking with HADDOCK3') \
                    .replace('[caprieval_0]', '[caprieval]') \
                    .replace('[caprieval_1]', '[caprieval]') \
                    .replace('[caprieval_2]', '[caprieval]') \
                    .replace('[caprieval_3]', '[caprieval]') \
                    .replace('[caprieval_4]', '[caprieval]')
                                                                                                         
    ## Write the file out again
    with open(output_file, 'w') as configfile:
      configfile.write(cfgdata)


## Read the experiments file
experiments = pd.read_excel("Experiments_updated.xlsx")
for index, experiment in experiments.iterrows():
    experiment_id = experiment['Combined']
    print(f"Preparing experiment: {experiment_id}")

    protein1 = experiment['Sample_x']
    protein2 = experiment['Protein_y']

    ## Make experiment folders
    print(f"\tMaking experiment folders...")
    output_path = "2_Docking/" + experiment_id + "/"
    os.makedirs(output_path, exist_ok=True)

    ## Generate AIR file for each experiment using active and passive residues
    start = experiment['Start']
    end = experiment['End']
    active1 = list(range(start, min(start + 140, end + 1)))
    passive1 = []
    active2 = list(range(1001, 1100))
    passive2 = []
    print(f"\tGenerating AIR files...")
    write_ambig_air_file(active1, passive1,
                        active2, passive2,
                        segid1='A', segid2='B',
                        output_file = f"{output_path}ambig.tbl")

    ## Generate config file for each experiment
    print(f"\tGenerating config file...")
    protein1_path = glob.glob(f"1_Folding/{protein1}/{protein1}*_rank_001_*.pdb")
    if protein1_path:
        protein1_path = protein1_path[0]
    else:
        raise FileNotFoundError(f"No file found for pattern 1_Folding/{protein1}/{protein1}*_rank_001_*.pdb")
    create_config(
            protein1 = protein1_path,
            protein2 = f"input/{protein2}_modified.pdb",
            output = output_path)

    print(f'\tDone preparing the experiment files for {experiment_id}!')