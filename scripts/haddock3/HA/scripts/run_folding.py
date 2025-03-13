#!/bin/python3

import pandas as pd
import os
import subprocess

experiments = pd.read_excel('Experiments.xlsx', sheet_name='Viral Proteins')

for index, row in experiments.iterrows():
    accession = row['accession']
    sequence = row['aa_sequence']
    with open(f"input/{accession}.fasta", "w") as fasta_file:
        fasta_file.write(f">{accession}\n{sequence}\n")

input_folder = '/scratch/sguirale/Influenza_SA/input/'
for fasta_file in os.listdir(input_folder):
    if fasta_file.endswith('.fasta'):
        fasta_path = os.path.join(input_folder, fasta_file)
        subprocess.run(['sbatch', 'scripts/2_run_alphafold.sh', fasta_path])
