#!/bin/python3

import pandas as pd
import os
import subprocess

protein1 = pd.read_excel('Experiments.xlsx', sheet_name='Protein')
# protein2 = pd.read_excel('Experiments.xlsx', sheet_name='Protein2')

for index, row in protein1.iterrows():
    accession = row['Sample']
    sequence = row['Sequence']
    with open(f"input/{accession}.fasta", "w") as fasta_file:
        fasta_file.write(f">{accession}\n{sequence}\n")

for index, row in protein2.iterrows():
    accession = row['Protein']
    sequence = row['Sequence']
    with open(f"input/{accession}.fasta", "w") as fasta_file:
        fasta_file.write(f">{accession}\n{sequence}\n")

input_folder = '/scratch/sguirale/Influenza_SA/PB2/input/'
for fasta_file in os.listdir(input_folder):
    if fasta_file.endswith('.fasta'):
        fasta_path = os.path.join(input_folder, fasta_file)
        accession = os.path.splitext(fasta_file)[0]
        if accession in protein1['Sample'].values or accession in protein2['Protein'].values:
            subprocess.run(['sbatch', 'scripts/alphafold.sh', fasta_path])
