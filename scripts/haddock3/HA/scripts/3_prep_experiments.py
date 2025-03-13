#!/bin/python3

## Importing libraries
import pandas as pd
import glob
import os

## Read the experiments file
experiments = pd.read_excel('Experiments.xlsx', sheet_name='Viral Proteins')

# Find the {accession}*rank_001*.pdb file path and add it as a new column
def find_pdb_file(accession):
    files = glob.glob(f'1_Folding/{accession}/{accession}*rank_001*.pdb')
    return files[0] if files else None

experiments['protein_pdb'] = experiments['accession'].apply(find_pdb_file)

# Get list of PDB files
pdb_files = [f for f in os.listdir('input') if f.endswith(".pdb")]

# Create an expanded list of pairwise combinations
data = []
for _, row in experiments.iterrows():
    accession = row['accession']
    aa_sequence = row['aa_sequence']
    protein_pdb = row['protein_pdb']
    for pdb_file in pdb_files:
        pdb_basename = os.path.splitext(pdb_file)[0]  # Remove .pdb extension
        pdb_path = os.path.join('input', pdb_file)
        accession_pdb_combination = f"{accession}_{pdb_basename}"
        
        data.append([accession_pdb_combination, accession, pdb_basename, aa_sequence, protein_pdb, pdb_path])

# Create a new DataFrame
output_df = pd.DataFrame(data, columns=["combination", "accession", "glycan", "aa_sequence", "protein_pdb", "glycan_pdb"])

# Save the output to an Excel file
output_df.to_excel('Experiments_updated.xlsx', index=False)
