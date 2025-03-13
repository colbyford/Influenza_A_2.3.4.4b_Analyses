#!/bin/python3

import os
import pymol
from pymol import cmd
import pandas as pd
import re

Protein2 = pd.read_excel('Experiments_updated.xlsx')
input_dir = '/scratch/sguirale/Influenza_SA/PB2/1_Folding'

def process_pdb(file_path, out_path):
    cmd.load(file_path, 'protein')
    cmd.alter('chain A', 'chain="B"')
    cmd.alter('all', 'resi=str(int(resi)+1000)')
    cmd.save(out_path)
    cmd.delete('all')

pymol.finish_launching(['pymol', '-cq'])

for protein in Protein2['Protein_y']:
    protein_dir = os.path.join(input_dir, protein)
    if os.path.isdir(protein_dir):
        for filename in os.listdir(protein_dir):
            if re.match(rf'{protein}.*_rank_001_.*\.pdb', os.path.join(protein_dir, filename)):
                input_file = os.path.join(protein_dir, filename)
                output_file = f'input/{protein}_modified.pdb'
                process_pdb(input_file, output_file)
