#!/bin/python3

import pandas as pd
import os

## Number of ligands to use
# if len(sys.argv) != 3:
#     print("Usage: python3 5_haddock3.py <start_number> <end_number>")
#     sys.exit(1)

# start_number = int(sys.argv[1])
# end_number = int(sys.argv[2])

## Import ZINC smiles dataframe
excel_file = "Experiments_updated.xlsx"
df = pd.read_excel(excel_file)
# selected_df = df.iloc[start_number:end_number + 1]

for index, row in df.iterrows():
    cmd = f"sbatch scripts/run_haddock3.sh 2_Docking/{row['combination']}/config.cfg"
    os.system(cmd)