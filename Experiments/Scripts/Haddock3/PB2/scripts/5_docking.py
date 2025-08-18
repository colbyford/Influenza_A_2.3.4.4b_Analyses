#!/bin/python3

import pandas as pd
import os

excel_file = "Experiments_updated.xlsx"
df = pd.read_excel(excel_file)

for index, row in df.iterrows():
    cmd = f"sbatch scripts/run_haddock3.sh 2_Docking/{row['Combined']}/config.cfg"
    os.system(cmd)
