#!/usr/bin/python3

import pandas as pd
import os

## Read in experiments data
experiments = pd.read_excel('Experiments_updated.xlsx')

## Gather Docking Metrics
# Define the columns to extract from capri_ss.tsv
columns_to_extract = ['model', 'score', 'bsa', 'desolv', 'elec', 'total', 'vdw']

# Initialize a list to store the extracted data
extracted_data = []

# Iterate over each row in the experiments DataFrame
for index, row in experiments.iterrows():
    # Get the value from the first column
    value = row.iloc[0]
    
    # Construct the path to the capri_ss.tsv file
    capri_file_path = f'2_Docking/{value}/output/11_caprieval/capri_ss.tsv'
    
    # Check if the file exists
    if os.path.exists(capri_file_path):
        # Read the capri_ss.tsv file
        capri_df = pd.read_csv(capri_file_path, sep='\t')
        
        # Find the row with a 1 in the caprieval_rank column
        capri_row = capri_df[capri_df['caprieval_rank'] == 1]
        
        # If such a row exists, extract the required columns
        if not capri_row.empty:
            extracted_row = capri_row[columns_to_extract].copy()

            # Replace '..' with '3_Docking/{value}/output' in the 'model' column
            extracted_row['model'] = extracted_row['model'].str.replace('..', f'2_Docking/{value}/output')
            
            # Append the extracted row to the list
            extracted_data.append(extracted_row.values.tolist()[0])
        else:
            extracted_data.append([None] * len(columns_to_extract))
    else:
        extracted_data.append([None] * len(columns_to_extract))

# Convert the extracted data to a DataFrame
extracted_df = pd.DataFrame(extracted_data, columns=columns_to_extract)

# Concatenate the extracted data with the original experiments DataFrame
experiments = pd.concat([experiments, extracted_df], axis=1)

## Run PRODIGY
for index, row in experiments.iterrows():
    # Skip rows where the model column is empty
    if pd.isna(row['model']):
        continue

    ## Create full path to best PDB file
    pdb_path = str(row['model'])

    ## Uncompress the PDB file
    # Create the 'Top Models' directory if it doesn't exist
    os.makedirs('Top_Models', exist_ok=True)

    # Get the value from the 'combinations' column
    combination_value = row['Combined']

    # Construct the new path in the 'Top Models' directory
    new_pdb_path = f'Top_Models/{combination_value}.pdb.gz'

    # Copy the PDB file to the new location
    os.system(f'cp {pdb_path}.gz {new_pdb_path}')

    # Uncompress the new PDB file
    uncompressed_pdb_path = new_pdb_path.replace('.gz', '')
    os.system(f'gunzip {new_pdb_path}')
    
    ## Run PRODIGY and parse stdout
    prodigy_output = os.popen(f'prodigy {uncompressed_pdb_path}').read()
    prodigy_output_lines = prodigy_output.split('\n')
    predicted_binding_affinity = float(prodigy_output_lines[-3].split(':')[1].replace(' ', ''))
    predicted_dissociation_constant = float(prodigy_output_lines[-2].split(':')[1].replace(' ', ''))

    experiments.loc[index, 'prodigy_deltaG_kcalpermol'] = predicted_binding_affinity
    experiments.loc[index, 'prodigy_dg_score'] = predicted_dissociation_constant

# Save the updated DataFrame to a new CSV file
experiments.to_excel('Experiments_final.xlsx', index=False)