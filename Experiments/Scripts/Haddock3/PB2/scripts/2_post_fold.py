#!/bin/python3

## Importing libraries
import pandas as pd
import glob
import os

## Read the experiments file
Protein1 = pd.read_excel('Experiments.xlsx', sheet_name='Protein')
Protein2 = pd.read_excel('Experiments.xlsx', sheet_name='Protein2')

# Create a new dataframe by pairing each sample in Protein1 with each protein in Protein2
combined_df = pd.merge(Protein1.assign(key=1), Protein2.assign(key=1), on='key').drop('key', axis=1)

# Combine the first column values of both dataframes
combined_df['Combined'] = combined_df.apply(lambda row: f"{row[0]}_{row[len(Protein1.columns) + 1]}", axis=1)

# Reorder columns to have 'Combined' as the first column
cols = combined_df.columns.tolist()
cols = [cols[-1]] + cols[:-1]
combined_df = combined_df[cols]

# Save the combined dataframe
combined_df.to_excel('Experiments_updated.xlsx', index=False)