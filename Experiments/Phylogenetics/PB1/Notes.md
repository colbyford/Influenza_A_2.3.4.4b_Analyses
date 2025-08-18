# How to run the entire pipeline

## Collect Data
Both DNA sequences and metadata should be downloaded from databases such as NCBI and GISAID.
We queried for 'H5N1 clade 2.3.4.4b'

## Combine Sequences
Run:
`python3 combine_fasta.py <input_directory>` # This script will also remove any duplicate accession IDs

## Deduplication and N Removal
Run:
`python3 clean-up.py <input.fasta> <output.fasta>`

## Alignment
Run MAFFT:
`sbatch mafft.sh <input.fasta> <output.fasta>`

Manual adjustments need to be done in Aliview after alignment.
- Trim sequences to reading frame
- Remove sequences that create gaps
- Output should be in NEXUS format with gaps replaced as '?'

## Get Metadata
If you have collected sequences from NCBI, you can collect the metadata using Entrez in the script below:
`python3 get-metadata-ncbi.py --years-only <input.fasta> <output.csv>`

This will not collect metadata fro GISAID. That must be done manually.
`python3 get-metadata-gisaid.py --id-col <protein name> --output_csv <ouput.csv> <input.fasta> <metadata.xls>`

## Parsimony Tree
Modify the config file `tnt_config.sh`
Config file will need the following modified:
- input file path
- random seed (or leave as standard '*')
- maxram (depending on system availability)

Run TNT:
`sbatch slurm_tnt.sh`

## Visualize with Mesquite along with year metadata
Metadata should be a CSV of two columns: "Accession ID" and "Year (collection)"

- Upload tnt tree to Mesquite and make a new project
- Move root of tree to the bottom (can click options to change direction of root)
- Import Nexus file first
- Import TNT tree using "linked tree"
- Make new empty matrix = 1
- Edit state names and input unique states (years) 
- Save copy of Nexus
- Make sed script to replace year with state number (sed -f script.txt metadata.csv) # Make sure it is tab delimited
- Replace metadata columns in copy of Nexus file
- Reopen in Mesquite
- Analysis:Tree
- Trace Character History
- Parsimony Ancestral States
- Standard Categorical Data
- Remove clades that are older than 2021 - Select taxa:Tree:Alter/Transform Tree:cut selected taxa

## Prepare for RAxML
 Following steps are to prepare tree and alignment from mesquite for RAxML

 - Step 1
    - File > Export #While on tree window
    - Simple Newick/Phylip Treefile
    - Save as .nwk

- Step 2
    - Make copy of mesquite nexus file
    - Modify copy to remove all blocks except 'character-matrix'
    - Run following script on nexus copy:
    `python3 prep4raxml.py <input.nex> <input.nwk> <output.phy>`

## Run RAxML
Take Phylip file from mesquite:
`sbatch raxml.sh <input.phy>`

Remove pruned taxa from Metadata:
`python3 prune_metadata.py <input.phy> <metadata.csv> <output.csv>`

## Run HyPhy

