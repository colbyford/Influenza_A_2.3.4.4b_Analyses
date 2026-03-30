#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH -t 1-
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G

## Modules
module load mafft

## Variables
input=$1
output="Alignments/$(basename $input .fasta).mafft.fasta"
mkdir -p Alignments

## Run MAFFT
mafft --auto --codon --reorder --thread 8 $input > $output
