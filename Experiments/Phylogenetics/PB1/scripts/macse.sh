#!/bin/bash
#SBATCH --job-name=macse
#SBATCH -t 10-
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

## Variables
input=$1
output_NT="Alignments/$(basename $input .fasta).macse_NT.fasta"
output_AA="Alignments/$(basename $input .fasta).macse_AA.fasta"
mkdir -p Alignments

## Run MACSE
macse -prog alignSequences -seq $input -out_NT $output_NT -out_AA $output_AA