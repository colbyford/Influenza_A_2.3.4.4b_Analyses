#!/bin/bash
#SBATCH --time=3-
#SBATCH --mem=32G
#SBATCH -n 4
#SBATCH -c 16

module load singularity

## Get config file for experiment
config=$1

## Run HADDOCK
singularity exec references/singularity_images/haddock3.sif haddock3 $config