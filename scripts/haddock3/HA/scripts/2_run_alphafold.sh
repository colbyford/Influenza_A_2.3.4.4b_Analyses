#!/bin/bash
#SBATCH --job-name="afrun"
#SBATCH --partition=GPU
#SBATCH --gres=gpu:1
#SBATCH --constraint=FP32
#SBATCH --time=10-
#SBATCH --mem=32G

module load singularity
module load cuda

base=$(basename "$1" .fasta)

singularity run --nv \
  -B /users/sguirale/alphafold2/alphafold_weights:/cache -B $(pwd):/work \
  /users/sguirale/alphafold2/colabfold_1.5.5-cuda11.8.0.sif \
  colabfold_batch $1 1_Folding/${base}
