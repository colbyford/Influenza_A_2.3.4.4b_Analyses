#!/bin/bash
#SBATCH --job-name=hyphy
#SBATCH --output=hyphy.out
#SBATCH --error=hyphy.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=5-

module load hyphy

expect hyphy_expect.sh
