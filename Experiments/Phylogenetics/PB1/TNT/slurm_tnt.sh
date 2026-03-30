#!/usr/bin/bash
#SBATCH --job-name=tnt
#SBATCH -t 3-
#SBATCH -n 16
#SBATCH --mem=64g

/users/sguirale/gTNT/tnt BGROUND proc tnt_config.txt
