# Running scripts separately in HPC

* Parse Input - sbatch -J "parse" -t 1- --mem=8G --wrap "python3 scripts/parse_input.py --csv input.csv"

* AlphaFold -  sbatch -J "afrun" -p GPU --gres=gpu:1 --constraint=FP32 -t 3- --mem=32G --wrap "python3 scripts/run_alphafold.py"

* Ligand Prep - sbatch -J "ligand_prep" -p GPU --gres=gpu:1 --constraint=FP32 -t 3- --mem=32G --wrap "python3 scripts/prepare_ligand.py"

* Docking Prep - sbatch -J "docking_prep" -t 1- --mem=16G --wrap "python3 scripts/prep_docking_jobs.py"

* Docking Run - sbatch -J "docking" -t 1- --mem=8G --wrap "python3 scripts/run_haddock3.py"

* Results - sbatch -J "results" -t 1- --mem=8G --wrap "python3 scripts/collect_haddock_metrics.py"
