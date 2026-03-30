#!/bin/bash
#SBATCH --job-name="af3_h5"
#SBATCH --partition=GPU
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --gres=gpu:A100:4
#SBATCH --time=1:00:00

module load singularity

## Base Directories
export AF3_RESOURCES_DIR=/apps/pkg/singularity/alphafold/3.0.0
export AF3_HPC_DATA_DIR=/projects/datasets/alphafold/3.0.0

## Image and Code
export AF3_IMAGE=${AF3_RESOURCES_DIR}/alphafold3.sif
export AF3_CODE_DIR=${AF3_RESOURCES_DIR}/code

## Weights and Databases
export AF3_MODEL_PARAMETERS_DIR=${AF3_HPC_DATA_DIR}/weights
export AF3_DATABASES_DIR=${AF3_HPC_DATA_DIR}/databases


## Experiment Resources
export AF3_EXPERIMENT_DIR=/users/cford38/Influenza_A_2.3.4.4b_Analyses/data/experiments/af3
export AF3_INPUT_DIR=${AF3_EXPERIMENT_DIR}/${1}
export AF3_OUTPUT_DIR=${AF3_EXPERIMENT_DIR}/${1}

singularity exec \
     --nv \
     --bind $AF3_INPUT_DIR:/root/af_input \
     --bind $AF3_OUTPUT_DIR:/root/af_output \
     --bind $AF3_MODEL_PARAMETERS_DIR:/root/models \
     --bind $AF3_DATABASES_DIR:/root/public_databases \
     $AF3_IMAGE \
     python ${AF3_CODE_DIR}/run_alphafold.py \
     --json_path=/root/af_input/alphafold_input.json \
     --model_dir=/root/models \
     --db_dir=/root/public_databases \
     --output_dir=/root/af_output