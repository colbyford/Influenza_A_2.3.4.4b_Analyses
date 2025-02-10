# OpenMM Docker Image

## Pull Docker Image

```bash
docker pull cford38/openmm:cuda12.5.0
```

## Run Docker Container

```bash
docker run -it --gpus all -v .:/mnt --name openmm cford38/openmm:cuda12.5.0
```

## Run OpenMM Simulation

```bash
experiment_path="/mnt/data/experiments/af3/PQ591824__4K63_I/pq591824__4k63_i"

python /mnt/scripts/af3/run_openmm.py --input_pdb_path ${experiment_path}/pq591824__4k63_i_model_protein.pdb --input_sdf_path ${experiment_path}/pq591824__4k63_i_model_mol_withH.sdf --output_pdb_path ${experiment_path}/pq591824__4k63_i_model_protein_minimized.pdb

```