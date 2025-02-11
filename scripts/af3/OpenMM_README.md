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

- Update CUDA Toolkit: `conda install -c conda-forge cudatoolkit=11.6`

```bash
experiment_path="/mnt/data/experiments/af3/EPI1846961__4K63_I/epi1846961__4k63_i/"

## Note that the PDB file has both the ligand (with hydrogens) and the protein (without hydrogens)
python /mnt/scripts/af3/run_openmm.py --input_pdb_path ${experiment_path}epi1846961__4k63_i_model_complex_molwithH.pdb --input_sdf_path ${experiment_path}epi1846961__4k63_i_model_mol_withH.sdf --output_pdb_path ${experiment_path}epi1846961__4k63_i_model_complex_minimized.pdb

```