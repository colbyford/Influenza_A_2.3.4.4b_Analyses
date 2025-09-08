## Amber relaxation with OpenMM
# python .\amber_relax.py . .\minimized

## Orbital MD simulation

docker run -v .\\md:/workspace --gpus all -it --rm --name orb_models -it cford38/orb_models:latest /bin/bash


python OrbMD.py --input_file /workspace/data/PQ809560_pocket2__baloxavir_out.xyz --output_file /workspace/data/PQ809560_pocket2__baloxavir_out.xyz