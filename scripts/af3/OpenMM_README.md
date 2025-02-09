# OpenMM Docker Image

## Pull Docker Image

```bash
docker pull cford38/openmm:cuda12.5.0
```

## Run Docker Container

```bash
docker run -it --gpus all -v .:/mnt --name openmm cford38/openmm:cuda12.5.0
```