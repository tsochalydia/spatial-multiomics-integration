# The container image (`integration.sif`)

Every step of the pipeline runs inside one Apptainer image, `integration.sif`. Snakemake
itself runs outside it, in `snakemake_env`. For each step it calls a Python interpreter
inside the image:

```
apptainer exec [--nv] -B <folders> integration.sif /opt/conda/envs/<env>/bin/python workflow/scripts/<step>.py …
```

`workflow/Snakefile` builds this command. The folders to bind are derived from the
config, and `--nv` (GPU access) is used for step 3 only.

---

## What is inside

A Debian-based `mambaorg/micromamba:2.8` image with three conda environments. Each is
created from a pinned environment file in `envs/`:

| Env in the image | Environment file | Used by | Key packages (versions in the current image) |
|------------------|------------------|---------|----------------------------------------------|
| `spatial` | `envs/spatial_env.yml` | step 1a | Python 3.11.15, spatialdata 0.7.3, spatialdata-io 0.7.1, anndata 0.12.19 |
| `int_retrieval_env` | `envs/int_retrieval_env.yml` | steps 1b, 2 | Python 3.14.3, spatialdata 0.7.2, tifffile 2026.3.3, anndata 0.12.10, scanpy 1.12 |
| `scvi_env` | `envs/scvi_env.yml` | step 3 (GPU) | Python 3.13.13, scvi-tools 1.4.2, torch 2.11.0+cu130 (CUDA 13.0), mudata 0.3.4, anndata 0.12.10 |

`envs/snakemake_env.yml` is **not** in the image; it's the env you install on the login
node to run Snakemake.

- **No env is activated on the `PATH`.** The workflow calls each interpreter by its full
  path (`/opt/conda/envs/<env>/bin/python`). These paths are listed in
  `workflow/Snakefile` (`PY = {...}`).
- **Readable by any user:** `/opt/conda` is readable by everyone, so the image works
  under any cluster user ID.
- **The entrypoint is only for `docker run`.** `container/entrypoint.sh` makes sure
  `HOME` is writable and loads micromamba's shell hook. `apptainer exec`, which is how
  the workflow runs, does not use it.
- **Caches go to the job's temporary folder.** The workflow points the numba,
  matplotlib and cache folders to the job's `$TMPDIR`, so a read-only home directory
  can't break scanpy.

---

## Current image

`apptainer inspect integration.sif` shows:

- **built:** 14 July 2026, with Singularity 4.5.0
- **built from:** a Docker archive (`bootstrap: docker-archive`, from `integration.tar`),
  i.e. an image built with Docker, then saved to a `.tar` file and converted
- **size:** 13 GB

---

## Rebuilding the image

Only needed if you can't get the existing `integration.sif`. This is how the current
image was made.

1. **On a computer with Docker and internet access**, from the
   repository root:
   ```bash
   docker build -f container/Dockerfile -t integration .
   docker save integration:latest -o integration.tar
   ```
   The `.dockerignore` makes Docker read only `envs/` and `container/`, not the large
   image files or results.
2. **Copy `integration.tar` to Euler and convert it.** Euler has Apptainer (the new name
   of Singularity; `singularity build` works the same way):
   ```bash
   apptainer build integration.sif docker-archive://integration.tar
   ```
3. Set `container.sif` in `config/config.yaml` to the new file, and check the points below.

### Before relying on a rebuilt image, check

- **Versions:** the key package versions match the table above (for example,
  `apptainer exec integration.sif /opt/conda/envs/scvi_env/bin/python -c "import torch; print(torch.__version__)"`).
- **torch CUDA build:** `envs/scvi_env.yml` pins `torch==2.11.0` via pip but does not name
  a PyPI index. The current image has the **CUDA 13** build (`+cu130`); check that a
  rebuild installs the same.
- **Duplicate pins in `envs/scvi_env.yml`:** `anndata` and `mudata` appear twice, once
  from conda (anndata 0.12.14, mudata 0.3.8) and once from pip (anndata 0.12.10, mudata
  0.3.4). The current image has the **pip** versions. The file also lists NVIDIA
  libraries for both CUDA 12 and CUDA 13.

---

## GPU drivers

torch in `scvi_env` is built for CUDA 13, so step 3 needs a GPU node whose NVIDIA driver
supports CUDA 13. If step 3 fails with a CUDA or driver error, set its partition in
`profiles/slurm/config.yaml` to one with newer drivers.
