# Xenium (RNA) + CODEX (protein) integration workflow

Snakemake pipeline that concatenates Xenium transcriptomics across TMA slides,
retrieves per-cell CODEX protein intensities, QCs and merges the two modalities,
and integrates them with **totalVI** to produce Leiden clusters for annotation.

All compute runs **inside `integration.sif`** via `apptainer exec`. Snakemake
itself runs from the tiny `snakemake_env` and submits **one SLURM job per step**.

```
 raw Xenium runs                                    CODEX tiffs
        │                                                │
    [step1a]  Xenium concatenation                       │
        │                                                │
        ├──▶ per-slide *.zarr (cell label images) ───────┤
        │                                                ▼
        │                                             [step1b]  per-cell CODEX intensities
        │                                                │
        ▼                                                ▼
 xenium_concatenated.h5ad                   ID_<slide>_intensity.parquet
        │                                                │
        └───────────────────────┬────────────────────────┘
                                ▼
                             [step2]  QC + merge RNA and protein
                                │
                                ▼
                         final_adata.h5ad
                                │
                                ▼
                             [step3]  totalVI integration (GPU)
                                │
                                ▼
          totalVI model, Leiden clusters, UMAPs, DE table
```

## Layout

| File | Purpose |
|------|---------|
| `Snakefile` | the 4 rules + `apptainer exec` wiring |
| `config.yaml` | **all** paths, parameters, SLURM resources — edit this, not the scripts |
| `scripts/step1a…step3…py` | the four step scripts (argparse + config, logic unchanged from the tested versions) |
| `profiles/slurm/config.yaml` | SLURM executor profile (one sbatch per rule) |
| `run.sh` | launcher (`dry` / `dag` / `run` / `unlock`) |
| `logs/` | per-step python stdout/stderr |

## Container env → step mapping (inside the .sif)

| Step | conda env in .sif | why |
|------|-------------------|-----|
| 1a | `spatial` | needs `spatialdata_io` |
| 1b | `int_retrieval_env` | `spatialdata` + `tifffile` |
| 2  | `int_retrieval_env` | `anndata` + `scanpy` |
| 3  | `scvi_env` (GPU, `--nv`) | `scvi` + `torch` cu130 |

## Run it

```bash
cd /cluster/project/moor/lydia/snakemake_workflow

# 1. Dry run — prints the plan, submits nothing, runs nothing:
./run.sh dry

# 2. Real run — submits one job per step, in order. Keep this process alive
#    (use tmux/screen) while it orchestrates:
tmux new -s smk
./run.sh run
```

Re-running only redoes steps whose inputs changed. To force a clean redo of a
step, delete its outputs (e.g. `rm .../snakemake_testing_qc/final_adata.h5ad`)
and run again. If a run is killed mid-job, `./run.sh unlock` clears the lock.

## Configuring

- **Which slides:** list the slide IDs under `slides:` in `config.yaml`. A slide ID
  is the number in the Xenium run folder name, e.g. `0056777` in
  `output-XETG00404__0056777__...`. The list does not apply to every step:
  - **Step 1a ignores it** and processes every Xenium run found under
    `xenium_raw`. To process only some slides, point `xenium_raw` to a folder
    that contains only those runs.
  - **Step 1b** only computes CODEX intensities for the listed slides.
  - **Step 2** keeps only cells that have CODEX intensities, so unlisted slides
    are dropped here. It reads **every** `ID_<slide>_intensity.parquet` in
    `parquet_dir`, so remove parquets left over from earlier runs with other slides.
- **QC / totalVI parameters:** the `qc:` and `integration:` sections.
- **SLURM resources:** the `resources:` and `slurm:` sections (per-step mem_mb,
  runtime in minutes, partition).
- **GPU step:** runs on `gpu.24h` with a 24 h limit. The container's torch is a
  CUDA-13 build; if a `gpu.24h` node has an older driver and torch complains
  about the CUDA version, set `gpu_partition: cuda13pr.24h` in `config.yaml`.
- **Account:** default association is used. If sbatch on your allocation demands
  `-A`, uncomment `slurm_account:` in `profiles/slurm/config.yaml` (e.g. `es_anmoor`).
