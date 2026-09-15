# Xenium (RNA) + CODEX (protein) integration workflow

A Snakemake pipeline for spatial multi-omics on tissue slides. It:

1. concatenates the Xenium cell × gene tables of several slides,
2. measures per-cell CODEX protein intensities on the Xenium cell segmentation,
3. quality-controls both modalities and merges them per cell,
4. integrates RNA and protein with **totalVI** and clusters the cells (Leiden), with
   UMAPs, differential-expression tables and marker plots for annotation.

All computation runs **inside one Apptainer image** (`integration.sif`). Snakemake runs
outside it, in a small conda env, and submits **one SLURM job per step**.

```
 raw Xenium runs                                    registered CODEX tiffs
        │                                                │
    [step1a]  Xenium concatenation                       │
        │                                                │
        ├──▶ per-slide *.zarr (cell label images) ───────┤
        │                                                ▼
        │                                             [step1b]  per-cell CODEX intensities
        ▼                                                ▼
 xenium_concatenated.h5ad                   ID_<slide>_intensity.parquet
        └───────────────────────┬────────────────────────┘
                                ▼
                             [step2]  QC + merge RNA and protein
                                ▼
                         final_adata.h5ad
                                ▼
                             [step3]  totalVI integration + clustering (GPU)
                                ▼
          totalVI model, Leiden clusters, UMAPs, DE tables, marker plots
```

New to the project? [docs/restructure_explained.md](docs/restructure_explained.md)
explains in plain words how the repository is organised and why.

---

## Repository layout

| Path | What it is |
|------|------------|
| `run.sh` | the launcher, the only thing you run (see [Running](#running)) |
| `config/config.example.yaml` | template for **your run**: input data, results folder, parameters, job sizes. Copy it to `config/config.yaml`. |
| `config/README.md` | every config option explained |
| `profiles/slurm/config.yaml.example` | template for **your cluster**: SLURM partitions, account. Copy it to `profiles/slurm/config.yaml`. |
| `workflow/Snakefile` | the four rules and how each runs inside the container |
| `workflow/scripts/` | the four step scripts (`step1a_…` to `step3_…`) |
| `envs/` | `snakemake_env.yml` (the env you install) and the three pinned envs inside the image |
| `container/` | `Dockerfile` + `entrypoint.sh` to build the image; see `container/README.md` |
| `docs/` | background documentation |

`config/config.yaml`, `profiles/slurm/config.yaml`, results, logs, `.snakemake/` and the
image files are git-ignored.

---

## Requirements

- A **SLURM** cluster with **Apptainer** (Singularity), and a GPU partition for step 3.
  The container's torch is a CUDA 13 build.
- **micromamba** (or conda/mamba) to create the Snakemake env.
- The image **`integration.sif`**. Ask the maintainers for the existing one, or build it
  (see `container/README.md`).

---

## Quick start

```bash
# 1. The env that runs Snakemake (once)
micromamba create -n snakemake_env -f envs/snakemake_env.yml

# 2. Your run settings: input folders, results_dir, path to integration.sif, slides
cp config/config.example.yaml config/config.yaml
$EDITOR config/config.yaml

# 3. Your cluster settings: partition names, optional account
cp profiles/slurm/config.yaml.example profiles/slurm/config.yaml
$EDITOR profiles/slurm/config.yaml

# 4. Rehearse: prints the plan, runs nothing
./run.sh dry

# 5. Run: the orchestrator runs as a small SLURM job, so you can log out
./run.sh submit
squeue -u $USER          # watch the step jobs
```

Paths in the config can be absolute, or relative to the repository root.

---

## Input data

Point `paths.xenium_raw` and `paths.codex_base` in `config/config.yaml` to:

**Xenium:** a folder with one raw Xenium output folder per slide (each contains
`experiment.xenium`). The slide ID is read from the folder name: the first block of
digits between double underscores, e.g. `0056777` in
`output-XETG00404__0056777__Region_1__20250612__144008`. Standard Xenium folder names
already look like this. The pipeline expects **one run per slide**.

**CODEX:** a folder with one sub-folder per slide. Each sub-folder holds one image per
channel. The images must be **registered to the Xenium image**: same height and width
as the Xenium cell label image. Two naming rules apply:

- **Folder name starts with the slide ID**, optionally after `ID_`: e.g. `0056777_tifs`
  or `ID_0056777__Region_1_scale0_tif`. A sub-folder named differently
  (e.g. `slide_0056777`) is ignored without a warning, and step 1b then fails because
  that slide has no CODEX data.
- **File name** `…_<CHANNEL>_REGISTERED_….tif` or `.tiff`, e.g.
  `morphology_focus_0000.ome.tif_8_aSMA_REGISTERED_scale0.tiff`. The channel name is
  the part just before `_REGISTERED_`, and file names are split at underscores, so
  **channel names must not contain `_`**: `HLA_DR` would be read as `DR` (use `HLA-DR`).
  Files without `_REGISTERED_` in the name are ignored.

**Which slides:** list them under `slides:` in the config.

- Slide IDs must be **digits only** (the scripts recognise only digits), written in
  quotes (see `config/README.md`).
- Steps 1a and 1b process only the listed slides.
- Step 1b **stops with an error** if `codex_base` contains a CODEX folder for a slide
  that has no Xenium data in this run. Keep only the folders of the listed slides there.
- Step 2 reads **every** `ID_<slide>_intensity.parquet` in the step 1b output folder.
  Use a fresh `results_dir` when you change the slide list.

---

## Running

| Command | What it does |
|---------|--------------|
| `./run.sh dry` | Dry run: prints the plan (steps, files, commands, resources); runs nothing |
| `./run.sh submit` | Runs the pipeline; the orchestrator is a small SLURM job (48 h limit) |
| `./run.sh run` | Same, but the orchestrator runs in your shell. Keep it open, e.g. in `tmux`. |
| `./run.sh unlock` | Removes a stale lock after a killed run |
| `./run.sh dag` | Writes `dag.svg` (or `dag.dot`) of the step graph |

- **Another config:** `./run.sh dry --configfile config/my_other_run.yaml`. Only that
  file is read.
- **Other Snakemake options** are passed on, e.g. `./run.sh dry --forceall`.
- **Finding Snakemake:** `run.sh` uses `$SNAKEMAKE` if set, else `snakemake` on the
  `PATH`, else `micromamba run -n snakemake_env snakemake`.
- **Partition or account for the orchestrator job:** set `SBATCH_PARTITION` /
  `SBATCH_ACCOUNT`.

**Re-running.** `run.sh` uses `--rerun-triggers mtime`: a step is redone only if one of
its outputs is missing or an input file is newer than its outputs. This protects
finished results. It also means that **changing a parameter in the config does not by
itself redo anything.** After changing, for example, a `qc:` value, force that step and
everything after it:

```bash
./run.sh submit --forcerun step2_data_qc      # redoes step 2 and then step 3
```

Rule names: `step1a_xenium_concatenation`, `step1b_codex_intensities`, `step2_data_qc`,
`step3_data_integration`.

> Running plain `snakemake` **without** `--rerun-triggers mtime` can decide to redo
> finished steps after a config change. Use `run.sh`.

---

## Outputs

Everything goes to `paths.results_dir`:

| Folder | Contents |
|--------|----------|
| `step1a/` | `xenium_concatenated.h5ad` (raw counts in `.X` and `layers["counts"]`, index `<cell_id>_<slide_ID>`), `zarr_dir/<run>.zarr` per slide |
| `step1b/parquet_dir/` | `ID_<slide>_intensity.parquet`: mean intensity per cell (`cell_uid`) and channel |
| `step2/` | `final_adata.h5ad` (RNA + protein in `obsm["protein_expression"]`), `filtered_adata.h5ad`, the CODEX tables before/after filtering, `tracking_counts.csv` (cells kept at each QC stage), `cell_id_x_group/` (per slide: kept / removed by Xenium QC / removed by CODEX QC, for Xenium Explorer) |
| `step3/` | `model.pt` (trained totalVI), `mdata_*.h5mu` (latent, denoised, Leiden, dendrogram), `elbo_training.png`, UMAPs (`umap_clusters*.png`, `forprob.png`, `denoised_protein.png`), `differential_expression.csv`, `DE_table_*.csv`, `dotplot_*.png`, `matrixplot_*.png` |
| `logs/` | one log per step (the script's output), plus `orchestrator_<jobid>.log` from `./run.sh submit` |

---

## Troubleshooting

- **Step 3 fails with a CUDA/driver error.** The node's GPU driver is too old for the
  CUDA 13 torch. In `profiles/slurm/config.yaml` set the step 3 partition to one with
  CUDA 13 drivers (on ETH Euler: `cuda13pr.24h`).
- **sbatch asks for an account.** Uncomment `slurm_account` in
  `profiles/slurm/config.yaml`; set `SBATCH_ACCOUNT` for `./run.sh submit`.
- **A job ran out of memory or time.** Raise `mem_mb_per_cpu`, `runtime` (minutes) or
  `threads` for that step under `resources:` in `config/config.yaml`. Total memory =
  `mem_mb_per_cpu × threads`.
- **"LockException" after a killed run.** Run `./run.sh unlock`.
- **A script can't find a file that exists.** The container only sees the input folders,
  `results_dir`, the scripts and the config folder. If your data contains symlinks to
  other places, add those folders under `container.extra_binds`.
- **Where is the error?** Look in `<results_dir>/logs/<step>.log`. SLURM's own job logs
  are kept only for failed jobs, under `.snakemake/slurm_logs/`.
