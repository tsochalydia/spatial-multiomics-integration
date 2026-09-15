# =============================================================================
# Xenium (RNA) + CODEX (protein) integration workflow
#
#   step1a  raw Xenium runs        -> per-slide *.zarr  +  concatenated .h5ad
#   step1b  zarr + CODEX tiffs     -> one intensity parquet per slide
#   step2   concat .h5ad + parquets-> final_adata.h5ad (both modalities merged)
#   step3   final_adata.h5ad       -> totalVI model, Leiden clusters, plots  [GPU]
#
# All compute runs INSIDE integration.sif via `apptainer exec`; Snakemake only
# orchestrates and submits one SLURM job per rule. Parameters live in
# config/config.yaml (or the file given with --configfile).
# =============================================================================

import os

# config/config.yaml is read only when no --configfile is given, so a different
# config (e.g. config/config.test.yaml) never inherits settings from it.
DEFAULT_CONFIG = "config/config.yaml"
if not workflow.config_settings.configfiles:
    configfile: DEFAULT_CONFIG

# Paths in the config may be absolute or relative to the directory Snakemake is
# run from (run.sh runs it from the repo root). They are made absolute here so
# the scripts and apptainer never depend on the current directory.
def abspath(path):
    return os.path.abspath(os.path.expanduser(path))

P      = {key: abspath(value) for key, value in config["paths"].items()}
SIF    = abspath(config["container"]["sif"])
SLURM  = config["slurm"]
RES    = config["resources"]
SLIDES = config["slides"]

SCRIPTS = os.path.join(workflow.basedir, "scripts")

# Every step writes into its own subfolder of paths.results_dir.
RESULTS         = P["results_dir"]
ZARR_DIR        = f"{RESULTS}/step1a/zarr_dir"                   # per-slide *.zarr
CONCAT_H5AD     = f"{RESULTS}/step1a/xenium_concatenated.h5ad"
PARQUET_DIR     = f"{RESULTS}/step1b/parquet_dir"                # ID_<slide>_intensity.parquet
QC_DIR          = f"{RESULTS}/step2"                             # final_adata.h5ad, ...
INTEGRATION_DIR = f"{RESULTS}/step3"                             # totalVI model, clusters, plots
LOGS            = f"{RESULTS}/logs"                              # one log per step

# Config file the step scripts read: the one passed with --configfile (e.g. the
# test config), otherwise config/config.yaml.
_CLI_CFGS = list(workflow.config_settings.configfiles)
CFG = os.path.abspath(_CLI_CFGS[-1] if _CLI_CFGS else DEFAULT_CONFIG)

# Python of each conda env inside integration.sif (fixed by micromamba.dockerfile).
PY = {
    "spatial":       "/opt/conda/envs/spatial/bin/python",            # step 1a
    "int_retrieval": "/opt/conda/envs/int_retrieval_env/bin/python",  # steps 1b + 2
    "scvi":          "/opt/conda/envs/scvi_env/bin/python",           # step 3 (GPU)
}

# Host folders the container must see: the inputs, results_dir, the scripts, the
# config file and any container.extra_binds. A folder already inside another one
# is dropped. $TMPDIR is added in apptainer() below.
def bind_list(paths):
    keep = []
    for p in sorted(set(paths)):
        if not any(p.startswith(k + os.sep) for k in keep):
            keep.append(p)
    return keep

BINDS = bind_list(
    [P["xenium_raw"], P["codex_base"], RESULTS, SCRIPTS, os.path.dirname(CFG)]
    + [abspath(b) for b in config["container"].get("extra_binds", [])]
)

# Per-slide CODEX parquets that step 1b must produce (explicit DAG targets).
PARQUETS = expand(f"{PARQUET_DIR}/ID_{{slide}}_intensity.parquet", slide=SLIDES)

# Step 3 plot/table names encode config values — built here the SAME way the
# script builds them, so they stay in sync when you change the thresholds.
IC = config["integration"]
_DE_STEM   = (f"top_{IC['top_rna']}_rna_top_{IC['top_pro']}_proteins_"
             f"log2BF_rna_{IC['bf_rna']}_protein_{IC['bf_protein']}")
DE_TABLE   = f"DE_table_{_DE_STEM}.csv"
DOTPLOT    = f"dotplot_{_DE_STEM}.png"
MATRIXPLOT = f"matrixplot_protein_log2BF_{IC['bf_protein']}.png"


# -- apptainer exec prefix ----------------------------------------------------
# Returns the full "mkdir ... && apptainer exec ... <env python>" command used as
# {params.exec} in each rule. numba/matplotlib caches are redirected into the
# job's local $TMPDIR so a non-writable HOME can't break scanpy's numba cache.
# nv=True exposes the GPU (needed by step 3 / scvi_env).
def apptainer(env_key, nv=False):
    binds = ",".join(BINDS)
    nvflag = "--nv " if nv else ""
    return (
        'mkdir -p "$TMPDIR/cache/numba" "$TMPDIR/cache/mpl" "$TMPDIR/cache/xdg" && '
        f'apptainer exec {nvflag}'
        f'-B {binds} -B "$TMPDIR" '
        '--env NUMBA_CACHE_DIR="$TMPDIR/cache/numba",'
        'MPLCONFIGDIR="$TMPDIR/cache/mpl",'
        'XDG_CACHE_HOME="$TMPDIR/cache/xdg" '
        f'{SIF} {PY[env_key]}'
    )


# =============================================================================
rule all:
    input:
        f"{INTEGRATION_DIR}/{MATRIXPLOT}",   # last file step 3 writes -> whole pipeline done


# =============================================================================
# Step 1a — Xenium concatenation                                    (env: spatial)
# =============================================================================
rule step1a_xenium_concatenation:
    input:
        xenium_raw=P["xenium_raw"],
    output:
        zarr=directory(ZARR_DIR),
        h5ad=CONCAT_H5AD,
    params:
        exec=apptainer("spatial"),
    log:
        f"{LOGS}/step1a_xenium_concatenation.log",
    threads: RES["step1a"]["threads"]
    resources:
        mem_mb_per_cpu=RES["step1a"]["mem_mb_per_cpu"],
        runtime=RES["step1a"]["runtime"],
        slurm_partition=SLURM["cpu_partition"],
    shell:
        r"""
        {params.exec} {SCRIPTS}/step1a_xenium_concatenation.py \
            --config-yaml {CFG} \
            --xenium-raw {input.xenium_raw} \
            --zarr-dir {output.zarr} \
            --output-h5ad {output.h5ad} \
            > {log} 2>&1
        """


# =============================================================================
# Step 1b — CODEX per-cell intensities                    (env: int_retrieval_env)
# =============================================================================
rule step1b_codex_intensities:
    input:
        zarr=rules.step1a_xenium_concatenation.output.zarr,
        codex_base=P["codex_base"],
    output:
        parquets=PARQUETS,
    params:
        exec=apptainer("int_retrieval"),
        out_dir=PARQUET_DIR,
    log:
        f"{LOGS}/step1b_codex_intensities.log",
    threads: RES["step1b"]["threads"]
    resources:
        mem_mb_per_cpu=RES["step1b"]["mem_mb_per_cpu"],
        runtime=RES["step1b"]["runtime"],
        slurm_partition=SLURM["cpu_partition"],
    shell:
        r"""
        {params.exec} {SCRIPTS}/step1b_codex_intensities.py \
            --config-yaml {CFG} \
            --xenium-zarr-base {input.zarr} \
            --codex-base {input.codex_base} \
            --output-dir {params.out_dir} \
            > {log} 2>&1
        """


# =============================================================================
# Step 2 — QC + modality merge                             (env: int_retrieval_env)
# =============================================================================
rule step2_data_qc:
    input:
        h5ad=rules.step1a_xenium_concatenation.output.h5ad,
        parquets=PARQUETS,
    output:
        final=f"{QC_DIR}/final_adata.h5ad",
        tracking=f"{QC_DIR}/tracking_counts.csv",
    params:
        exec=apptainer("int_retrieval"),
        out_dir=QC_DIR,
        parquet_dir=PARQUET_DIR,
    log:
        f"{LOGS}/step2_data_qc.log",
    threads: RES["step2"]["threads"]
    resources:
        mem_mb_per_cpu=RES["step2"]["mem_mb_per_cpu"],
        runtime=RES["step2"]["runtime"],
        slurm_partition=SLURM["cpu_partition"],
    shell:
        r"""
        {params.exec} {SCRIPTS}/step2_data_qc.py \
            --config-yaml {CFG} \
            --input-h5ad {input.h5ad} \
            --parquet-dir {params.parquet_dir} \
            --output-dir {params.out_dir} \
            > {log} 2>&1
        """


# =============================================================================
# Step 3 — totalVI integration + Leiden clustering            (env: scvi_env, GPU)
# =============================================================================
rule step3_data_integration:
    input:
        final=rules.step2_data_qc.output.final,
    output:
        mdata=f"{INTEGRATION_DIR}/mdata_leiden_dendrogram.h5mu",
        de=f"{INTEGRATION_DIR}/differential_expression.csv",
        de_table=f"{INTEGRATION_DIR}/{DE_TABLE}",
        elbo=f"{INTEGRATION_DIR}/elbo_training.png",
        umap=f"{INTEGRATION_DIR}/umap_clusters.png",
        umap_on_data=f"{INTEGRATION_DIR}/umap_clusters_on_data.png",
        forprob=f"{INTEGRATION_DIR}/forprob.png",
        denoised_protein=f"{INTEGRATION_DIR}/denoised_protein.png",
        dotplot=f"{INTEGRATION_DIR}/{DOTPLOT}",
        matrixplot=f"{INTEGRATION_DIR}/{MATRIXPLOT}",
    params:
        exec=apptainer("scvi", nv=True),
        save_dir=INTEGRATION_DIR,
    log:
        f"{LOGS}/step3_data_integration.log",
    threads: RES["step3"]["threads"]
    resources:
        mem_mb_per_cpu=RES["step3"]["mem_mb_per_cpu"],
        runtime=RES["step3"]["runtime"],
        slurm_partition=SLURM["gpu_partition"],
        # The SLURM executor plugin builds --gpus itself from this resource;
        # passing --gpus through slurm_extra is rejected by its validator.
        gpu=SLURM["gpus"],
    shell:
        r"""
        {params.exec} {SCRIPTS}/step3_data_integration.py \
            --config-yaml {CFG} \
            --input-h5ad {input.final} \
            --save-dir {params.save_dir} \
            > {log} 2>&1
        """
