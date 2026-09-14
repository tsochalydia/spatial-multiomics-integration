#!/usr/bin/env bash
#SBATCH --job-name=smk_orchestrator
#SBATCH --output=logs/orchestrator_%j.log
#SBATCH --time=48:00:00          # must outlive the WHOLE pipeline (all 4 steps + queue waits)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G         # the orchestrator only waits + submits; the real work is in the step jobs
#SBATCH --partition=normal.120h
# ---------------------------------------------------------------------------
# Submit the whole workflow as ONE small orchestrator job:
#
#     sbatch submit.sh
#
# This job runs Snakemake, which then submits step1a, 1b, 2, 3 as their own
# SLURM jobs (with their own mem/GPU) in dependency order. Watch them with
# `squeue -u $USER`; per-step output is in logs/step*.log.
#
# (Equivalent to running ./run.sh run on the login node, just wrapped in sbatch.)
# ---------------------------------------------------------------------------
set -euo pipefail

MM=/cluster/project/moor/lydia/bin/micromamba
export MAMBA_ROOT_PREFIX=/cluster/project/moor/lydia/micromamba

# Absolute path: under sbatch, $0 is SLURM's spool copy, so dirname "$0" is wrong.
WF=/cluster/project/moor/lydia/snakemake_workflow
cd "$WF"                         # workflow dir (Snakefile + config.yaml + profiles/)

$MM run -n snakemake_env snakemake \
    --workflow-profile "$WF/profiles/slurm" \
    --rerun-triggers mtime
