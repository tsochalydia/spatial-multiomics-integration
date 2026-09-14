#!/usr/bin/env bash
# Launch the Xenium+CODEX integration workflow.
#
#   ./run.sh dry     # dry run: print the plan (no jobs submitted, nothing runs)
#   ./run.sh dag     # write dag.svg of the rule graph
#   ./run.sh run     # submit one SLURM job per step, in order (login-node orchestrator)
#   ./run.sh unlock  # release a stale .snakemake lock after a killed run
#
# Snakemake runs from snakemake_env and submits each step into integration.sif.
# Keep this shell (or a tmux/screen) alive while `run` is going — it's the
# orchestrator that tracks the jobs. For a long run, wrap it in tmux.
set -euo pipefail

MM=/cluster/project/moor/lydia/bin/micromamba
export MAMBA_ROOT_PREFIX=/cluster/project/moor/lydia/micromamba

cd "$(dirname "$0")"                       # workflow dir (holds Snakefile + config.yaml)
PROFILE=profiles/slurm
MODE="${1:-dry}"

case "$MODE" in
  dry)
    $MM run -n snakemake_env snakemake -n -p --rerun-triggers mtime
    ;;
  dag)
    $MM run -n snakemake_env snakemake --dag | \
      $MM run -n snakemake_env dot -Tsvg > dag.svg 2>/dev/null || \
      echo "graphviz 'dot' not in snakemake_env; run './run.sh dry' instead"
    echo "wrote dag.svg"
    ;;
  run)
    $MM run -n snakemake_env snakemake --workflow-profile "$PROFILE" --rerun-triggers mtime
    ;;
  unlock)
    $MM run -n snakemake_env snakemake --unlock
    ;;
  *)
    echo "usage: $0 {dry|dag|run|unlock}"; exit 1
    ;;
esac
