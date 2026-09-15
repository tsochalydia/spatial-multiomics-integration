#!/usr/bin/env bash
# Launch the Xenium+CODEX integration workflow.
#
#   ./run.sh dry     # dry run: print the plan (no jobs submitted, nothing runs)
#   ./run.sh run     # submit one SLURM job per step, in order; this shell is the
#                    #   orchestrator, so keep it alive (tmux/screen) until done
#   ./run.sh submit  # same as run, but the orchestrator itself runs as a small
#                    #   SLURM job, so you can log out
#   ./run.sh dag     # write dag.svg (or dag.dot) of the step graph
#   ./run.sh unlock  # release a stale .snakemake lock after a killed run
#
# Config: config.yaml by default; another one with --configfile, e.g.
#   ./run.sh dry --configfile test_data/config.test.yaml
# Any other options are passed on to Snakemake, e.g.  ./run.sh dry --forceall
#
# Snakemake is taken from $SNAKEMAKE if set, else `snakemake` on PATH (activated
# snakemake_env), else `micromamba run -n snakemake_env snakemake`.
# dry/run/submit use --rerun-triggers mtime: a step reruns only when its inputs are
# newer than its outputs, so editing a config comment or path does not redo steps.
# submit writes the orchestrator log to <results_dir>/logs/orchestrator_<jobid>.log;
# set SBATCH_PARTITION / SBATCH_ACCOUNT if your cluster needs them.
set -euo pipefail

cd "$(dirname "$(readlink -f "$0")")"      # repo root (holds Snakefile + config.yaml)
REPO=$PWD
PROFILE=profiles/slurm
MODE="${1:-dry}"
shift || true

# -- how to call Snakemake ----------------------------------------------------
if [ -n "${SNAKEMAKE:-}" ]; then
    read -r -a SMK <<< "$SNAKEMAKE"
elif command -v snakemake > /dev/null; then
    SMK=(snakemake)
elif MM="${MAMBA_EXE:-$(command -v micromamba || true)}" && [ -n "$MM" ]; then
    SMK=("$MM" run -n snakemake_env snakemake)
else
    echo "Snakemake not found: activate snakemake_env or set SNAKEMAKE=/path/to/snakemake" >&2
    exit 1
fi

# -- config file in use (for the orchestrator log location) --------------------
CONFIG=config.yaml
ARGS=("$@")
for i in "${!ARGS[@]}"; do
    case "${ARGS[$i]}" in
        --configfile)   CONFIG="${ARGS[$((i + 1))]:-$CONFIG}" ;;
        --configfile=*) CONFIG="${ARGS[$i]#*=}" ;;
    esac
done

results_dir() {   # paths.results_dir from a config file
    sed -n -E "s/^[[:space:]]+results_dir:[[:space:]]*[\"']?([^\"'#[:space:]]+).*/\1/p" "$1" | head -1
}

# -----------------------------------------------------------------------------
case "$MODE" in
    dry)
        "${SMK[@]}" -n -p --rerun-triggers mtime "$@"
        ;;
    run)
        "${SMK[@]}" --workflow-profile "$PROFILE" --rerun-triggers mtime "$@"
        ;;
    submit)
        RESULTS=$(results_dir "$CONFIG")
        [ -n "$RESULTS" ] || { echo "no paths.results_dir in $CONFIG" >&2; exit 1; }
        LOGDIR=$(readlink -m "$RESULTS/logs")
        mkdir -p "$LOGDIR"
        # The orchestrator only waits and submits; the steps get their own jobs.
        # --time must outlive the whole pipeline (all steps + queue waits).
        JOB=$(sbatch --parsable \
            --job-name=smk_orchestrator \
            --output="$LOGDIR/orchestrator_%j.log" \
            --time=48:00:00 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G \
            --wrap="$(printf '%q ' "$REPO/run.sh" run "$@")")
        echo "submitted orchestrator job $JOB (log: $LOGDIR/orchestrator_$JOB.log)"
        echo "watch the step jobs with: squeue -u $USER"
        ;;
    dag)
        if command -v dot > /dev/null; then
            "${SMK[@]}" --dag "$@" | dot -Tsvg > dag.svg
            echo "wrote dag.svg"
        else
            "${SMK[@]}" --dag "$@" > dag.dot
            echo "graphviz 'dot' not found; wrote dag.dot instead"
        fi
        ;;
    unlock)
        "${SMK[@]}" --unlock "$@"
        ;;
    *)
        echo "usage: $0 {dry|run|submit|dag|unlock} [--configfile FILE] [snakemake options]" >&2
        exit 1
        ;;
esac
