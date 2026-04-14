#!/bin/bash
#SBATCH --job-name=pgs_meta
#SBATCH --cpus-per-task=64
#SBATCH --mem=490G
#SBATCH --time=14-0
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# ──────────────────────────────────────────────────────────────────────────────
# run_pipeline.sh — SLURM submission script for PGS meta-analysis pipeline
#
# Usage:
#   sbatch run_pipeline.sh --outdir allcases
#   sbatch run_pipeline.sh --outdir idhwt_cases  --idh-subtype wt
#   sbatch run_pipeline.sh --outdir idhmt_cases  --idh-subtype mt
#   sbatch run_pipeline.sh --outdir idhmt_codel  --idh-subtype mt --pq-subtype codel
#   sbatch run_pipeline.sh --outdir idhmt_intact --idh-subtype mt --pq-subtype intact
#   sbatch run_pipeline.sh --outdir test_run     --test
#
# --outdir is required. All output (logs, results, plots) goes there.
# All Python stdout/stderr is captured to <outdir>/pipeline.out.txt so the
# SLURM output files are intentionally empty.
#
# This script can be called from any directory — subscripts are located
# relative to the script's own directory, not the working directory.
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Locate the directory this script lives in ────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Parse --outdir from args, collect remainder for Python ───────────────────
OUTDIR=""
PYTHON_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        *)
            PYTHON_ARGS+=("$1")
            shift
            ;;
    esac
done

if [[ -z "$OUTDIR" ]]; then
    echo "ERROR: --outdir is required" >&2
    echo "Usage: sbatch run_pipeline.sh --outdir <directory> [pipeline options]" >&2
    exit 1
fi

# ── Create output directory ───────────────────────────────────────────────────
mkdir -p "$OUTDIR"

# ── Redirect all subsequent stdout/stderr to a log file in outdir ────────────
exec > "${OUTDIR}/pipeline.out.txt" 2>&1

# ── Job info ──────────────────────────────────────────────────────────────────
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "CPUs:   $SLURM_CPUS_PER_TASK"
echo "Outdir: $OUTDIR"
echo "Start:  $(date)"

# ── Environment activation ───────────────────────────────────────────────────
# Uncomment and edit ONE of the following:
# conda activate pgs_env
# source /path/to/venv/bin/activate
# module load R/4.3 python/3.10

module load r

# ── Run pipeline ─────────────────────────────────────────────────────────────
python3 "$SCRIPT_DIR/04_run_pipeline.py" \
    --n-jobs ${SLURM_CPUS_PER_TASK:-16} \
    --verbose \
    --outdir "$OUTDIR" \
    "${PYTHON_ARGS[@]}"

echo "End:    $(date)"
echo "Done."
