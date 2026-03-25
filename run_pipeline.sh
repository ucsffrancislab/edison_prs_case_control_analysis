#!/bin/bash
#SBATCH --job-name=pgs_meta
#SBATCH --cpus-per-task=64
#SBATCH --mem=490G
#SBATCH --time=14-0
#SBATCH --output=slurm_%j.out.txt
#SBATCH --error=slurm_%j.err.txt

# ──────────────────────────────────────────────────────────────────────────────
# run_pipeline.sh — SLURM submission script for PGS meta-analysis pipeline
#
# Usage:
#   sbatch run_pipeline.sh            # full run
#   sbatch run_pipeline.sh --test     # test mode (50 models)
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "CPUs:   $SLURM_CPUS_PER_TASK"
echo "Start:  $(date)"

# ── Environment activation ───────────────────────────────────────────────────
# Uncomment and edit ONE of the following:
# conda activate pgs_env
# source /path/to/venv/bin/activate
# module load R/4.3 python/3.10

module load r

# ── Run pipeline ─────────────────────────────────────────────────────────────
python3 04_run_pipeline.py \
    --n-jobs ${SLURM_CPUS_PER_TASK:-16} \
    --verbose \
    "$@"

echo "End:    $(date)"
echo "Done."
