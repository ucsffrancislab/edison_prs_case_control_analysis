#!/bin/bash
#SBATCH --job-name=pgs_meta
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=120G
#SBATCH --time=14-0
#SBATCH --output=/dev/null
#SBATCH --error=/dev/null

# ──────────────────────────────────────────────────────────────────────────────
# run_pipeline.sh — SLURM submission script for PGS meta-analysis pipeline
#
# Usage:
#   sbatch ~/github/ucsffrancislab/edison_prs_case_control_analysis/run_pipeline.sh --outdir allcases
#   sbatch ~/github/ucsffrancislab/edison_prs_case_control_analysis/run_pipeline.sh --outdir idhwt_cases  --idh-subtype wt
#   sbatch ~/github/ucsffrancislab/edison_prs_case_control_analysis/run_pipeline.sh --outdir idhmt_cases  --idh-subtype mt
#   sbatch ~/github/ucsffrancislab/edison_prs_case_control_analysis/run_pipeline.sh --outdir idhmt_codel  --idh-subtype mt --pq-subtype codel
#   sbatch ~/github/ucsffrancislab/edison_prs_case_control_analysis/run_pipeline.sh --outdir idhmt_intact --idh-subtype mt --pq-subtype intact
#   sbatch ~/github/ucsffrancislab/edison_prs_case_control_analysis/run_pipeline.sh --outdir test_run     --test
#
# --outdir is required. All output (logs, results, plots) goes there.
# All Python stdout/stderr is captured to <outdir>/pipeline.out.txt so the
# SLURM output files are intentionally empty.
#
# This script can be called from any directory — subscripts are located
# via scontrol (under SLURM) or dirname "$0" (interactive/local runs).
# ──────────────────────────────────────────────────────────────────────────────

set -euo pipefail

# ── Locate the directory this script lives in ────────────────────────────────
# Under SLURM, sbatch copies the script to a spool directory so dirname "$0"
# points to the wrong place. Use scontrol to recover the original path.
# Outside SLURM (interactive/local), dirname "$0" works fine.
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    PIPELINE_DIR=$(dirname "$(scontrol show job "$SLURM_JOB_ID" \
        | awk '/Command=/{sub(/.*Command=/, ""); print $1}')")
else
    PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
fi


#DATA_DIR="${1:-input}"
#OUTDIR="${2:-results}"
#EXTRA_ARGS="${@:3}"

# ── Parse --outdir from args, collect remainder for Python ───────────────────
DATA_DIR="input"
OUTDIR="results"
PYTHON_ARGS=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --datadir)
            DATA_DIR="$2"
            shift 2
            ;;
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
echo "Job ID:       ${SLURM_JOB_ID:-local}"
echo "Node:         $(hostname)"
echo "CPUs:         ${SLURM_CPUS_PER_TASK:-$(nproc)}"
echo "Outdir:       $OUTDIR"
echo "Pipeline dir: $PIPELINE_DIR"
echo "Start:        $(date)"

# ── Environment activation ───────────────────────────────────────────────────
# Uncomment and edit ONE of the following:
# conda activate pgs_env
# source /path/to/venv/bin/activate
# module load R/4.3 python/3.10

module load r

# ── Run pipeline ─────────────────────────────────────────────────────────────
python3 "$PIPELINE_DIR/04_run_pipeline.py" \
    --n-jobs ${SLURM_CPUS_PER_TASK:-16} \
    --verbose \
    --outdir "$OUTDIR" \
    "${PYTHON_ARGS[@]}"

echo "End:    $(date)"
echo "Done."
