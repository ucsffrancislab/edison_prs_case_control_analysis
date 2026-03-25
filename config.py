#!/usr/bin/env python3
"""
config.py — Central configuration for PGS case/control logistic regression meta-analysis pipeline.

Edit BASE_DIR and OUTPUT_DIR to point at your data and desired output locations.
Everything else can remain as-is for standard runs.
"""

from pathlib import Path

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE_DIR = Path("input")          # Root directory containing cohort data files
OUTPUT_DIR = Path("results")  # Pipeline outputs (tables, logs)
PLOT_DIR = Path("plots")      # Plot outputs

# ── Cohort definitions ─────────────────────────────────────────────────────────
# Each cohort maps to:
#   - covariates_file : filename of covariates CSV (relative to BASE_DIR)
#   - scores_file     : filename of z-scored PGS CSV (relative to BASE_DIR)
#   - extra_covariates: list of additional covariates beyond BASE_COVARIATES
#   - sample_id_cov   : column name for sample ID in covariates file
#   - sample_id_scores: column name for sample ID in scores file
COHORTS = {
    "cidr": {
        "covariates_file": "cidr-covariates.csv",
        "scores_file": "cidr.scores.z-scores.txt.gz",
        #"extra_covariates": ["source"],       # two recruitment sources (CIDR, MDSAML) # cases and controls separate could be issue
        "extra_covariates": [],       # two recruitment sources (CIDR, MDSAML) # cases and controls separate could be issue
        "sample_id_cov": "IID",
        "sample_id_scores": "sample",
    },
    "onco": {
        "covariates_file": "onco-covariates.csv",
        "scores_file": "onco.scores.z-scores.txt.gz",
        "extra_covariates": ["source"],       # two recruitment sources (AGS, Mayo)
        "sample_id_cov": "IID",
        "sample_id_scores": "sample",
    },
    "i370": {
        "covariates_file": "i370-covariates.csv",
        "scores_file": "i370.scores.z-scores.txt.gz",
        "extra_covariates": [],               # single source, no extra covariates
        "sample_id_cov": "IID",
        "sample_id_scores": "sample",
    },
    "tcga": {
        "covariates_file": "tcga-covariates.csv",
        "scores_file": "tcga.scores.z-scores.txt.gz",
        "extra_covariates": ["source"],       # two sources (TCGA, WTCCC)
        "sample_id_cov": "IID",
        "sample_id_scores": "sample",
    },
}

# ── Covariates ─────────────────────────────────────────────────────────────────
BASE_COVARIATES = ["age", "sex", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8"]

# ── QC thresholds ──────────────────────────────────────────────────────────────
MIN_CASES = 10          # Minimum number of cases to fit a model
MIN_CONTROLS = 10       # Minimum number of controls to fit a model
MISSINGNESS_THRESHOLD = 0.20  # Drop covariates missing > this fraction

# ── Multiple-testing correction ────────────────────────────────────────────────
FDR_ALPHA = 0.05

# ── Plotting ───────────────────────────────────────────────────────────────────
TOP_N_FOREST = 20       # Number of top models to create individual forest plots for
TOP_N_SUMMARY = 50      # Number of top models for summary forest plot

# ── Subtype filtering ──────────────────────────────────────────────────────────
# Column names in the covariates files for subtype classification.
# Set to None if your cohort files use different column names — override via CLI.
IDH_COLUMN = "idh"      # Expected values: "wt", "mt"  (or None to skip)
PQ_COLUMN  = "pq"       # Expected values: "codel", "intact"  (or None to skip)

# ── Test mode ──────────────────────────────────────────────────────────────────
TEST_N_MODELS = 50      # Number of PGS models to run in --test mode

# ── Logging ────────────────────────────────────────────────────────────────────
VERBOSE = True
LOG_LEVEL = "INFO"      # DEBUG / INFO / WARNING

