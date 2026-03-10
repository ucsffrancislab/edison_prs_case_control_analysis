#!/usr/bin/env python3
"""
utils.py — Shared utilities for PGS case/control logistic regression meta-analysis pipeline.

Provides:
    - setup_logging()                : configure logging to console + file
    - load_cohort_data()             : load and merge PGS scores with covariates
    - drop_high_missingness_covariates() : remove covariates above missingness threshold
    - drop_zero_variance_covariates(): remove constant covariates
    - check_min_samples()            : verify minimum case/control counts
"""

import logging
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# Logging
# ─────────────────────────────────────────────────────────────────────────────

def setup_logging(verbose: bool = True, log_level: str = "INFO",
                  output_dir: Path = Path("results"), log_filename: str = "pipeline.log"):
    """Configure Python logging with timestamps; log to both console and file.

    Parameters
    ----------
    verbose : bool
        If True, set console to the specified log_level; otherwise WARNING.
    log_level : str
        One of DEBUG, INFO, WARNING, ERROR, CRITICAL.
    output_dir : Path
        Directory for the log file.
    log_filename : str
        Name of the log file.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / log_filename

    level = getattr(logging, log_level.upper(), logging.INFO)
    console_level = level if verbose else logging.WARNING

    # Root logger
    root = logging.getLogger()
    root.setLevel(logging.DEBUG)  # capture everything; handlers filter

    # Remove existing handlers (avoid duplicates on re-import)
    root.handlers.clear()

    fmt = logging.Formatter("[%(asctime)s] %(levelname)-8s %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(console_level)
    ch.setFormatter(fmt)
    root.addHandler(ch)

    # File handler
    fh = logging.FileHandler(log_path, mode="a")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(fmt)
    root.addHandler(fh)

    logging.info("Logging initialised — console=%s  file=%s", console_level, log_path)
    return root


# ─────────────────────────────────────────────────────────────────────────────
# Data loading
# ─────────────────────────────────────────────────────────────────────────────

def load_cohort_data(cohort_name: str, cohort_cfg: dict, base_dir: Path,
                     pgs_ids: list = None):
    """Load and merge PGS z-scores with covariates for one cohort.

    Steps:
        1. Read covariates CSV, filter to exclude==0
        2. Read PGS scores (only requested columns if pgs_ids given)
        3. Merge on sample ID
        4. Encode categorical covariates (sex as 0/1, source as 0/1 if present)
        5. Log shape and missingness

    Parameters
    ----------
    cohort_name : str
        Label for this cohort (used in logging).
    cohort_cfg : dict
        Entry from config.COHORTS (must have covariates_file, scores_file, etc.).
    base_dir : Path
        Root data directory.
    pgs_ids : list, optional
        Subset of PGS model IDs to load.  If None, load all.

    Returns
    -------
    df : pd.DataFrame
        Merged dataframe ready for modelling.
    all_pgs_ids : list
        List of PGS model column names present in df.
    """
    logger = logging.getLogger(__name__)

    # --- 1. Covariates -------------------------------------------------------
    cov_path = base_dir / cohort_cfg["covariates_file"]
    logger.info("[%s] Loading covariates from %s", cohort_name, cov_path)
    cov = pd.read_csv(cov_path)

    # Filter to exclude == 0 FIRST
    n_before = len(cov)
    cov = cov[cov["exclude"] == 0].copy()
    n_excluded = n_before - len(cov)
    logger.info("[%s] Excluded %d samples (exclude==1); %d remain",
                cohort_name, n_excluded, len(cov))

    # Rename sample ID column to 'sample_id' for merging
    sid_cov = cohort_cfg["sample_id_cov"]
    cov = cov.rename(columns={sid_cov: "sample_id"})

    # --- 2. PGS scores -------------------------------------------------------
    scores_path = base_dir / cohort_cfg["scores_file"]
    logger.info("[%s] Loading PGS scores from %s", cohort_name, scores_path)

    sid_scores = cohort_cfg["sample_id_scores"]

    if pgs_ids is not None:
        usecols = [sid_scores] + list(pgs_ids)
        scores = pd.read_csv(scores_path, usecols=usecols)
    else:
        scores = pd.read_csv(scores_path)

    scores = scores.rename(columns={sid_scores: "sample_id"})
    all_pgs_ids = [c for c in scores.columns if c != "sample_id"]
    logger.info("[%s] Loaded %d PGS models for %d samples",
                cohort_name, len(all_pgs_ids), len(scores))

    # --- 3. Merge -------------------------------------------------------------
    df = cov.merge(scores, on="sample_id", how="inner")
    logger.info("[%s] After merge: %d samples, %d columns",
                cohort_name, len(df), len(df.columns))

    # --- 4. Encode categoricals -----------------------------------------------
    # Sex: F->0, M->1
    if "sex" in df.columns:
        df["sex"] = df["sex"].map({"F": 0, "M": 1}).astype(float)

    # Source: encode as 0/1 if present and has 2 levels
    if "source" in df.columns:
        unique_sources = df["source"].dropna().unique()
        if len(unique_sources) == 2:
            source_map = {unique_sources[0]: 0, unique_sources[1]: 1}
            logger.info("[%s] Encoding source: %s", cohort_name, source_map)
            df["source"] = df["source"].map(source_map).astype(float)
        elif len(unique_sources) == 1:
            logger.info("[%s] Source has only 1 level — will be dropped as zero-variance",
                        cohort_name)

    # --- 5. Missingness summary -----------------------------------------------
    from config import BASE_COVARIATES
    extra_covs = cohort_cfg.get("extra_covariates", [])
    all_covs = BASE_COVARIATES + extra_covs
    for c in all_covs:
        if c in df.columns:
            n_miss = df[c].isna().sum()
            pct = n_miss / len(df) * 100
            if n_miss > 0:
                logger.info("[%s] Covariate '%s': %d missing (%.1f%%)",
                            cohort_name, c, n_miss, pct)

    logger.info("[%s] Cases=%d, Controls=%d",
                cohort_name, (df["case"] == 1).sum(), (df["case"] == 0).sum())

    return df, all_pgs_ids


# ─────────────────────────────────────────────────────────────────────────────
# Covariate QC
# ─────────────────────────────────────────────────────────────────────────────

def drop_high_missingness_covariates(df: pd.DataFrame, covariates: list,
                                     threshold: float = 0.20,
                                     cohort_name: str = ""):
    """Remove covariates where >threshold fraction of values are missing.

    Parameters
    ----------
    df : pd.DataFrame
        Data (only rows that will be used in analysis).
    covariates : list of str
        Candidate covariate names.
    threshold : float
        Maximum allowed fraction of missing values (default 0.20).
    cohort_name : str
        Label for logging.

    Returns
    -------
    kept : list of str
        Covariates that pass the missingness filter.
    """
    logger = logging.getLogger(__name__)
    logger.info("[%s] Checking all covariates for high missingness", cohort_name)
    kept = []
    #logger.info(df.head())
    for c in covariates:
        if c not in df.columns:
            logger.debug("[%s] Covariate '%s' not in dataframe — skipping", cohort_name, c)
            continue
        frac_miss = df[c].isna().mean()
        if frac_miss > threshold:
            logger.info("[%s] Dropping covariate '%s' — %.1f%% missing (threshold %.0f%%)",
                        cohort_name, c, frac_miss * 100, threshold * 100)
        else:
            kept.append(c)
    return kept


def drop_zero_variance_covariates(df: pd.DataFrame, covariates: list,
                                  cohort_name: str = ""):
    """Remove covariates with zero variance in the provided data.

    Parameters
    ----------
    df : pd.DataFrame
        Data subset (e.g., after dropping missing PGS/outcome).
    covariates : list of str
        Candidate covariate names.
    cohort_name : str
        Label for logging.

    Returns
    -------
    kept : list of str
        Covariates with non-zero variance.
    """
    logger = logging.getLogger(__name__)
    logger.info("[%s] Checking all covariates for variance", cohort_name)
    kept = []
    for c in covariates:
        if c not in df.columns:
            continue
        vals = df[c].dropna()
        logger.debug("[%s] Covariate '%s' len:%d; nuniq:%d", cohort_name, c, len(vals), vals.nunique())
        if len(vals) == 0 or vals.nunique() <= 1:
            logger.info("[%s] Dropping zero-variance covariate '%s'; len(%d); nuniq(%d)", cohort_name, c, len(vals), vals.nunique() )
        else:
            kept.append(c)
    return kept


# ─────────────────────────────────────────────────────────────────────────────
# Sample size check
# ─────────────────────────────────────────────────────────────────────────────

def check_min_samples(n_cases: int, n_controls: int,
                      min_cases: int = 10, min_controls: int = 10) -> bool:
    """Return True if both case and control counts meet minimums."""
    return n_cases >= min_cases and n_controls >= min_controls

def drop_outcome_separated_covariates(df: pd.DataFrame, covariates: list,
                                       outcome_col: str = "case",
                                       cohort_name: str = ""):
    """Remove covariates that perfectly separate the outcome (causing perfect separation).

    A covariate perfectly separates the outcome if each level of the covariate
    maps to only one level of the outcome. This causes logistic regression to fail
    with a singular matrix.

    Parameters
    ----------
    df : pd.DataFrame
    covariates : list of str
    outcome_col : str
    cohort_name : str

    Returns
    -------
    kept : list of str
    """
    logger = logging.getLogger(__name__)
    logger.info("[%s] Checking all covariates for outcome separation", cohort_name)
    kept = []
    for c in covariates:
        if c not in df.columns:
            continue
        # For each unique value of the covariate, check if outcome is constant
        vals = df[[c, outcome_col]].dropna()
        if len(vals) == 0:
            continue

        # Only check categorical-like covariates (few unique values)
        # Continuous variables with many unique values are not subject to this check
        n_unique = vals[c].nunique()
        if n_unique > 20:
            kept.append(c)
            continue

        # Check if covariate perfectly predicts outcome
        groups = vals.groupby(c)[outcome_col].nunique()
        if (groups == 1).all() and groups.shape[0] > 1:
            logger.warning("[%s] Dropping covariate '%s' — perfectly separates outcome "
                          "(each level maps to a single outcome value)", cohort_name, c)
        else:
            kept.append(c)
    return kept
