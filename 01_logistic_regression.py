#!/usr/bin/env python3
"""
01_logistic_regression.py — Per-cohort logistic regression of case/control status on PGS.

For each cohort × PGS model:
    case ~ PGS + age + sex + PC1..PC8 [+ source if applicable]

Results are written to results/{cohort}_logistic_results.tsv.

Usage:
    python 01_logistic_regression.py                       # full run
    python 01_logistic_regression.py --test                # TEST_N_MODELS only
    python 01_logistic_regression.py --cohort onco --test  # single cohort, test mode
    python 01_logistic_regression.py --n-jobs 16           # parallel threads
"""

import argparse
import logging
import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from joblib import Parallel, delayed

# Local imports
from config import (BASE_COVARIATES, BASE_DIR, COHORTS, FDR_ALPHA, LOG_LEVEL,
                    MIN_CASES, MIN_CONTROLS, MISSINGNESS_THRESHOLD,
                    OUTPUT_DIR, TEST_N_MODELS, VERBOSE)
from utils import (check_min_samples, drop_high_missingness_covariates,
                   drop_outcome_separated_covariates,
                   drop_zero_variance_covariates, load_cohort_data,
                   setup_logging)

warnings.filterwarnings("ignore", category=FutureWarning)


# ─────────────────────────────────────────────────────────────────────────────
# Single-model fit
# ─────────────────────────────────────────────────────────────────────────────

def fit_one_model(pgs_id: str, df: pd.DataFrame, covariates: list,
                  cohort_name: str, missingness_threshold: float,
                  min_cases: int, min_controls: int):
    """Fit logistic regression for one PGS model in one cohort.

    Returns a dict with results or failure info.
    """
    result = {
        "pgs_id": pgs_id,
        "cohort": cohort_name,
        "log_or": np.nan,
        "se": np.nan,
        "z": np.nan,
        "pvalue": np.nan,
        "n_cases": 0,
        "n_controls": 0,
        "n_total": 0,
        "converged": False,
        "covariate_list": "",
        "status": "failed",
        "error_msg": "",
    }

    logger = logging.getLogger(__name__)
    logger.info("pgs: %s", pgs_id)
    try:
        # Subset to outcome + PGS + covariates columns
        cols_needed = ["case", pgs_id] + [c for c in covariates if c in df.columns]
        sub = df[cols_needed].copy()
        #logger.info("sub size: %s", sub.shape)

        # Drop rows missing outcome or PGS
        sub = sub.dropna(subset=["case", pgs_id])

        # Drop high-missingness covariates (within this subset)
        usable_covs = drop_high_missingness_covariates(
            sub, covariates, threshold=missingness_threshold, cohort_name=cohort_name
        )
        logger.info("pgs: %s ; cohort: %s  ; covariates remaining: %s",pgs_id, cohort_name, usable_covs)

        # Drop rows with missing covariates (after removing high-miss covariates)
        sub = sub.dropna(subset=[pgs_id, "case"] + usable_covs)

        # Drop zero-variance covariates
        usable_covs = drop_zero_variance_covariates(sub, usable_covs, cohort_name=cohort_name)
        logger.info("pgs: %s ; cohort: %s  ; covariates remaining: %s",pgs_id, cohort_name, usable_covs)

        # Drop covariates that perfectly separate the outcome
        usable_covs = drop_outcome_separated_covariates(sub, usable_covs,
                                                         outcome_col="case",
                                                         cohort_name=cohort_name)
        logger.info("pgs: %s ; cohort: %s  ; covariates remaining: %s",pgs_id, cohort_name, usable_covs)

        # Sample counts
        n_cases = int((sub["case"] == 1).sum())
        n_controls = int((sub["case"] == 0).sum())
        n_total = len(sub)
        result["n_cases"] = n_cases
        result["n_controls"] = n_controls
        result["n_total"] = n_total

        if not check_min_samples(n_cases, n_controls, min_cases, min_controls):
            result["status"] = "skipped_low_n"
            result["error_msg"] = f"cases={n_cases}, controls={n_controls}"
            return result

        # Build design matrix: PGS + covariates + intercept
        X_cols = [pgs_id] + usable_covs
        X = sub[X_cols].astype(float)
        X = sm.add_constant(X)
        y = sub["case"].astype(float)

        # Fit logistic regression
        model = sm.Logit(y, X)
        fit = model.fit(disp=0, maxiter=100, method="newton")

        result["log_or"] = fit.params[pgs_id]
        result["se"] = fit.bse[pgs_id]
        result["z"] = fit.tvalues[pgs_id]
        result["pvalue"] = fit.pvalues[pgs_id]
        result["converged"] = bool(fit.mle_retvals.get("converged", False))
        result["covariate_list"] = ";".join(usable_covs)
        result["status"] = "success"

    except Exception as e:
        result["status"] = "failed"
        result["error_msg"] = str(e)[:200]

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Cohort-level runner
# ─────────────────────────────────────────────────────────────────────────────

def run_cohort(cohort_name: str, cohort_cfg: dict, pgs_ids: list,
               n_jobs: int = 1, test_mode: bool = False,
               idh_subtype: str = None, pq_subtype: str = None,
               idh_column: str = None, pq_column: str = None,
               outdir: Path = None):
    """Run logistic regression across all PGS models for one cohort.

    Parameters
    ----------
    cohort_name : str
    cohort_cfg : dict
    pgs_ids : list
        PGS model IDs to evaluate (already subsetted in test mode).
    n_jobs : int
    test_mode : bool
    idh_subtype : str, optional
        Restrict cases to this IDH value (e.g. "wt" or "mt").
    pq_subtype : str, optional
        Restrict cases to this 1p19q value (e.g. "codel" or "intact").
    idh_column : str, optional
        Override the IDH column name in the covariates file.
    pq_column : str, optional
        Override the 1p19q column name in the covariates file.

    Returns
    -------
    results_df : pd.DataFrame
    """
    logger = logging.getLogger(__name__)
    logger.info("=" * 70)
    logger.info("Starting cohort: %s  (%d models, n_jobs=%d)",
                cohort_name, len(pgs_ids), n_jobs)

    # Log active subtype filters
    if idh_subtype or pq_subtype:
        logger.info("[%s] Subtype filters — IDH: %s  1p19q: %s",
                    cohort_name, idh_subtype or "all", pq_subtype or "all")

    # Load data (only the PGS columns we need)
    df, available_pgs = load_cohort_data(cohort_name, cohort_cfg, BASE_DIR,
                                          pgs_ids=pgs_ids,
                                          idh_subtype=idh_subtype,
                                          pq_subtype=pq_subtype,
                                          idh_column=idh_column,
                                          pq_column=pq_column)

    # Determine which PGS IDs are actually available
    pgs_to_run = [p for p in pgs_ids if p in df.columns]
    logger.info("[%s] %d / %d requested PGS models found in data",
                cohort_name, len(pgs_to_run), len(pgs_ids))

    # Build full covariate list for this cohort
    covariates = list(BASE_COVARIATES) + cohort_cfg.get("extra_covariates", [])

    # Pre-filter: drop cohort-level high-missingness covariates once
    # (per-model check will refine further)
    cohort_covs = drop_high_missingness_covariates(
        df, covariates, threshold=MISSINGNESS_THRESHOLD, cohort_name=cohort_name
    )
    # Drop covariates that perfectly separate outcome at cohort level
    cohort_covs = drop_outcome_separated_covariates(
        df, cohort_covs, outcome_col="case", cohort_name=cohort_name
    )
    logger.info("[%s] Cohort-level covariates after missingness + separation filter: %s",
                cohort_name, cohort_covs)

    # Progress logging interval
    log_every = 10 if test_mode else 500

    def _fit_with_progress(i, pid):
        if (i + 1) % log_every == 0:
            logger.info("[%s] Progress: %d / %d models", cohort_name, i + 1, len(pgs_to_run))
        return fit_one_model(pid, df, cohort_covs, cohort_name,
                             MISSINGNESS_THRESHOLD, MIN_CASES, MIN_CONTROLS)

    # Parallel execution
    t0 = time.time()
    results = Parallel(n_jobs=n_jobs, prefer="threads")(
        delayed(_fit_with_progress)(i, pid)
        for i, pid in enumerate(pgs_to_run)
    )
    elapsed = time.time() - t0

    results_df = pd.DataFrame(results)

    # Summary
    n_success = (results_df["status"] == "success").sum()
    n_failed = (results_df["status"] == "failed").sum()
    n_skipped = (results_df["status"] == "skipped_low_n").sum()
    logger.info("[%s] Done in %.1fs — %d success, %d failed, %d skipped (low N)",
                cohort_name, elapsed, n_success, n_failed, n_skipped)

    # Save
    save_dir = outdir if outdir is not None else OUTPUT_DIR
    save_dir.mkdir(parents=True, exist_ok=True)
    subtype_tag = _build_subtype_tag(idh_subtype, pq_subtype)
    out_path = save_dir / f"{cohort_name}{subtype_tag}_logistic_results.tsv"
    results_df.to_csv(out_path, sep="\t", index=False)
    logger.info("[%s] Results saved to %s", cohort_name, out_path)

    return results_df


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _build_subtype_tag(idh_subtype: str = None, pq_subtype: str = None) -> str:
    """Build a filename-safe tag string from active subtype filters.

    Examples
    --------
    ("wt", None)        -> "_IDHwt"
    ("mt", "codel")     -> "_IDHmt_codel"
    (None, None)        -> ""
    """
    parts = []
    if idh_subtype:
        parts.append(f"IDH{idh_subtype}")
    if pq_subtype:
        parts.append(pq_subtype)
    return ("_" + "_".join(parts)) if parts else ""


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Per-cohort PGS logistic regression")
    parser.add_argument("--test", action="store_true",
                        help=f"Run only {TEST_N_MODELS} models for quick validation")
    parser.add_argument("--verbose", action="store_true", default=VERBOSE)
    parser.add_argument("--cohort", type=str, default=None,
                        help="Run a single cohort (e.g., onco)")
    parser.add_argument("--n-jobs", type=int,
                        default=int(os.environ.get("SLURM_CPUS_PER_TASK", 1)),
                        help="Number of parallel threads")
    parser.add_argument("--idh-subtype", type=str, default=None,
                        metavar="VALUE",
                        help="Restrict cases to this IDH value (e.g. wt or mt). "
                             "Controls are always retained.")
    parser.add_argument("--pq-subtype", type=str, default=None,
                        metavar="VALUE",
                        help="Restrict cases to this 1p19q value (e.g. codel or intact). "
                             "Controls are always retained.")
    parser.add_argument("--idh-column", type=str, default=None,
                        metavar="COLNAME",
                        help="Override the IDH column name in covariates files "
                             f"(default: config.IDH_COLUMN)")
    parser.add_argument("--pq-column", type=str, default=None,
                        metavar="COLNAME",
                        help="Override the 1p19q column name in covariates files "
                             f"(default: config.PQ_COLUMN)")
    parser.add_argument("--outdir", type=str, default=None,
                        metavar="DIR",
                        help="Output directory for results and logs. "
                             "Created if it does not exist. Overrides config OUTPUT_DIR.")
    args = parser.parse_args()

    # Resolve output directory
    outdir = Path(args.outdir) if args.outdir else OUTPUT_DIR
    outdir.mkdir(parents=True, exist_ok=True)

    # Logging
    log_tag = "test" if args.test else "full"
    subtype_tag = _build_subtype_tag(args.idh_subtype, args.pq_subtype)
    setup_logging(args.verbose, LOG_LEVEL, outdir,
                  log_filename=f"01_logistic{subtype_tag}_{log_tag}.log")
    logger = logging.getLogger(__name__)
    logger.info("Pipeline step 01 — logistic regression (mode=%s)", log_tag)
    if args.idh_subtype or args.pq_subtype:
        logger.info("Subtype filters — IDH: %s  1p19q: %s",
                    args.idh_subtype or "all", args.pq_subtype or "all")

    # Determine which cohorts to run
    cohorts_to_run = {args.cohort: COHORTS[args.cohort]} if args.cohort else COHORTS
    logger.debug("cohorts (%s)", cohorts_to_run)

    # Discover PGS model IDs from the first cohort's scores file
    first_cohort_cfg = list(COHORTS.values())[0]
    scores_path = BASE_DIR / first_cohort_cfg["scores_file"]
    logger.info("Reading PGS model list from %s", scores_path)
    header = pd.read_csv(scores_path, nrows=0).columns.tolist()
    sid = first_cohort_cfg["sample_id_scores"]
    all_pgs_ids = [c for c in header if c != sid]
    logger.info("Total PGS models available: %d", len(all_pgs_ids))

    if args.test:
        pgs_ids = all_pgs_ids[:TEST_N_MODELS]
        logger.info("TEST MODE: limiting to first %d models", TEST_N_MODELS)
    else:
        pgs_ids = all_pgs_ids

    # Run each cohort
    all_results = []
    for cname, ccfg in cohorts_to_run.items():
        logger.debug("running cohort (%s,%s)", cname, ccfg)
        res = run_cohort(cname, ccfg, pgs_ids,
                         n_jobs=args.n_jobs,
                         test_mode=args.test,
                         idh_subtype=args.idh_subtype,
                         pq_subtype=args.pq_subtype,
                         idh_column=args.idh_column,
                         pq_column=args.pq_column,
                         outdir=outdir)
        all_results.append(res)

    # Combined summary
    combined = pd.concat(all_results, ignore_index=True)
    logger.info("=" * 70)
    logger.info("OVERALL SUMMARY: %d total model-cohort pairs", len(combined))
    logger.info("  Success:  %d", (combined["status"] == "success").sum())
    logger.info("  Failed:   %d", (combined["status"] == "failed").sum())
    logger.info("  Skipped:  %d", (combined["status"] == "skipped_low_n").sum())

    return combined


if __name__ == "__main__":
    main()
