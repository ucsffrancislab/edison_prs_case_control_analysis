#!/usr/bin/env python3
"""
04_run_pipeline.py — Master runner for PGS case/control meta-analysis pipeline.

Calls steps 01–03 in order, passing through --test and --verbose flags.
Logs start/end time and wall-clock duration for each step.

Usage:
    python 04_run_pipeline.py                # full run
    python 04_run_pipeline.py --test         # test mode (50 models)
    python 04_run_pipeline.py --test --verbose
"""

import argparse
import logging
import subprocess
import sys
import time
from pathlib import Path

from config import LOG_LEVEL, OUTPUT_DIR, PLOT_DIR, VERBOSE
from utils import setup_logging


def run_step(step_name: str, cmd: list, logger: logging.Logger):
    """Execute a pipeline step as a subprocess.

    Parameters
    ----------
    step_name : str
        Human-readable step label.
    cmd : list
        Command + arguments.
    logger : logging.Logger
    """
    logger.info("─" * 60)
    logger.info("STEP: %s", step_name)
    logger.info("CMD:  %s", " ".join(cmd))
    t0 = time.time()

    proc = subprocess.run(cmd, capture_output=True, text=True)

    elapsed = time.time() - t0
    logger.info("  stdout:\n%s", proc.stdout[-3000:] if len(proc.stdout) > 3000 else proc.stdout)
    if proc.stderr:
        logger.info("  stderr:\n%s", proc.stderr[-3000:] if len(proc.stderr) > 3000 else proc.stderr)

    if proc.returncode != 0:
        logger.error("STEP %s FAILED (exit code %d) after %.1fs",
                      step_name, proc.returncode, elapsed)
        logger.error("stderr:\n%s", proc.stderr)
        sys.exit(proc.returncode)
    else:
        logger.info("STEP %s completed in %.1fs", step_name, elapsed)


def main():
    parser = argparse.ArgumentParser(description="PGS meta-analysis pipeline runner")
    parser.add_argument("--test", action="store_true", help="Test mode")
    parser.add_argument("--verbose", action="store_true", default=VERBOSE)
    parser.add_argument("--n-jobs", type=int, default=1)
    args = parser.parse_args()

    log_tag = "test" if args.test else "full"
    setup_logging(args.verbose, LOG_LEVEL, OUTPUT_DIR,
                  log_filename=f"04_pipeline_{log_tag}.log")
    logger = logging.getLogger(__name__)

    logger.info("=" * 70)
    logger.info("PGS CASE/CONTROL META-ANALYSIS PIPELINE")
    logger.info("Mode: %s", "TEST" if args.test else "FULL")
    logger.info("=" * 70)

    t_start = time.time()

    # Step 01: Logistic regression
    cmd_01 = [sys.executable, "01_logistic_regression.py",
              "--n-jobs", str(args.n_jobs)]
    if args.test:
        cmd_01.append("--test")
    if args.verbose:
        cmd_01.append("--verbose")
    run_step("01_logistic_regression", cmd_01, logger)

    # Step 02: Meta-analysis (R)
    cmd_02 = ["Rscript", "02_meta_analysis.R",
              "--results-dir", str(OUTPUT_DIR)]
    run_step("02_meta_analysis", cmd_02, logger)

    # Step 03: Plots (R)
    cmd_03 = ["Rscript", "03_plots.R",
              "--results-dir", str(OUTPUT_DIR),
              "--plot-dir", str(PLOT_DIR)]
    run_step("03_plots", cmd_03, logger)

    total = time.time() - t_start
    logger.info("=" * 70)
    logger.info("PIPELINE COMPLETE — total wall-clock: %.1fs (%.1f min)", total, total / 60)
    logger.info("=" * 70)


if __name__ == "__main__":
    main()
