#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# 02_meta_analysis.R — Random-effects meta-analysis across cohorts
#
# For each PGS model:
#   - If ≥2 cohorts have results → rma(yi, sei, method="REML")
#   - If 1 cohort  → pass-through single-cohort estimate
#   - If 0 cohorts → skip
#
# Applies Benjamini–Hochberg FDR correction across all models.
#
# Usage:
#   Rscript 02_meta_analysis.R [--results-dir results]
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(metafor)
  library(data.table)
})

# ── Parse args ───────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
results_dir <- "results"
if ("--results-dir" %in% args) {
  idx <- which(args == "--results-dir")
  results_dir <- args[idx + 1]
}

cat(sprintf("[META] Results directory: %s\n", results_dir))

# ── Read per-cohort results ──────────────────────────────────────────────────
files <- list.files(results_dir, pattern = "_logistic_results\\.tsv$",
                    full.names = TRUE)
cat(sprintf("[META] Found %d cohort result files:\n", length(files)))
for (f in files) cat(sprintf("       %s\n", f))

dat_list <- lapply(files, fread, sep = "\t")
dat <- rbindlist(dat_list)

# Keep only successful fits
dat <- dat[status == "success"]
cat(sprintf("[META] Total successful model-cohort pairs: %d\n", nrow(dat)))

# ── Meta-analysis per PGS model ─────────────────────────────────────────────
pgs_ids <- unique(dat$pgs_id)
cat(sprintf("[META] Unique PGS models with ≥1 successful fit: %d\n", length(pgs_ids)))

meta_results <- list()

for (pid in pgs_ids) {
  sub <- dat[pgs_id == pid]
  n_cohorts <- nrow(sub)

  row <- list(
    pgs_id        = pid,
    pooled_log_or = NA_real_,
    pooled_se     = NA_real_,
    ci_lower      = NA_real_,
    ci_upper      = NA_real_,
    meta_pvalue   = NA_real_,
    i2            = NA_real_,
    q_stat        = NA_real_,
    q_pvalue      = NA_real_,
    tau2          = NA_real_,
    n_cohorts     = n_cohorts,
    single_cohort = (n_cohorts == 1)
  )

  if (n_cohorts >= 2) {
    tryCatch({
      fit <- rma(yi = sub$log_or, sei = sub$se, method = "REML")
      row$pooled_log_or <- as.numeric(fit$beta)
      row$pooled_se     <- fit$se
      row$ci_lower      <- fit$ci.lb
      row$ci_upper      <- fit$ci.ub
      row$meta_pvalue   <- fit$pval
      row$i2            <- fit$I2
      row$q_stat        <- fit$QE
      row$q_pvalue      <- fit$QEp
      row$tau2          <- fit$tau2
    }, error = function(e) {
      cat(sprintf("[META] WARNING: rma failed for %s: %s\n", pid, e$message))
    })
  } else {
    # Single-cohort pass-through
    row$pooled_log_or <- sub$log_or[1]
    row$pooled_se     <- sub$se[1]
    row$ci_lower      <- sub$log_or[1] - 1.96 * sub$se[1]
    row$ci_upper      <- sub$log_or[1] + 1.96 * sub$se[1]
    row$meta_pvalue   <- sub$pvalue[1]
    row$i2            <- 0
    row$q_stat        <- 0
    row$q_pvalue      <- 1
    row$tau2          <- 0
  }

  meta_results[[length(meta_results) + 1]] <- row
}

meta_dt <- rbindlist(meta_results)

# ── FDR correction ───────────────────────────────────────────────────────────
meta_dt[, fdr_qvalue := p.adjust(meta_pvalue, method = "BH")]

# ── Save ─────────────────────────────────────────────────────────────────────
out_path <- file.path(results_dir, "meta_analysis_results.tsv")
fwrite(meta_dt, out_path, sep = "\t")
cat(sprintf("[META] Saved meta-analysis results to %s\n", out_path))

# ── Summary ──────────────────────────────────────────────────────────────────
cat(sprintf("[META] SUMMARY:\n"))
cat(sprintf("  Total models meta-analysed:  %d\n", nrow(meta_dt)))
cat(sprintf("  Models with ≥2 cohorts:      %d\n", sum(!meta_dt$single_cohort)))
cat(sprintf("  Single-cohort models:        %d\n", sum(meta_dt$single_cohort)))
n_sig <- sum(meta_dt$fdr_qvalue < 0.05, na.rm = TRUE)
cat(sprintf("  Significant at FDR < 0.05:   %d\n", n_sig))
n_nom <- sum(meta_dt$meta_pvalue < 0.05, na.rm = TRUE)
cat(sprintf("  Nominally significant (p<0.05): %d\n", n_nom))
