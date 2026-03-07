#!/usr/bin/env Rscript
# ──────────────────────────────────────────────────────────────────────────────
# 03_plots.R — Visualizations for PGS case/control meta-analysis
#
# Generates:
#   1. Volcano plot  (plots/volcano_plot.pdf + .png)
#   2. Funnel plot   (plots/funnel_plot_all.pdf)
#   3. Individual forest plots for top N models  (plots/forest_{pgs_id}.pdf)
#   4. Summary forest plot of top 50 pooled estimates
#
# Usage:
#   Rscript 03_plots.R [--results-dir results] [--plot-dir plots]
#                       [--top-n-forest 20] [--top-n-summary 50]
# ──────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(metafor)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# ── Parse args ───────────────────────────────────────────────────────────────
args <- commandArgs(trailingOnly = TRUE)
results_dir   <- "results"
plot_dir      <- "plots"
top_n_forest  <- 20
top_n_summary <- 50

parse_arg <- function(flag, default) {
  if (flag %in% args) {
    idx <- which(args == flag)
    return(args[idx + 1])
  }
  return(default)
}

results_dir   <- parse_arg("--results-dir", results_dir)
plot_dir      <- parse_arg("--plot-dir", plot_dir)
top_n_forest  <- as.integer(parse_arg("--top-n-forest", top_n_forest))
top_n_summary <- as.integer(parse_arg("--top-n-summary", top_n_summary))

dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("[PLOTS] Reading from: %s   Writing to: %s\n", results_dir, plot_dir))

# ── Load data ────────────────────────────────────────────────────────────────
meta <- fread(file.path(results_dir, "meta_analysis_results.tsv"))
cat(sprintf("[PLOTS] Loaded %d meta-analysis results\n", nrow(meta)))

# Also load per-cohort results for forest plots
cohort_files <- list.files(results_dir, pattern = "_logistic_results\\.tsv$",
                           full.names = TRUE)
cohort_dat <- rbindlist(lapply(cohort_files, fread, sep = "\t"))
cohort_dat <- cohort_dat[status == "success"]

# ──────────────────────────────────────────────────────────────────────────────
# 1. Volcano plot
# ──────────────────────────────────────────────────────────────────────────────
cat("[PLOTS] Generating volcano plot...\n")

meta[, neg_log10_p := -log10(meta_pvalue)]
meta[, sig_category := ifelse(fdr_qvalue < 0.05, "FDR < 0.05",
                         ifelse(meta_pvalue < 0.05, "Nominal p < 0.05", "NS"))]
meta[, sig_category := factor(sig_category, levels = c("FDR < 0.05", "Nominal p < 0.05", "NS"))]

# Label top 20 by p-value
meta[, rank_p := rank(meta_pvalue, ties.method = "first")]
meta[, label := ifelse(rank_p <= 20, pgs_id, "")]

volcano_colors <- c("FDR < 0.05" = "red", "Nominal p < 0.05" = "orange", "NS" = "grey60")

p_volcano <- ggplot(meta, aes(x = pooled_log_or, y = neg_log10_p, color = sig_category)) +
  geom_point(alpha = 0.7, size = 1.5) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 30,
                  show.legend = FALSE) +
  scale_color_manual(values = volcano_colors, name = "Significance") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  labs(x = "Pooled log(OR)", y = expression(-log[10](p)),
       title = "Volcano plot: PGS association with glioma case/control status") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(plot_dir, "volcano_plot.pdf"), p_volcano, width = 10, height = 7)
ggsave(file.path(plot_dir, "volcano_plot.png"), p_volcano, width = 10, height = 7, dpi = 300)
cat("[PLOTS] Saved volcano_plot.pdf and .png\n")

# ──────────────────────────────────────────────────────────────────────────────
# 2. Funnel plot (all models with ≥2 cohorts)
# ──────────────────────────────────────────────────────────────────────────────
cat("[PLOTS] Generating funnel plot...\n")

multi_cohort <- meta[single_cohort == FALSE & !is.na(pooled_log_or)]

if (nrow(multi_cohort) > 0) {
  pdf(file.path(plot_dir, "funnel_plot_all.pdf"), width = 8, height = 6)
  funnel(x = multi_cohort$pooled_log_or, sei = multi_cohort$pooled_se,
         xlab = "Pooled log(OR)", main = "Funnel plot: All multi-cohort PGS models")
  dev.off()
  cat("[PLOTS] Saved funnel_plot_all.pdf\n")
} else {
  cat("[PLOTS] WARNING: No multi-cohort models for funnel plot\n")
}

# ──────────────────────────────────────────────────────────────────────────────
# 3. Individual forest plots for top N models
# ──────────────────────────────────────────────────────────────────────────────
cat(sprintf("[PLOTS] Generating forest plots for top %d models...\n", top_n_forest))

meta_sorted <- meta[order(meta_pvalue)]
top_models <- head(meta_sorted$pgs_id, top_n_forest)

for (pid in top_models) {
  sub <- cohort_dat[pgs_id == pid]
  if (nrow(sub) == 0) next

  tryCatch({
    safe_pid <- gsub("[^A-Za-z0-9_]", "_", pid)
    pdf_path <- file.path(plot_dir, sprintf("forest_%s.pdf", safe_pid))
    pdf(pdf_path, width = 8, height = max(4, 2 + nrow(sub)))

    if (nrow(sub) >= 2) {
      fit <- rma(yi = sub$log_or, sei = sub$se, method = "REML",
                 slab = sub$cohort)
      forest(fit, xlab = "log(OR)", main = sprintf("Forest plot: %s", pid),
             header = TRUE)
    } else {
      forest(x = sub$log_or, sei = sub$se, slab = sub$cohort,
             xlab = "log(OR)", main = sprintf("Forest plot: %s (single cohort)", pid),
             header = TRUE)
    }

    dev.off()
    cat(sprintf("[PLOTS] Saved %s\n", pdf_path))
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    cat(sprintf("[PLOTS] WARNING: Forest plot failed for %s: %s\n", pid, e$message))
  })
}

# ──────────────────────────────────────────────────────────────────────────────
# 4. Summary forest plot — top 50 pooled estimates
# ──────────────────────────────────────────────────────────────────────────────
cat("[PLOTS] Generating summary forest plot (top 50)...\n")

top_summary <- head(meta_sorted, min(top_n_summary, nrow(meta_sorted)))
top_summary <- top_summary[!is.na(pooled_log_or)]

if (nrow(top_summary) > 0) {
  # Compute OR and CI for display
  top_summary[, or_display := exp(pooled_log_or)]
  top_summary[, or_lower := exp(ci_lower)]
  top_summary[, or_upper := exp(ci_upper)]
  top_summary[, fdr_label := ifelse(fdr_qvalue < 0.05, "*", "")]

  # Order by p-value (best at top)
  top_summary[, pgs_label := paste0(pgs_id, fdr_label)]
  top_summary[, pgs_label := factor(pgs_label, levels = rev(pgs_label))]

  p_summary <- ggplot(top_summary, aes(x = pooled_log_or, y = pgs_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, color = "grey30") +
    geom_point(aes(color = sig_category), size = 2.5) +
    scale_color_manual(values = volcano_colors, name = "Significance", drop = FALSE) +
    labs(x = "Pooled log(OR) with 95% CI",
         y = NULL,
         title = sprintf("Top %d PGS models by meta-analysis p-value", nrow(top_summary)),
         subtitle = "* = FDR < 0.05") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 7),
          legend.position = "bottom")

  ggsave(file.path(plot_dir, "forest_summary_top50.pdf"), p_summary,
         width = 10, height = max(8, nrow(top_summary) * 0.25))
  cat("[PLOTS] Saved forest_summary_top50.pdf\n")
} else {
  cat("[PLOTS] WARNING: No results for summary forest plot\n")
}

cat("[PLOTS] All plotting complete.\n")
