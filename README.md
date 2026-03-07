
# PGS Case/Control Logistic Regression Meta-Analysis Pipeline

## Overview

This pipeline performs a case/control logistic regression meta-analysis across multiple glioma cohorts, testing the association of ~5,000 polygenic scores (PGS) with glioma risk.

**Architecture:** per-cohort logistic regression → result aggregation → random-effects meta-analysis → visualization.

## Files

| File | Description |
|------|-------------|
| `config.py` | Central configuration: paths, cohort definitions, thresholds |
| `utils.py` | Shared utilities: logging, data loading, covariate QC |
| `01_logistic_regression.py` | Per-cohort logistic regression (parallelised with joblib) |
| `02_meta_analysis.R` | Random-effects meta-analysis with `metafor` (REML) |
| `03_plots.R` | Volcano, funnel, and forest plots |
| `04_run_pipeline.py` | Master runner that calls steps 01–03 in sequence |
| `run_pipeline.sh` | SLURM submission script |

## Quick Start

### Test mode (50 models, fast validation)

```bash
# Local (single-threaded)
python 04_run_pipeline.py --test --verbose

# Local (multi-threaded)
python 04_run_pipeline.py --test --verbose --n-jobs 8

# SLURM
sbatch run_pipeline.sh --test
```

### Full run (~5,100 models)

```bash
# Local
python 04_run_pipeline.py --verbose --n-jobs 16

# SLURM
sbatch run_pipeline.sh
```

### Run individual steps

```bash
# Step 1: Logistic regression (all cohorts)
python 01_logistic_regression.py --test --verbose --n-jobs 4

# Step 1: Single cohort
python 01_logistic_regression.py --cohort onco --test --verbose

# Step 2: Meta-analysis
Rscript 02_meta_analysis.R --results-dir results

# Step 3: Plots
Rscript 03_plots.R --results-dir results --plot-dir plots
```

## Data Layout

The pipeline expects the following files (paths configured in `config.py`):

```
{BASE_DIR}/
    onco-covariates.csv
    onco.scores.z-scores.txt.gz
    i370-covariates.csv
    i370.scores.z-scores.txt.gz
    tcga-covariates.csv
    tcga.scores.z-scores.txt.gz
```

### Covariates files

Each CSV must contain: `IID` (sample ID), `case` (0/1), `exclude` (0/1), `age`, `sex` (M/F), `PC1`–`PC8`, and optionally `source`.

### PGS scores files

Comma-separated (gzipped). First column `sample` (sample ID), remaining columns are PGS model IDs.

## Output

```
results/
    onco_logistic_results.tsv    # per-cohort logistic regression results
    i370_logistic_results.tsv
    tcga_logistic_results.tsv
    meta_analysis_results.tsv    # pooled meta-analysis with FDR correction
    *.log                        # log files

plots/
    volcano_plot.pdf / .png      # volcano plot of all PGS models
    funnel_plot_all.pdf           # funnel plot (multi-cohort models only)
    forest_{PGS_ID}.pdf           # individual forest plots (top N)
    forest_summary_top50.pdf      # summary forest of top 50 models
```

## Cohort-specific notes

- **onco**: Includes `source` (AGS/Mayo) as an additional covariate. Age has ~47% missingness and is automatically dropped by the >20% missingness rule.
- **i370**: Single source (AGS), all covariates complete.
- **tcga**: Includes `source` (TCGA/WTCCC) as covariate. WTCCC controls lack age data (~89% missing overall), so age is automatically dropped.
- All cohorts: samples with `exclude == 1` are removed before any analysis.
- Tumor-specific fields (grade, idh, pq, treated) are NOT used as covariates.

## Configuration

Edit `config.py` to change:
- `BASE_DIR` / `OUTPUT_DIR` / `PLOT_DIR` — file paths
- `COHORTS` — add or remove cohorts
- `MIN_CASES` / `MIN_CONTROLS` — minimum sample sizes
- `TEST_N_MODELS` — number of models in test mode
- `FDR_ALPHA` — significance threshold

## Requirements

**Python:** pandas, numpy, statsmodels, joblib  
**R:** metafor, data.table, ggplot2, ggrepel
