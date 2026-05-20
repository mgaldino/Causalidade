#' ==============================================================================
#' Script: 99_run_all.R
#' Projeto: Obra Transparente DiD Analysis
#' Description: Master script - runs the complete replication pipeline
#'
#' Usage:
#'   source(here::here("code", "99_run_all.R"))
#'
#' Or from the command line:
#'   Rscript code/99_run_all.R
#'
#' Pipeline:
#'   Step 1. Download SIMEC snapshots from GitHub (00a) - if missing
#'   Step 2. DiD analysis with 6 periods + Event Study (01) - MAIN SPECIFICATION
#'   Step 3. Tables and figures for the paper (05)
#'   Step 4. Robustness: political covariates (06) - ideology, margin
#'   Step 5. Robustness: Wild Cluster Bootstrap + Permutation Test (07)
#'
#' Outputs:
#'   - data/processed/did_panel_6periods.rds
#'   - data/processed/did_results_6periods.rds
#'   - output/*.rds (for the paper)
#'   - output/tables/*.csv, *.tex
#'   - output/figures/*.png
#'
#' Author: Manoel Galdino
#' Date: 2026-02-10
#' ==============================================================================

# Setup ------------------------------------------------------------------------
library(here)

# Record start time
.pipeline_start_time <- Sys.time()

cat("\n")
cat("################################################################################\n")
cat("#                                                                              #\n")
cat("#     OBRA TRANSPARENTE - DiD ANALYSIS REPLICATION PIPELINE                    #\n")
cat("#                                                                              #\n")
cat("################################################################################\n")
cat("\n")
cat("Start:", format(.pipeline_start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Working directory:", here(), "\n")
cat("\n")

# Verify working directory
if (basename(here()) != "did-obra-transparente") {
  stop("ERROR: Wrong working directory. Run from the project root.")
}

# =============================================================================
# CHECK REQUIRED DATA FILES
# =============================================================================

cat("--- Checking required data files ---\n\n")

# Essential data files (all in data/raw/)
required_data <- c(
  "data/raw/simec_2017_05.Rdata",
  "data/raw/simec_2019_08.Rdata",
  "data/raw/treated_works.xlsx",
  "data/raw/escolas_brasil_2015.xlsx",
  "data/raw/simec_2023_10.csv",
  "data/raw/municipal_covariates.rds"
)

# SIMEC snapshot CSVs (can be downloaded if missing)
required_snapshots <- c(
  "data/raw/simec_snapshots/obras_08032018.csv",
  "data/raw/simec_snapshots/obras_upload28092018.csv"
)

# Political covariates for robustness (Step 4)
required_political <- c(
  "data/raw/elections_mayor_results.csv",
  "data/raw/party_ideology.csv"
)

# Data source descriptions for error messages
data_sources <- c(
  "data/raw/simec_2017_05.Rdata" = "SIMEC snapshot May 2017 (from Transparencia Brasil)",
  "data/raw/simec_2019_08.Rdata" = "SIMEC snapshot Aug 2019 (from Transparencia Brasil)",
  "data/raw/treated_works.xlsx" = "List of treated works (from project records)",
  "data/raw/escolas_brasil_2015.xlsx" = "School construction data 2015 (from FNDE/SIMEC)",
  "data/raw/simec_2023_10.csv" = "SIMEC snapshot Oct 2023 (from FNDE/SIMEC)",
  "data/raw/municipal_covariates.rds" = "Municipal covariates (from IBGE/Atlas Brasil)"
)

# Check essential data files
missing_data <- required_data[!file.exists(here(required_data))]
if (length(missing_data) > 0) {
  missing_details <- paste(
    sprintf("  - %s\n    Source: %s", missing_data, data_sources[missing_data]),
    collapse = "\n"
  )
  stop("ERROR: Essential data files missing:\n", missing_details,
       "\n\nSee REPLICATION.md for details on required data files.")
}
cat("Essential data files: OK\n")

# Check snapshots (offer download if missing)
missing_snapshots <- required_snapshots[!file.exists(here(required_snapshots))]
if (length(missing_snapshots) > 0) {
  cat("\nSIMEC snapshots missing. Downloading...\n\n")
  source(here("code", "00a_download_data.R"))

  still_missing <- required_snapshots[!file.exists(here(required_snapshots))]
  if (length(still_missing) > 0) {
    stop("ERROR: Could not download:\n  ", paste(still_missing, collapse = "\n  "))
  }
}
cat("SIMEC snapshots: OK\n")

# Check political covariates
missing_political <- required_political[!file.exists(here(required_political))]
if (length(missing_political) > 0) {
  cat("WARNING: Political covariate files missing. Robustness test (Step 4) will be skipped.\n")
}
cat("\n")

# =============================================================================
# STEP 1: MAIN DiD ANALYSIS (6 periods + Event Study)
# =============================================================================

cat("================================================================================\n")
cat("STEP 1: DiD ANALYSIS (6 PERIODS + EVENT STUDY)\n")
cat("================================================================================\n\n")

source(here("code", "01_did_analysis_6periods.R"))

# =============================================================================
# STEP 2: TABLES AND FIGURES FOR THE PAPER
# =============================================================================

cat("\n================================================================================\n")
cat("STEP 2: TABLES AND FIGURES FOR THE PAPER\n")
cat("================================================================================\n\n")

source(here("code", "05_paper_tables_figures.R"))

# =============================================================================
# STEP 3: ROBUSTNESS - POLITICAL COVARIATES
# =============================================================================

cat("\n================================================================================\n")
cat("STEP 3: ROBUSTNESS - POLITICAL COVARIATES\n")
cat("================================================================================\n\n")

if (all(file.exists(here(required_political)))) {
  source(here("code", "06_robustness_covariates.R"))
} else {
  cat("SKIPPED: Political covariate files not found.\n")
  cat("Missing:", paste(basename(missing_political), collapse = ", "), "\n")
}

# =============================================================================
# STEP 4: ROBUSTNESS - WILD CLUSTER BOOTSTRAP + PERMUTATION TEST
# =============================================================================

cat("\n================================================================================\n")
cat("STEP 4: ROBUSTNESS - WILD CLUSTER BOOTSTRAP + PERMUTATION TEST\n")
cat("================================================================================\n\n")

source(here("code", "07_robustness_wild_bootstrap.R"))

# =============================================================================
# FINAL SUMMARY
# =============================================================================

.pipeline_end_time <- Sys.time()
elapsed <- difftime(.pipeline_end_time, .pipeline_start_time, units = "secs")

cat("\n")
cat("################################################################################\n")
cat("#                                                                              #\n")
cat("#     PIPELINE COMPLETE                                                        #\n")
cat("#                                                                              #\n")
cat("################################################################################\n")
cat("\n")
cat("End:", format(.pipeline_end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total time:", sprintf("%.1f", elapsed), "seconds\n")
cat("\n")

# List generated outputs
cat("OUTPUTS GENERATED:\n")
cat("==================\n\n")

cat("Processed data (data/processed/):\n")
processed_files <- list.files(here("data", "processed"), pattern = "\\.rds$")
for (f in processed_files) {
  size_kb <- round(file.size(here("data", "processed", f)) / 1024, 0)
  cat(sprintf("  - %s (%d KB)\n", f, size_kb))
}

cat("\nPaper inputs (output/):\n")
if (dir.exists(here("output"))) {
  rds_files <- list.files(here("output"), pattern = "\\.rds$|.RDS$")
  for (f in rds_files) {
    cat(sprintf("  - %s\n", f))
  }
}

cat("\nTables (output/tables/):\n")
if (dir.exists(here("output", "tables"))) {
  table_files <- list.files(here("output", "tables"), pattern = "\\.(csv|tex|html)$")
  for (f in table_files) {
    cat(sprintf("  - %s\n", f))
  }
}

cat("\nFigures (output/figures/):\n")
if (dir.exists(here("output", "figures"))) {
  fig_files <- list.files(here("output", "figures"), pattern = "\\.png$")
  for (f in fig_files) {
    cat(sprintf("  - %s\n", f))
  }
}

cat("\n")
cat("================================================================================\n")
cat("KEY RESULTS:\n")
cat("================================================================================\n")

# Load results for summary
results <- readRDS(here("data", "processed", "did_results_6periods.rds"))

cat("\n1. STATIC DiD:\n")
cat(sprintf("   ATT = %.3f (SE = %.3f, p = %.3f)\n",
            coef(results$model_twfe)["treat_post"],
            fixest::se(results$model_twfe)["treat_post"],
            fixest::pvalue(results$model_twfe)["treat_post"]))

cat("\n2. EVENT STUDY (coefficients relative to t=-1):\n")
for (i in 1:nrow(results$es_coefs)) {
  r <- results$es_coefs[i, ]
  if (r$periodo_rel == -1) {
    cat(sprintf("   t=%+d: 0.000 (reference)\n", r$periodo_rel))
  } else {
    cat(sprintf("   t=%+d: %+.3f %s\n", r$periodo_rel, r$coef, r$sig))
  }
}

# Show robustness results if available
robustness_file <- here("output", "tables", "robustness_covariates.csv")
if (file.exists(robustness_file)) {
  cat("\n3. ROBUSTNESS (POLITICAL COVARIATES):\n")
  rob <- read.csv(robustness_file)
  for (i in 1:nrow(rob)) {
    sig <- ifelse(rob$pvalue[i] < 0.01, "***",
                  ifelse(rob$pvalue[i] < 0.05, "**",
                         ifelse(rob$pvalue[i] < 0.1, "*", "")))
    cat(sprintf("   %-15s: %+.3f (%.3f) %s [N=%d]\n",
                rob$model[i], rob$att[i], rob$se[i], sig, rob$n_obs[i]))
  }
}

cat("\n================================================================================\n")
cat("Replication completed successfully!\n")
cat("================================================================================\n\n")
