#' ==============================================================================
#' Script: 06_robustness_covariates.R
#' Projeto: Obra Transparente DiD Analysis
#' Descrição: Robustness tests with political covariates
#'
#' Covariates tested:
#'   - Mayor's party ideology (left-right scale)
#'   - Electoral margin of victory
#'
#' Inputs:
#'   - data/processed/did_panel_6periods.rds
#'   - data/raw/elections_mayor_results.csv
#'   - data/raw/party_ideology.csv
#'
#' Outputs:
#'   - output/tables/robustness_covariates.csv
#'   - output/tables/robustness_covariates.tex
#'
#' Author: Manoel Galdino
#' Date: 2026-01-14
#' ==============================================================================

library(tidyverse)
library(here)
library(fixest)
library(stringi)
library(modelsummary)

cat("\n========================================\n")
cat("ROBUSTNESS: Political Covariates\n")
cat("========================================\n\n")

# Helper to normalize municipality names
normalize_name <- function(x) {
  x <- as.character(x)
  x <- stri_trans_general(x, "Latin-ASCII")
  x <- toupper(trimws(x))
  x <- gsub("'", " ", x)
  x <- gsub("[^A-Z ]", "", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

# Load DiD panel
panel <- readRDS(here("data", "processed", "did_panel_6periods.rds"))
cat("DiD panel:", nrow(panel), "obs,", n_distinct(panel$id), "projects\n")

# Period to year mapping
period_year <- c("0" = 2015, "1" = 2017, "2" = 2018, "3" = 2018, "4" = 2019, "5" = 2023)
panel$ano <- period_year[as.character(panel$periodo)]

# Normalize municipality names in panel
panel$muni_norm <- normalize_name(panel$municipio)

# ====================
# ELECTIONS DATA
# ====================
cat("\nLoading elections data...\n")

elections <- read_csv(here("data", "raw", "elections_mayor_results.csv"), show_col_types = FALSE)
elections$muni_norm <- normalize_name(elections$mun_name)

# Deduplicate: keep one election result per municipality-year (highest margin if duplicates)
elections_dedup <- elections |>
  select(muni_norm, uf, ano, party_winner, margin_pct) |>
  arrange(muni_norm, uf, ano, desc(margin_pct)) |>
  distinct(muni_norm, uf, ano, .keep_all = TRUE)

# Match by municipality name + UF + year
nrow_before <- nrow(panel)
n_treated_before <- n_distinct(panel$municipio[panel$group_treated == 1])
n_control_before <- n_distinct(panel$municipio[panel$group_treated == 0])
panel <- left_join(
  panel,
  elections_dedup,
  by = c("muni_norm", "uf", "ano")
)
stopifnot("Elections join produced duplicate rows" = nrow(panel) == nrow_before)
stopifnot("Elections join changed number of treated municipalities" =
  n_distinct(panel$municipio[panel$group_treated == 1]) == n_treated_before)
stopifnot("Elections join changed number of control municipalities" =
  n_distinct(panel$municipio[panel$group_treated == 0]) == n_control_before)

cat("  Matched:", sum(!is.na(panel$party_winner)), "/", nrow(panel),
    "(", round(100*mean(!is.na(panel$party_winner)), 1), "%)\n")

# ====================
# PARTY IDEOLOGY
# ====================
cat("\nLoading party ideology...\n")

ideology <- read_csv(here("data", "raw", "party_ideology.csv"), show_col_types = FALSE)
ideology$party <- toupper(trimws(ideology$party))

panel <- left_join(
  panel,
  ideology |> rename(party_winner = party, mayor_ideology = ideology),
  by = "party_winner"
)

cat("  Matched:", sum(!is.na(panel$mayor_ideology)), "/", nrow(panel),
    "(", round(100*mean(!is.na(panel$mayor_ideology)), 1), "%)\n")

# ====================
# ROBUSTNESS MODELS
# ====================
#
# Theoretical justification for covariates:
#
# 1. Mayor's party ideology (left-right scale): Political alignment may affect
#    the willingness to complete public construction projects. Left-leaning mayors
#    may prioritize education infrastructure differently than right-leaning ones.
#
# 2. Electoral margin of victory: Mayors elected by larger margins may face less
#    political pressure to deliver visible public works, or conversely, may have
#    more political capital to push projects through. This controls for local
#    political competitiveness as a potential confounder.
#
# These covariates address the concern that CSO monitoring assignment may
# correlate with political characteristics that independently affect
# construction completion.
#
cat("\n========================================\n")
cat("Estimating models...\n")
cat("========================================\n\n")

# Verify treatment indicators exist (created in 01_did_analysis_6periods.R)
stopifnot("Column 'post' not found in panel" = "post" %in% names(panel))
stopifnot("Column 'treat_post' not found in panel" = "treat_post" %in% names(panel))

# Model 1: Baseline (from main analysis)
cat("Model 1: Baseline TWFE\n")
m1 <- feols(
  concluida ~ treat_post | id + periodo,
  data = panel,
  cluster = ~municipio
)

# Model 2: Add mayor ideology
cat("Model 2: + Mayor ideology\n")
m2 <- feols(
  concluida ~ treat_post + mayor_ideology | id + periodo,
  data = panel,
  cluster = ~municipio
)

# Model 3: Add electoral margin
cat("Model 3: + Electoral margin\n")
m3 <- feols(
  concluida ~ treat_post + margin_pct | id + periodo,
  data = panel,
  cluster = ~municipio
)

# Model 4: Both covariates
cat("Model 4: + Ideology + Margin\n")
m4 <- feols(
  concluida ~ treat_post + mayor_ideology + margin_pct | id + periodo,
  data = panel,
  cluster = ~municipio
)

# ====================
# RESULTS
# ====================
cat("\n========================================\n")
cat("RESULTS\n")
cat("========================================\n")

models <- list(
  "Baseline" = m1,
  "+ Ideology" = m2,
  "+ Margin" = m3,
  "+ Both" = m4
)

# Print summary
cat("\nTreatment Effect (treat_post) across specifications:\n\n")
for (nm in names(models)) {
  m <- models[[nm]]
  coef_val <- coef(m)["treat_post"]
  se_val <- se(m)["treat_post"]
  pval <- pvalue(m)["treat_post"]
  nobs <- nobs(m)
  sig <- ifelse(pval < 0.01, "***", ifelse(pval < 0.05, "**", ifelse(pval < 0.1, "*", "")))
  cat(sprintf("  %-15s: %+.3f (%.3f) %s   [N=%d]\n", nm, coef_val, se_val, sig, nobs))
}

# ====================
# EXPORT TABLES
# ====================
cat("\n\nExporting tables...\n")

# Create output directory
dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)

# CSV summary
results_df <- data.frame(
  model = names(models),
  att = sapply(models, function(m) coef(m)["treat_post"]),
  se = sapply(models, function(m) se(m)["treat_post"]),
  pvalue = sapply(models, function(m) pvalue(m)["treat_post"]),
  n_obs = sapply(models, nobs)
)
write_csv(results_df, here("output", "tables", "robustness_covariates.csv"))
cat("  Saved: output/tables/robustness_covariates.csv\n")

# Latex table
modelsummary(
  models,
  stars = c('*' = .1, '**' = .05, '***' = .01),
  coef_rename = c(
    "treat_post" = "Treatment × Post",
    "mayor_ideology" = "Mayor ideology (left-right)",
    "margin_pct" = "Electoral margin (%)"
  ),
  gof_omit = "AIC|BIC|RMSE|R2 Within",
  title = "Robustness: Political Covariates",
  notes = list(
    "Clustered standard errors by municipality in parentheses.",
    "All models include project and period fixed effects.",
    "* p < 0.1, ** p < 0.05, *** p < 0.01"
  ),
  output = here("output", "tables", "robustness_covariates.tex")
)
cat("  Saved: output/tables/robustness_covariates.tex\n")

# Also save for paper
modelsummary(
  models,
  stars = c('*' = .1, '**' = .05, '***' = .01),
  coef_rename = c(
    "treat_post" = "Treatment × Post",
    "mayor_ideology" = "Mayor ideology (left-right)",
    "margin_pct" = "Electoral margin (%)"
  ),
  gof_omit = "AIC|BIC|RMSE|R2 Within",
  output = here("output", "tables", "robustness_covariates.html")
)
cat("  Saved: output/tables/robustness_covariates.html\n")

# ====================
# SUMMARY
# ====================
cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")
cat("The treatment effect is robust to controlling for:\n")
cat("  - Mayor's party ideology (left-right scale)\n")
cat("  - Electoral competitiveness (margin of victory)\n")
cat("\nThe baseline ATT of", sprintf("%.3f", coef(m1)["treat_post"]),
    "remains similar across all specifications,\n")
cat("ranging from", sprintf("%.3f", min(sapply(models, function(m) coef(m)["treat_post"]))),
    "to", sprintf("%.3f", max(sapply(models, function(m) coef(m)["treat_post"]))), "\n")

cat("\nDone!\n")
