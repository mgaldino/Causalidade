#' ==============================================================================
#' Script: 07_robustness_wild_bootstrap.R
#' Projeto: Obra Transparente DiD Analysis
#' Descrição: Few-cluster inference robustness checks
#'
#' Motivation:
#'   With only 21 treated municipalities, cluster-robust standard errors may be
#'   unreliable (Cameron, Gelbach, and Miller, 2008). This script implements:
#'     1. Wild Cluster Bootstrap via fwildclusterboot (Roodman et al., 2019)
#'     2. Permutation / Randomization Inference test
#'
#' Inputs:
#'   - data/processed/did_panel_6periods.rds
#'
#' Outputs:
#'   - output/tables/wild_bootstrap_results.csv
#'   - output/wild_bootstrap.rds
#'
#' Author: Manoel Galdino
#' Date: 2026-02-11
#' ==============================================================================

library(tidyverse)
library(here)
library(fixest)
library(fwildclusterboot)

set.seed(20260210)

cat("\n========================================\n")
cat("ROBUSTNESS: Few-Cluster Inference\n")
cat("Cameron, Gelbach, and Miller (2008)\n")
cat("========================================\n\n")

# Load data
panel <- readRDS(here("data", "processed", "did_panel_6periods.rds"))

# Verify treatment variables exist (created in 01_did_analysis_6periods.R)
stopifnot("Column 'post' not found in panel" = "post" %in% names(panel))
stopifnot("Column 'treat_post' not found in panel" = "treat_post" %in% names(panel))
stopifnot("Column 'rel_time' not found in panel" = "rel_time" %in% names(panel))

municipalities <- panel |>
  select(municipio, group_treated) |>
  distinct()

n_treated <- sum(municipalities$group_treated == 1)
n_total <- nrow(municipalities)

cat("Panel:", nrow(panel), "obs,", n_total, "municipalities\n")
cat("  Treated:", n_treated, "municipalities\n")
cat("  Control:", n_total - n_treated, "municipalities\n\n")

# =============================================================================
# METHOD 1: WILD CLUSTER BOOTSTRAP (fwildclusterboot)
# =============================================================================

cat("=== METHOD 1: Wild Cluster Bootstrap ===\n")
cat("(Project-level data, Rademacher weights, fwildclusterboot)\n\n")

# Ensure event study dummies exist
panel <- panel |>
  mutate(
    treat_m3 = group_treated * (rel_time == -3),
    treat_m2 = group_treated * (rel_time == -2),
    treat_0  = group_treated * (rel_time == 0),
    treat_p1 = group_treated * (rel_time == 1),
    treat_p2 = group_treated * (rel_time == 2)
  )

# --- Static DiD ---
cat("Static DiD:\n")
fit_static <- feols(concluida ~ treat_post | id + periodo,
                    data = panel, cluster = ~municipio)

boot_static <- boottest(fit_static, param = "treat_post",
                        clustid = c("municipio"), B = 9999,
                        type = "rademacher")

boot_static_tidy <- tidy(boot_static)

cat(sprintf("  Coefficient: %.4f\n", coef(fit_static)["treat_post"]))
cat(sprintf("  Cluster SE:  %.4f\n", se(fit_static)["treat_post"]))
cat(sprintf("  p(cluster):  %.4f\n", pvalue(fit_static)["treat_post"]))
cat(sprintf("  p(wild):     %.4f\n", boot_static_tidy$p.value))
cat(sprintf("  Wild CI 95%%: [%.4f, %.4f]\n\n",
            boot_static_tidy$conf.low, boot_static_tidy$conf.high))

# --- Event Study ---
cat("Event Study: Wild bootstrap for each coefficient...\n\n")

fit_es <- feols(
  concluida ~ treat_m3 + treat_m2 + treat_0 + treat_p1 + treat_p2 | id + periodo,
  data = panel, cluster = ~municipio
)

es_coefs <- c("treat_m3", "treat_m2", "treat_0", "treat_p1", "treat_p2")
es_labels <- c("t=-3 (2015)", "t=-2 (2017)", "t=0 (Set/18)", "t=+1 (2019)", "t=+2 (2023)")
es_results <- list()

for (i in seq_along(es_coefs)) {
  coef_name <- es_coefs[i]
  cat(sprintf("  %s (%s)... ", coef_name, es_labels[i]))

  boot_result <- boottest(fit_es, param = coef_name,
                          clustid = c("municipio"), B = 9999,
                          type = "rademacher")

  es_results[[coef_name]] <- tidy(boot_result)
  cat(sprintf("coef=%+.4f, p(cluster)=%.4f, p(wild)=%.4f\n",
              coef(fit_es)[coef_name],
              pvalue(fit_es)[coef_name],
              es_results[[coef_name]]$p.value))
}

# =============================================================================
# METHOD 2: PERMUTATION / RANDOMIZATION INFERENCE
# =============================================================================

cat("\n\n=== METHOD 2: Randomization Inference ===\n")
cat("(Project-level data with original specification)\n\n")

cat("Randomly reassigning treatment to municipalities (keeping N_treated =", n_treated, ")...\n")

# Observed ATT from project-level model (main specification)
att_observed <- coef(fit_static)["treat_post"]
cat(sprintf("  Observed ATT: %.4f\n", att_observed))

B_perm <- 4999
att_perm <- numeric(B_perm)

for (b in seq_len(B_perm)) {
  perm_treated <- sample(municipalities$municipio, size = n_treated, replace = FALSE)

  panel$group_perm <- as.numeric(panel$municipio %in% perm_treated)
  panel$treat_post_perm <- panel$group_perm * panel$post

  fit_perm <- tryCatch(
    feols(concluida ~ treat_post_perm | id + periodo, data = panel, cluster = ~municipio),
    error = function(e) NULL
  )

  if (!is.null(fit_perm) && "treat_post_perm" %in% names(coef(fit_perm))) {
    att_perm[b] <- coef(fit_perm)["treat_post_perm"]
  } else {
    att_perm[b] <- NA
  }

  if (b %% 1000 == 0) cat(sprintf("  ... %d/%d permutations done\n", b, B_perm))
}

att_perm_valid <- att_perm[!is.na(att_perm)]
p_perm <- mean(abs(att_perm_valid) >= abs(att_observed))

cat(sprintf("\n  Permutation p-value (two-sided): %.4f\n", p_perm))
cat(sprintf("  Effective permutations: %d\n", length(att_perm_valid)))
cat(sprintf("  Observed ATT rank: %d/%d\n",
            sum(abs(att_perm_valid) >= abs(att_observed)),
            length(att_perm_valid)))

# Clean up temp columns
panel$group_perm <- NULL
panel$treat_post_perm <- NULL

# =============================================================================
# SUMMARY TABLE
# =============================================================================

cat("\n\n========================================\n")
cat("SUMMARY: Few-Cluster Inference Results\n")
cat("========================================\n\n")

cat(sprintf("%-15s | %8s | %8s | %10s | %10s\n",
            "Coefficient", "Estimate", "SE", "p(cluster)", "p(wild)"))
cat(sprintf("%-15s-|-%8s-|-%8s-|-%10s-|-%10s\n",
            "---------------", "--------", "--------", "----------", "----------"))

cat(sprintf("%-15s | %+8.4f | %8.4f | %10.4f | %10.4f\n",
            "Static ATT",
            coef(fit_static)["treat_post"],
            se(fit_static)["treat_post"],
            pvalue(fit_static)["treat_post"],
            boot_static_tidy$p.value))

for (i in seq_along(es_coefs)) {
  r <- es_results[[es_coefs[i]]]
  cat(sprintf("%-15s | %+8.4f | %8.4f | %10.4f | %10.4f\n",
              es_labels[i],
              coef(fit_es)[es_coefs[i]],
              se(fit_es)[es_coefs[i]],
              pvalue(fit_es)[es_coefs[i]],
              r$p.value))
}

cat(sprintf("\nPermutation test (static ATT, project-level): p = %.4f\n", p_perm))

# =============================================================================
# SAVE RESULTS
# =============================================================================

results_df <- data.frame(
  coefficient = c("Static ATT", es_labels),
  estimate = c(coef(fit_static)["treat_post"],
               sapply(es_coefs, function(x) coef(fit_es)[x])),
  se_cluster = c(se(fit_static)["treat_post"],
                 sapply(es_coefs, function(x) se(fit_es)[x])),
  p_cluster = c(pvalue(fit_static)["treat_post"],
                sapply(es_coefs, function(x) pvalue(fit_es)[x])),
  p_wild_bootstrap = c(boot_static_tidy$p.value,
                       sapply(es_results, function(r) r$p.value)),
  p_permutation = c(p_perm, rep(NA, length(es_coefs)))
)

dir.create(here("output", "tables"), recursive = TRUE, showWarnings = FALSE)
write_csv(results_df, here("output", "tables", "wild_bootstrap_results.csv"))
cat("\nSaved: output/tables/wild_bootstrap_results.csv\n")

saveRDS(
  list(
    fit_static = fit_static,
    boot_static = boot_static,
    fit_es = fit_es,
    event_study_boot = es_results,
    permutation = list(att_observed = att_observed, att_perm = att_perm_valid, p_perm = p_perm),
    summary = results_df
  ),
  here("output", "wild_bootstrap.rds")
)
cat("Saved: output/wild_bootstrap.rds\n")

cat("\nDone!\n")
sessionInfo()
