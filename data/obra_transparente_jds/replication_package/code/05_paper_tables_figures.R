#' ==============================================================================
#' Script: 05_paper_tables_figures.R
#' Descricao: Gera todas as tabelas e figuras necessarias para o paper
#'
#' Inputs:
#'   - data/processed/did_panel_6periods.rds (do script 01)
#'   - data/processed/did_results_6periods.rds (do script 01)
#'   - data/raw/municipal_covariates.rds (covariáveis municipais)
#'
#' Outputs:
#'   - output/tab_01.rds (Table 1: Summary Statistics)
#'   - output/table2.rds (Table 2: Completion Rates)
#'   - output/did_static.rds (Table 3: DiD Static)
#'   - output/did_event_study.rds (Table 4: Event Study)
#'   - output/obra_transparente.RDS (Data for figures in paper)
#'   - output/tables/*.tex and *.csv versions
#'   - output/figures/*.png
#'
#' Autor: Manoel Galdino
#' Data: 2026-01-13
#' ==============================================================================

# Setup ----
library(tidyverse)
library(here)
library(fixest)
library(modelsummary)
library(kableExtra)

cat("\n")
cat("================================================================\n")
cat("  GERANDO TABELAS E FIGURAS PARA O PAPER\n")
cat("================================================================\n")
cat("Data/Hora:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Criar diretorios necessarios ----
dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(here("output", "figures"), showWarnings = FALSE, recursive = TRUE)

# Carregar dados ----
cat("Carregando dados processados...\n")

if (!file.exists(here("data", "processed", "did_panel_6periods.rds"))) {
  stop("ERRO: Arquivo did_panel_6periods.rds nao encontrado. Execute 01_did_analysis_6periods.R primeiro.")
}

did_panel <- readRDS(here("data", "processed", "did_panel_6periods.rds"))
did_results <- readRDS(here("data", "processed", "did_results_6periods.rds"))

cat("  Painel carregado:", nrow(did_panel), "observacoes\n\n")

# =============================================================================
# TABLE 1: SUMMARY STATISTICS
# =============================================================================

cat("--- Table 1: Summary Statistics ---\n")

# Load municipal covariates (required for Table 1)
if (!file.exists(here("data", "raw", "municipal_covariates.rds"))) {
  stop("ERRO: Arquivo data/raw/municipal_covariates.rds nao encontrado. Necessario para Table 1.")
}

cat("  Usando covariáveis municipais existentes...\n")
data_ot_desc <- readRDS(here("data", "raw", "municipal_covariates.rds"))

# Criar tabela de estatísticas descritivas
cols_stats <- c("idhm", "ext_pobres", "rpc", "ext_prop_pobres",
                "p10_ricos_desigualdade", "freq_escola_15_1",
                "analfabetismo", "populacao")

# Filtrar colunas existentes
cols_stats <- cols_stats[cols_stats %in% names(data_ot_desc)]

label_map <- c(
  idhm = "Human Development Index (HDI)",
  ext_pobres = "Extreme Poverty Rate (%)",
  rpc = "Per Capita Income (R$/month)",
  ext_prop_pobres = "Poverty Rate (%)",
  p10_ricos_desigualdade = "Income Share of Top 10% (%)",
  freq_escola_15_1 = "School Enrollment (Age 15-17) (%)",
  analfabetismo = "Illiteracy Rate (%)",
  populacao = "Population (thousands)"
)

df_stats <- data_ot_desc |>
  group_by(treatment) |>
  summarise(
    across(
      all_of(cols_stats),
      list(
        Mean = ~mean(.x, na.rm = TRUE),
        SD = ~sd(.x, na.rm = TRUE),
        Min = ~min(.x, na.rm = TRUE),
        Max = ~max(.x, na.rm = TRUE)
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  ) |>
  pivot_longer(
    -treatment,
    names_to = c("Variable", "Stat"),
    names_pattern = "(.*)_(.*)"
  ) |>
  pivot_wider(
    names_from = treatment,
    values_from = value,
    names_prefix = "trt"
  ) |>
  rename(Control = trt0, Treated = trt1) |>
  mutate(Variable = recode(Variable, !!!label_map))

descr <- df_stats |>
  pivot_wider(
    names_from = Stat,
    values_from = c(Control, Treated)
  )

# Salvar Table 1
saveRDS(descr, here("output", "tab_01.rds"))

# Também salvar em output/tables/
write_csv(descr, here("output", "tables", "table1_summary_stats.csv"))

cat("  Table 1 salva em output/tab_01.rds\n")

# =============================================================================
# TABLE 2: COMPLETION RATES (kable format for paper)
# =============================================================================

cat("\n--- Table 2: Completion Rates ---\n")

# Calcular taxas de conclusão por período e grupo
completion_rates <- did_panel |>
  group_by(periodo, group_treated) |>
  summarise(
    n = n(),
    concluidas = sum(concluida),
    pct = mean(concluida),
    .groups = "drop"
  ) |>
  mutate(group = ifelse(group_treated == 1, "Treatment", "Control")) |>
  select(periodo, group, pct) |>
  pivot_wider(names_from = group, values_from = pct) |>
  mutate(
    Diff = Treatment - Control,
    DiD = Diff - lag(Diff)
  ) |>
  # Now round for display (Diff and DiD computed from exact values)
  mutate(
    Control = round(Control * 100, 1),
    Treatment = round(Treatment * 100, 1),
    Diff = round(Diff * 100, 1),
    DiD = round(DiD * 100, 1)
  )

# Criar objeto kable para o paper
table2_kable <- completion_rates |>
  kable(
    format = "latex",
    booktabs = TRUE,
    col.names = c("Period", "Control", "Treatment", "Diff", "DiD"),
    caption = "Construction Completion Rates by Period and Treatment Status"
  ) |>
  kable_styling(latex_options = c("hold_position"))

# Salvar Table 2
saveRDS(table2_kable, here("output", "table2.rds"))
write_csv(completion_rates, here("output", "tables", "table2_completion_rates.csv"))

cat("  Table 2 salva em output/table2.rds\n")

# =============================================================================
# DADOS PARA FIGURAS (obra_transparente.RDS)
# =============================================================================

cat("\n--- Preparando dados para figuras ---\n")

# Criar arquivo compatível com o código do paper
obra_transparente <- did_panel |>
  select(id, municipio, uf, concluida, group_treated, periodo) |>
  mutate(
    # Manter compatibilidade com código antigo do paper
    indicator_muni_ot = group_treated
  )

saveRDS(obra_transparente, here("output", "obra_transparente.RDS"))
cat("  Dados salvos em output/obra_transparente.RDS\n")

# =============================================================================
# TABLE 3: DiD STATIC (modelo fixest)
# =============================================================================

cat("\n--- Table 3: DiD Static ---\n")

# Usar modelo do script 01
model_twfe <- did_results$model_twfe

# Salvar modelo para o paper
saveRDS(model_twfe, here("output", "did_static.rds"))

# Também criar tabela formatada
msummary_did <- modelsummary(
  list("TWFE" = model_twfe),
  stars = c('*' = 0.1, '**' = 0.05, '***' = 0.01),
  gof_map = c("nobs", "r.squared"),
  coef_rename = c("treat_post" = "Treatment × Post"),
  output = "data.frame"
)

write_csv(msummary_did, here("output", "tables", "table3_did_static.csv"))

cat("  Table 3 salva em output/did_static.rds\n")

# =============================================================================
# TABLE 4: EVENT STUDY (modelo fixest)
# =============================================================================

cat("\n--- Table 4: Event Study ---\n")

# Usar modelo do script 01
model_es <- did_results$model_es

# Salvar modelo para o paper
saveRDS(model_es, here("output", "did_event_study.rds"))

# Também criar tabela formatada
es_coefs <- did_results$es_coefs
write_csv(es_coefs, here("output", "tables", "table4_event_study.csv"))

cat("  Table 4 salva em output/did_event_study.rds\n")

# =============================================================================
# FIGURAS PARA O PAPER
# =============================================================================

cat("\n--- Gerando Figuras ---\n")

# Figure 1: Completion Trends (Parallel Trends)
cat("  Figure 1: Completion Trends...\n")

# Period labels for x-axis
period_labels <- c("0" = "2015", "1" = "2017", "2" = "Mar\n2018",
                   "3" = "Sep\n2018", "4" = "2019", "5" = "2023")

fig_trends <- did_panel |>
  mutate(group = ifelse(group_treated == 1, "Treatment", "Control")) |>
  group_by(periodo, group) |>
  summarise(taxa_conclusao = mean(concluida), .groups = "drop") |>
  ggplot(aes(x = periodo, y = taxa_conclusao, color = group)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Control" = "#E41A1C", "Treatment" = "#377EB8")) +
  scale_x_continuous(breaks = 0:5, labels = period_labels) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
  labs(
    title = "Construction Completion Rates Over Time",
    subtitle = "Treatment begins at Period 3 (September 2018)",
    x = "Period",
    y = "Completion Rate",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(here("output", "figures", "fig1_completion_trends.png"),
       fig_trends, width = 8, height = 6, dpi = 300)

# Figure 2: Event Study Plot
cat("  Figure 2: Event Study Plot...\n")

es_plot_data <- did_results$es_coefs |>
  mutate(
    ci_lower = ifelse(is.na(se), NA_real_, coef - 1.96 * se),
    ci_upper = ifelse(is.na(se), NA_real_, coef + 1.96 * se)
  )

fig_event_study <- ggplot(es_plot_data, aes(x = periodo_rel, y = coef)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "gray50") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2, color = "#377EB8") +
  geom_point(size = 3, color = "#377EB8") +
  scale_x_continuous(breaks = -3:2, labels = c("-3", "-2", "-1\n(ref)", "0", "+1", "+2")) +
  labs(
    title = "Event Study: Dynamic Treatment Effects",
    subtitle = "Reference period: t = -1 (March 2018). 6 periods: 2015-2023.",
    x = "Periods Relative to Treatment",
    y = "Treatment Effect (Completion Probability)",
    caption = "Note: Error bars show 95% confidence intervals. Clustered SE by municipality."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.caption = element_text(hjust = 0, color = "gray50")
  )

ggsave(here("output", "figures", "fig2_event_study.png"),
       fig_event_study, width = 8, height = 6, dpi = 300)

# Figure 3: Event Study with different style
cat("  Figure 3: Event Study (alternative style)...\n")

fig_es_alt <- ggplot(es_plot_data, aes(x = periodo_rel, y = coef)) +
  geom_hline(yintercept = 0, color = "gray30") +
  geom_vline(xintercept = -0.5, linetype = "dashed", color = "red", alpha = 0.5) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), fill = "#377EB8", alpha = 0.2) +
  geom_line(color = "#377EB8", linewidth = 1) +
  geom_point(size = 4, color = "#377EB8") +
  annotate("text", x = -0.3, y = max(es_plot_data$ci_upper, na.rm = TRUE) * 0.9,
           label = "Treatment\nstarts", hjust = 0, size = 3, color = "red") +
  scale_x_continuous(breaks = -3:2) +
  labs(
    title = "Dynamic Treatment Effects of CSO Monitoring",
    subtitle = "6 periods (2015-2023), 3 pre-treatment periods",
    x = "Time Relative to Treatment",
    y = "Estimated Treatment Effect",
    caption = "Shaded area: 95% CI. Reference: t = -1."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

ggsave(here("output", "figures", "fig3_event_study_alt.png"),
       fig_es_alt, width = 8, height = 6, dpi = 300)

# =============================================================================
# LEAVE-ONE-OUT ROBUSTNESS (drop each treated municipality)
# =============================================================================

cat("\n--- Leave-One-Out Robustness ---\n")

treated_munis <- did_panel |>
  filter(group_treated == 1) |>
  pull(municipio) |>
  unique()

loo_results <- map_dfr(treated_munis, function(m) {
  df_loo <- did_panel |> filter(municipio != m)
  fit_loo <- feols(concluida ~ treat_post | id + periodo,
                   data = df_loo, cluster = ~municipio)
  tibble(
    dropped = m,
    att = coef(fit_loo)["treat_post"],
    se = se(fit_loo)["treat_post"],
    pvalue = pvalue(fit_loo)["treat_post"]
  )
})

cat(sprintf("  Full sample ATT: %.4f\n", coef(model_twfe)["treat_post"]))
cat(sprintf("  LOO range: [%.4f, %.4f]\n", min(loo_results$att), max(loo_results$att)))
cat(sprintf("  All significant at 10%%: %s\n", all(loo_results$pvalue < 0.10)))
cat(sprintf("  All significant at 5%%:  %s\n\n", all(loo_results$pvalue < 0.05)))

saveRDS(loo_results, here("output", "leave_one_out.rds"))
write_csv(loo_results, here("output", "tables", "leave_one_out.csv"))

cat("\n")
cat("================================================================\n")
cat("  TABELAS E FIGURAS GERADAS COM SUCESSO\n")
cat("================================================================\n\n")

cat("Arquivos gerados:\n")
cat("\nPara o paper (output/):\n")
cat("  - tab_01.rds (Table 1: Summary Statistics)\n")
cat("  - table2.rds (Table 2: Completion Rates - kable)\n")
cat("  - did_static.rds (Table 3: DiD Static Model)\n")
cat("  - did_event_study.rds (Table 4: Event Study Model)\n")
cat("  - obra_transparente.RDS (Data for figures)\n")
cat("\nPara output/tables/:\n")
cat("  - table1_summary_stats.csv\n")
cat("  - table2_completion_rates.csv\n")
cat("  - table3_did_static.csv\n")
cat("  - table4_event_study.csv\n")
cat("\nPara output/figures/:\n")
cat("  - fig1_completion_trends.png\n")
cat("  - fig2_event_study.png\n")
cat("  - fig3_event_study_alt.png\n")

cat("\n")
sessionInfo()
