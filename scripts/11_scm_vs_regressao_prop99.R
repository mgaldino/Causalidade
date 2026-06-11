# Comparando controle sintetico e regressao no exemplo da Proposicao 99
#
# Objetivo:
# Mostrar, de forma aplicada, a afirmacao de Abadie e L'Hour (2021):
# uma regressao que prediz o contrafactual da unidade tratada tambem usa uma
# combinacao linear dos controles. Com intercepto, os pesos implicitos somam 1,
# mas podem ser negativos ou muito grandes, diferentemente dos pesos convexos
# do SCM.
#
# Rode a partir da raiz do repositorio:
# Rscript --vanilla scripts/11_scm_vs_regressao_prop99.R

options(scipen = 999)
set.seed(123)

required_packages <- c("dplyr", "tidyr", "ggplot2", "readr", "tibble", "tidysynth")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Instale os pacotes faltantes antes de rodar o script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)

out_dir <- file.path("outputs", "scm_vs_regressao_prop99")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

treated_state <- "California"
treatment_year <- 1988
pre_years <- 1970:1988
post_years <- 1989:2000
post_year_for_regression <- 2000

data(smoking, package = "tidysynth")
smoking_df <- tibble::as_tibble(smoking)

message("Dados carregados: ", nrow(smoking_df), " observacoes, ", ncol(smoking_df), " variaveis.")
message("Unidades: ", dplyr::n_distinct(smoking_df$state), "; anos: ",
        min(smoking_df$year), "-", max(smoking_df$year), ".")

missing_summary <- smoking_df |>
  dplyr::summarise(
    dplyr::across(
      dplyr::everything(),
      ~ sum(is.na(.x)),
      .names = "n_missing_{.col}"
    )
  ) |>
  tidyr::pivot_longer(
    cols = dplyr::everything(),
    names_to = "variable",
    values_to = "n_missing"
  )

readr::write_csv(
  missing_summary,
  file.path(out_dir, "00_diagnostico_missing.csv")
)

# ---------------------------------------------------------------------------
# 1. SCM canonico com tidysynth, seguindo o capitulo do curso
# ---------------------------------------------------------------------------

scm_fit <- smoking_df |>
  tidysynth::synthetic_control(
    outcome = cigsale,
    unit = state,
    time = year,
    i_unit = treated_state,
    i_time = treatment_year,
    generate_placebos = FALSE
  ) |>
  tidysynth::generate_predictor(
    time_window = 1980:1988,
    ln_income = mean(lnincome, na.rm = TRUE),
    ret_price = mean(retprice, na.rm = TRUE),
    youth = mean(age15to24, na.rm = TRUE)
  ) |>
  tidysynth::generate_predictor(
    time_window = 1984:1988,
    beer_sales = mean(beer, na.rm = TRUE)
  ) |>
  tidysynth::generate_predictor(
    time_window = 1975,
    cigsale_1975 = cigsale
  ) |>
  tidysynth::generate_predictor(
    time_window = 1980,
    cigsale_1980 = cigsale
  ) |>
  tidysynth::generate_predictor(
    time_window = 1988,
    cigsale_1988 = cigsale
  ) |>
  tidysynth::generate_weights(
    optimization_window = pre_years,
    margin_ipop = 0.02,
    sigf_ipop = 7,
    bound_ipop = 6
  ) |>
  tidysynth::generate_control()

scm_weights <- tidysynth::grab_unit_weights(scm_fit) |>
  dplyr::rename(state = unit, scm_weight = weight)

scm_series <- tidysynth::grab_synthetic_control(scm_fit) |>
  dplyr::rename(year = time_unit, actual = real_y, scm_counterfactual = synth_y)

# ---------------------------------------------------------------------------
# 2. Regressao de outcome: o estimador que parece inocente
# ---------------------------------------------------------------------------
#
# Para cada estado do donor pool, construimos preditores pre-tratamento.
# Depois estimamos, apenas nos controles:
#
#   cigsale_2000 ~ cigsale_1970 + ... + cigsale_1988 + covariaveis pre
#
# e predizemos cigsale_2000 para a California.
#
# Com intercepto, a predicao pode ser escrita como:
#
#   yhat_CA,2000 = soma_j h_j * y_j,2000
#
# onde h_j sao os pesos implicitos da regressao:
#
#   h = x_CA' (X'X)^(-1) X'
#
# O ponto central: esses pesos nao aparecem no summary(lm), nao sao esparsos
# e podem ser negativos. Ainda assim, eles definem o contrafactual.

pre_outcome_predictors <- smoking_df |>
  dplyr::filter(year %in% pre_years) |>
  dplyr::select(state, year, cigsale) |>
  tidyr::pivot_wider(
    names_from = year,
    values_from = cigsale,
    names_prefix = "cigsale_"
  )

covariate_predictors <- smoking_df |>
  dplyr::filter(year %in% 1980:1988) |>
  dplyr::group_by(state) |>
  dplyr::summarise(
    ln_income = mean(lnincome, na.rm = TRUE),
    ret_price = mean(retprice, na.rm = TRUE),
    youth = mean(age15to24, na.rm = TRUE),
    .groups = "drop"
  )

beer_predictor <- smoking_df |>
  dplyr::filter(year %in% 1984:1988) |>
  dplyr::group_by(state) |>
  dplyr::summarise(
    beer_sales = mean(beer, na.rm = TRUE),
    .groups = "drop"
  )

post_outcome <- smoking_df |>
  dplyr::filter(year == post_year_for_regression) |>
  dplyr::select(state, y_post = cigsale)

regression_data <- pre_outcome_predictors |>
  dplyr::left_join(covariate_predictors, by = "state") |>
  dplyr::left_join(beer_predictor, by = "state") |>
  dplyr::left_join(post_outcome, by = "state")

if (anyNA(regression_data)) {
  stop("A base da regressao contem missing values. Inspecione regression_data.", call. = FALSE)
}

donor_regression_data <- regression_data |>
  dplyr::filter(state != treated_state)

treated_regression_data <- regression_data |>
  dplyr::filter(state == treated_state)

predictor_cols <- setdiff(names(regression_data), c("state", "y_post"))
regression_formula <- stats::reformulate(predictor_cols, response = "y_post")
x_formula <- stats::reformulate(predictor_cols)

ols_fit <- stats::lm(regression_formula, data = donor_regression_data)

X_donors <- stats::model.matrix(x_formula, donor_regression_data)
x_treated <- stats::model.matrix(x_formula, treated_regression_data)

if (qr(X_donors)$rank < ncol(X_donors)) {
  stop("A matriz X da regressao nao tem posto completo; reduza os preditores.", call. = FALSE)
}

regression_implicit_weights <- as.numeric(
  x_treated %*% solve(crossprod(X_donors), t(X_donors))
)

regression_weights <- tibble::tibble(
  state = donor_regression_data$state,
  regression_weight = regression_implicit_weights
)

regression_counterfactual_2000 <- sum(
  regression_weights$regression_weight * donor_regression_data$y_post
)

lm_predict_2000 <- as.numeric(
  stats::predict(ols_fit, newdata = treated_regression_data)
)

if (!isTRUE(all.equal(regression_counterfactual_2000, lm_predict_2000, tolerance = 1e-8))) {
  stop("A predicao por pesos implicitos nao bate com predict(lm).", call. = FALSE)
}

actual_2000 <- treated_regression_data$y_post

# Como h independe de y, o mesmo vetor de pesos gera a serie contrafactual
# que a regressao produziria se estimasse uma regressao separada para cada ano.
regression_series <- smoking_df |>
  dplyr::filter(state != treated_state) |>
  dplyr::select(state, year, cigsale) |>
  dplyr::left_join(regression_weights, by = "state") |>
  dplyr::group_by(year) |>
  dplyr::summarise(
    regression_counterfactual = sum(regression_weight * cigsale),
    .groups = "drop"
  )

comparison_series <- scm_series |>
  dplyr::left_join(regression_series, by = "year") |>
  dplyr::mutate(
    scm_gap = actual - scm_counterfactual,
    regression_gap = actual - regression_counterfactual
  )

weights_comparison <- scm_weights |>
  dplyr::full_join(regression_weights, by = "state") |>
  dplyr::mutate(
    scm_weight = dplyr::coalesce(scm_weight, 0),
    regression_weight = dplyr::coalesce(regression_weight, 0),
    abs_regression_weight = abs(regression_weight)
  ) |>
  dplyr::arrange(dplyr::desc(abs_regression_weight))

weight_summary <- tibble::tibble(
  method = c("SCM", "Regressao OLS"),
  n_weights = c(nrow(weights_comparison), nrow(weights_comparison)),
  sum_weights = c(sum(weights_comparison$scm_weight), sum(weights_comparison$regression_weight)),
  min_weight = c(min(weights_comparison$scm_weight), min(weights_comparison$regression_weight)),
  max_weight = c(max(weights_comparison$scm_weight), max(weights_comparison$regression_weight)),
  n_negative = c(sum(weights_comparison$scm_weight < -1e-10),
                 sum(weights_comparison$regression_weight < -1e-10)),
  n_greater_than_one = c(sum(weights_comparison$scm_weight > 1 + 1e-10),
                         sum(weights_comparison$regression_weight > 1 + 1e-10)),
  n_effectively_nonzero = c(sum(abs(weights_comparison$scm_weight) > 0.001),
                            sum(abs(weights_comparison$regression_weight) > 0.001)),
  l1_norm = c(sum(abs(weights_comparison$scm_weight)),
              sum(abs(weights_comparison$regression_weight)))
)

counterfactual_summary <- comparison_series |>
  dplyr::filter(year == post_year_for_regression) |>
  dplyr::transmute(
    year,
    actual_california = actual,
    scm_counterfactual,
    regression_counterfactual,
    scm_gap = actual - scm_counterfactual,
    regression_gap = actual - regression_counterfactual
  )

readr::write_csv(weights_comparison, file.path(out_dir, "01_pesos_scm_vs_regressao.csv"))
readr::write_csv(weight_summary, file.path(out_dir, "02_resumo_pesos.csv"))
readr::write_csv(counterfactual_summary, file.path(out_dir, "03_contrafactual_2000.csv"))
readr::write_csv(comparison_series, file.path(out_dir, "04_series_contrafactuais.csv"))
saveRDS(ols_fit, file.path(out_dir, "modelo_regressao_ols.rds"))

# ---------------------------------------------------------------------------
# 3. Figuras
# ---------------------------------------------------------------------------

weights_long <- weights_comparison |>
  dplyr::select(state, scm_weight, regression_weight) |>
  tidyr::pivot_longer(
    cols = c(scm_weight, regression_weight),
    names_to = "method",
    values_to = "weight"
  ) |>
  dplyr::mutate(
    method = dplyr::recode(
      method,
      scm_weight = "SCM: pesos convexos e explicitos",
      regression_weight = "Regressao: pesos implicitos"
    ),
    state = factor(state, levels = weights_comparison$state)
  )

figure_weights <- ggplot(weights_long, aes(x = state, y = weight, fill = method)) +
  geom_col(show.legend = FALSE) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey30") +
  coord_flip() +
  facet_wrap(~ method, ncol = 2, scales = "free_x") +
  labs(
    title = "Figura 1. Pesos do controle sintético e pesos implícitos da regressão",
    subtitle = "A regressão usa pesos que somam 1, mas muitos são negativos; o SCM restringe os pesos ao simplex.",
    x = NULL,
    y = "Peso"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "none",
    panel.grid.major.y = element_blank()
  )

ggsave(
  filename = file.path(out_dir, "figura_1_pesos_scm_vs_regressao.png"),
  plot = figure_weights,
  width = 11,
  height = 8,
  dpi = 300
)

series_long <- comparison_series |>
  dplyr::select(year, actual, scm_counterfactual, regression_counterfactual) |>
  tidyr::pivot_longer(
    cols = c(actual, scm_counterfactual, regression_counterfactual),
    names_to = "series",
    values_to = "cigsale"
  ) |>
  dplyr::mutate(
    series = dplyr::recode(
      series,
      actual = "California observada",
      scm_counterfactual = "California sintética (SCM)",
      regression_counterfactual = "Contrafactual por regressão"
    )
  )

figure_series <- ggplot(series_long, aes(x = year, y = cigsale, colour = series, linetype = series)) +
  geom_line(linewidth = 0.9) +
  geom_vline(xintercept = treatment_year, linewidth = 0.5, linetype = "dotted") +
  scale_colour_manual(
    values = c(
      "California observada" = "#b2182b",
      "California sintética (SCM)" = "#2166ac",
      "Contrafactual por regressão" = "#4d4d4d"
    )
  ) +
  labs(
    title = "Figura 2. Contrafactual da Califórnia: SCM versus regressão",
    subtitle = "Linha vertical: Proposição 99 em 1988.",
    x = "Ano",
    y = "Vendas per capita de cigarros",
    colour = NULL,
    linetype = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(out_dir, "figura_2_contrafactual_scm_vs_regressao.png"),
  plot = figure_series,
  width = 10,
  height = 6,
  dpi = 300
)

gap_long <- comparison_series |>
  dplyr::filter(year %in% c(pre_years, post_years)) |>
  dplyr::select(year, scm_gap, regression_gap) |>
  tidyr::pivot_longer(
    cols = c(scm_gap, regression_gap),
    names_to = "method",
    values_to = "gap"
  ) |>
  dplyr::mutate(
    method = dplyr::recode(
      method,
      scm_gap = "SCM",
      regression_gap = "Regressão"
    )
  )

figure_gap <- ggplot(gap_long, aes(x = year, y = gap, colour = method)) +
  geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey40") +
  geom_vline(xintercept = treatment_year, linewidth = 0.5, linetype = "dotted") +
  geom_line(linewidth = 0.9) +
  labs(
    title = "Figura 3. Gap estimado: observado menos contrafactual",
    subtitle = "Gaps positivos indicam que a Califórnia observada está acima do contrafactual; negativos, abaixo.",
    x = "Ano",
    y = "Gap em vendas per capita",
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(out_dir, "figura_3_gap_scm_vs_regressao.png"),
  plot = figure_gap,
  width = 10,
  height = 6,
  dpi = 300
)

# ---------------------------------------------------------------------------
# 4. Mini-simulacao: regressao pode produzir peso maior que 1
# ---------------------------------------------------------------------------
#
# No exemplo da California, os pesos da regressao ficam negativos, mas nenhum
# passa de 1 nesta especificacao. O caso abaixo mostra a geometria da
# extrapolacao: se a unidade tratada esta fora do suporte dos controles em x,
# a regressao com intercepto prediz por uma combinacao linear que soma 1, mas
# usa um peso negativo e outro maior que 1.

toy_controls <- tibble::tibble(
  unit = c("Controle A", "Controle B"),
  x = c(0, 1),
  y = c(10, 20)
)

toy_treated <- tibble::tibble(x = 2)

X_toy <- stats::model.matrix(~ x, toy_controls)
x_toy <- stats::model.matrix(~ x, toy_treated)

toy_weights <- as.numeric(x_toy %*% solve(crossprod(X_toy), t(X_toy)))

toy_output <- tibble::tibble(
  unit = toy_controls$unit,
  regression_weight = toy_weights,
  y = toy_controls$y,
  contribution = regression_weight * y
)

readr::write_csv(toy_output, file.path(out_dir, "05_mini_simulacao_pesos_maiores_que_um.csv"))

# ---------------------------------------------------------------------------
# 5. Resumo no console
# ---------------------------------------------------------------------------

message("\nResumo dos pesos:")
print(weight_summary)

message("\nContrafactual em ", post_year_for_regression, ":")
print(counterfactual_summary)

message("\nMaiores pesos implicitos da regressao em valor absoluto:")
print(
  weights_comparison |>
    dplyr::select(state, scm_weight, regression_weight, abs_regression_weight) |>
    dplyr::slice_head(n = 10)
)

message("\nMini-simulacao de extrapolacao:")
print(toy_output)
message("Soma dos pesos da mini-simulacao: ", round(sum(toy_output$regression_weight), 6))

message("\nArquivos salvos em: ", normalizePath(out_dir))
