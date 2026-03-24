# CLAUDE.md — Curso de Inferência Causal

## Projeto
Livro-texto "Curso de Inferência Causal" (bookdown, português).

## Convenções
- Notação: $D_i$, $Y_{it}(0)$, $Y_{it}(1)$, $\tau^{ATT}$, $\alpha_i + \lambda_t$ para TWFE
- Referências: `references.bib` + citações `@key` no texto (pandoc-citeproc). Caps antigos ainda têm listas manuais — migrar ao editar.
- Cap 08 (DiD) é referência de estilo
- Render: `rmarkdown::render("filename.Rmd", output_format="html_document")`
- R: usar `/usr/local/bin/Rscript` com script files (não `-e` inline)

## TODOs

### Cap 11 — Controle Sintético e Estimação Contrafactual

- [ ] **Ler e incorporar Augmented SCM (Ben-Michael, Feller & Rothstein, 2021)**
  - Paper baixado: `ben-michael_augmented_scm_2021.pdf` (no repo)
  - Adicionar 1-2 parágrafos entre SCM clássico (seção 3) e framework contrafactual (seção 4)
  - Explicar bias correction e conceito de dupla robustez como precursor do TROP
  - Referência: Ben-Michael, E., Feller, A., & Rothstein, J. (2021). The augmented synthetic control method. *JASA*, 116(536), 1789–1803.

- [ ] **Incorporar Chernozhukov, Wüthrich & Zhu (2021) mais adequadamente**
  - Atualmente mencionado apenas de passagem (1 frase)
  - Expandir: explicar inferência conformal como alternativa à permutação de unidades
  - Lógica: permutação de períodos (não unidades) → validade sem aleatorização
  - Referência: Chernozhukov, V., Wüthrich, K., & Zhu, Y. (2021). An exact and robust conformal inference method for counterfactual and synthetic controls. *JASA*, 116(536), 1849–1864.

### Cap 08 — DiD: Covariáveis (REVISAR — ler os papers antes de finalizar)

- [ ] **Ler Caetano, Callaway, Payne & Sant'Anna (2024) — "DiD with Time-Varying Covariates"**
  - Paper no repo: baixar de arxiv 2202.02903v3
  - 5 razões pelas quais TWFE + time-varying covariates falha
  - Propõe estimandos doubly-robust e estratégias de imputação
  - Verificar se a seção adicionada ao cap 08 está alinhada com o paper

- [ ] **Ler Lin & Zhang (2022, Economics Letters) — "Interpreting coefficients in dynamic TWFE with time-varying covariates"**
  - Paper no repo: `1-s2.0-S0165176522001823-main.pdf`
  - Foco: "covariate effect bias" — viés que persiste mesmo com efeitos homogêneos e 2 períodos
  - Complementa Caetano et al. — entender a relação entre os dois papers

- [ ] **Ler Abadie (2005) — "Semiparametric DiD Estimators" (RES)**
  - Foundational para PTA condicional com IPW
  - Citado nos slides do Borusyak mas não explicado no cap 08

- [ ] **Verificar de Chaisemartin & D'Haultfœuille (2023) sobre unit-specific trends**
  - Slides Borusyak (D4, slides 6-8): estimadores apropriados lidam com state-specific trends
  - Exemplo Friedberg (1998) vs Wolfers (2006)
