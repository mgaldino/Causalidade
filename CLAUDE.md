# CLAUDE.md — Curso de Inferência Causal

## Projeto
Livro-texto "Curso de Inferência Causal" (bookdown, português).

## Convenções
- Notação: $D_i$, $Y_{it}(0)$, $Y_{it}(1)$, $\tau^{ATT}$, $\alpha_i + \lambda_t$ para TWFE
- Referências: inline autor-ano no texto, lista manual ao final de cada capítulo (sem .bib)
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
