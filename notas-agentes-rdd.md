# Notas dos agentes sobre o capítulo de RDD

Data da revisão: 14 de maio de 2026  
Capítulo revisado: `07-RDD.Rmd`  
Base de referência: commit `545fe65`, após renderização completa do livro em `docs/`

## Síntese executiva

O capítulo ficou mais forte depois da revisão: agora separa melhor o desenho sharp do fuzzy, conecta fuzzy RD à aula de variáveis instrumentais, apresenta os frameworks de continuidade e aleatorização local, e explicita que bandwidths automáticos de `rdrobust` pertencem ao framework de continuidade. Essa distinção é essencial para evitar uma interpretação equivocada de "experimento local" em janelas grandes.

Os agentes convergiram em quatro mensagens centrais:

1. A identificação principal em RD vem da continuidade dos resultados potenciais no cutoff, não de continuidade literal da running variable.
2. Fuzzy RD deve ser ensinado como IV local: forma reduzida dividida pela primeira etapa, com interpretação LATE para compliers no cutoff.
3. Inferência moderna deve priorizar `rdrobust`, bandwidth MSE/CER-optimal e intervalos robust bias-corrected.
4. Diagnósticos ajudam a avaliar plausibilidade, mas não provam identificação; covariáveis e ML aumentam precisão, não salvam desenhos mal identificados.

## Agente independente: RDD e econometria

### Parecer técnico

O agente técnico avaliou que a versão inicial do capítulo cobria os tópicos importantes, mas misturava identificação, estimação, inferência e diagnóstico. A recomendação principal foi reorganizar o núcleo do capítulo para que os alunos consigam distinguir:

- o estimando causal;
- a suposição de identificação;
- o estimador local;
- a inferência;
- os diagnósticos de plausibilidade.

### Pontos críticos identificados

- A condição central não é "continuidade de X", mas continuidade dos resultados potenciais condicionais no cutoff.
- A condição de Lee deve ser formulada como ausência de manipulação precisa da running variable.
- Covariáveis balanceadas são diagnóstico, não prova de validade.
- A fórmula sharp RD deve ser apresentada como diferença entre limites laterais.
- A running variable deve ser centralizada em `X_i - c` quando se escreve uma regressão linear local.
- RD não é "falta de overlap" no sentido de matching; é um problema de estimação de fronteira.
- O trade-off correto é: janelas menores reduzem viés e aumentam variância; janelas maiores aumentam precisão, mas dependem mais de suavidade e forma funcional.
- Fuzzy RD precisava de uma apresentação mais clara como Wald local.
- O resultado não precisa saltar no cutoff; quem precisa saltar é o tratamento ou a probabilidade de tratamento.
- Bins são ferramenta visual, não método de estimação.
- `rdrobust` deve ser ensinado com destaque para inferência robusta com correção de viés.

### Status após implementação

Implementado no capítulo:

- Fórmula sharp RD por limites laterais.
- Fuzzy RD como Wald local e IV local.
- Centralização da running variable em `X_i - c`.
- Correção do trade-off viés-variância.
- Separação entre frameworks de continuidade e aleatorização local.
- Alerta sobre bandwidth MSE/CER-optimal não justificar automaticamente local randomization.
- Uso de `rdrobust` com MSE/CER-optimal e interpretação de intervalos robust bias-corrected.
- Nota sobre manipulação e testes de equivalência.
- Nota sobre covariáveis como precisão/diagnóstico, não identificação.
- Seções curtas sobre running variable discreta, clusters, RD dinâmico, erro de medida e RD espacial.

Ainda pode melhorar em próxima rodada:

- Dar mais exemplos aplicados reais de fuzzy RD em ciência política.
- Expandir inferência com clusters quando houver aplicações com municípios, escolas, eleições repetidas ou distritos.
- Reduzir densidade de algumas simulações para melhorar ritmo de aula.
- Melhorar rótulos de eixos em alguns gráficos gerados automaticamente por `rdplot`.

## Agente independente: revisão de literatura 2022-2026

### Parecer bibliográfico

O agente de lit review recomendou organizar os desenvolvimentos recentes em três blocos:

- Fundamentos modernos de RD: continuidade versus aleatorização local.
- Inferência e diagnósticos: fuzzy RD, manipulação, covariáveis, discrete scores, clusters.
- Extensões: RD dinâmico, boundary/geographic RD, múltiplos cutoffs/scores e erro de medida.

### Trabalhos prioritários incorporados

- Cattaneo e Titiunik (2022) e Cattaneo, Idrobo e Titiunik (2024): usados para organizar a aula em dois frameworks, continuidade e aleatorização local.
- Noack e Rothe (2024): incorporado na seção de fuzzy RD como alerta sobre inferência bias-aware, primeira etapa fraca, score discreto e donut designs.
- Hartman (2021) e Fitzgerald (2025): incorporados na seção de McCrary/rddensity para deixar claro que falhar em rejeitar manipulação não prova ausência de manipulação.
- Kreiß e Rothe (2023) e Noack, Olma e Rothe (2025): incorporados na seção de covariáveis, com a mensagem de que covariáveis/ML servem para precisão e diagnóstico, não para salvar identificação.
- Hsu e Shen (2024): incorporado na seção de RD dinâmico ou repetido.
- Eckles, Ignatiadis, Wager e Wu (2025): incorporado na seção de erro de medida e noise-induced randomization.
- Cattaneo, Titiunik e Yu (2025): incorporado na seção de múltiplos cutoffs, boundary/geographic RD e múltiplas dimensões.

### Trabalhos recomendados para uma próxima versão

- Pei, Lee, Card e Weber (2022), sobre ordem local polinomial em RD.
- Fusejima, Ishihara e Sawada (2025), sobre testes unificados para RD e o problema de múltiplos diagnósticos.
- Noack, Olma e Rothe (2026 WP), sobre inferência em RD com dados clusterizados.
- Litschwartz (2022 WP), como alerta sobre local randomization quando a running variable é score de teste.
- Cattaneo, Keele e Titiunik (2023), especialmente para aplicações com running variables discretas em saúde e políticas públicas.

## Notas da revisão pedagógica

### Avaliação

Nota pedagógica: B+/A-

O capítulo agora tem uma sequência didática mais defensável para uma aula posterior a IV:

1. Sharp RD como salto determinístico no tratamento.
2. Fuzzy RD como IV local.
3. Identificação por continuidade dos resultados potenciais.
4. Estimação local com `rdrobust`.
5. Bandwidth MSE/CER-optimal e inferência robust bias-corrected.
6. Diagnósticos de manipulação, covariáveis e placebos.
7. Casos difíceis e extensões recentes.

### Pontos fortes pedagógicos

- A ponte com IV no fuzzy RD está clara e ajuda a aula a se conectar com o capítulo anterior.
- A distinção entre continuidade e aleatorização local evita uma simplificação comum em aulas de RD.
- O alerta sobre janelas grandes em close elections está bem colocado depois da discussão de MSE/CER-optimal.
- O capítulo está mais aplicado, com exemplos em R ao longo da exposição.
- O checklist final ajuda os alunos a transferirem a aula para leitura e avaliação de papers.

### Riscos pedagógicos remanescentes

- A sequência de simulações e gráficos pode ficar densa em uma aula só.
- Algumas extensões recentes aparecem como notas curtas; isso é adequado para uma aula aplicada, mas exigiria mais desenvolvimento em uma aula avançada.
- O professor deve enfatizar oralmente que diagnósticos são evidência auxiliar, não testes definitivos da suposição causal.

## Notas da auditoria visual

### Avaliação

Nota visual: B

A renderização final não apresenta warnings visíveis no HTML do capítulo. Os warnings que apareciam em gráficos e chamadas de pacotes foram removidos ou neutralizados. O DAG de PCRD foi substituído por uma versão mais legível.

### Pontos fortes visuais

- O gráfico de sensibilidade a bandwidth e critério de seleção é mais informativo do que uma sequência de outputs separados.
- O DAG de PCRD ficou legível no HTML, sem labels cortados.
- O capítulo separa melhor gráficos, outputs e interpretação textual.
- A saída renderizada ficou utilizável para leitura online.

### Problemas visuais remanescentes

- Alguns gráficos produzidos por `rdplot` ainda carregam rótulos genéricos de eixo, como "X axis" e "Y axis".
- Alguns outputs longos de pacotes econométricos continuam visualmente densos para leitura em HTML.
- Se o capítulo for usado como material principal de aula, vale transformar alguns outputs em tabelas resumidas ou gráficos próprios.

## Referências completas dos trabalhos recentes

Cattaneo, M. D., & Titiunik, R. (2022). Regression discontinuity designs. *Annual Review of Economics*, 14, 821-851. https://doi.org/10.1146/annurev-economics-051520-021409

Cattaneo, M. D., Idrobo, N., & Titiunik, R. (2024). *A Practical Introduction to Regression Discontinuity Designs: Extensions*. Cambridge University Press. https://doi.org/10.1017/9781009441896

Noack, C., & Rothe, C. (2024). Bias-aware inference in fuzzy regression discontinuity designs. *Econometrica*, 92(3), 687-711. https://doi.org/10.3982/ECTA19466

Fitzgerald, J. (2025). Manipulation tests in regression discontinuity design: The need for equivalence testing. *MetaArXiv*. https://doi.org/10.31219/osf.io/2dgrp_v1

Hartman, E. (2021). Equivalence testing for regression discontinuity designs. *Political Analysis*, 29(4), 505-521. https://doi.org/10.1017/pan.2020.43

Kreiß, A., & Rothe, C. (2023). Inference in regression discontinuity designs with high-dimensional covariates. *The Econometrics Journal*, 26(2), 105-123. https://doi.org/10.1093/ectj/utac029

Noack, C., Olma, T., & Rothe, C. (2025). Flexible covariate adjustments in regression discontinuity designs. Working paper, revised April 2025. https://arxiv.org/abs/2107.07942

Hsu, Y.-C., & Shen, S. (2024). Dynamic regression discontinuity under treatment effect heterogeneity. *Quantitative Economics*, 15(4), 1035-1064. https://doi.org/10.3982/QE2150

Eckles, D., Ignatiadis, N., Wager, S., & Wu, H. (2025). Noise-induced randomization in regression discontinuity designs. *Biometrika*, 112(2), asaf003. https://doi.org/10.1093/biomet/asaf003

Cattaneo, M. D., Titiunik, R., & Yu, R. R. (2025a). Estimation and inference in boundary discontinuity designs. arXiv:2505.05670. https://arxiv.org/abs/2505.05670

Cattaneo, M. D., Titiunik, R., & Yu, R. R. (2025b). rd2d: Causal inference in boundary discontinuity designs. arXiv:2505.07989. https://arxiv.org/abs/2505.07989

## Próximas melhorias recomendadas

1. Adicionar uma aplicação empírica curta de fuzzy RD em R, preferencialmente com primeira etapa visual e Wald local.
2. Customizar rótulos de gráficos `rdplot` para remover labels genéricos.
3. Incluir uma subseção curta sobre discrete running variables com um exemplo em R.
4. Expandir a discussão de clusters quando houver aplicações com unidades administrativas.
5. Transformar os outputs longos de `rdrobust` em tabelas didáticas para leitura em HTML.
