# Revisao de Literatura: Shift-Share Instruments (Bartik Instruments)

**Data**: 2026-02-16 (atualizado)
**Objetivo**: Textos pedagogicos, exposicao matematica, aplicacoes em CP e RI (por cientistas politicos em journals de politica), e implementacao em R

---

## 1. Visao geral do campo

Os instrumentos shift-share (tambem chamados de instrumentos Bartik) combinam duas fontes de variacao para gerar variacao exogena em uma variavel endogena: (i) *shares* -- parcelas pre-determinadas que capturam a exposicao diferencial de unidades locais a choques agregados, e (ii) *shifts* -- choques agregados (nacionais ou internacionais) que afetam diferentes setores/categorias de forma heterogenea. A interacao entre shares e shifts gera uma medida de exposicao local ao choque que, sob certas condicoes, pode ser usada como variavel instrumental.

O campo passou por uma revolucao metodologica entre 2018-2025, com tres contribuicoes centrais que clarificaram as condicoes de identificacao: Goldsmith-Pinkham, Sorkin & Swift (2020), que mostram que a identificacao vem da exogeneidade das *shares*; Borusyak, Hull & Jaravel (2022), que mostram um caminho alternativo onde a identificacao vem da exogeneidade dos *shifts*; e Adao, Kolesar & Morales (2019), que resolvem problemas de inferencia. Borusyak, Hull & Jaravel (2025) sintetizam tudo em um guia pratico publicado no *Journal of Economic Perspectives*.

---

## 2. Exposicao matematica detalhada

### 2.1 O instrumento shift-share: construcao

Considere $N$ regioes (indexadas por $i$) e $K$ setores/industrias (indexadas por $k$). Para cada regiao $i$, observamos a parcela (*share*) do setor $k$ no emprego (ou outra variavel) local em um periodo base $t_0$:

$$w_{ik} = \frac{L_{ik,t_0}}{L_{i,t_0}}$$

onde $L_{ik,t_0}$ e o emprego no setor $k$ na regiao $i$ no periodo base, e $L_{i,t_0}$ e o emprego total na regiao $i$. Note que $\sum_k w_{ik} = 1$.

Observamos tambem um choque (*shift*) agregado no nivel do setor $k$, denotado $g_k$. No exemplo classico de Bartik (1991), $g_k$ e a taxa de crescimento nacional do emprego no setor $k$:

$$g_k = \frac{\Delta L_k^{Nacional}}{L_{k,t_0}^{Nacional}}$$

O **instrumento shift-share** (ou instrumento Bartik) para a regiao $i$ e definido como:

$$B_i = \sum_{k=1}^{K} w_{ik} \cdot g_k$$

Este instrumento calcula a media ponderada dos choques setoriais, usando as shares pre-determinadas como pesos. A intuicao e que regioes com composicao setorial diferente serao diferencialmente afetadas pelos mesmos choques agregados.

### 2.2 O modelo estrutural e o IV

O objetivo tipico e estimar o efeito causal de uma variavel endogena $X_i$ sobre uma variavel de resultado $Y_i$. O modelo estrutural e:

$$Y_i = \alpha + \beta X_i + \varepsilon_i$$

onde $X_i$ e endogena (por exemplo, crescimento do emprego local, que e tanto causa quanto consequencia de outros fatores locais). OLS e inconsistente porque $\text{Cov}(X_i, \varepsilon_i) \neq 0$.

**Primeiro estagio (first stage)**: O instrumento $B_i$ preve $X_i$:

$$X_i = \gamma_0 + \gamma_1 B_i + \nu_i$$

A condicao de relevancia exige $\gamma_1 \neq 0$: as shares combinadas com os shifts devem efetivamente prever a variavel endogena.

**Segundo estagio (second stage)**: Substitui-se $X_i$ por $\hat{X}_i$:

$$Y_i = \alpha + \beta \hat{X}_i + \eta_i$$

A condicao de exclusao exige $\text{Cov}(B_i, \varepsilon_i) = 0$: o instrumento so afeta $Y_i$ atraves de $X_i$. E aqui que a recente literatura diverge sobre *de onde* vem essa exogeneidade.

### 2.3 Duas estrategias de identificacao

#### Estrategia 1: Exogeneidade das shares (GPSS 2020)

Goldsmith-Pinkham, Sorkin & Swift (2020) demonstram um resultado de equivalencia fundamental: o estimador IV usando $B_i = \sum_k w_{ik} g_k$ como instrumento e **numericamente identico** a um estimador GMM que usa cada share $w_{ik}$ (para $k = 1, ..., K$) como instrumento separado, com os shifts $g_k$ servindo apenas como pesos na combinacao otima:

$$\hat{\beta}^{IV} = \sum_k \hat{\alpha}_k \hat{\beta}_k$$

onde $\hat{\beta}_k$ e a estimativa just-identified usando apenas $w_{ik}$ como instrumento (para o setor $k$), e $\hat{\alpha}_k$ e o **peso de Rotemberg** que mede a influencia relativa do setor $k$ na estimativa global. Os pesos de Rotemberg satisfazem $\sum_k \hat{\alpha}_k = 1$ mas podem ser **negativos**, o que e um sinal de alerta.

**Hipotese de identificacao**: As shares iniciais $w_{ik}$ sao exogenas -- isto e, nao-correlacionadas com fatores nao-observados que afetam $Y_i$. Em termos formais:

$$\mathbb{E}[w_{ik} \cdot \varepsilon_i] = 0 \quad \forall k$$

Os shifts $g_k$ nao precisam ser exogenos; eles apenas afetam a ponderacao e a relevancia do instrumento.

**Diagnosticos**: Calcular os Rotemberg weights $\hat{\alpha}_k$ para identificar quais setores mais influenciam a estimativa. Testar balanceamento das shares mais influentes contra covariaveis pre-determinadas.

#### Estrategia 2: Exogeneidade dos shifts (BHJ 2022)

Borusyak, Hull & Jaravel (2022) oferecem um framework alternativo. O resultado de equivalencia deles e diferente: a regressao IV no nivel regional e numericamente equivalente a uma regressao IV no **nivel do choque** (setor), onde a variavel de resultado e o tratamento sao agregados ao nivel do choque usando as shares como pesos:

$$\bar{Y}_k = \sum_i s_{ik} Y_i, \quad \bar{X}_k = \sum_i s_{ik} X_i$$

onde $s_{ik} = w_{ik} / \sum_i w_{ik}$ normaliza os pesos. A regressao no nivel do choque:

$$\bar{Y}_k = \alpha + \beta \bar{X}_k + \bar{\varepsilon}_k$$

instrumentada por $g_k$, produz **exatamente a mesma estimativa** $\hat{\beta}$ que a regressao regional com $B_i$.

**Hipotese de identificacao**: Os choques $g_k$ sao quasi-aleatoriamente atribuidos -- isto e, nao-correlacionados com a media ponderada (por shares) dos fatores nao-observados regionais:

$$\mathbb{E}[g_k \cdot \bar{\varepsilon}_k] = 0 \quad \forall k$$

As shares $w_{ik}$ podem ser endogenas. A identificacao vem inteiramente da aleatoriedade dos choques, dado que ha muitos choques independentes ($K \to \infty$).

**Diagnosticos**: Testes de pre-tendencias e balanceamento no nivel do choque. Verificar se $g_k$ correlaciona-se com $\bar{\varepsilon}_k$ pre-tratamento.

### 2.4 O problema de inferencia (AKM 2019)

Adao, Kolesar & Morales (2019) identificam um problema critico: erros-padrao convencionais (incluindo clustered por regiao ou estado) sao **severamente subdimensionados** em regressoes shift-share. Em exercicios de placebo, testes ao nivel nominal de 5% rejeitam a hipotese nula verdadeira ate **55%** das vezes.

O problema surge porque residuos $\hat{\varepsilon}_i$ sao correlacionados entre regioes com composicao setorial similar, independentemente de proximidade geografica. Duas regioes $i$ e $j$ com shares parecidos ($w_{ik} \approx w_{jk}$) compartilham exposicao aos mesmos choques setoriais, induzindo correlacao em $\hat{\varepsilon}_i$ e $\hat{\varepsilon}_j$.

A variancia correta do estimador IV deve levar em conta essa estrutura. AKM propoem:

$$\hat{V}_{AKM} = \left(\frac{\partial \hat{\beta}}{\partial g}\right)' \hat{\Sigma}_g \left(\frac{\partial \hat{\beta}}{\partial g}\right)$$

onde $\hat{\Sigma}_g$ captura a correlacao entre os residuos no nivel dos choques. Na pratica, os intervalos de confianca corrigidos por AKM sao substancialmente mais largos.

**Alternativa BHJ**: Em vez de corrigir os erros-padrao no nivel regional, rodar a regressao diretamente no nivel do choque (via `ssaggregate`) e usar erros-padrao convencionais (heteroscedasticidade-robustos) no nivel do choque. Se os choques sao independentes, a inferencia convencional e valida.

### 2.5 Resumo das condicoes

| Condicao | Exogenous Shares (GPSS) | Exogenous Shifts (BHJ) |
|----------|------------------------|----------------------|
| **Identificacao** | $\mathbb{E}[w_{ik} \cdot \varepsilon_i] = 0$ | $\mathbb{E}[g_k \cdot \bar{\varepsilon}_k] = 0$ |
| **Relevancia** | $\gamma_1 \neq 0$ (first stage forte) | $\gamma_1 \neq 0$ |
| **Inferencia** | Rotemberg weights + clustering usual | AKM SEs ou regressao no nivel do choque |
| **N. de choques** | Pode funcionar com poucos | Precisa de muitos ($K$ grande) |
| **Shares endogenos?** | Nao | Sim |
| **Shifts endogenos?** | Sim (so afetam pesos) | Nao |

---

## 3. Textos pedagogicos e metodologicos fundamentais

### 3.1 Trabalhos seminais

| Autor(es) | Ano | Titulo | Journal/Editora | Contribuicao |
|-----------|-----|--------|-----------------|--------------|
| Bartik | 1991 | *Who Benefits from State and Local Economic Development Policies?* | W.E. Upjohn Institute | Trabalho original. Shares = emprego industrial local; Shifts = crescimento nacional por industria. |
| Card | 2001 | "Immigrant Inflows, Native Outflows..." | *J. Labor Economics*, 19(1): 22-64 | Instrumento de "enclave": shares = distribuicao geografica historica de imigrantes; shifts = fluxos nacionais por pais de origem. |
| Autor, Dorn & Hanson | 2013 | "The China Syndrome" | *AER*, 103(6): 2121-2168 | "China Shock". Shares = emprego industrial local; Shifts = importacoes chinesas por industria, instrumentadas por exportacoes chinesas para outros paises. |
| Goldsmith-Pinkham, Sorkin & Swift | 2020 | "Bartik Instruments: What, When, Why, and How" | *AER*, 110(8): 2586-2624 | Equivalencia com GMM de shares individuais. Identificacao via **exogeneidade das shares**. Rotemberg weights. |
| Adao, Kolesar & Morales | 2019 | "Shift-Share Designs: Theory and Inference" | *QJE*, 134(4): 1949-2010 | SEs convencionais rejeitam em 55% quando nominal e 5%. Propoe variancia corrigida (AKM). |
| Borusyak, Hull & Jaravel | 2022 | "Quasi-Experimental Shift-Share Research Designs" | *Rev. Econ. Studies*, 89(1): 181-213 | Identificacao via **exogeneidade dos shifts**. Equivalencia com regressao no nivel do choque. |
| Borusyak & Hull | 2023 | "Nonrandom Exposure to Exogenous Shocks" | *Econometrica*, 91(6): 2155-2185 | Generaliza para "formula instruments". Instrumento recentrado. |
| Jaeger, Ruist & Stuhler | 2018 | "Shift-Share Instruments and the Impact of Immigration" | NBER WP 24285 | Vies dinamico quando distribuicao de imigrantes e estavel no tempo. |

### 3.2 Guias praticos e pedagogicos

| Autor(es) | Ano | Titulo | Fonte | Descricao |
|-----------|-----|--------|-------|-----------|
| **Borusyak, Hull & Jaravel** | **2025** | **"A Practical Guide to Shift-Share Instruments"** | ***J. Econ. Perspectives*, 39(1): 181-204** | **Melhor recurso pedagogico.** Checklists praticos para as duas estrategias. |
| Borusyak, Hull & Jaravel | 2025 | "Design-Based Identification with Formula Instruments" | *Econometrics Journal*, 28(1): 83-108 | Revisao tecnica: shift-share dentro da classe de "formula instruments". |
| Breuer | 2022 | "Bartik Instruments: An Applied Introduction" | *J. Financial Reporting*, 7(1): 49-67 | Introducao acessivel com simulacoes. [GitHub](https://github.com/mb4468/bartik_intro) |
| Hull (Mixtape Sessions) | -- | Shift-Share IV Workshop | [GitHub](https://github.com/Mixtape-Sessions/Shift-Share) | Workshop com slides, exercicios em Stata/R e videos. Gratuito. |
| Tilburg Science Hub | -- | "The Shift-Share Instrumental Variable" | [Tutorial](https://tilburgsciencehub.com/topics/analyze/causal-inference/instrumental-variables/shiftshare/) | Tutorial conciso (~6 min). |
| World Bank (Dev Impact) | -- | "Rethinking Identification Under the Bartik Shift-Share Instrument" | [Blog](https://blogs.worldbank.org/en/impactevaluations/rethinking-identification-under-bartik-shift-share-instrument) | Post acessivel sobre implicacoes de GPSS e Jaeger et al. |
| UC Berkeley D-Lab | -- | "A Practical Guide to Shift-Share Instruments" | [Blog](https://dlab.berkeley.edu/news/practical-guide-shift-share-instruments-and-what-i-learned-replicating-china-shock) | Walkthrough replicando o China Shock com framework BHJ. |

---

## 4. Implementacao em R

### 4.1 Construcao do instrumento

O instrumento $B_i = \sum_k w_{ik} g_k$ e uma simples multiplicacao matricial:

```r
# shares: matriz N x K de exposure shares (cada linha soma 1)
# shocks: vetor K x 1 de choques setoriais
bartik_instrument <- shares %*% shocks
```

### 4.2 Estimacao IV

Diversos pacotes R permitem estimar 2SLS:

| Pacote | Funcao | Vantagens |
|--------|--------|-----------|
| `fixest` | `feols(y ~ controls \| FEs \| endogenous ~ instrument)` | Rapido, suporta FE de alta dimensao, SEs clusterizados |
| `ivreg` | `ivreg(y ~ endogenous + controls \| instrument + controls)` | Classico, sintaxe familiar |
| `estimatr` | `iv_robust()` | SEs robustos/clusterizados built-in |

Exemplo com `fixest`:

```r
library(fixest)

# Construir instrumento
dados$bartik <- as.numeric(shares %*% shocks)

# IV com efeitos fixos de estado
modelo <- feols(y ~ controles | estado | x_endogeno ~ bartik, data = dados)
summary(modelo)
```

### 4.3 Pacotes especializados para shift-share

#### `ShiftShareSE` (Kolesar) -- **disponivel no CRAN**

Implementa os erros-padrao AKM (Adao, Kolesar & Morales, 2019):

```r
install.packages("ShiftShareSE")
library(ShiftShareSE)

# ivreg_ss: IV com SEs corrigidos para shift-share
resultado <- ivreg_ss(y ~ x_endogeno | bartik,
                       X = shares_matrix,    # matriz N x K de shares
                       data = dados)
summary(resultado)
```

Este pacote e essencial porque, como demonstrado por AKM, erros-padrao convencionais (mesmo clusterizados por geografia) sao severamente subdimensionados em designs shift-share.

#### `ssaggregate` (Kyle Butts) -- **GitHub apenas**

Implementa a agregacao ao nivel do choque de Borusyak, Hull & Jaravel (2022):

```r
# install: remotes::install_github("kylebutts/ssaggregate")
library(ssaggregate)

# Agregar dados regionais para o nivel do choque
dados_choque <- ssaggregate(
  data = dados_long,     # formato longo: regiao x setor
  shares = "share",
  n = "regiao",
  s = "setor",
  t = "ano",
  y = c("y", "x_endogeno")
)

# Agora rodar IV convencional no nivel do choque
feols(y ~ 1 | x_endogeno ~ g_k, data = dados_choque)
```

### 4.4 Comparacao R vs. Stata

| Funcionalidade | Stata | R | Status R |
|---------------|-------|---|----------|
| Agregacao BHJ | `ssaggregate` (Borusyak et al.) | `ssaggregate` (Kyle Butts) | GitHub apenas, nao no CRAN |
| SEs AKM | `ShiftShareSE` (ado file) | `ShiftShareSE` (Kolesar) | **No CRAN** |
| Rotemberg weights | `bartik_weight` (GPSS) | **Nenhum pacote** | Gap principal |
| IV basico | `ivregress`, `ivreg2` | `fixest`, `ivreg`, `estimatr` | Excelente suporte |

### 4.5 A lacuna dos Rotemberg weights no R

O diagnostico de Rotemberg weights -- central na abordagem GPSS (2020) -- **nao possui um pacote R dedicado**. Em Stata, o comando `bartik_weight` de Goldsmith-Pinkham, Sorkin & Swift calcula automaticamente esses pesos. Em R, o pesquisador precisa:

1. Adaptar o codigo de replicacao disponivel no GitHub dos autores ([github.com/paulgp/bartik-weight](https://github.com/paulgp/bartik-weight)), que inclui scripts R nao-empacotados;
2. Ou implementar manualmente o calculo, que envolve rodar $K$ regressoes just-identified (uma por setor) e computar os pesos.

Essa lacuna tem consequencias praticas: pesquisadores que usam R e adotam a estrategia de exogenous shares (GPSS) nao tem uma ferramenta plug-and-play para o diagnostico mais importante dessa abordagem. Isso pode:
- **Desencorajar** o uso de diagnosticos adequados, levando a resultados menos robustos;
- **Enviesar a escolha metodologica** em favor da abordagem BHJ (exogenous shifts), que tem melhor suporte em R via `ssaggregate`, independentemente de qual estrategia e mais adequada ao problema;
- **Criar uma barreira de entrada** para pesquisadores de CP e RI que tipicamente usam R (nao Stata), dificultando a adocao da metodologia.

Para o futuro, seria util que um pacote como `bartikweights` fosse desenvolvido para R e publicado no CRAN, consolidando a implementacao dos Rotemberg weights com interface padronizada.

---

## 5. Aplicacoes em Ciencia Politica

*Nota: Esta secao inclui apenas papers publicados em journals de ciencia politica (APSR, AJPS, JOP, CPS, Electoral Studies, BJPS, IO, ISQ) por cientistas politicos ou equipes com cientistas politicos como autores principais.*

### 5.1 Comercio e voto nos EUA

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Baccini & Weymouth | 2021 | "Gone For Good: Deindustrialization, White Voter Backlash, and US Presidential Voting" | *APSR*, 115(2): 550-567 | Share de trabalhadores por industria manufatureira por condado | Importacoes chinesas por industria (2000-2015) | Brancos votam mais em republicanos em condados com mais demissoes na manufatura; negros votam mais em democratas. Desindustrializacao ameaca status do grupo dominante. |
| Feigenbaum & Hall | 2015 | "How Legislators Respond to Localized Economic Shocks" | *J. Politics*, 77(4) | Shares de emprego industrial por distrito congressual | Crescimento de importacoes chinesas por industria | Choques de importacao fazem legisladores votar de forma mais protecionista em projetos de comercio, sem efeito em outras votacoes nem em reeleicao. |
| Margalit | 2011 | "Costly Jobs: Trade-related Layoffs, Government Compensation, and Voting in U.S. Elections" | *APSR*, 105(1): 166-188 | Exposicao local a demissoes relacionadas a comercio (via Trade Adjustment Assistance) | Choques de comercio que geram demissoes | Eleitores sao mais sensiveis a perda de empregos por competicao estrangeira do que por outros fatores. Compensacao governamental (TAA) mitiga o backlash politico. |

### 5.2 Comercio e nacionalismo/populismo na Europa

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Colantone & Stanig | 2018 | "Global Competition and Brexit" | *APSR*, 112(2): 201-218 | Shares de emprego regional (NUTS-2/3) por industria no UK | Importacoes chinesas por industria (instrumentadas por importacoes chinesas nos EUA) | Voto Leave sistematicamente maior em regioes mais expostas a competicao chinesa. Um desvio-padrao no choque = ~2pp de Leave. |
| Colantone & Stanig | 2018 | "The Trade Origins of Economic Nationalism" | *AJPS*, 62(4): 936-953 | Shares de emprego regional (NUTS-2) por industria em 15 paises europeus | Importacoes chinesas por industria (instrumentadas por importacoes nos EUA) | Maior exposicao ao choque chines aumenta apoio a partidos nacionalistas e de extrema-direita. Shift geral para a direita. 1988-2007. |
| Hays, Lim & Spoon | 2019 | "The Path from Trade to Right-Wing Populism in Europe" | *Electoral Studies* | Shares de emprego regional por industria manufatureira | Crescimento de importacoes de paises de baixo salario por industria | Regioes europeias mais expostas a competicao comercial tiveram aumento no voto em partidos populistas de direita, via inseguranca economica. |
| Rommel & Walter | 2022 | "The Electoral Consequences of Offshoring" | *Comparative Political Studies* | Shares de emprego regional em setores suscetiveis a offshoring | Mudancas nacionais/globais na intensidade de offshoring por industria | Regioes mais expostas a offshoring tiveram aumento no apoio a partidos de direita radical e reducao no apoio a partidos mainstream. |

### 5.3 Valores autoritarios e choques comerciais

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Ballard-Rosa, Malik, Rickard & Scheve | 2021 | "The Economic Origins of Authoritarian Values: Evidence from Local Trade Shocks in the United Kingdom" | *Comparative Political Studies*, 54(13): 2321-2353 | Shares de emprego industrial local no UK | Crescimento de importacoes chinesas por industria | Regioes britanicas mais afetadas por importacoes chinesas tem valores mais autoritarios. Mecanismo de frustracao-agressao. Venceu o David A. Lake Award no IPES 2017. |
| Ballard-Rosa, Jensen & Scheve | 2022 | "Economic Decline, Social Identity, and Authoritarian Values in the United States" | *Intl Studies Quarterly*, 66(1) | Shares de emprego industrial em commuting zones | Crescimento de importacoes chinesas por industria | Individuos em regioes diversas com competicao chinesa intensa tem valores mais autoritarios. Declinio economico afeta identidade social de grupos dominantes. |

### 5.4 Comercio e welfare state historico

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Scheve & Serlin | 2023 | "The German Trade Shock and the Rise of the Neo-Welfare State in Early Twentieth-Century Britain" | *APSR*, 117(2): 557-574 | Shares de emprego industrial por constituency britanica (c. 1880) | Crescimento de importacoes alemas por industria (1880-1910) | Importacoes alemas pioraram mercado de trabalho e mudaram crencas sobre deservingness dos pobres. Contribuiu para o surgimento do welfare state. |

### 5.5 Automacao e preferencias redistributivas

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Thewissen & Rueda | 2019 | "Automation and the Welfare State: Technological Change as a Determinant of Redistribution Preferences" | *Comparative Political Studies* | Shares de emprego em ocupacoes suscetiveis a automacao | Mudancas nacionais na intensidade de automacao por industria | Trabalhadores mais expostos a automacao expressam preferencias mais fortes por redistribuicao. Efeito mediado por contexto institucional. |

### 5.6 Comercio e politica no Brasil

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Campello & Urdinez | 2021 | "Voter and Legislator Responses to Localized Trade Shocks from China in Brazil" | *Comparative Political Studies*, 54(7): 1131-1162 | Shares de emprego/producao municipal por industria | Crescimento de importacoes/exportacoes chinesas por setor | Respostas assimetricas: regioes afetadas por importacoes veem China como risco, mas regioes beneficiadas por exportacoes *nao* veem como oportunidade. |

---

## 6. Aplicacoes em Relacoes Internacionais

*Nota: Inclui papers em journals de RI (IO, ISQ, J. Contemporary China) e CPS quando o tema e explicitamente de RI, por autores com formacao em ciencia politica/RI. Papers de economistas em journals de economia sobre temas de RI sao listados separadamente como "Referencias complementares de economia".*

### 6.1 Comercio e politica externa/voto presidencial

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Jensen, Quinn & Weymouth | 2017 | "Winners and Losers in International Trade: The Effects on U.S. Presidential Voting" | *International Organization*, 71(3): 423-457 | Shares de emprego por industria tradable (bens e servicos) | Mudancas nacionais em penetracao de importacoes e crescimento de exportacoes | Importacoes reduzem voto no incumbente; exportacoes aumentam. Swing states com manufatura de baixa qualificacao foram particularmente vulneraveis. |
| Broz, Frieden & Weymouth | 2021 | "Populism in Place: The Economic Geography of the Globalization Backlash" | *International Organization* | Shares de emprego industrial local | Mudancas em competicao de importacao e exposicao a automacao (robos) por industria | Regioes simultaneamente expostas a importacoes e automacao tiveram o maior aumento no voto populista (Trump 2016). Interacao dos dois choques mais poderosa que cada um isoladamente. |
| Kuk, Seligsohn & Zhang | 2018 | "From Tiananmen to Outsourcing: The Effect of Rising Import Competition on Congressional Voting Towards China" | *J. Contemporary China*, 27(109) | Shares de emprego manufatureiro por distrito | Crescimento de importacoes chinesas por industria | Apos 2003 (pos-OMC), legisladores de distritos afetados tornaram-se sistematicamente mais hostis a China em votos de politica externa. |
| Owen | 2017 | "Exposure to Offshoring and the Politics of Trade Liberalization: Debate and Votes on DR-CAFTA" | *International Studies Quarterly* | Shares de emprego em industrias com alto potencial de offshoring por distrito | Tendencias nacionais de offshoring por industria | Congressistas de distritos mais expostos a offshoring votaram mais contra o DR-CAFTA, mostrando como vulnerabilidade economica local molda comportamento legislativo em politica comercial. |

### 6.2 Choques comerciais e atitudes/valores

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Ballard-Rosa, Jensen & Scheve | 2022 | "Economic Decline, Social Identity, and Authoritarian Values in the United States" | *Intl Studies Quarterly*, 66(1) | Shares de emprego industrial em commuting zones | Crescimento de importacoes chinesas | Competicao chinesa em regioes diversas gera valores mais autoritarios. Declinio economico ameaca identidade social de grupos dominantes. |
| Campello & Urdinez | 2021 | "Voter and Legislator Responses to Localized Trade Shocks from China in Brazil" | *Comparative Political Studies*, 54(7) | Shares municipais por industria | Importacoes/exportacoes chinesas por setor | Perdas de comercio moldam percepcoes de politica externa mais que ganhos. Aplicacao rara fora dos EUA/Europa. |

### 6.3 Referencias complementares de economia sobre temas de RI

Estes papers sao de economistas em journals de economia, mas abordam temas centrais de RI e sao amplamente citados na disciplina:

| Autor(es) | Ano | Titulo | Journal | Tema RI |
|-----------|-----|--------|---------|---------|
| Dube & Vargas | 2013 | "Commodity Price Shocks and Civil Conflict: Evidence from Colombia" | *Rev. Econ. Studies* | Precos de commodities e conflito civil. Dois mecanismos opostos: custo de oportunidade (cafe) vs. rapacidade (petroleo). |
| Gallea | 2023 | "Weapons and War: The Effect of Arms Transfers on Internal Conflict" | *J. Dev. Economics* | Transferencias de armas e conflito na Africa. Shares = participacao historica de cada fornecedor; Shifts = guerras dos fornecedores fora da Africa. |
| Fitchett & Wesselbaum | 2022 | "Does Aid Drive Migration?" | *Intl Migration Review* | Ajuda internacional e migracao. Aumento de 10% na ajuda aumenta migracao em ~2%. |

---

## 7. Padroes recorrentes na construcao do instrumento

### 7.1 Choques comerciais (o mais comum em CP/RI)
- **Shares** = participacao industrial local no emprego
- **Shifts** = crescimento de importacoes (geralmente da China) por industria, instrumentado por importacoes em terceiros paises
- **Uso**: Polarizacao, nacionalismo, populismo, Brexit, voto extrema-direita, politica externa

### 7.2 Offshoring e automacao (crescente)
- **Shares** = emprego local em ocupacoes/industrias suscetiveis a offshoring ou automacao
- **Shifts** = tendencias nacionais de offshoring ou adocao de robos por industria
- **Uso**: Preferencias por redistribuicao, voto populista, comportamento legislativo

### 7.3 Imigracao (instrumento de enclave tipo Card)
- **Shares** = padroes historicos de assentamento de imigrantes por pais de origem
- **Shifts** = fluxos nacionais de imigracao por pais, impulsionados por eventos exogenos
- **Uso**: Voto anti-imigracao, welfare state, backlash politico

---

## 8. Gaps identificados

### Gaps na Ciencia Politica
- **Pouca aplicacao no Brasil e America Latina**: Alem de Campello & Urdinez (2021), ha pouquissimas aplicacoes em CP na regiao
- **Politica subnacional**: Maioria dos estudos foca em eleicoes nacionais; aplicacoes para eleicoes municipais ou estaduais sao escassas
- **Decomposicao de voto por partido**: O capitulo 15 ja aponta nessa direcao -- contribuicao pedagogica original
- **Politica fiscal subnacional**: Choques de receita municipal (royalties, transferencias federais) como instrumento shift-share

### Gaps em Relacoes Internacionais
- **Sancoes internacionais**: Nenhum estudo usando shift-share para medir exposicao diferencial a sancoes
- **Cooperacao internacional**: Ausencia de aplicacoes em estudos sobre cooperacao, OIs ou regimes
- **China Shock fora dos EUA/Europa**: Poucos estudos (apenas Campello & Urdinez para o Brasil)
- **Conflito**: Literatura dominada por economistas; poucos papers em journals de RI

### Gaps metodologicos/software
- **Rotemberg weights em R**: Principal lacuna do ecossistema R (ver secao 4.5)
- **`ssaggregate` nao esta no CRAN**: Barreira para adocao
- **Inferencia com poucos choques**: Muitos contextos de CP/RI tem poucos choques relevantes
- **Dados em painel + shift-share**: Pouca orientacao pedagogica sobre a combinacao

---

## 9. Sugestoes para o capitulo 15

1. **Referencia pedagogica central**: Borusyak, Hull & Jaravel (2025, JEP)
2. **Exemplo motivador de CP**: Colantone & Stanig (2018, APSR) sobre Brexit
3. **Exemplo motivador de RI**: Jensen, Quinn & Weymouth (2017, IO)
4. **Incluir derivacao matematica** da secao 2 deste documento, adaptada ao exemplo de voto por partido do capitulo
5. **Implementacao em R** com `fixest` para IV e `ShiftShareSE` para SEs corretos
6. **Discutir as duas estrategias de identificacao** com exemplos concretos

---

## 10. Referencias-chave (formato APSA)

Adao, Rodrigo, Michal Kolesar, and Eduardo Morales. 2019. "Shift-Share Designs: Theory and Inference." *Quarterly Journal of Economics* 134(4): 1949-2010.

Baccini, Leonardo, and Stephen Weymouth. 2021. "Gone For Good: Deindustrialization, White Voter Backlash, and US Presidential Voting." *American Political Science Review* 115(2): 550-567.

Ballard-Rosa, Cameron, Amalie Jensen, and Kenneth Scheve. 2022. "Economic Decline, Social Identity, and Authoritarian Values in the United States." *International Studies Quarterly* 66(1): sqab027.

Ballard-Rosa, Cameron, Mashail Malik, Stephanie Rickard, and Kenneth Scheve. 2021. "The Economic Origins of Authoritarian Values: Evidence From Local Trade Shocks in the United Kingdom." *Comparative Political Studies* 54(13): 2321-2353.

Bartik, Timothy J. 1991. *Who Benefits from State and Local Economic Development Policies?* Kalamazoo, MI: W.E. Upjohn Institute for Employment Research.

Borusyak, Kirill, and Peter Hull. 2023. "Nonrandom Exposure to Exogenous Shocks." *Econometrica* 91(6): 2155-2185.

Borusyak, Kirill, Peter Hull, and Xavier Jaravel. 2022. "Quasi-Experimental Shift-Share Research Designs." *Review of Economic Studies* 89(1): 181-213.

Borusyak, Kirill, Peter Hull, and Xavier Jaravel. 2025. "A Practical Guide to Shift-Share Instruments." *Journal of Economic Perspectives* 39(1): 181-204.

Breuer, Matthias. 2022. "Bartik Instruments: An Applied Introduction." *Journal of Financial Reporting* 7(1): 49-67.

Broz, J. Lawrence, Jeffry Frieden, and Stephen Weymouth. 2021. "Populism in Place: The Economic Geography of the Globalization Backlash." *International Organization*.

Campello, Daniela, and Francisco Urdinez. 2021. "Voter and Legislator Responses to Localized Trade Shocks from China in Brazil." *Comparative Political Studies* 54(7): 1131-1162.

Card, David. 2001. "Immigrant Inflows, Native Outflows, and the Local Labor Market Impacts of Higher Immigration." *Journal of Labor Economics* 19(1): 22-64.

Colantone, Italo, and Piero Stanig. 2018a. "Global Competition and Brexit." *American Political Science Review* 112(2): 201-218.

Colantone, Italo, and Piero Stanig. 2018b. "The Trade Origins of Economic Nationalism." *American Journal of Political Science* 62(4): 936-953.

Dube, Oeindrila, and Juan Vargas. 2013. "Commodity Price Shocks and Civil Conflict: Evidence from Colombia." *Review of Economic Studies* 80(4): 1384-1421.

Feigenbaum, James, and Andrew Hall. 2015. "How Legislators Respond to Localized Economic Shocks." *Journal of Politics* 77(4).

Fitchett, Hamish, and Dennis Wesselbaum. 2022. "Does Aid Drive Migration?" *International Migration Review* 56(3).

Gallea, Quentin. 2023. "Weapons and War: The Effect of Arms Transfers on Internal Conflict." *Journal of Development Economics* 160.

Goldsmith-Pinkham, Paul, Isaac Sorkin, and Henry Swift. 2020. "Bartik Instruments: What, When, Why, and How." *American Economic Review* 110(8): 2586-2624.

Hays, Jude, Junghyun Lim, and Jae-Jae Spoon. 2019. "The Path from Trade to Right-Wing Populism in Europe." *Electoral Studies*.

Jaeger, David A., Joakim Ruist, and Jan Stuhler. 2018. "Shift-Share Instruments and the Impact of Immigration." NBER Working Paper 24285.

Jensen, J. Bradford, Dennis Quinn, and Stephen Weymouth. 2017. "Winners and Losers in International Trade." *International Organization* 71(3): 423-457.

Kuk, John Seungmin, Deborah Seligsohn, and Jiakun Jack Zhang. 2018. "From Tiananmen to Outsourcing." *Journal of Contemporary China* 27(109).

Margalit, Yotam. 2011. "Costly Jobs: Trade-related Layoffs, Government Compensation, and Voting in U.S. Elections." *American Political Science Review* 105(1): 166-188.

Owen, Erica. 2017. "Exposure to Offshoring and the Politics of Trade Liberalization." *International Studies Quarterly*.

Rommel, Tobias, and Stefanie Walter. 2022. "The Electoral Consequences of Offshoring." *Comparative Political Studies*.

Scheve, Kenneth, and Theo Serlin. 2023. "The German Trade Shock and the Rise of the Neo-Welfare State in Early Twentieth-Century Britain." *American Political Science Review* 117(2): 557-574.

Thewissen, Stefan, and David Rueda. 2019. "Automation and the Welfare State." *Comparative Political Studies*.
