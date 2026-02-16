# Revisao de Literatura: Shift-Share Instruments (Bartik Instruments)

**Data**: 2026-02-16
**Objetivo**: Mapear textos pedagogicos/explicativos e aplicacoes em Ciencia Politica e Relacoes Internacionais

---

## 1. Visao geral do campo

Os instrumentos shift-share (tambem chamados de instrumentos Bartik) combinam duas fontes de variacao para gerar variacao exogena em uma variavel endogena: (i) *shares* -- parcelas pre-determinadas que capturam a exposicao diferencial de unidades locais a choques agregados, e (ii) *shifts* -- choques agregados (nacionais ou internacionais) que afetam diferentes setores/categorias de forma heterogenea. A interacao entre shares e shifts gera uma medida de exposicao local ao choque que, sob certas condicoes, pode ser usada como variavel instrumental.

O campo passou por uma revolucao metodologica entre 2018-2025, com tres contribuicoes centrais que clarificaram as condicoes de identificacao: Goldsmith-Pinkham, Sorkin & Swift (2020), que mostram que a identificacao vem da exogeneidade das *shares*; Borusyak, Hull & Jaravel (2022), que mostram um caminho alternativo onde a identificacao vem da exogeneidade dos *shifts*; e Adao, Kolesar & Morales (2019), que resolvem problemas de inferencia. Borusyak, Hull & Jaravel (2025) sintetizam tudo em um guia pratico publicado no *Journal of Economic Perspectives*.

---

## 2. Textos pedagogicos e metodologicos fundamentais

### 2.1 Trabalhos seminais

| Autor(es) | Ano | Titulo | Journal/Editora | Contribuicao |
|-----------|-----|--------|-----------------|--------------|
| Bartik | 1991 | *Who Benefits from State and Local Economic Development Policies?* | W.E. Upjohn Institute | Trabalho original que introduziu o instrumento. Combina shares de emprego industrial local com shifts de crescimento nacional por industria para prever crescimento local de emprego. |
| Card | 2001 | "Immigrant Inflows, Native Outflows, and the Local Labor Market Impacts of Higher Immigration" | *J. Labor Economics*, 19(1): 22-64 | Aplicacao influente na economia da imigracao. Cria o instrumento de "enclave" combinando distribuicao geografica historica de imigrantes por pais de origem (shares) com fluxos nacionais de imigracao (shifts). |
| Autor, Dorn & Hanson | 2013 | "The China Syndrome: Local Labor Market Effects of Import Competition in the United States" | *AER*, 103(6): 2121-2168 | A aplicacao moderna mais famosa ("China Shock"). Shares = participacao industrial local no emprego; Shifts = crescimento de importacoes chinesas por industria, instrumentado por exportacoes chinesas para outros paises. |
| Goldsmith-Pinkham, Sorkin & Swift | 2020 | "Bartik Instruments: What, When, Why, and How" | *AER*, 110(8): 2586-2624 | Demonstra que o IV shift-share e numericamente equivalente a um estimador GMM ponderado que usa as shares individuais como instrumentos separados. A identificacao vem da **exogeneidade das shares**. Introduz os **Rotemberg weights** como ferramenta diagnostica. |
| Adao, Kolesar & Morales | 2019 | "Shift-Share Designs: Theory and Inference" | *QJE*, 134(4): 1949-2010 | Mostra que erros-padrao convencionais sao severamente subdimensionados em regressoes shift-share (rejeicao de 55% quando o nominal e 5%). Desenvolve estimadores de variancia robustos a correlacao entre regioes com composicao setorial similar. |
| Borusyak, Hull & Jaravel | 2022 | "Quasi-Experimental Shift-Share Research Designs" | *Rev. Econ. Studies*, 89(1): 181-213 | Framework alternativo onde a identificacao vem da **exogeneidade dos shifts**. Resultado de equivalencia: a regressao IV no nivel regional pode ser obtida identicamente de uma regressao IV no nivel do choque. Disponibiliza o comando `ssaggregate` para Stata e R. |
| Borusyak & Hull | 2023 | "Nonrandom Exposure to Exogenous Shocks" | *Econometrica*, 91(6): 2155-2185 | Generaliza o framework para "formula instruments" -- instrumentos que combinam choques exogenos com variaveis pre-determinadas (potencialmente endogenas) via forma funcional conhecida. Introduz o **instrumento recentrado**. |
| Jaeger, Ruist & Stuhler | 2018 | "Shift-Share Instruments and the Impact of Immigration" | NBER WP 24285 | Identifica vies dinamico: quando a distribuicao geografica de imigrantes e estavel no tempo, o instrumento confunde respostas de curto e longo prazo. Propoe procedimento de "instrumentacao multipla". |

### 2.2 Guias praticos e pedagogicos

| Autor(es) | Ano | Titulo | Fonte | Descricao |
|-----------|-----|--------|-------|-----------|
| **Borusyak, Hull & Jaravel** | **2025** | **"A Practical Guide to Shift-Share Instruments"** | ***J. Econ. Perspectives*, 39(1): 181-204** | **O melhor recurso pedagogico disponivel.** Sintetiza as duas estrategias de identificacao (exogenous shares vs. exogenous shifts) com checklists praticos. Cobre construcao, verificacao de hipoteses, inferencia e shares incompletos. |
| Borusyak, Hull & Jaravel | 2025 | "Design-Based Identification with Formula Instruments: A Review" | *The Econometrics Journal*, 28(1): 83-108 | Revisao mais tecnica que situa shift-share dentro da classe mais ampla de "formula instruments". |
| Breuer | 2022 | "Bartik Instruments: An Applied Introduction" | *J. Financial Reporting*, 7(1): 49-67 | Introducao acessivel com codigo de simulacao. Explica que instrumentos Bartik isolam variacao do tratamento decorrente do impacto diferencial de choques comuns em unidades com exposicoes pre-determinadas distintas. [GitHub com codigo](https://github.com/mb4468/bartik_intro) |
| Ferri | 2022 | "Novel Shift-Share Instruments and Their Applications" | Boston College WP 1053 | Argumenta que a essencia do shift-share e decompor uma variavel endogena em identidade contabil, preservar o componente mais exogeno e neutralizar o mais endogeno. |
| Hull (Mixtape Sessions) | -- | Shift-Share IV Workshop | [GitHub](https://github.com/Mixtape-Sessions/Shift-Share) | Workshop de meio dia com slides, exercicios em Stata/R e videos. Materiais abertos e gratuitos. |
| Schiff | 2021-2023 | Lecture Notes: Graduate Urban Economics | [PDF](https://www.nathanschiff.com/webdocs/grad_urban/grad_urban_2021/bartik_goldsmith_etal_HANDOUT_2021.pdf) | Notas de aula estruturadas sobre GPSS (2020), incluindo exercicios e codigo Stata. |
| Tilburg Science Hub | -- | "The Shift-Share Instrumental Variable" | [Tutorial online](https://tilburgsciencehub.com/topics/analyze/causal-inference/instrumental-variables/shiftshare/) | Tutorial conciso (~6 min) cobrindo conceito, construcao e condicoes de validade. |
| World Bank (Dev Impact Blog) | -- | "Rethinking Identification Under the Bartik Shift-Share Instrument" | [Blog](https://blogs.worldbank.org/en/impactevaluations/rethinking-identification-under-bartik-shift-share-instrument) | Post acessivel discutindo implicacoes de GPSS (2020) e Jaeger et al. (2018) para economistas do desenvolvimento. |
| UC Berkeley D-Lab | -- | "A Practical Guide to Shift-Share Instruments" | [Blog](https://dlab.berkeley.edu/news/practical-guide-shift-share-instruments-and-what-i-learned-replicating-china-shock) | Walkthrough pratico replicando o China Shock com o framework BHJ. Cobre shares incompletos e erros-padrao corretos. |

### 2.3 As duas estrategias de identificacao

| | **Exogenous Shares** (GPSS 2020) | **Exogenous Shocks** (BHJ 2022) |
|---|---|---|
| **Hipotese-chave** | Shares iniciais sao as-good-as-randomly assigned | Shifts (choques) sao as-good-as-randomly assigned |
| **Shares podem ser endogenos?** | Nao | Sim |
| **Shifts podem ser endogenos?** | Sim (so afetam pesos) | Nao |
| **N. de choques necessario** | Pode funcionar com poucos | Precisa de muitos (assintotica de amostra grande) |
| **Inferencia** | Rotemberg weights + clustering usual | Erros-padrao AKM (2019) ou regressao no nivel do choque |
| **Diagnosticos** | Rotemberg weights, testes de balanceamento nas shares-chave | Pre-tendencias em agregados no nivel do choque |
| **Aplicacoes tipicas** | Elasticidades de oferta de trabalho, especializacao historica exogena | Imigracao (enclave), comercio (China Shock) |

### 2.4 Software

| Pacote | Plataforma | Descricao |
|--------|-----------|-----------|
| `ssaggregate` | [Stata](https://ideas.repec.org/c/boc/bocode/s458526.html), [R](https://github.com/kylebutts/ssaggregate) | Implementa resultado de equivalencia BHJ (2022). Constroi agregados no nivel do choque. |
| `ShiftShareSE` | [Stata](https://github.com/zhangxiang0822/ShiftShareSEStata) | Implementa estimador de variancia de Adao, Kolesar & Morales (2019). |

---

## 3. Aplicacoes em Ciencia Politica

### 3.1 Comercio e polarizacao politica nos EUA

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Autor, Dorn, Hanson & Majlesi | 2020 | "Importing Political Polarization? The Electoral Consequences of Rising Trade Exposure" | *AER*, 110(10): 3139-3183 | Shares de emprego industrial local em commuting zones | Crescimento de importacoes chinesas por industria (instrumentado por importacoes chinesas em 8 paises de alta renda) | Condados majoritariamente brancos expostos ao comercio migraram para conservadores GOP; condados majoritariamente minoritarios migraram para democratas liberais. Efeito se estende a eleicao presidencial de 2016. |
| Baccini & Weymouth | 2021 | "Gone For Good: Deindustrialization, White Voter Backlash, and US Presidential Voting" | *APSR*, 115(2): 550-567 | Share de trabalhadores por industria manufatureira em cada condado | Importacoes chinesas por industria (2000-2015) | Brancos votam mais em republicanos em condados com mais demissoes na manufatura; negros em localidades afetadas votam mais em democratas. Desindustrializacao ameaca status do grupo dominante. |
| Feigenbaum & Hall | 2015 | "How Legislators Respond to Localized Economic Shocks" | *J. Politics*, 77(4) | Shares de emprego industrial por distrito congressual | Crescimento de importacoes chinesas por industria | Choques de importacao fazem legisladores votar de forma mais protecionista em projetos de comercio, mas sem efeito em outras votacoes. Sem efeito em reeleicao ou controle partidario. |
| Che, Lu, Pierce, Schott & Tao | 2022 | "Did Trade Liberalization with China Influence US Elections?" | *J. Intl Economics*, 139 | Shares de emprego industrial por condado (1990) | "NTR gap" -- diferenca entre tarifas Smoot-Hawley e tarifas NTR aplicadas por industria | Condados mais expostos a liberalizacao com China migraram para democratas apos 2000, mas efeito diminuiu apos Tea Party em 2010. |

### 3.2 Comercio e nacionalismo na Europa

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Colantone & Stanig | 2018 | "Global Competition and Brexit" | *APSR*, 112(2): 201-218 | Shares de emprego regional (NUTS-2/3) por industria no UK | Importacoes chinesas por industria (instrumentadas por importacoes chinesas nos EUA) | Voto Leave sistematicamente maior em regioes mais expostas a competicao chinesa. Um desvio-padrao no choque = ~2pp a mais de Leave. |
| Colantone & Stanig | 2018 | "The Trade Origins of Economic Nationalism" | *AJPS*, 62(4): 936-953 | Shares de emprego regional (NUTS-2) por industria em 15 paises europeus | Importacoes chinesas por industria (instrumentadas por importacoes chinesas nos EUA) | Maior exposicao ao choque chines aumenta apoio a partidos nacionalistas, shift geral para a direita e apoio a partidos de extrema-direita. Cobre 1988-2007. |
| Dippel, Gold & Heblich | 2022 | "Globalization and Its (Dis-)Content" | *Economic Journal*, 132(641): 199-217 | Shares de emprego industrial por condado alemao | Crescimento de importacoes/exportacoes com China e Leste Europeu por industria | Apenas partidos de extrema-direita respondem significativamente. Competicao de importacao aumenta voto neles; acesso a exportacao diminui. 2/3 do efeito mediado por emprego manufatureiro. |
| Fetzer | 2019 | "Did Austerity Cause Brexit?" | *AER*, 109(11): 3849-3886 | Sensibilidade regional a tipos especificos de gasto em welfare (composicao local de beneficiarios) | Episodios narrativos de consolidacao fiscal nacional (cortes em housing benefit, disability benefits etc.) | Reformas de welfare por austeridade aumentaram apoio ao UKIP em 3.5-11.9pp. Sem austeridade, referendo teria resultado em Remain. |

### 3.3 Imigracao e politica

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Tabellini | 2020 | "Gifts of the Immigrants, Woes of the Natives" | *Rev. Econ. Studies*, 87(1): 454-486 | Share de individuos de cada pais europeu em cada condado dos EUA em 1900 | Fluxos de imigracao por pais de origem (1910-1930), impulsionados por WWI e Immigration Acts | Imigracao gerou reacao politica hostil: eleicao de congressistas mais conservadores e apoio a cotas. Mas aumentou emprego e producao industrial. Backlash cresce com distancia cultural. |
| Halla, Wagner & Zweimuller | 2017 | "Immigration and Voting for the Far Right" | *JEEA*, 15(6): 1341-1378 | Padroes historicos de assentamento de guest workers (Turquia e Yugoslavia) em comunidades austriacas (1960s-70s) | Fluxos nacionais de imigracao para Austria (especialmente 1990s, guerras iugoslavas) | Imigracao aumenta significativamente voto no FPO (extrema-direita). Efeito driven por imigrantes de baixa/media qualificacao; alta qualificacao nao tem efeito. |
| Dustmann, Vasiljeva & Piil Damm | 2019 | "Refugee Migration and Electoral Outcomes" | *Rev. Econ. Studies*, 86(5): 2035-2091 | Alocacao quase-aleatoria de refugiados pelo governo dinamarques | Variacao em coortes de refugiados alocados por ciclo eleitoral | Mais refugiados = mais voto em partidos anti-imigracao, exceto em municipios muito urbanos (efeito invertido). Preocupacoes fiscais dominam. |

### 3.4 Austeridade e extremismo

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Gabriel, Klein & Pessoa | 2026 | "The Political Costs of Austerity" | *Rev. Econ. Statistics*, 108(1): 145 | Sensibilidades regionais a gastos publicos nacionais | Episodios narrativos de consolidacao fiscal nacional | Consolidacao fiscal aumenta voto em partidos extremos (~3pp por 1% de reducao em gasto regional), reduz turnout e aumenta fragmentacao politica. Cobre varios paises europeus. |

### 3.5 Comercio e welfare state historico

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Scheve & Serlin | 2023 | "The German Trade Shock and the Rise of the Neo-Welfare State in Early Twentieth-Century Britain" | *APSR*, 117(2): 557-574 | Shares de emprego industrial por constituency britanica (c. 1880) | Crescimento de importacoes alemas por industria (1880-1910) | Importacoes alemas pioraram mercado de trabalho e mudaram crencas sobre deservingness dos pobres. Menos apoio a conservadores e surgimento de termos como "unemployment" em vez de "vagrancy". Contribuiu para o surgimento do welfare state. |

---

## 4. Aplicacoes em Relacoes Internacionais

### 4.1 Comercio e atitudes de politica externa

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Jensen, Quinn & Weymouth | 2017 | "Winners and Losers in International Trade: The Effects on U.S. Presidential Voting" | ***International Organization***, 71(3): 423-457 | Shares de emprego por industria tradable (bens e servicos) no nivel de condado | Mudancas nacionais em penetracao de importacoes e crescimento de exportacoes por industria | Importacoes reduzem voto no incumbente; exportacoes aumentam. Swing states com manufatura de baixa qualificacao foram particularmente vulneraveis. |
| Ballard-Rosa, Jensen & Scheve | 2022 | "Economic Decline, Social Identity, and Authoritarian Values in the United States" | ***Intl Studies Quarterly***, 66(1): sqab027 | Shares de emprego industrial local em commuting zones | Crescimento de importacoes chinesas por industria (metodologia Autor et al. 2013) | Individuos em regioes diversas enfrentando competicao chinesa intensa tem valores mais autoritarios. Declinio economico afeta identidade social de grupos historicamente dominantes. |
| Ballard-Rosa, Malik, Rickard & Scheve | 2021 | "The Economic Origins of Authoritarian Values: Evidence from Local Trade Shocks in the United Kingdom" | *Comparative Political Studies*, 54(13): 2321-2353 | Shares de emprego industrial local no UK | Crescimento de importacoes chinesas por industria | Regioes britanicas mais afetadas por importacoes chinesas tem valores significativamente mais autoritarios. Mecanismo de frustracao-agressao. Venceu o David A. Lake Award no IPES 2017. |
| Kuk, Seligsohn & Zhang | 2018 | "From Tiananmen to Outsourcing: The Effect of Rising Import Competition on Congressional Voting Towards China" | *J. Contemporary China*, 27(109) | Shares de emprego manufatureiro por distrito | Crescimento de importacoes chinesas por industria | Apos 2003 (pos-acessao a OMC), legisladores de distritos afetados por importacoes chinesas tornaram-se sistematicamente mais hostis a China em votos de politica externa. Congresso agora mais hostil que apos Tiananmen. |
| Campello & Urdinez | 2021 | "Voter and Legislator Responses to Localized Trade Shocks from China in Brazil" | *Comparative Political Studies*, 54(7): 1131-1162 | Shares de emprego/producao municipal por industria | Crescimento de importacoes/exportacoes chinesas por setor | Respostas assimetricas: regioes afetadas por importacoes veem China como risco, mas regioes beneficiadas por exportacoes *nao* veem como oportunidade. Perdas moldam percepcoes de politica externa mais que ganhos. |

### 4.2 Comercio, anti-europeismo e integracao internacional

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Lechler | 2019 | "Employment Shocks and Anti-EU Sentiment" | *European J. Political Economy*, 59: 266-295 | Shares de emprego industrial regional em ~260 regioes NUTS-II europeias | Taxas de crescimento de emprego setorial em toda a Europa | Choques negativos de emprego aumentam sentimento anti-UE e voto eurosceptico em eleicoes para o Parlamento Europeu. Efeitos mais fortes para desempregados e trabalhadores de baixa qualificacao. |

### 4.3 Choques de commodities e conflito armado

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Dube & Vargas | 2013 | "Commodity Price Shocks and Civil Conflict: Evidence from Colombia" | *Rev. Econ. Studies*, 80(4): 1384-1421 | Intensidade de producao municipal para cada commodity (cafe/petroleo) | Mudancas exogenas em precos internacionais de commodities | Queda no preco do cafe aumenta violencia (canal de custo de oportunidade -- salarios menores barateiam recrutamento); aumento no preco do petroleo aumenta violencia (canal de rapacidade -- rendas maiores tornam captura de recursos mais atraente). |
| Gallea | 2023 | "Weapons and War: The Effect of Arms Transfers on Internal Conflict" | *J. Development Economics*, 160 | Porcentagem historica de armas que cada pais africano importava de cada fornecedor | Se cada pais fornecedor esta em guerra *fora* da Africa em dado ano (choque de oferta) | Importacoes de armas na Africa aumentam conflito interno. Aumento de 10% em importacoes de armas eleva risco de conflito em 0.16pp (8% da baseline). Quando um fornecedor entra em guerra, suas exportacoes caem, reduzindo conflito nos receptores. |
| Nobauer | WP | "Trade and Conflict in Myanmar: A Reverse China Shock" | Working paper | Shares de producao/exportacao de bens minerais por township em Myanmar | Crescimento nacional de exportacoes por commodity | Exposicao a exportacoes de mineracao associada a maior conflito violento, predominantemente em townships de minorias etnicas. |

### 4.4 Ajuda internacional e migracao

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Fitchett & Wesselbaum | 2022 | "Does Aid Drive Migration? Evidence from a Shift-Share Instrument" | *International Migration Review*, 56(3) | Padroes historicos bilaterais de alocacao de ajuda (quanto cada doador da a cada receptor como share do total) | Mudancas agregadas nos orcamentos de ajuda dos doadores ao longo do tempo | Relacao positiva: aumento de 10% na ajuda aumenta migracao em ~2%. Forma de U invertido. Contradiz objetivo politico de usar ajuda para reduzir migracao do Sul Global. |

### 4.5 Gasto militar e desigualdade

| Autor(es) | Ano | Titulo | Journal | Shares | Shifts | Achado principal |
|-----------|-----|--------|---------|--------|--------|-----------------|
| Albalate, Bel & Elias | 2020 | "The Effect of Military Spending on Income Inequality: Evidence from NATO Countries" | *Empirical Economics*, 58(3) | Composicao historica de gasto de defesa por pais (por tipo de despesa militar) | Tendencias agregadas de gasto de defesa na OTAN/global | Gasto de defesa causa *reducao* na desigualdade de renda, com elasticidade de curto prazo de -0.075. Conscricao militar tem efeito redistributivo. 14 paises da OTAN, 1977-2007. |

---

## 5. Padroes recorrentes na construcao do instrumento

A literatura revela tres padroes dominantes de construcao de instrumentos shift-share:

### 5.1 Choques comerciais (o mais comum)
- **Shares** = participacao industrial local no emprego
- **Shifts** = crescimento de importacoes (geralmente da China) por industria, instrumentado por importacoes em terceiros paises
- **Uso**: Efeitos politicos de globalizacao, polarizacao, nacionalismo, Brexit, voto de extrema-direita

### 5.2 Imigracao (instrumento de enclave tipo Card)
- **Shares** = padroes historicos de assentamento de imigrantes por pais de origem
- **Shifts** = fluxos nacionais de imigracao por pais, impulsionados por eventos exogenos (guerras, politicas migratoriass)
- **Uso**: Efeitos da imigracao sobre mercado de trabalho e voto

### 5.3 Politica fiscal/austeridade
- **Shares** = sensibilidade regional a tipos especificos de gasto publico (composicao local de beneficiarios)
- **Shifts** = episodios narrativos de consolidacao fiscal nacional
- **Uso**: Efeitos politicos de austeridade

---

## 6. Gaps identificados

### Gaps na Ciencia Politica
- **Pouca aplicacao no Brasil e America Latina**: Alem de Campello & Urdinez (2021) sobre o Brazil e a China, ha pouquissimas aplicacoes shift-share em CP na regiao
- **Politica subnacional**: A maioria dos estudos foca em eleicoes nacionais ou europeias; aplicacoes para eleicoes municipais ou estaduais sao escassas
- **Decomposicao de voto por partido**: O proprio capitulo 15 ja aponta nessa direcao com o exemplo de voto por estado/partido -- essa e uma contribuicao pedagogica original
- **Politica fiscal subnacional**: Choques de receita municipal (royalties, transferencias federais) como instrumento shift-share para estudar politica local

### Gaps em Relacoes Internacionais
- **Conflito e comercio alem de commodities**: Dube & Vargas (2013) e o principal exemplo; poucos estudos usam shift-share para estudar conflito fora do contexto de commodities
- **Sancoes internacionais**: Nenhum estudo encontrado que use shift-share para medir exposicao diferencial a sancoes
- **Cooperacao internacional**: Ausencia quase total de aplicacoes em estudos sobre cooperacao, organizacoes internacionais ou regimes
- **China alem dos EUA e Europa**: Poucos estudos sobre o "China Shock" politico em paises em desenvolvimento alem do Brasil
- **Shift-share para poder/capacidades**: A decomposicao de poder em componentes (militar, economico, cultural) nao e adequada (como notado no cap. 15), mas volume de comercio por parceiro/industria e um candidato viavel

### Gaps metodologicos
- **Inferencia com poucos choques**: A maioria das aplicacoes usa o framework de muitos choques (BHJ 2022), mas muitos contextos de CP/RI tem poucos choques relevantes
- **Dados em painel com shift-share**: Pouca orientacao pedagogica sobre como combinar shift-share com metodos de painel (FE, DiD)
- **Validacao de shares vs. shifts**: Poucos trabalhos aplicados em CP/RI discutem explicitamente qual das duas estrategias de identificacao estao usando

---

## 7. Sugestoes para o capitulo 15

1. **Referencia pedagogica central**: Borusyak, Hull & Jaravel (2025, JEP) e a referencia obrigatoria como guia pratico
2. **Exemplo motivador de CP**: Colantone & Stanig (2018, APSR) sobre Brexit e uma excelente ilustracao para alunos
3. **Exemplo motivador de RI**: Jensen, Quinn & Weymouth (2017, IO) e publicado em journal de RI e usa shift-share de forma clara
4. **Conflito**: Dube & Vargas (2013) e um exemplo elegante de dois mecanismos opostos
5. **A decomposicao de voto proposta no capitulo e original** e pode ser desenvolvida como contribuicao pedagogica propria
6. **Incluir discussao das duas estrategias de identificacao** (exogenous shares vs. exogenous shifts) com exemplos concretos de cada uma

---

## 8. Referencias-chave (formato APSA)

Adao, Rodrigo, Michal Kolesar, and Eduardo Morales. 2019. "Shift-Share Designs: Theory and Inference." *Quarterly Journal of Economics* 134(4): 1949-2010.

Albalate, Daniel, Germa Bel, and Ferran Elias. 2020. "The Effect of Military Spending on Income Inequality: Evidence from NATO Countries." *Empirical Economics* 58(3).

Autor, David H., David Dorn, and Gordon H. Hanson. 2013. "The China Syndrome: Local Labor Market Effects of Import Competition in the United States." *American Economic Review* 103(6): 2121-2168.

Autor, David H., David Dorn, Gordon H. Hanson, and Kaveh Majlesi. 2020. "Importing Political Polarization? The Electoral Consequences of Rising Trade Exposure." *American Economic Review* 110(10): 3139-3183.

Baccini, Leonardo, and Stephen Weymouth. 2021. "Gone For Good: Deindustrialization, White Voter Backlash, and US Presidential Voting." *American Political Science Review* 115(2): 550-567.

Ballard-Rosa, Cameron, Amalie Jensen, and Kenneth Scheve. 2022. "Economic Decline, Social Identity, and Authoritarian Values in the United States." *International Studies Quarterly* 66(1): sqab027.

Ballard-Rosa, Cameron, Mashail Malik, Stephanie Rickard, and Kenneth Scheve. 2021. "The Economic Origins of Authoritarian Values: Evidence From Local Trade Shocks in the United Kingdom." *Comparative Political Studies* 54(13): 2321-2353.

Bartik, Timothy J. 1991. *Who Benefits from State and Local Economic Development Policies?* Kalamazoo, MI: W.E. Upjohn Institute for Employment Research.

Borusyak, Kirill, and Peter Hull. 2023. "Nonrandom Exposure to Exogenous Shocks." *Econometrica* 91(6): 2155-2185.

Borusyak, Kirill, Peter Hull, and Xavier Jaravel. 2022. "Quasi-Experimental Shift-Share Research Designs." *Review of Economic Studies* 89(1): 181-213.

Borusyak, Kirill, Peter Hull, and Xavier Jaravel. 2025. "A Practical Guide to Shift-Share Instruments." *Journal of Economic Perspectives* 39(1): 181-204.

Borusyak, Kirill, Peter Hull, and Xavier Jaravel. 2025. "Design-Based Identification with Formula Instruments: A Review." *The Econometrics Journal* 28(1): 83-108.

Breuer, Matthias. 2022. "Bartik Instruments: An Applied Introduction." *Journal of Financial Reporting* 7(1): 49-67.

Campello, Daniela, and Francisco Urdinez. 2021. "Voter and Legislator Responses to Localized Trade Shocks from China in Brazil." *Comparative Political Studies* 54(7): 1131-1162.

Card, David. 2001. "Immigrant Inflows, Native Outflows, and the Local Labor Market Impacts of Higher Immigration." *Journal of Labor Economics* 19(1): 22-64.

Che, Yi, Yi Lu, Justin Pierce, Peter Schott, and Zhigang Tao. 2022. "Did Trade Liberalization with China Influence US Elections?" *Journal of International Economics* 139: 103652.

Colantone, Italo, and Piero Stanig. 2018a. "Global Competition and Brexit." *American Political Science Review* 112(2): 201-218.

Colantone, Italo, and Piero Stanig. 2018b. "The Trade Origins of Economic Nationalism: Import Competition and Voting Behavior in Western Europe." *American Journal of Political Science* 62(4): 936-953.

Dippel, Christian, Robert Gold, and Stephan Heblich. 2022. "Globalization and Its (Dis-)Content: Trade Shocks and Voting Behavior." *The Economic Journal* 132(641): 199-217.

Dube, Oeindrila, and Juan Vargas. 2013. "Commodity Price Shocks and Civil Conflict: Evidence from Colombia." *Review of Economic Studies* 80(4): 1384-1421.

Dustmann, Christian, Kristine Vasiljeva, and Anna Piil Damm. 2019. "Refugee Migration and Electoral Outcomes." *Review of Economic Studies* 86(5): 2035-2091.

Feigenbaum, James, and Andrew Hall. 2015. "How Legislators Respond to Localized Economic Shocks: Evidence from Chinese Import Competition." *Journal of Politics* 77(4).

Ferri, Benjamin. 2022. "Novel Shift-Share Instruments and Their Applications." Boston College Working Papers in Economics, No. 1053.

Fetzer, Thiemo. 2019. "Did Austerity Cause Brexit?" *American Economic Review* 109(11): 3849-3886.

Fitchett, Hamish, and Dennis Wesselbaum. 2022. "Does Aid Drive Migration? Evidence from a Shift-Share Instrument." *International Migration Review* 56(3).

Gabriel, Ricardo, Mathias Klein, and Ana Sofia Pessoa. 2026. "The Political Costs of Austerity." *Review of Economics and Statistics* 108(1): 145.

Gallea, Quentin. 2023. "Weapons and War: The Effect of Arms Transfers on Internal Conflict." *Journal of Development Economics* 160.

Goldsmith-Pinkham, Paul, Isaac Sorkin, and Henry Swift. 2020. "Bartik Instruments: What, When, Why, and How." *American Economic Review* 110(8): 2586-2624.

Halla, Martin, Alexander Wagner, and Josef Zweimuller. 2017. "Immigration and Voting for the Far Right." *Journal of the European Economic Association* 15(6): 1341-1378.

Jaeger, David A., Joakim Ruist, and Jan Stuhler. 2018. "Shift-Share Instruments and the Impact of Immigration." NBER Working Paper 24285.

Jensen, J. Bradford, Dennis Quinn, and Stephen Weymouth. 2017. "Winners and Losers in International Trade: The Effects on U.S. Presidential Voting." *International Organization* 71(3): 423-457.

Kuk, John Seungmin, Deborah Seligsohn, and Jiakun Jack Zhang. 2018. "From Tiananmen to Outsourcing: The Effect of Rising Import Competition on Congressional Voting Towards China." *Journal of Contemporary China* 27(109).

Lechler, Marie. 2019. "Employment Shocks and Anti-EU Sentiment." *European Journal of Political Economy* 59: 266-295.

Malgouyres, Clement. 2017. "Trade Shocks and Far-Right Voting: Evidence from French Presidential Elections." European University Institute Working Paper RSCAS 2017/21.

Margalit, Yotam. 2011. "Costly Jobs: Trade-related Layoffs, Government Compensation, and Voting in U.S. Elections." *American Political Science Review* 105(1): 166-188.

Scheve, Kenneth, and Theo Serlin. 2023. "The German Trade Shock and the Rise of the Neo-Welfare State in Early Twentieth-Century Britain." *American Political Science Review* 117(2): 557-574.

Tabellini, Marco. 2020. "Gifts of the Immigrants, Woes of the Natives: Lessons from the Age of Mass Migration." *Review of Economic Studies* 87(1): 454-486.
