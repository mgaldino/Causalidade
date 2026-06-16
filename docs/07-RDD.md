# Desenho de Regressão Descontínua

## Roteiro da aula

Na aula de hoje, iremos aprender como desenhos de regressão descontínua usam regras com pontos de corte para identificar efeitos causais locais.

A aula terá uma ênfase aplicada. Primeiro, vamos construir a intuição do desenho, distinguir *sharp RD* de *fuzzy RD* e conectar fuzzy RD à aula anterior de variáveis instrumentais. Depois, veremos como estimar RD em R com `rdrobust`, como interpretar bandwidths automáticos e como checar a plausibilidade do desenho. Ao final, discutiremos casos difíceis: eleições apertadas, variáveis de ordenação discretas, múltiplos cutoffs, RD espacial e diagnósticos de manipulação.

## Características-chave da RDD

A Regressão Descontínua (RDD) é caracterizada por uma variável de ordenação $X_i$ que determina, ou altera fortemente, a chance de receber tratamento em um ponto de corte conhecido $c$. Por convenção, $X_i$ é chamada de *running variable*, *assignment variable* ou *forcing variable*. Como no restante do livro, o tratamento é denotado por $D_i$, com $D_i=1$ para tratados e $D_i=0$ para controles.

### Determinação do Tratamento

Em um desenho RDD *sharp*, uma unidade é tratada se $X_i \geq c$ e não tratada se $X_i < c$. Assim, $D_i$ é uma função determinística de $X_i$: $D_i = 1(X_i \geq c)$. A *running variable* determina completamente quem recebe tratamento.

### Observação e Corte

É essencial observar $X$ e conhecer o ponto de corte, ou limiar, $c$.

O desenho não exige que a *running variable* seja perfeitamente contínua em sentido literal. A análise precisa de suporte suficiente dos dois lados do cutoff para estimar limites laterais, além de continuidade dos resultados potenciais condicionais ao redor de $c$. *Running variables* discretas não invalidam automaticamente um RD, mas mudam o tipo de evidência disponível: com poucos valores distintos perto do cutoff, a extrapolação local fica mais forte, os testes de densidade perdem interpretação simples e a inferência precisa ser mais cautelosa.

A suposição de identificação é que os resultados potenciais $Y_i(0)$ e $Y_i(1)$ variem suavemente ao redor do ponto de corte. Essa suposição não é diretamente testável, pois nunca observamos os dois resultados potenciais da mesma unidade. Lee (2008) enfatiza uma condição suficiente: as unidades não podem manipular a *running variable* com precisão suficiente para escolher em qual lado do cutoff ficar. Isso implica que covariáveis de pré-tratamento deveriam se comportar suavemente ao redor do cutoff. Essa implicação é diagnosticável nas variáveis observadas, mas não prova a suposição para variáveis não observadas.

### Estimativa dos Efeitos do Tratamento

No sharp RD, o estimando causal é um ATE local no cutoff:

$$
\tau_{SRD} = E[Y_i(1)-Y_i(0)\mid X_i=c].
$$

Esse estimando é local porque se refere às unidades no ponto de corte, não à população inteira. Sob continuidade dos resultados potenciais, podemos identificá-lo pelo salto no resultado observado em $c$. Escreverei os limites laterais como $x \to c^-$ para aproximação pela esquerda e $x \to c^+$ para aproximação pela direita:

Parte da literatura chama esse objeto de LATE porque ele é local em $X_i=c$. Neste capítulo, usarei "ATE local no cutoff" para o sharp RD. A escolha evita confundir esse estimando com o LATE de variáveis instrumentais, que reaparece no fuzzy RD como efeito para compliers.

$$
\tau_{SRD} = \lim_{x \to c^+} E[Y_i|X_i=x] - \lim_{x \to c^-} E[Y_i|X_i=x].
$$

No desenho *sharp*, essa comparação é equivalente a comparar $\lim_{x \to c^-} E[Y_i | X_i = x, D_i=0]$ com $\lim_{x \to c^+} E[Y_i | X_i = x, D_i=1]$, porque à direita de $c$ todos recebem tratamento e à esquerda ninguém recebe. Sob continuidade dos resultados potenciais:

- $\lim_{x \to c^-} E[Y_i | X_i = x] \approx E[Y_i(0) | X_i = c]$
- $\lim_{x \to c^+} E[Y_i | X_i = x] \approx E[Y_i(1) | X_i = c]$

Se fôssemos usar regressão linear, o modelo seria:

$$
Y_i = \alpha + \tau D_i + \beta_1(X_i-c) + \beta_2D_i(X_i-c) + \epsilon_i,
$$

em que $D_i = 1(X_i \geq c)$. Ao centralizar a *running variable* em $c$, $\tau$ passa a ser diretamente o salto estimado no cutoff.

## Fuzzy RD como IV local

Depois de entender o sharp RD, a extensão para fuzzy RD é direta. Pode acontecer de o ponto de corte não determinar perfeitamente quem recebe tratamento, mas apenas aumentar a probabilidade de tratamento. Nesse caso, a regra de elegibilidade funciona como um instrumento local.

Defina $Z_i = 1(X_i \geq c)$. Em um fuzzy RD, $Z_i$ afeta o tratamento efetivamente recebido $D_i$, mas nem todos cumprem a regra. A estimativa é um Wald local:

$$
\tau_{FRD} =
\frac{\lim_{x \to c^+} E[Y_i|X_i=x] - \lim_{x \to c^-} E[Y_i|X_i=x]}
{\lim_{x \to c^+} E[D_i|X_i=x] - \lim_{x \to c^-} E[D_i|X_i=x]}.
$$

\begin{table}

\caption{(\#tab:tabela-ponte-fuzzy-iv)Ponte entre variáveis instrumentais e fuzzy RD.}
\centering
\begin{tabular}[t]{l|l|l}
\hline
Objeto & Na aula de IV & No fuzzy RD\\
\hline
$Z_i$ & Instrumento & Elegibilidade gerada pelo cutoff\\
\hline
$D_i$ & Tratamento recebido & Tratamento efetivamente recebido\\
\hline
Numerador & Forma reduzida: efeito de Z sobre Y & Salto local no resultado em c\\
\hline
Denominador & Primeira etapa: efeito de Z sobre D & Salto local na probabilidade de tratamento em c\\
\hline
Estimando & LATE para compliers & LATE dos compliers no cutoff\\
\hline
\end{tabular}
\end{table}

A interpretação é um LATE no cutoff: o efeito médio local para as unidades cujo tratamento muda por causa da regra de elegibilidade. Como em IV, a interpretação exige relevância da primeira etapa, monotonicidade, exclusão e continuidade dos resultados potenciais e do tratamento potencial ao redor do ponto de corte.

O fuzzy RD merece atenção especial porque a inferência usual pode se comportar mal em três situações: primeira etapa fraca, *running variable* discreta e *donut designs*. Um *donut design* exclui uma pequena faixa muito próxima ao cutoff, por exemplo $X_i \in (c-\varepsilon, c+\varepsilon)$, para evitar observações potencialmente manipuladas, arredondadas ou contaminadas exatamente no limiar.

Noack e Rothe (2024) propõem inferência *bias-aware* para fuzzy RD. A intuição é simples: no fuzzy RD, estimamos uma razão entre dois saltos. Se o salto no tratamento é pequeno ou estimado com muita incerteza, a razão pode ficar instável. A inferência precisa reconhecer essa instabilidade, em vez de tratar a primeira etapa como se fosse forte e precisamente estimada.

Exemplo em ciência política: uma regra eleitoral pode definir o número mínimo de votos para obter cadeiras, mas migração partidária, coligações ou decisões administrativas podem fazer com que a regra altere fortemente, mas não determine perfeitamente, o tratamento substantivo de interesse.

## Identificação: continuidade e manipulação

A suposição de continuidade faz o desenho funcionar. Em uma aplicação de *close elections*, por exemplo, candidatos que vencem por margem mínima recebem o tratamento de incumbência. Habilidades, carisma e recursos de campanha podem afetar tanto a vitória quanto resultados futuros. O RD exige que essas características não mudem de forma descontínua exatamente no cutoff de 50%. O que salta no cutoff é o status eleitoral, não a qualidade latente dos candidatos.

A condição de Lee (2008) fornece uma forma prática de pensar esse problema. As unidades não podem manipular $X_i$ com precisão suficiente para escolher de que lado do cutoff ficar. Se conseguem, a comparação local deixa de parecer plausível: unidades logo acima e logo abaixo de $c$ podem diferir em características observadas e não observadas.

Parte dessa suposição deixa rastros nos dados. Podemos procurar acúmulo de observações perto do cutoff e testar se covariáveis pré-tratamento saltam em $c$. Esses diagnósticos são úteis, mas não encerram o argumento. A validade do desenho também depende de conhecimento institucional sobre quem poderia manipular $X_i$, com que precisão e em que direção.

Tendo estabelecido as suposições de identificação, passamos agora à questão prática: como estimar o efeito do tratamento em um RDD?

## Dois frameworks para RD

A literatura recente organiza a análise de RD em dois frameworks principais: o framework de continuidade e o framework de aleatorização local [Cattaneo e Titiunik 2022; Cattaneo, Idrobo e Titiunik 2024].

No framework de continuidade, o estimando é a diferença entre limites laterais no cutoff. A validade depende da suavidade dos resultados potenciais ao redor de $c$. Esse é o framework padrão por trás de métodos como regressão local linear, bandwidths MSE/CER-optimal e inferência robusta à correção de viés.

No framework de aleatorização local, escolhemos uma janela pequena ao redor do cutoff e tratamos a atribuição ao tratamento dentro dessa janela como se fosse aproximadamente aleatória. Essa abordagem pode ser útil quando há poucos valores discretos da *running variable* perto do cutoff ou quando temos forte conhecimento institucional de que pequenas diferenças no score são essencialmente acidentais. Essa interpretação exige uma janela substantivamente defensável e diagnósticos de balanceamento. Ela não decorre automaticamente do bandwidth escolhido por `rdrobust`.

## Estimação em RDD

RDD não tem sobreposição no sentido usado em matching: em um *sharp RD*, unidades abaixo do cutoff não recebem tratamento e unidades acima recebem. Isso não é um defeito do desenho; é exatamente a regra que gera identificação. A dificuldade de estimação é outra: precisamos estimar duas funções condicionais nos limites laterais do cutoff.

Por isso, RD é um problema de estimação de fronteira. Quanto mais usamos observações distantes do cutoff, mais dependemos de suavidade e forma funcional. Quanto mais restringimos a amostra a observações muito próximas, menos viés introduzimos, mas maior fica a incerteza estatística.

O método padrão moderno é regressão polinomial local, em geral local linear, estimada separadamente dos dois lados do cutoff e ponderando mais as observações próximas de $c$. Na prática, usamos `rdrobust`, que implementa seleção de bandwidth, estimação local polinomial e inferência robusta com correção de viés.

A identificação dos efeitos do tratamento ocorre no limite, à medida que $X_i \rightarrow c$. Quanto mais usarmos observações distantes de $c$ em $X$, mais dependeremos de extrapolação e das suposições sobre a forma funcional.

### Bandwidth e inferência

Janelas menores usam observações mais próximas de $c$. Elas reduzem a dependência de forma funcional e o viés potencial, mas também reduzem o tamanho amostral efetivo e aumentam a variância.

Janelas maiores usam mais observações. Elas aumentam a precisão estatística, mas dependem mais de suavidade e de escolhas de forma funcional.

A ideia é restringir a estimativa a uma janela ao redor de $X_i = c$, que pode ter tamanhos diferentes à esquerda ou à direita. Estes métodos buscam equilibrar a precisão das estimativas minimizando viés e variância conforme a proximidade do ponto de corte $c$.

No pacote `rdrobust`, a largura de banda pode ser escolhida automaticamente. Duas escolhas aparecem com frequência:

1. MSE-optimal bandwidth, em que MSE significa *Mean Squared Error* ou erro quadrático médio: escolhe a largura de banda para minimizar o erro quadrático médio do estimador. Em termos práticos, busca um bom equilíbrio entre viés e variância para a estimativa pontual.

2. CER-optimal bandwidth, em que CER significa *Coverage Error Rate* ou taxa de erro de cobertura: escolhe a largura de banda para melhorar a cobertura dos intervalos de confiança, especialmente quando usamos inferência com correção robusta de viés (*robust bias-corrected inference*).

Essas larguras de banda automáticas pertencem ao framework de continuidade, não ao framework de aleatorização local. Se `rdrobust` escolhe, por exemplo, uma janela de 10 pontos percentuais para cada lado do cutoff em uma aplicação com eleições apertadas, isso não significa que todas as eleições vencidas ou perdidas por até 10 pontos percentuais possam ser interpretadas como se fossem aleatórias. Uma janela desse tamanho pode ser útil para estimar os limites laterais de $E[Y|X=x]$ sob continuidade/suavidade dos resultados potenciais, mas ela dificilmente sustenta a interpretação literal de um experimento local.

Portanto, a interpretação correta depende do framework. No framework de continuidade, usamos observações próximas ao cutoff, ponderadas pela distância, para estimar limites laterais. No framework de aleatorização local, precisamos justificar substantivamente uma janela em que a atribuição ao tratamento seja plausivelmente "como se aleatória"; essa janela deve ser defendida com conhecimento institucional e diagnósticos de balanceamento, não simplesmente herdada do bandwidth automático de `rdrobust`.

## Regras arbitrárias

RDDs aparecem quando uma regra institucional usa um limiar para mudar acesso, elegibilidade ou intensidade de tratamento. Programas de transferência de renda podem depender de renda familiar; aprovação no ensino superior pode depender de uma nota mínima; políticas ambientais podem mudar quando a propriedade ultrapassa certo tamanho; regras eleitorais podem depender de população municipal, margem de vitória ou número de votos. O ponto comum é que a regra cria uma mudança discreta em $D_i$ quando $X_i$ cruza $c$.

## Simulação


``` r
set.seed(123)
N <- 1000
X <- runif(N, -5, 5)
Y0 <- rnorm(n = N, mean = X, sd = 1)       # resultado potencial sob D = 0
Y1 <- rnorm(n = N, mean = X + 2, sd = 1)   # resultado potencial sob D = 1
D <- as.integer(X >= 0)
Y <- Y1 * D + Y0 * (1 - D)
```




\begin{figure}
\includegraphics[alt={Gráfico em degrau mostrando D igual a zero à esquerda do cutoff e igual a um à direita.}]{07-RDD_files/figure-latex/plot-d-assignment-1} \caption{Regra sharp RD: o tratamento D muda mecanicamente no cutoff X = 0.}(\#fig:plot-d-assignment)
\end{figure}

Começamos olhando $Y_i(0)$, o resultado que cada unidade teria sem tratamento.


\begin{figure}
\includegraphics[alt={Dispersão de Y(0) contra X sem salto visível no cutoff.}]{07-RDD_files/figure-latex/plot-po-y0-1} \caption{Resultado potencial sob controle: Y(0) varia suavemente ao redor do cutoff.}(\#fig:plot-po-y0)
\end{figure}

Agora olhamos $Y_i(1)$, o resultado que cada unidade teria sob tratamento.


\begin{figure}
\includegraphics[alt={Dispersão de Y(1) contra X sem salto abrupto no cutoff.}]{07-RDD_files/figure-latex/plot-po-y1-1} \caption{Resultado potencial sob tratamento: Y(1) também é suave no cutoff.}(\#fig:plot-po-y1)
\end{figure}

Colocar $Y_i(0)$ e $Y_i(1)$ no mesmo gráfico mostra a lógica do desenho: os resultados potenciais são suaves, mas só observamos um deles para cada unidade.


\begin{figure}
\includegraphics[alt={Dispersão dos resultados potenciais Y(0) e Y(1) por lado do cutoff.}]{07-RDD_files/figure-latex/plot-po-y1-Y0-1} \caption{Resultados potenciais simulados por status de tratamento observado.}(\#fig:plot-po-y1-Y0-1)
\end{figure}
\begin{figure}
\includegraphics[alt={Curvas suavizadas de Y(0) e Y(1) contra X, ambas contínuas no cutoff.}]{07-RDD_files/figure-latex/plot-po-y1-Y0-2} \caption{Funções suavizadas de Y(0) e Y(1): o efeito de RD vem do salto no tratamento, não de uma quebra nos resultados potenciais.}(\#fig:plot-po-y1-Y0-2)
\end{figure}

O resultado observado combina esses dois resultados potenciais de acordo com a regra de tratamento.
\begin{figure}
\includegraphics[alt={Dispersão do resultado observado contra X com mudança de nível no cutoff.}]{07-RDD_files/figure-latex/plot-observed-1} \caption{Resultado observado: o salto em Y aparece porque o tratamento muda no cutoff.}(\#fig:plot-observed)
\end{figure}

A simulação acima ilustra os elementos fundamentais de um RDD. Vejamos agora em que condições as estimativas de RDD são válidas.

## Quando o RDD funciona?

A suposição chave para o RDD não é que o resultado salte no cutoff. O resultado pode não saltar se o efeito causal for zero. O que precisa saltar é o tratamento, ou a probabilidade de tratamento no caso fuzzy. Para que o salto no resultado possa ser interpretado causalmente, os resultados potenciais e covariáveis de pré-tratamento devem variar suavemente ao redor do cutoff. Vamos ver o que isso significa comparando cinco gráficos: alguns em que a estimativa de RD é válida, ainda que com diferentes graus de validade externa, e um em que a descontinuidade em outra variável torna o desenho inválido.



\begin{figure}
\includegraphics[alt={Curvas suavizadas de Y(0) e Y(1) sem quebra no cutoff.}]{07-RDD_files/figure-latex/caso-rdd-suave-1} \caption{Caso 1: resultados potenciais suaves e salto claro no tratamento.}(\#fig:caso-rdd-suave)
\end{figure}

\begin{figure}
\includegraphics[alt={Curvas de resultados potenciais com heterogeneidade e continuidade em zero.}]{07-RDD_files/figure-latex/caso-rdd-heterogeneo-1} \caption{Caso 2: efeitos heterogêneos, mas resultados potenciais ainda suaves no cutoff.}(\#fig:caso-rdd-heterogeneo)
\end{figure}

\begin{figure}
\includegraphics[alt={Curvas de resultados potenciais com heterogeneidade forte e continuidade no cutoff.}]{07-RDD_files/figure-latex/caso-rdd-heterogeneo-forte-1} \caption{Caso 3: heterogeneidade mais forte, com continuidade preservada no cutoff.}(\#fig:caso-rdd-heterogeneo-forte)
\end{figure}

\begin{figure}
\includegraphics[alt={Curvas de Y(0) e Y(1) que se encontram no cutoff.}]{07-RDD_files/figure-latex/caso-rdd-efeito-nulo-1} \caption{Caso 4: nenhum efeito exatamente no cutoff, apesar de as funções diferirem longe dele.}(\#fig:caso-rdd-efeito-nulo)
\end{figure}

\begin{figure}
\includegraphics[alt={Curva de Y(0) com salto artificial ao redor do cutoff, invalidando o desenho.}]{07-RDD_files/figure-latex/caso-rdd-invalido-1} \caption{Caso 5: desenho inválido, pois Y(0) apresenta uma descontinuidade no cutoff.}(\#fig:caso-rdd-invalido)
\end{figure}

Os quatro primeiros casos preservam a continuidade dos resultados potenciais no cutoff. Eles diferem em heterogeneidade e validade externa, mas a comparação local continua interpretável. O quinto caso quebra a lógica do desenho: mesmo sem tratamento, $Y_i(0)$ mudaria abruptamente em $c$.




## Dados brutos versus bins

\begin{figure}
\includegraphics[alt={Gráfico de médias de Y por bins igualmente espaçados da running variable.}]{07-RDD_files/figure-latex/plot-binscatter-1} \caption{Binscatter com bins igualmente espaçados: os pontos mostram médias locais de Y.}(\#fig:plot-binscatter)
\end{figure}

Como escolher os bins?

1. Espaçamentos iguais ou quantis?
2. Quantos bins?

No exemplo, escolhi espaçamento igual e 20 bins. Podemos usar quantis.

\begin{figure}
\includegraphics[alt={Gráfico de médias de Y por bins definidos por quantis da running variable.}]{07-RDD_files/figure-latex/plot-binscatter-quantis-1} \caption{Binscatter com bins por quantis: cada bin contém número semelhante de observações.}(\#fig:plot-binscatter-quantis)
\end{figure}

Esses gráficos servem para visualização, não para estimação. A escolha dos bins afeta como a evidência aparece no gráfico, mas a estimativa principal deve vir da regressão local com bandwidth explícito. Bins por quantis ajudam a evitar bins vazios e tornam os pontos visualmente mais estáveis, mas também podem esconder variações na densidade da *running variable*. Por isso, o gráfico de RD deve ser acompanhado por diagnósticos de densidade.

Sobre o número de bins, Cattaneo et al. (2020) discutem duas abordagens. A primeira minimiza o IMSE, *Integrated Mean Squared Error*, isto é, o erro quadrático médio integrado ao longo da função estimada. A segunda busca imitar a variância da nuvem de pontos original, o que costuma produzir mais bins.

O pacote `rdplot` implementa essas escolhas automaticamente.

\begin{figure}
\includegraphics[alt={RD plot com pontos por bins e ajuste local nos dois lados do cutoff.}]{07-RDD_files/figure-latex/plot-binscatter-quantis-rdplot-qs-1} \caption{RD plot com seleção automática de bins por quantis.}(\#fig:plot-binscatter-quantis-rdplot-qs)
\end{figure}

\begin{figure}
\includegraphics[alt={RD plot com bins definidos para aproximar a variância visual dos dados brutos.}]{07-RDD_files/figure-latex/plot-binscatter-quantis-rdplot-qsmv-1} \caption{RD plot com seleção de bins mimicking-variance.}(\#fig:plot-binscatter-quantis-rdplot-qsmv)
\end{figure}

## Testes de permutação para balanceamento

Para checar balanceamento, podemos usar testes de permutação em covariáveis pré-tratamento. No exemplo abaixo, `w_pre` representa uma covariável observável medida antes do tratamento. O teste pergunta se sua distribuição muda de forma anormal no cutoff.


``` r
library(RATest)
df_perm <- df_u
resultado <- RDperm(
  W = "w_pre",
  z = "x",
  data = df_perm,
  cutoff = 0
)

df_perm %>%
  mutate(lado = ifelse(x < 0, "À esquerda do cutoff", "À direita do cutoff")) %>%
  ggplot(aes(x = w_pre, colour = lado, linetype = lado)) +
  stat_ecdf(linewidth = 0.9) +
  scale_colour_manual(values = c("À esquerda do cutoff" = "#0072B2", "À direita do cutoff" = "#D55E00")) +
  scale_linetype_manual(values = c("À esquerda do cutoff" = "solid", "À direita do cutoff" = "longdash")) +
  labs(
    x = "Covariável pré-tratamento",
    y = "Distribuição acumulada empírica",
    colour = "Lado",
    linetype = "Lado"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")
```

\begin{figure}
\includegraphics[alt={Funções de distribuição acumulada empíricas de uma covariável pré-tratamento à esquerda e à direita do cutoff.}]{07-RDD_files/figure-latex/permutation-test-1} \caption{Teste de permutação para uma covariável pré-tratamento observável.}(\#fig:permutation-test)
\end{figure}

Canay & Kamat (2018) utilizaram esse teste para revisitar o trabalho de Lee (2008) e descobriram que havia problema de balanceamento.

Caughey and Sekhon (2011) na political analysis mostraram que de fato havia problemas de balanceamento no estudo de Lee (2008).

Testes de permutação são mais naturais no framework de aleatorização local. Eles exigem que o pesquisador escolha uma janela em que a atribuição ao tratamento seja plausivelmente "como se aleatória". Essa janela não deve ser simplesmente o bandwidth MSE/CER-optimal de `rdrobust`: bandwidths automáticos estimam limites laterais sob continuidade, enquanto a aleatorização local é uma afirmação mais forte sobre o mecanismo de atribuição dentro de uma janela.



Do artigo na *Political Analysis*:

\begin{figure}
\includegraphics[width=0.9\linewidth,alt={Trecho do artigo de Caughey e Sekhon sobre problemas de balanceamento.}]{imagens/quote} \caption{Caughey e Sekhon (2011): alerta sobre balanceamento em eleições apertadas.}(\#fig:quote-1)
\end{figure}
\begin{figure}
\includegraphics[width=0.9\linewidth,alt={Segundo trecho do artigo de Caughey e Sekhon sobre close elections.}]{imagens/quote2} \caption{Caughey e Sekhon (2011): implicação para a interpretação de close elections.}(\#fig:quote-2)
\end{figure}

Houve um debate na ciência política sobre isso. Erikson & Rader (2017) e Cuesta & Imai (2016) argumentam que o RDD é identificado. Até onde eu sei, cientistas políticos não revisitaram a controvérsia com os novos métodos desenvolvidos pelos economistas.

De todo modo, as evidências de De Magalhães et al. (2025) sugerem que a recomendação que estou adotando no curso de quais práticas usar são as melhores e mais robustas.



## Teste de McCrary e densidade

Um dos principais desafios à identificação causal em RDDs é a possibilidade de manipulação por parte dos agentes sobre ficar acima ou abaixo do ponto de corte. A lógica esperada é que se o tratamento é desejável, indivíduos tentarão receber o tratamento, levando a um gap justamente abaixo do ponto de corte. Se o tratamento é indesejável (efeitos negativos), indivíduos vão evitar o tratamento, levando a um gap justamente acima do ponto de corte.

O exemplo mais evidente para nós cientistas políticos é a aprovação de um projeto de lei no legislativo. Nós sabemos que os legisladores agem estrategicamente retirando propostas que não vão ser aprovadas ou postergando a votação, até terem a maioria, ainda que por margem mínima. Nesse caso, a aplicação de RDD produzirá estimativas viesadas. McCrary, em um artigo de 2008, argumentou que tais casos apareceriam como descontinuidade na densidade da *running variable* ao redor do ponto de corte. Eis o gráfico feito por McCrary em seu estudo original:

\begin{figure}
\includegraphics[width=1\linewidth,alt={Histograma e densidades laterais indicando possível descontinuidade na densidade da running variable.}]{imagens/manipulation test} \caption{Figura clássica do teste de manipulação em McCrary (2008).}(\#fig:McCay-test)
\end{figure}

Para formalizar essa ideia, McCrary estima os limites da densidade pela esquerda e pela direita e avalia se a diferença (do logaritmo) das estimativas é estatisticamente significante diferente de zero. Portanto, rejeitar a hipótese nula é encontrar evidências de que há manipulação. Cattaneo, Jansson e Ma (2018, 2020) introduziram uma versão alternativa do teste, com espírito similar.

Na prática, de um ponto de vista retórico, o que pesquisadores querem é falhar em rejeitar a nula. Como o teste pode ter baixo poder, ausência de evidência não quer dizer evidência de ausência. Hartman (2021) e Fitzgerald (2025) defendem uma leitura mais exigente: em vez de apenas perguntar se rejeitamos manipulação, deveríamos perguntar se os dados permitem descartar manipulação substantivamente relevante. Isso leva a testes de equivalência, em que o pesquisador define uma margem de manipulação considerada pequena o bastante para não comprometer o desenho.

\begin{figure}
\includegraphics[alt={Densidades estimadas da variável de ordenação à esquerda e à direita do cutoff.}]{07-RDD_files/figure-latex/McCay-test1-1} \caption{Diagnóstico visual de densidade: não há acúmulo evidente de observações no cutoff.}(\#fig:McCay-test1)
\end{figure}

Uma implementação mais moderna é o teste de densidade de Cattaneo, Jansson e Ma. A tabela resume o diagnóstico em vez de despejar todo o output do pacote.

\begin{table}

\caption{(\#tab:Cattaneo-test)Teste de densidade de Cattaneo, Jansson e Ma para os dados simulados.}
\centering
\begin{tabular}[t]{l|r}
\hline
Medida & Valor\\
\hline
Diferença estimada na densidade & 0.011\\
\hline
Estatística t & 0.279\\
\hline
p-valor & 0.780\\
\hline
\end{tabular}
\end{table}

Na Tabela \@ref(tab:Cattaneo-test), a hipótese nula é continuidade da densidade de $X_i$ no cutoff. A diferença estimada na densidade é 0.011, isto é, praticamente zero neste exemplo. A estatística $t$ é 0.279 e o p-valor é 0.78. Como o p-valor é alto, não rejeitamos a hipótese nula de continuidade da densidade. A interpretação substantiva é que, nesses dados simulados, não há evidência de acúmulo ou buraco de observações exatamente no cutoff. Isso não prova ausência de manipulação; apenas indica que o diagnóstico não encontrou uma descontinuidade detectável na densidade.

Essa é uma área ativa de pesquisa. Para a aula, o ponto é simples: `rddensity` é um diagnóstico útil, mas não substitui conhecimento institucional sobre quem poderia manipular $X$, com que precisão e em que direção.


## Robustez e inferência em R

A análise aplicada deve separar três objetos: o gráfico, a estimativa e a inferência. O gráfico ajuda a avaliar a plausibilidade da descontinuidade. A estimativa vem da regressão local. A inferência moderna deve usar os intervalos robustos com correção de viés de `rdrobust`, não apenas o erro padrão convencional.


``` r
library(rdrobust)

df <- df_aux

# Estimativa sharp RD com bandwidth MSE-optimal.
modelo_mse <- rdrobust(y = df$y, x = df$x, c = 0, p = 1, bwselect = "mserd")

# Mesmo desenho, mas com bandwidth CER-optimal.
modelo_cer <- rdrobust(y = df$y, x = df$x, c = 0, p = 1, bwselect = "cerrd")
```

No output de `rdrobust`, a linha Conventional usa a estimativa e o erro padrão convencionais. A linha Bias-Corrected corrige a estimativa pontual para viés. A linha Robust combina a estimativa corrigida com erro padrão apropriado para inferência robusta. Para reportar resultados principais, o padrão moderno é dar destaque ao intervalo robusto. A tabela abaixo mostra apenas os elementos que interessam para a aula.

\begin{table}

\caption{(\#tab:tabela-rdrobust)Estimativas sharp RD com inferência robust bias-corrected.}
\centering
\begin{tabular}[t]{l|r|r|r|r|r}
\hline
Especificação & Estimativa & IC 95\% inferior & IC 95\% superior & h esquerda & h direita\\
\hline
MSE-optimal & 2.411 & 1.785 & 3.037 & 1.172 & 1.172\\
\hline
CER-optimal & 2.515 & 1.821 & 3.208 & 0.830 & 0.830\\
\hline
\end{tabular}
\end{table}

\begin{figure}
\includegraphics[alt={Gráfico de pontos com intervalos de confiança para especificações MSE, CER e bandwidths fixos.}]{07-RDD_files/figure-latex/plot-efeitos-1} \caption{Sensibilidade da estimativa RD à escolha de bandwidth.}(\#fig:plot-efeitos)
\end{figure}

Esse gráfico é mais informativo do que uma sequência de gráficos separados porque coloca lado a lado a estimativa principal e alternativas razoáveis. Se a conclusão substantiva muda drasticamente com pequenas mudanças no bandwidth, o desenho precisa ser apresentado com mais cautela.

### Covariáveis em RD

Covariáveis de pré-tratamento podem melhorar precisão e ajudar nos diagnósticos, mas não salvam a identificação se a continuidade dos resultados potenciais for implausível. A literatura recente sobre covariáveis e aprendizado de máquina em RD reforça esse ponto: covariáveis devem ser vistas como ajuste de precisão e ferramenta de diagnóstico, não como substituto do desenho [Kreiß e Rothe 2023; Noack, Olma e Rothe 2025].


``` r
set.seed(1234)
df$w_pre <- 0.3 * df$x + rnorm(nrow(df))

modelo_sem_cov <- rdrobust(y = df$y, x = df$x, c = 0, p = 1)
modelo_com_cov <- rdrobust(
  y = df$y,
  x = df$x,
  c = 0,
  p = 1,
  covs = as.matrix(df$w_pre)
)
```

\begin{table}

\caption{(\#tab:tabela-rd-covs)Covariáveis pré-tratamento podem mudar precisão, mas não substituem a identificação.}
\centering
\begin{tabular}[t]{l|r|r|r|r|r}
\hline
especificacao & Estimativa & ic\_95\_inf & ic\_95\_sup & h\_esquerda & h\_direita\\
\hline
Sem covariável & 2.411 & 1.785 & 3.037 & 1.172 & 1.172\\
\hline
Com covariável pré-tratamento & 2.432 & 1.790 & 3.075 & 1.077 & 1.077\\
\hline
\end{tabular}
\end{table}

### Fuzzy RD em R

No fuzzy RD, o argumento `fuzzy` recebe a variável de tratamento efetivamente recebido. O pacote estima a razão entre o salto no resultado e o salto na probabilidade de tratamento.


``` r
set.seed(4321)
prob_tratamento <- ifelse(df$x >= 0, 0.8, 0.2)
d_fuzzy <- rbinom(nrow(df), size = 1, prob = prob_tratamento)
y_fuzzy <- df$y0 + 2 * d_fuzzy + rnorm(nrow(df), sd = 0.5)

modelo_fuzzy <- rdrobust(
  y = y_fuzzy,
  x = df$x,
  c = 0,
  fuzzy = d_fuzzy
)

modelo_rf_fuzzy <- rdrobust(y = y_fuzzy, x = df$x, c = 0)
modelo_primeira_fuzzy <- rdrobust(y = d_fuzzy, x = df$x, c = 0)
```

\begin{table}

\caption{(\#tab:tabela-fuzzy-rd-r)Fuzzy RD em R: forma reduzida, primeira etapa e razão de Wald local.}
\centering
\begin{tabular}[t]{l|r|r|r}
\hline
Componente & Estimativa & IC 95\% inferior & IC 95\% superior\\
\hline
Forma reduzida: salto em Y & 1.394 & 0.713 & 2.075\\
\hline
Primeira etapa: salto em D & 0.611 & 0.400 & 0.822\\
\hline
Wald local: salto em Y / salto em D & 2.167 & 1.221 & 3.114\\
\hline
\end{tabular}
\end{table}

A Tabela \@ref(tab:tabela-fuzzy-rd-r) deve ser lida em três passos. A primeira linha é a forma reduzida: a elegibilidade gerada pelo cutoff aumenta $Y$ em 1.394 unidades. A segunda linha é a primeira etapa: a elegibilidade aumenta a probabilidade de receber tratamento em 0.611 pontos de probabilidade. A terceira linha divide esses dois saltos. Portanto, o LATE não são as três estimativas. O LATE é apenas a razão de Wald local, 2.167, interpretada como o efeito médio local para os compliers no cutoff.

## Testes placebo

Testando descontinuidade em covariáveis predeterminadas: covariáveis que não devem ser afetadas pelo tratamento não devem apresentar salto no ponto de corte. Esse é um diagnóstico de plausibilidade, não uma prova da suposição de continuidade.

Testando descontinuidades em outros pontos: verificar a existência de descontinuidades em pontos arbitrários ao longo da variável de ordenação. Se encontramos vários "efeitos" longe do cutoff real, isso sugere que a especificação está capturando forma funcional ou outras descontinuidades, não necessariamente o tratamento.

Uso de VDs placebos: se uma variável dependente que não deveria ser afetada pelo tratamento apresentar descontinuidade significativa, isso levanta dúvidas sobre a validade do desenho RD.

Avaliação de sensibilidade às covariáveis: as estimativas de RD não devem ser altamente sensíveis à inclusão ou exclusão de covariáveis.

Ao mesmo tempo, muitos testes diagnósticos criam outro problema: se testarmos muitas covariáveis, muitos cutoffs placebo e muitos outcomes placebo, alguns resultados "significativos" podem aparecer por acaso. Por isso, diagnósticos devem ser planejados e interpretados como um conjunto de evidências. A literatura recente propõe testes mais unificados para RD justamente para reduzir esse problema de múltiplas comparações.

## Extensões e casos difíceis

Algumas *running variables* são discretas: idade em anos, notas arredondadas, população municipal em números inteiros ou margens eleitorais reportadas com poucas casas decimais. O fato de $X_i$ ser discreta não invalida automaticamente um RD. O problema aparece quando há poucos valores de $X_i$ perto do cutoff. Nesse caso, a estimativa depende mais de extrapolação local, e testes de densidade ou balanceamento podem ter poder limitado.

Muitos RDDs em ciência política e políticas públicas usam dados agrupados: municípios, escolas, distritos, eleições repetidas ou unidades administrativas. Clusterizar erros padrão pode ser necessário, mas não resolve o problema por si só. A inferência depende de quantos clusters existem perto do cutoff e de como $X_i$ varia dentro e entre clusters.

Outras extensões exigem definir o estimando com mais cuidado. Regras com cutoff podem operar repetidamente no tempo, como eleições em vários ciclos, políticas que renovam elegibilidade ou limites fiscais anuais. Hsu e Shen (2024) mostram que, nesses casos, precisamos distinguir efeito do primeiro tratamento, tratamento acumulado e exposição dinâmica.

Também há desenhos em que a *running variable* observada mede com erro uma variável latente, como testes escolares, índices de risco ou scores administrativos. Eckles et al. (2025) mostram que, quando conhecemos a estrutura desse ruído, a incerteza no score pode induzir uma forma de aleatorização útil para identificação. Essa estratégia não substitui o RD padrão; ela define outro caminho quando o erro de medida é parte conhecida do desenho.

Por fim, desenhos com múltiplos cutoffs ou fronteiras geográficas não se reduzem facilmente a um único ponto $c$. Em RD espacial, unidades de lados diferentes de uma fronteira administrativa podem receber políticas diferentes. Cattaneo, Titiunik e Yu (2025) mostram que tratar a distância até a fronteira como uma *running variable* unidimensional pode introduzir viés quando a fronteira tem quinas ou irregularidades. Quando possível, a análise deve usar a informação bidimensional da localização.

## PCRD

Marshall (2024) na AJPS introduz a nomenclatura do desenho de pesquisa Politician characteristic regression discontinuity (PCRD). Basicamente, o argumento é que RDD não permite identificar efeito de características de políticos (como gênero, profissão, raça, ideologia, alinhamento com governo federal etc.)

"In contrast, the treatment in PCRD designs — which instead seek to estimate the LATE of an elected politician characteristic — is defined by possessing (or not) predetermined characteristic X, conditional on narrowly winning an election. (...) restricting attention to close elections entails conditioning on candidate vote shares that may be affected by X. (...) [It] generally introduce bias — even when X is independent of other predetermined variables and the weak continuity assumption underpinning standard RD designs holds." (p. 495)

Basicamente, Marshall está dizendo que nesses casos, close election é um collider, e isso abre as portas para vieses de variáveis que causem $y$ e se a eleição é apertada.

A recomendação não é resolver o problema adicionando controles. O pesquisador precisa escolher entre três caminhos. O primeiro é defender uma hipótese forte: ou a característica $X_i$ não afeta a votação, ou os diferenciais compensatórios que tornam a eleição apertada não afetam o resultado de interesse. Por exemplo, ao estudar o efeito de eleger mulheres, seria preciso argumentar que gênero não afeta votos naquela eleição, ou que diferenças de experiência, competência, partido e recursos entre candidatas e candidatos em disputas apertadas não afetam o outcome. Em muitas aplicações, essa defesa é difícil.

O segundo caminho é tratar os testes de covariáveis de candidatos como diagnóstico de viés, não como validação automática. Se mulheres eleitas em disputas apertadas são mais experientes que homens eleitos em disputas apertadas, essa diferença ajuda a entender o que o estimador está comparando. Se não encontramos diferença observável, ainda não provamos que não há diferenciais compensatórios não observados.

O terceiro caminho é mudar o estimando. Em vez de afirmar que o PCRD identifica o efeito "all else equal" de gênero, raça, profissão ou partido, o pesquisador pode interpretar a estimativa como um efeito composto: eleger uma pessoa com a característica $X_i$ junto com o conjunto de características que permite que esse tipo de candidato vença uma eleição apertada. Essa interpretação pode ser útil, mas é menos limpa teoricamente e tem validade externa limitada.

\begin{figure}
\includegraphics[alt={DAG em que gênero e competência afetam a probabilidade de eleição apertada e também o resultado, tornando eleição apertada um collider.}]{07-RDD_files/figure-latex/DAG-1} \caption{PCRD: condicionar em eleições apertadas pode abrir um caminho de collider.}(\#fig:DAG)
\end{figure}

## Checklist para um paper de RDD

1. Descrever a regra institucional que gera o cutoff e quem poderia manipulá-la.

2. Definir se o desenho é sharp ou fuzzy. Em fuzzy RD, reportar a primeira etapa e interpretar o estimando como LATE local.

3. Apresentar gráfico RD com bins apropriados, curva local e cutoff visível.

4. Reportar estimativas com regressão local linear, kernel e bandwidth explícitos.

5. Usar inferência robusta com correção de viés, com bandwidth MSE/CER-optimal ou justificativa clara para bandwidth escolhido.

6. Separar a interpretação de continuidade da interpretação de aleatorização local. Bandwidth automático de `rdrobust` não define uma janela "como se aleatória".

7. Testar balanceamento de covariáveis de pré-tratamento.

8. Avaliar densidade da *running variable* no cutoff com `rddensity`/McCrary, lembrando que falhar em rejeitar manipulação não prova ausência de manipulação.

9. Fazer testes placebo em cutoffs e outcomes substantivamente justificados, evitando uma pescaria de especificações.

10. Examinar sensibilidade a bandwidths, polinômios locais, covariáveis e *donut RD* (exclusão de uma faixa muito próxima ao cutoff) quando houver risco de manipulação muito perto do cutoff.

11. Discutir casos difíceis: running variable discreta, clusters, múltiplos cutoffs, RD espacial, eleições apertadas ou tratamento dinâmico.

12. Estabelecer a validade da estratégia antes de discutir apenas resultados.

13. Em PCRD, não interpretar automaticamente a estimativa como efeito isolado de uma característica do político. Explicitar se a hipótese forte de Marshall é plausível; caso contrário, apresentar a estimativa como efeito composto ou como exercício de bounds/sensibilidade.

## Referências

Canay, I. A., & Kamat, V. (2018). Approximate permutation tests and induced order statistics in the regression discontinuity design. The Review of Economic Studies, 85(3), 1577-1608.

Cattaneo, M. D., Idrobo, N., & Titiunik, R. (2024). *A practical introduction to regression discontinuity designs: Extensions*. Cambridge University Press. https://doi.org/10.1017/9781009441896

Cattaneo, M. D., Idrobo, N., & Titiunik, R. (2020). *A practical introduction to regression discontinuity designs: Foundations*. Cambridge University Press.

Cattaneo, M. D., Jansson, M., & Ma, X. (2020). Simple local polynomial density estimators. *Journal of the American Statistical Association*, 115(531), 1449-1455.

Cattaneo, M. D., & Titiunik, R. (2022). Regression discontinuity designs. *Annual Review of Economics*, 14, 821-851. https://doi.org/10.1146/annurev-economics-051520-021409

Cattaneo, M. D., Titiunik, R., & Yu, R. R. (2025a). Estimation and inference in boundary discontinuity designs. arXiv:2505.05670. https://arxiv.org/abs/2505.05670

Cattaneo, M. D., Titiunik, R., & Yu, R. R. (2025b). rd2d: Causal inference in boundary discontinuity designs. arXiv:2505.07989. https://arxiv.org/abs/2505.07989

Caughey, D., & Sekhon, J. S. (2011). Elections and the regression discontinuity design: Lessons from close U.S. House races, 1942-2008. *Political Analysis*, 19(4), 385-408.

De Magalhães, L., Hangartner, D., Hirvonen, S., Meriläinen, J., Ruiz, N. A., & Tukiainen, J. (2025). When Can We Trust Regression Discontinuity Design Estimates from Close Elections? Evidence from Experimental Benchmarks. Political Analysis, 1-8.

De la Cuesta, B., & Imai, K. (2016). Misunderstandings about the regression discontinuity design in the study of close elections. *Annual Review of Political Science*, 19(1), 375-396.

Eckles, D., Ignatiadis, N., Wager, S., & Wu, H. (2025). Noise-induced randomization in regression discontinuity designs. *Biometrika*, 112(2), asaf003. https://doi.org/10.1093/biomet/asaf003

Erikson, R. S., & Rader, K. (2017). Much ado about nothing: RDD and the incumbency advantage. *Political Analysis*, 25(2), 269-275.

Fitzgerald, J. (2025). Manipulation tests in regression discontinuity design: The need for equivalence testing. *MetaArXiv*. https://doi.org/10.31219/osf.io/2dgrp_v1

Gelman, A., & Imbens, G. (2019). Why high-order polynomials should not be used in regression discontinuity designs. Journal of Business & Economic Statistics, 37(3), 447-456.

Hartman, E. (2021). Equivalence testing for regression discontinuity designs. *Political Analysis*, 29(4), 505-521. https://doi.org/10.1017/pan.2020.43

Hsu, Y.-C., & Shen, S. (2024). Dynamic regression discontinuity under treatment effect heterogeneity. *Quantitative Economics*, 15(4), 1035-1064. https://doi.org/10.3982/QE2150

Kreiß, A., & Rothe, C. (2023). Inference in regression discontinuity designs with high-dimensional covariates. *The Econometrics Journal*, 26(2), 105-123. https://doi.org/10.1093/ectj/utac029

Lee, D. S. (2008). Randomized experiments from non-random selection in U.S. House elections. *Journal of Econometrics*, 142(2), 675-697.

Marshall, J. (2024). Can close election regression discontinuity designs identify effects of winning politician characteristics? *American Journal of Political Science*, 68(2), 494-510.

McCrary, J. (2008). Manipulation of the running variable in the regression discontinuity design: A density test. *Journal of Econometrics*, 142(2), 698-714.

Noack, C., Olma, T., & Rothe, C. (2025). Flexible covariate adjustments in regression discontinuity designs. Working paper, revised April 2025. https://arxiv.org/abs/2107.07942

Noack, C., & Rothe, C. (2024). Bias-aware inference in fuzzy regression discontinuity designs. *Econometrica*, 92(3), 687-711. https://doi.org/10.3982/ECTA19466

- Tutorial: https://congressdata.joshuamccrain.com/regression_discontinuity.html
