# Parecer — Parecerista 2 (Método & Inferência)

## Recomendação: R&R major

## Resumo do paper
O artigo estima os retornos eleitorais de curto e longo prazo da resposta governamental às enchentes do Elba de 2002 na Alemanha. Usando DiD em distritos eleitorais, os autores encontram ganho de cerca de 7 p.p. para o SPD em 2002, persistência parcial em 2005 e desaparecimento do efeito em 2009.

## Avaliação geral
O paper é metodologicamente forte: o choque é plausivelmente exógeno em timing, o desenho DiD é transparente, há placebo pré-tratamento, ajustes por redistritamento, testes com covariáveis, evidência de saliência via surveys e uma análise espacial útil. A principal fragilidade é que o estimando causal é mais amplo do que a linguagem do paper sugere: o desenho identifica melhor o efeito conjunto de estar em uma área atingida pela enchente e exposta à resposta política, não o efeito isolado do gasto de alívio ou da “gratidão” dos eleitores. Para um top journal, eu pediria revisão substancial de identificação, inferência espacial e escopo causal, mas vejo o núcleo empírico como promissor.

## Comentários maiores

1. **O tratamento combina desastre, saliência, dano e resposta governamental.**  
   A variável `Flooded` indica exposição a eventos como alerta, evacuação, rompimento ou estabilização de diques. Isso não mede diretamente a dose de ajuda pública recebida. O efeito estimado pode refletir gratidão pelo auxílio, mas também medo, solidariedade local, maior saliência de competência executiva, dano econômico, presença militar, cobertura de mídia ou rejeição ao challenger. Sugestão: reescrever o claim principal como efeito da “exposição à enchente e à resposta governamental” e tratar “voter gratitude” como interpretação mecanística, não como estimando diretamente identificado. Idealmente, incorporar dados administrativos de pagamentos por distrito/município ou intensidade de dano.

2. **A suposição de parallel trends é plausível, mas ainda subtestada.**  
   O placebo 1994-1998 é convincente e ajuda muito. Ainda assim, há apenas uma diferença pré-tratamento observada antes do choque principal. Como quase todos os tratados estão no Leste, a comparação nacional pode absorver tendências políticas regionais divergentes apenas parcialmente. Sugestão: incluir event-study com todos os pré-períodos disponíveis, mesmo que limitados; restringir controles a distritos do Leste, aos mesmos Länder, a bandas de distância do Elba, ou a distritos fluviais comparáveis; e reportar sensibilidade a tendências lineares regionais ou `Land × year FE`, quando houver variação suficiente.

3. **A inferência parece otimista diante da concentração espacial do tratamento.**  
   Os erros-padrão são clusterizados por distrito, o que é melhor que robustos simples, mas o tratamento é espacialmente concentrado em 29 distritos e choques eleitorais podem ser correlacionados dentro de regiões, bacias hidrográficas e Länder. Isso afeta especialmente o efeito de 2005, menor e substantivamente central. Sugestão: reportar inferência por randomization/permutation, placebo-in-space com pseudo-enchentes em rios/regiões comparáveis, wild cluster bootstrap ou métodos com dependência espacial tipo Conley. Também seria útil mostrar se os resultados sobrevivem a exclusão leave-one-state-out ou leave-one-cluster-out.

4. **O grupo de controle pode estar contaminado por spillovers.**  
   Os autores reconhecem que distritos próximos, mas não codificados como tratados, podem ter sido indiretamente afetados. Isso pode atenuar os efeitos, mas também muda o estimando: o contraste deixa de ser “atingidos vs não atingidos” e passa a depender da distância, cobertura de mídia e redes regionais. Sugestão: estimar efeitos por anéis de distância, excluir distritos próximos não tratados, ou modelar tratamento contínuo por distância/intensidade. Isso também ajudaria a separar exposição direta de saliência regional.

5. **O cálculo de “preço por voto” é interessante, mas causalmente frágil.**  
   Dividir o fundo total de €7,1 bilhões pelos votos adicionais assume implicitamente que todo o gasto relevante foi alocado aos distritos tratados e que a relação gasto-voto é causal e homogênea. Sem gasto observado por distrito, esse cálculo deve ser apresentado como back-of-the-envelope, não como estimativa causal de retorno marginal do spending. Sugestão: suavizar a linguagem ou obter dados de desembolso regional.

6. **O mecanismo de persuasão é sugestivo, não decisivo.**  
   As análises com Forsa são valiosas: mostram saliência da enchente, estabilidade de turnout e mudanças em intenção de voto entre antigos eleitores do SPD e CDU/CSU. Mas são dados agregados/repeated cross-section, não painel individual. Logo, não demonstram que os beneficiários diretos mudaram de voto por gratidão. Sugestão: explicitar que a evidência de mecanismo é consistente com persuasão, não uma mediação causal; reportar tamanhos amostrais efetivos para subgrupos em áreas atingidas; e evitar inferência forte sobre comportamento individual.

7. **A interpretação de longo prazo precisa de maior cautela.**  
   O efeito de 2005, cerca de 2 p.p., é substantivamente importante, mas mais vulnerável a choques políticos intervenientes: reformas Hartz, desgaste do SPD, mudanças econômicas locais e realinhamentos no Leste. A comparação 1998-2005 exige uma versão mais forte de parallel trends do que 1998-2002. Sugestão: apresentar 2005 como evidência de persistência residual sob pressupostos mais fortes, com análise de sensibilidade e controles/estratificações regionais adicionais.

8. **Múltiplas comparações e exploração de resultados secundários deveriam ser organizadas.**  
   O paper examina vários anos, partidos, turnout, intenção de voto, distância, issue salience e votos distritais. Isso é adequado, mas a significância do efeito de longo prazo deve ser lida nesse contexto. Sugestão: declarar outcomes primários e secundários; tratar mecanismo e heterogeneidade como exploratórios; e reportar intervalos de confiança com interpretação substantiva, não apenas significância.

## Comentários menores

1. A Tabela 1 é rica, mas muito densa. Um coefficient plot com modelos principais e intervalos facilitaria a leitura.
2. Faltam estatísticas descritivas claras comparando tratados e controles no pré-tratamento: voto SPD, desemprego, densidade, composição setorial, população, participação e região.
3. A definição de `Flooded` é ampla. Seria útil mostrar resultados separados para “fortemente afetados” versus “afetados leves” no texto principal, não apenas mencionar o apêndice.
4. A seção sobre redistritamento é cuidadosa. Ainda assim, eu gostaria de ver uma checagem de sensibilidade usando apenas unidades com fronteiras estáveis ou menor reponderação espacial.
5. A conclusão deveria distinguir com mais precisão “retornos à resposta a desastre” de “retornos a gasto público”. O primeiro é melhor identificado que o segundo.
6. As figuras de polling usam envelopes de 90%; para consistência com o restante do artigo, eu reportaria também 95% ou explicaria a escolha.

## Referências sugeridas

- Bertrand, Duflo & Mullainathan, “How Much Should We Trust Differences-in-Differences Estimates?”
- Abadie, “Semiparametric Difference-in-Differences Estimators.”
- Athey & Imbens, “Identification and Inference in Nonlinear Difference-in-Differences Models.”
- Conley, “GMM Estimation with Cross Sectional Dependence.”
- Cameron, Gelbach & Miller, “Bootstrap-Based Improvements for Inference with Clustered Errors.”
- Abadie, Diamond & Hainmueller, “Synthetic Control Methods for Comparative Case Studies.”
- Rosenbaum, *Design of Observational Studies*.
- Imbens & Wooldridge, “Recent Developments in the Econometrics of Program Evaluation.”
- Roth, “Pretest with Caution: Event-Study Estimates after Testing for Parallel Trends.”
- Rambachan & Roth, “A More Credible Approach to Parallel Trends.”
