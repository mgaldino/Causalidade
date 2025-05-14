
library(here)
library(tidyverse)
library(data.table)
library(readxl)
library(janitor)
library(readr)
library(stringr)
library(modelsummary)
library(fixest)
load(file=here("data", "obras_fim_seg_fase.Rdata" ))
load(file=here("data", "obras_inicio_projeto.Rdata" ))
load(here("data", "sit_obras_final.Rdata"))
list.files(here("data"))
simec23 <- fread(here("data", "simec 25-10-23 - simec.csv"), encoding = "UTF-8") %>%
  clean_names()
simec_14 <- read_excel(here("data", "Escolas Brasil.xlsx")) %>%
  clean_names()
sit_obras_final <- read_csv(here("data", "sit_obras_final.csv"), 
                            locale = locale(encoding = "ISO-8859-1"))

glimpse(simec_14)
glimpse(obras_fim_seg_fase)
glimpse(obras_inicio_projeto)
glimpse(simec23)
glimpse(sit_obras_final)

max(simec23$data_da_ultima_vistoria_do_estado_ou_municipio, na.rm=T)

sit_obras_final %>%
  group_by(projeto_obra_transparente) %>%
  summarise(n())

obras4 <- simec23 %>%
  select(id, municipio, uf, situacao) %>%
  mutate(data_coleta = as.Date("2023-11-01"),
         concluida = ifelse(grepl("Concluída", situacao), 1, 0),
         cancelada = ifelse(grepl("Cancelada", situacao), 1, 0)) 

obras3 <- sit_obras_final %>%
  select(id, municipio, uf, concluida, cancelada) %>%
  mutate(data_coleta = as.Date("2018-11-01")) 

ot <- sit_obras_final %>%
  select(id, municipio, uf, projeto_obra_transparente) %>%
  rename(treat = projeto_obra_transparente)

obras0 <- simec_14 %>%
  mutate(id = str_extract(nome_da_obra, "\\([0-9]+\\)"),
         id = gsub("\\(", "", id),
         id = gsub("\\)", "", id)) %>%
  select(id, municipio, uf, situacao_da_obra) %>%
  mutate(data_coleta = as.Date("2014-12-01"),
         concluida = ifelse(grepl("Concluída", situacao_da_obra), 1, 0),
         cancelada = ifelse(grepl("Cancelada", situacao_da_obra), 1, 0)) %>%
  select(-situacao_da_obra) %>%
  mutate(id = as.numeric(id))


obras1 <- obras_inicio_projeto %>%
  select(id, municipio, uf, situacao) %>%
  mutate(data_coleta = as.Date("2017-05-01"),
         concluida = ifelse(grepl("Concluída", situacao), 1, 0),
         cancelada = ifelse(grepl("Cancelada", situacao), 1, 0)) %>%
  select(-situacao)

obras2 <- obras_fim_seg_fase %>%
  select(id, municipio, uf, situacao) %>%
  mutate(data_coleta = as.Date("2018-12-01"),
         concluida = ifelse(grepl("Concluída", situacao), 1, 0),
         cancelada = ifelse(grepl("Cancelada", situacao), 1, 0)) %>%
  select(-situacao) %>%
  mutate(id = as.numeric(id))

obras0_final <- obras0 %>%
  mutate(excluir = concluida == 1)

obras1_final <- obras1 %>%
  mutate(excluir = concluida == 1)

obras2_final <- obras2 %>%
  mutate(excluir = FALSE)

obras3_final <- obras3 %>%
  mutate(excluir = FALSE)

obras4_final <- obras4 %>%
  mutate(excluir = FALSE)

obras <- bind_rows(obras0_final, obras1_final, obras2_final, obras3_final, obras4_final) %>%
  inner_join(select(ot, id, treat), by = join_by(id)) %>%
  inner_join(select(obras_inicio_projeto, id, percentual_de_execucao), by = join_by(id))


# install.packages("did")
library(did)
library(HonestDiD)

obras_final <- obras %>%
  mutate(periodo = ifelse(data_coleta < "2015-01-01", 1,
                          ifelse(data_coleta < "2018-01-01", 2,
                                 ifelse(data_coleta < "2018-12-01", 3,
                                        ifelse(data_coleta < "2019-12-01", 4, 5))))) %>%
  rename(percentual_de_execucao_maio_17 = percentual_de_execucao) %>%
  group_by(municipio) %>%
  mutate(time_treated = max(treat),
         time_treated1 = time_treated * 4,
         bin_time_treated = ifelse(time_treated1 > 0, 1,0),
         exclusao = max(excluir)) %>%
  ungroup() %>%
  filter(uf %in% c("SP", "PR", "SC", "RS", "MG")) %>%
  # filter(exclusao == 0) %>%
  mutate(pre = ifelse(periodo < 3, 1, 0),
         pos = ifelse(periodo > 2, 1, 0),
         post_treat = pos*bin_time_treated) %>%
  rename(group_treated = time_treated) %>%
  select(id, municipio, concluida, group_treated, periodo, time_treated1)

saveRDS(obras_final, file =  here("Dados", "obra_transparente.RDS"))
  

glimpse(obras_final)

obras_final2x2 <- obras_final %>%
  filter(periodo %in% c(2,5))

did <- lm(concluida ~ as.factor(time_treated) + factor(pos) +   as.factor(time_treated):factor(pos) ,data = obras_final2x2)
summary(did)

#Run the TWFE spec
twfe_results <- fixest::feols(concluida ~ post_treat| municipio + periodo,
                              cluster = "municipio",
                              data = obras_final)


msummary(twfe_results, stars = c('*' = .1, '**' = .05, '***' = .01))

twfe_dynamic <- fixest::feols(concluida ~ i(time, treatedgroup)| municipio + periodo,
                              cluster = "municipio",
                              data = obras_final)


msummary(twfe_dynamic, stars = c('*' = .1, '**' = .05, '***' = .01))
# aumento de 0.9 pontos percentuais na conclusão de obras.
obras_final %>%
  summarise(mean(concluida))

obras_final %>%
  group_by(bin_time_tretead) %>%
  summarise(mean(concluida))

library(ggplot2)

obras_final %>%
  group_by(periodo, group_treated) %>%
  summarise(taxa_conclusao = mean(concluida)) %>%
  ggplot(aes(y=taxa_conclusao, x= periodo, group =group_treated)) + geom_point() +
  geom_line()

est_did = feols(concluida ~ i(periodo, group_treated, ref=3) | id + periodo, obras_final)
summary(est_did)


# Passar de 51% para 52%. Efeito pequeno.

obras_evt <- obras_final %>%
  # supondo que já exista uma coluna time_treated (inteiro com 4, ou NA para controles)
  mutate(event_time = periodo - time_treated1) 

est_did = feols(concluida ~ i(event_time, group_treated, -1) | id + periodo, obras_evt)
summary(est_did)
            
example_attgt <- att_gt(yname = "concluida",
                        tname = "periodo",
                        idname = "id",
                        gname = "time_treated",
                        data = obras_final
)


ggdid(example_attgt)

summary(example_attgt)

set.seed(1814)

sp_list <- reset.sim(time.periods = 4, n = 40, ipw = TRUE, reg = TRUE)
time.periods <- 4

sp_list$te.e <- 1:time.periods

exemplo <- build_sim_dataset(sp_list, panel = TRUE)
View(exemplo)
# summarize the results


# We interact the variable 'period' with the variable 'treat'
data(base_did)
head(base_did)

base_did %>%
  group_by(period, treat) %>%
  summarise(taxa_conclusao = mean(y)) %>%
  ggplot(aes(y=taxa_conclusao, x= period, group =treat)) + geom_point() +
  geom_line()

est_did = feols(y ~ x1 + i(period, treat, 5) | id + period, base_did)

est_did = feols(concluida ~ i(periodo, time_treated, 5) | id + periodo, obras_final)
msummary(est_did, stars = c('*' = .1, '**' = .05, '***' = .01))
iplot(est_did)

# Using i() for factors
est_bis = feols(concluida ~ i(periodo, keep = 3:5) + i(periodo, time_treated, 5) | id, obras_final)
coefplot(est_bis, keep = "trea")
