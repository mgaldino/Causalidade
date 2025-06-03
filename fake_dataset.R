# Versão simplificada para download rápido
library(tidyverse)
library(httr)

# Função simples para baixar
download_fakebr_simple <- function(max_files = 500) {
  
  base_url <- "https://raw.githubusercontent.com/roneysco/Fake.br-Corpus/refs/heads/master/full_texts/"
  
  # Função para baixar uma categoria
  download_category <- function(folder, label, category_name) {
    
    cat(glue::glue("Baixando {category_name}...\n"))
    
    results <- 1:max_files %>%
      map_dfr(~ {
        url <- glue::glue("{base_url}{folder}/{.x}.txt")
        
        tryCatch({
          response <- GET(url, timeout(5))
          
          if (status_code(response) == 200) {
            text <- content(response, "text", encoding = "UTF-8") %>% str_trim()
            
            if (nchar(text) > 50) {
              return(tibble(
                text = text,
                label = label,
                category = category_name,
                file_id = .x
              ))
            }
          }
          return(NULL)
        }, error = function(e) NULL)
      }) %>%
      filter(!is.na(text))
    
    cat(glue::glue("✅ {nrow(results)} arquivos de {category_name}\n"))
    return(results)
  }
  
  # Baixa ambas categorias
  fake_data <- download_category("fake", 1, "fake")
  true_data <- download_category("true", 0, "true")
  
  # Combina e embaralha
  final_data <- bind_rows(fake_data, true_data) %>%
    slice_sample(prop = 1) %>%
    mutate(id = row_number())
  
  return(final_data)
}

# Executa
dataset <- download_fakebr_simple(max_files = 3000)
glimpse(dataset)
# Salva
library(here)
write_csv(dataset, here("data", "fakebr_dataset.csv"))

# Resumo
dataset %>%
  count(category) %>%
  mutate(prop = round(n/sum(n), 3))