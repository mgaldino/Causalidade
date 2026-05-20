#' ==============================================================================
#' Script: 00a_download_data.R
#' Descricao: Download dos dados brutos dos repositorios da Transparencia Brasil
#'
#' Este script clona os repositorios necessarios e copia os arquivos SIMEC
#' para a estrutura de dados do projeto.
#'
#' Fontes:
#'   - avaliacao_impacto_092019: Snapshots SIMEC 2017-2019
#'   - campanha_TDP_2018: Lista dos 21 municipios tratados (ot.Rdata)
#'   - original/data/: Snapshot SIMEC 2023 (local)
#'
#' Outputs:
#'   - data/raw/simec_snapshots/obras30052017.csv (periodo 1)
#'   - data/raw/simec_snapshots/obras_08032018.csv (periodo 2)
#'   - data/raw/simec_snapshots/obras_upload28092018.csv (periodo 3)
#'   - data/raw/simec_snapshots/obras2019_08_16.csv (periodo 4)
#'   - data/raw/simec_snapshots/simec 25-10-23 - simec.csv (periodo 5)
#'   - data/raw/ot.Rdata (lista de municipios tratados)
#'
#' Autor: Manoel Galdino
#' Data: 2026-01-05
#' ==============================================================================

# Setup ----
library(here)

cat("\n=== Download de Dados Brutos ===\n")
cat("Data/Hora:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Definir diretorios ----
temp_dir <- tempdir()
output_dir <- here("data", "raw", "simec_snapshots")

# Criar diretorio de saida se nao existir ----
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("Diretorio criado:", output_dir, "\n")
}

# URLs dos repositorios ----
repos <- list(
  avaliacao = "https://github.com/Transparencia-Brasil/avaliacao_impacto_092019.git",
  campanha = "https://github.com/Transparencia-Brasil/campanha_TDP_2018.git"
)

# Funcao para clonar repositorio ----
clone_repo <- function(url, dest_name) {
  dest_path <- file.path(temp_dir, dest_name)

  if (dir.exists(dest_path)) {
    cat("Repositorio ja existe:", dest_name, "\n")
    return(dest_path)
  }

  cat("Clonando", dest_name, "...\n")
  result <- system2(
    "git",
    args = c("clone", "--depth", "1", url, dest_path),
    stdout = TRUE,
    stderr = TRUE
  )

  if (!dir.exists(dest_path)) {
    stop("Falha ao clonar repositorio: ", dest_name)
  }

  cat("Clonado com sucesso:", dest_name, "\n")
  return(dest_path)
}

# Clonar repositorios ----
cat("\n--- Clonando repositorios ---\n")
path_avaliacao <- clone_repo(repos$avaliacao, "avaliacao_impacto_092019")
path_campanha <- clone_repo(repos$campanha, "campanha_TDP_2018")

# Mapear arquivos SIMEC ----
simec_files <- list(
  periodo_1 = list(
    source = file.path(path_avaliacao, "Bancos", "Planilhas SIMEC", "obras30052017.csv"),
    dest = "obras30052017.csv",
    date = "2017-05-30",
    desc = "Maio 2017 (baseline)"
  ),
  periodo_2 = list(
    source = file.path(path_avaliacao, "Bancos", "Planilhas SIMEC", "obras_08032018.csv"),
    dest = "obras_08032018.csv",
    date = "2018-03-08",
    desc = "Marco 2018"
  ),
  periodo_3 = list(
    source = file.path(path_avaliacao, "Bancos", "Planilhas SIMEC", "obras_upload28092018.csv"),
    dest = "obras_upload28092018.csv",
    date = "2018-09-28",
    desc = "Setembro 2018"
  ),
  periodo_4 = list(
    source = file.path(path_avaliacao, "Bancos", "Planilhas SIMEC", "obras2019_08_16.csv"),
    dest = "obras2019_08_16.csv",
    date = "2019-08-16",
    desc = "Agosto 2019"
  )
)

# Copiar arquivos SIMEC do GitHub ----
cat("\n--- Copiando snapshots SIMEC ---\n")

for (periodo in names(simec_files)) {
  file_info <- simec_files[[periodo]]
  dest_path <- file.path(output_dir, file_info$dest)

  if (file.exists(file_info$source)) {
    file.copy(file_info$source, dest_path, overwrite = TRUE)
    cat(sprintf("[OK] %s: %s (%s)\n", periodo, file_info$desc, file_info$date))
  } else {
    cat(sprintf("[ERRO] %s: Arquivo nao encontrado - %s\n", periodo, file_info$source))
  }
}

# Copiar periodo 5 (2023) da pasta original ----
cat("\n--- Copiando snapshot 2023 (local) ---\n")

simec_2023_source <- here("original", "data", "simec 25-10-23 - simec.csv")
simec_2023_dest <- file.path(output_dir, "simec 25-10-23 - simec.csv")

if (file.exists(simec_2023_source)) {
  file.copy(simec_2023_source, simec_2023_dest, overwrite = TRUE)
  cat("[OK] periodo_5: Outubro 2023 (2023-10-25)\n")
} else {
  cat("[ERRO] periodo_5: Arquivo nao encontrado -", simec_2023_source, "\n")
  cat("       Este arquivo deve estar em original/data/\n")
}

# Copiar ot.Rdata (lista de municipios tratados) ----
cat("\n--- Copiando lista de municipios tratados ---\n")

ot_source <- file.path(path_campanha, "ot.Rdata")
ot_dest <- here("data", "raw", "ot.Rdata")

if (file.exists(ot_source)) {
  file.copy(ot_source, ot_dest, overwrite = TRUE)
  cat("[OK] ot.Rdata copiado para data/raw/\n")
} else {
  cat("[ERRO] ot.Rdata nao encontrado em", ot_source, "\n")
}

# Verificar arquivos copiados ----
cat("\n--- Verificacao Final ---\n")

arquivos_esperados <- c(
  file.path(output_dir, "obras30052017.csv"),
  file.path(output_dir, "obras_08032018.csv"),
  file.path(output_dir, "obras_upload28092018.csv"),
  file.path(output_dir, "obras2019_08_16.csv"),
  file.path(output_dir, "simec 25-10-23 - simec.csv"),
  here("data", "raw", "ot.Rdata")
)

todos_presentes <- all(file.exists(arquivos_esperados))

if (todos_presentes) {
  cat("\nTodos os arquivos foram copiados com sucesso!\n")
  cat("\nArquivos em data/raw/simec_snapshots/:\n")
  for (f in list.files(output_dir)) {
    size_mb <- round(file.size(file.path(output_dir, f)) / 1024 / 1024, 1)
    cat(sprintf("  - %s (%.1f MB)\n", f, size_mb))
  }
} else {
  cat("\nATENCAO: Alguns arquivos estao faltando!\n")
  for (f in arquivos_esperados) {
    status <- ifelse(file.exists(f), "[OK]", "[FALTA]")
    cat(sprintf("  %s %s\n", status, basename(f)))
  }
}

# Log final ----
cat("\n=== Download Concluido ===\n")
cat("Data/Hora:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Session info ----
sessionInfo()
