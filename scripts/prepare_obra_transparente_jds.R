# Prepare the Obra Transparente data and estimates used in 08-DiD.Rmd.
#
# The raw data and replication scripts are stored under
# data/obra_transparente_jds/replication_package/. This wrapper keeps the course
# chapter reproducible without depending on the separate paper repository.

find_course_root <- function(start = getwd()) {
  path <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(path, "Causalidade.Rproj"))) {
      return(path)
    }

    parent <- dirname(path)
    if (identical(parent, path)) {
      stop("Could not find Causalidade.Rproj from ", start, call. = FALSE)
    }
    path <- parent
  }
}

invisible(Sys.setlocale("LC_CTYPE", "pt_BR.UTF-8"))
invisible(Sys.setlocale("LC_COLLATE", "pt_BR.UTF-8"))

run_replication_script <- function(script, package_root) {
  script_path <- file.path(package_root, "code", script)
  if (!file.exists(script_path)) {
    stop("Replication script not found: ", script_path, call. = FALSE)
  }

  script_env <- new.env(parent = globalenv())
  script_env$here <- function(...) {
    file.path(package_root, ...)
  }

  sys.source(script_path, envir = script_env)
  invisible(script_env)
}

course_root <- find_course_root()
package_root <- file.path(course_root, "data", "obra_transparente_jds", "replication_package")

if (!dir.exists(package_root)) {
  stop("Obra Transparente replication package not found: ", package_root, call. = FALSE)
}

dir.create(file.path(package_root, "data", "processed"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(package_root, "output"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(package_root, "output", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(package_root, "output", "figures"), recursive = TRUE, showWarnings = FALSE)

cat("Preparing Obra Transparente JDS objects in:\n")
cat("  ", package_root, "\n\n", sep = "")

run_replication_script("01_did_analysis_6periods.R", package_root)
run_replication_script("05_paper_tables_figures.R", package_root)

if (!file.exists(file.path(package_root, "output", "tables", "wild_bootstrap_results.csv"))) {
  message(
    "Wild bootstrap table not found. The chapter can be rendered without it, ",
    "but the inference table should be copied from the submitted replication ",
    "outputs or regenerated from code/07_robustness_wild_bootstrap.R in the ",
    "paper repository."
  )
}

cat("\nObra Transparente JDS preparation complete.\n")
