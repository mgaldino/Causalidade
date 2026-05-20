# Obra Transparente DiD Analysis

Replication package for:

> **"Can Civil Society Organizations Improve Public Service Delivery? Evidence from Construction Monitoring in Brazil"**
>
> Manoel Galdino, Bianca Vaz Mondo, Juliana Mari Sakai, and Natalia Paiva

## Overview

This repository contains the data and code to replicate the quantitative analysis of the **Obra Transparente** project's impact on public school and nursery construction completion rates in Brazil.

The Obra Transparente project was implemented by Transparencia Brasil from May 2017 to June 2019, engaging 21 local Civil Society Organizations (CSOs) to monitor federally funded construction works (ProInfancia program) in municipalities across South and Southeast Brazil.

## Key Findings

- **Static ATT**: CSO monitoring increased construction completion by 8.3 percentage points (p = 0.013)
- **Dynamic effects**: The treatment effect grows over time, reaching 18 p.p. by October 2023 (p = 0.039)
- **Parallel trends**: Pre-treatment coefficients are close to zero (joint F-test p > 0.90)
- **Robust inference**: Results hold under wild cluster bootstrap (p = 0.017)

## Quick Start

```r
# Install required packages
install.packages(c("tidyverse", "here", "fixest", "readxl", "janitor",
                    "data.table", "modelsummary", "kableExtra", "stringi"))

# Run complete replication pipeline
source(here::here("code", "99_run_all.R"))
```

See **[REPLICATION.md](REPLICATION.md)** for detailed instructions, expected results, and troubleshooting.

## Structure

```
code/                          # R scripts
  99_run_all.R                 # Master script (run this)
  01_did_analysis_6periods.R   # Main DiD + event study
  05_paper_tables_figures.R    # Tables and figures
  06_robustness_covariates.R   # Political covariates robustness
  07_robustness_wild_bootstrap.R # Wild cluster bootstrap
data/raw/                      # Raw data (included in repository)
output/                        # Generated tables, figures, and paper
```

## Data

All data required for replication are included in `data/raw/`. Sources include SIMEC/FNDE (construction project status), IPEA/IBGE (municipal characteristics), and TSE (election results). See [REPLICATION.md](REPLICATION.md) for full data documentation.

## Requirements

- R >= 4.0
- See [REPLICATION.md](REPLICATION.md) for the full list of R packages

## Citation

```bibtex
@article{galdino2025cso,
  title={Can Civil Society Organizations Improve Public Service Delivery? Evidence from Construction Monitoring in Brazil},
  author={Galdino, Manoel and Mondo, Bianca Vaz and Sakai, Juliana Mari and Paiva, Nat\'{a}lia},
  year={2025}
}
```

## License

MIT License. See `LICENSE` for details.

## Contact

Manoel Galdino - mgaldino@usp.br
