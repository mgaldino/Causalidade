# Replication Guide

## Obra Transparente DiD Analysis

**Paper**: "Can Civil Society Organizations Improve Public Service Delivery? Evidence from Construction Monitoring in Brazil" (Galdino, Mondo, Sakai, and Paiva)

**Last updated**: 2026-02-10

---

## Requirements

### Software

- **R** >= 4.0.0
- **RStudio** (recommended, but optional)

### R Packages

Install the required packages by running:

```r
install.packages(c(
  "tidyverse",
  "here",
  "fixest",
  "readxl",
  "janitor",
  "data.table",
  "modelsummary",
  "kableExtra",
  "stringi",
  "scales",
  "bookdown",
  "estimatr",
  "gt",
  "fwildclusterboot",
  "lubridate",
  "knitr",
  "rmarkdown"
))
```

| Package | Tested Version | Purpose |
|---------|----------------|---------|
| tidyverse | 2.0.0 | Data manipulation and visualization |
| here | 1.0.1 | Relative file paths |
| fixest | 0.12.1 | DiD and event study estimation |
| readxl | 1.4.3 | Reading Excel files |
| janitor | 2.2.0 | Data cleaning |
| data.table | 1.17.0 | Fast CSV reading |
| modelsummary | 2.3.0 | Regression tables |
| kableExtra | 1.4.0 | Table formatting |
| stringi | 1.8.4 | String normalization |
| scales | 1.3.0 | Axis formatting |
| bookdown | 0.40 | Manuscript compilation |
| estimatr | 1.0.4 | Robust standard errors |
| gt | 0.11.1 | Table formatting |
| fwildclusterboot | 0.14.3 | Wild cluster bootstrap inference |
| lubridate | 1.9.3 | Date manipulation |
| knitr | 1.48 | Dynamic report generation |
| rmarkdown | 2.28 | Manuscript rendering |

---

## Project Structure

```
did-obra-transparente/
├── REPLICATION.md              # This file
├── README.md                   # Project overview
├── code/
│   ├── 00a_download_data.R     # Download SIMEC snapshots from GitHub
│   ├── 01_did_analysis_6periods.R  # Main DiD analysis (6 periods)
│   ├── 05_paper_tables_figures.R   # Tables and figures for the paper
│   ├── 06_robustness_covariates.R  # Robustness: political covariates
│   ├── 07_robustness_wild_bootstrap.R # Robustness: wild cluster bootstrap
│   └── 99_run_all.R               # Master script (run this)
├── data/
│   ├── raw/                    # Raw data files (included in repository)
│   │   ├── escolas_brasil_2015.xlsx
│   │   ├── simec_2017_05.Rdata
│   │   ├── simec_2019_08.Rdata
│   │   ├── simec_2023_10.csv
│   │   ├── treated_works.xlsx
│   │   ├── municipal_covariates.rds
│   │   ├── elections_mayor_results.csv
│   │   ├── party_ideology.csv
│   │   └── simec_snapshots/    # Downloaded automatically
│   ├── processed/              # Generated panel data (created by scripts)
│   └── metadata/               # Data documentation
├── output/
│   ├── *.rds                   # Intermediate outputs for paper compilation
│   ├── tables/                 # Generated tables (.csv, .tex)
│   ├── figures/                # Generated figures (.png)
│   └── paper/                  # Manuscript (.Rmd, .pdf, .docx)
└── notes/                      # Additional documentation
```

---

## Data Sources

All data required for replication are included in `data/raw/`. The original sources are:

| File | Period | Source | Description |
|------|--------|--------|-------------|
| `escolas_brasil_2015.xlsx` | Aug 2015 | SIMEC/FNDE via TB | ProInfancia construction status |
| `simec_2017_05.Rdata` | May 2017 | SIMEC/FNDE via TB | Construction status at project start |
| `simec_snapshots/obras_08032018.csv` | Mar 2018 | GitHub (TB) | Auto-downloaded |
| `simec_snapshots/obras_upload28092018.csv` | Sep 2018 | GitHub (TB) | Auto-downloaded |
| `simec_2019_08.Rdata` | Aug 2019 | SIMEC/FNDE via TB | Construction status at project end |
| `simec_2023_10.csv` | Oct 2023 | SIMEC/FNDE via TB | Long-term follow-up |
| `treated_works.xlsx` | - | Transparencia Brasil | List of treated construction works |
| `municipal_covariates.rds` | - | IPEA/IBGE | Municipal socioeconomic characteristics |
| `elections_mayor_results.csv` | - | TSE | Mayoral election results (robustness) |
| `party_ideology.csv` | - | Literature | Party ideology scores (robustness) |

**GitHub repositories** (for automatic snapshot download):
- https://github.com/Transparencia-Brasil/avaliacao_impacto_092019
- https://github.com/Transparencia-Brasil/campanha_TDP_2018

---

## Step-by-Step Replication

### 1. Set Up the Environment

```bash
# Clone the repository
git clone [REPOSITORY_URL]
cd did-obra-transparente
```

```r
# Verify working directory in R
library(here)
here()  # Should show the project path
```

### 2. Run the Complete Pipeline

**Option A: From R/RStudio**
```r
source(here::here("code", "99_run_all.R"))
```

**Option B: From the command line**
```bash
Rscript code/99_run_all.R
```

The master script will:
1. Verify that all required data files are present
2. Download SIMEC snapshots from GitHub (if missing)
3. Build a balanced panel of 6 periods (2015-2023)
4. Estimate the static DiD and event study models
5. Generate all tables and figures
6. Run robustness checks (political covariates, wild cluster bootstrap)

### 3. Compile the Manuscript (Optional)

```r
rmarkdown::render(here::here("output", "paper", "paper_v2026.Rmd"))
```

### 4. Verify Outputs

After execution, the following files are generated:

**Processed data** (`data/processed/`):
- `did_panel_6periods.rds` - Balanced panel (~21,000 obs)
- `did_results_6periods.rds` - Complete results

**Paper inputs** (`output/`):
- `tab_01.rds` - Table 1: Summary Statistics
- `table2.rds` - Table 2: Completion Rates
- `did_static.rds` - Table 3: Static DiD Model
- `did_event_study.rds` - Table 4: Event Study Model
- `obra_transparente.RDS` - Data for paper figures
- `wild_bootstrap.rds` - Wild cluster bootstrap results

**Tables** (`output/tables/`):
- `table1_summary_stats.csv`
- `table2_completion_rates.csv`
- `table3_did_static.csv`
- `table4_event_study.csv`
- `robustness_covariates.csv` / `.tex` / `.html`
- `wild_bootstrap_results.csv`

**Figures** (`output/figures/`):
- `fig1_completion_trends.png`
- `fig2_event_study.png`
- `fig3_event_study_alt.png`

---

## Expected Results

### Data

- **3,494 construction projects** in 1,845 municipalities (balanced panel)
- **188 projects** in 21 treated municipalities
- 6 periods: Aug 2015, May 2017, Mar 2018, Sep 2018, Aug 2019, Oct 2023
- Treatment begins in Sep 2018 (period 3)

### Static DiD (TWFE)

| ATT | SE | 95% CI | p-value |
|:---:|:--:|:------:|:-------:|
| **+0.083** | 0.034 | [0.017, 0.150] | **0.013** |

### Event Study

| Relative Period | Coefficient | SE | p-value |
|:---------------:|:-----------:|:--:|:-------:|
| t = -3 | -0.013 | 0.062 | 0.837 |
| t = -2 | +0.002 | 0.014 | 0.871 |
| t = -1 | 0.000 | (ref) | - |
| t = 0 | +0.018 | 0.023 | 0.444 |
| t = +1 | +0.042 | 0.034 | 0.219 |
| t = +2 | **+0.180** | 0.088 | **0.039** |

### Robustness

- Wild cluster bootstrap (static ATT): p = 0.017
- Wild cluster bootstrap (t=+2): p = 0.050
- Political covariates: ATT stable across all specifications (0.082-0.084)

---

## Runtime

The complete pipeline runs in approximately **3-5 minutes** (most time is spent on the wild cluster bootstrap in Step 4).

---

## Troubleshooting

### Error: "Wrong working directory"

```r
setwd("/path/to/did-obra-transparente")
```

### Error: "Essential data files missing"

Verify that all files listed in `data/raw/` above are present. These should be included in the repository.

### Error: "Could not download snapshots"

Run the download script manually:
```r
source(here::here("code", "00a_download_data.R"))
```

Or download the CSV files manually from:
- https://github.com/Transparencia-Brasil/avaliacao_impacto_092019/tree/master/Bancos/Planilhas%20SIMEC

### Table 1 shows simplified version

The full Table 1 with municipal covariates requires `data/raw/municipal_covariates.rds`. If this file is missing, a simplified version based on the DiD panel is generated instead.

---

## Contact

- **Corresponding author**: Manoel Galdino (mgaldino@usp.br)
- **Transparencia Brasil**: https://www.transparencia.org.br/

---

## License

This code is released under the MIT License.
