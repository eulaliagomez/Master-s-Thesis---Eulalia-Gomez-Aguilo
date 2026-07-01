# MSc Thesis: Sensitivity Analysis of Bayesian Parameter Estimation
# in Age-Structured SEIR COVID-19 Models — Netherlands 2020–2021

**Author:** Eulàlia Gómez Aguiló
**Programme:** MSc Statistics and Data Science - Track: EMOS (European Master in Official Statistics)
**Institution:** Leiden University / Statistics Netherlands (CBS)
**Supervisor:** Prof. Dr. Frank P. Pijpers

---

## Overview

This repository contains all code, data, and outputs used in the analysis for my Master's thesis. The goal is to make the analysis fully reproducible from raw RIVM surveillance data through to the final sensitivity analysis figures.

The thesis fits a stochastic age-structured SEIR model to weekly COVID-19 case counts (young <60 / old ≥60) for four epidemic waves in the Netherlands using particle Markov chain Monte Carlo (PMCMC) within the POMP framework. A sensitivity analysis then applies three types of controlled data degradation — proportional under-reporting, stochastic under-reporting, and temporal aggregation — and measures the resulting change in the Bayesian posterior for each estimated parameter.

**Note:** Only final versions of scripts are included. Intermediate drafts and exploratory scripts have been excluded.

---

## Repository structure
```text
repo/
├── Analysis/
│ ├── Sensitivity_analysis/
│ └── Wave-by-wave_baseline_model/
│ ├── full_POMP_pipeline_V1/
│ └── full_POMP_pipeline_V2_k_fixed_age_groups/
├── Data/
│ ├── Degraded_data/
│ ├── Original_data/
│ └── Processed_data/
├── Data_preprocessing/
├── Figure_production/
└── Outcomes/
```


---

## Folder descriptions

### Analysis

Contains all modelling code: the wave-level POMP baseline models and the
sensitivity analysis pipeline.

#### Analysis / Sensitivity_analysis

Contains the full sensitivity analysis pipeline and per-wave runner scripts.

| File | Description |
|---|---|
| `sensitivity_pipeline.R` | Master pipeline. Defines `wave_config` (priors, parameters, file paths for all four waves), helper functions for building POMP objects and running PMCMC, and the three degradation mechanisms (systematic reduction, stochastic Binomial thinning, temporal aggregation). Source this file before running any wave-specific script. |
| `run_wave1_sensitivity.R` | Runs all degradation scenarios for Wave 1 and saves results to `Outcomes/results/wave1/`. Sources `sensitivity_pipeline.R`. |
| `run_wave2_sensitivity.R` | Same structure for Wave 2. |
| `run_wave3_sensitivity.R` | Same structure for Wave 3. |
| `run_wave4_sensitivity.R` | Same structure for Wave 4. |
| `sensitivity_plots_wave{1-4}.R` | Generates the six sensitivity summary figures per wave (posterior means and CIs, CI width ratios, mean shift %, theory vs observed, coverage heatmap, systematic vs stochastic comparison). |
| `sensitivity_prior_posterior_plots.R` | Generates prior vs posterior density overlay plots for all degradation scenarios. Requires baseline `.rds` files to already exist in `Outcomes/`. |
| `wave1_identifiability_check.R` | Identifiability experiment for Wave 1: runs PMCMC from three deliberately off-starting parameter values under flat (non-informative) priors to test whether the posterior is genuinely data-driven or prior-dominated. |

**How the pipeline works:**

1. Source `sensitivity_pipeline.R` to load the `wave_config` list and all
   helper functions into your R session.
2. Run the desired `run_waveN_sensitivity.R` script. This will:
   - Load the processed data CSVs from `Data/Processed_data/`.
   - Build the POMP model object for that wave.
   - Run the two-phase PMCMC baseline (if no saved `.rds` exists yet).
   - Apply each degradation scenario and run PMCMC on the modified data.
   - Save all results as `.rds` files in `Outcomes/results/waveN/`.
   - Skip any run whose `.rds` already exists (safe to re-run after interruption).
3. Run the corresponding `sensitivity_plots_waveN.R` to produce figures.

All `.rds` result files follow the naming convention
`waveN_<degradation_type>_<level>.rds`
(e.g. `wave1_systematic_f0.30.rds`, `wave2_aggregate_b4.rds`).

---

#### Analysis / Wave-by-wave_baseline_model / full_POMP_pipeline_V1

First model version: single-population (non-age-structured) SEIR.
Kept for reference and to document the model building process described in
Appendix C of the thesis. Not used in the final analysis.

#### Analysis / Wave-by-wave_baseline_model / full_POMP_pipeline_V2_k_fixed_age_groups

Final baseline model version: two-group age-structured SEIR with
age-specific overdispersion parameters. These are the scripts that produce
the baseline posteriors used as the reference in all sensitivity comparisons.

| File | Description |
|---|---|
| `covid19_age_wave1.R` | Wave 1 (24 Feb – 5 Jul 2020, ancestral strain). Three NPI phases. Shared overdispersion k=10. Parameters: β₁=0.520, β₂=0.120, β₃=0.180; ρ_y=0.060, ρ_o=0.300. |
| `covid19_age_wave2_v2.R` | Wave 2 (7 Sep 2020 – 17 Jan 2021, ancestral strain). Four NPI phases (double-humped incidence). Age-specific k_y=30, k_o=10. |
| `covid19_age_wave3_v2.R` | Wave 3 (1 Feb – 20 Jun 2021, Alpha variant). Three NPI phases. k_y=30, k_o=10. |
| `covid19_age_wave4_final.R` | Wave 4 (21 Jun – 3 Oct 2021, Delta variant). Three NPI phases. k_y=30, k_o=10. |

Each script is self-contained and follows the same structure:
1. Load processed data from `Data/Processed_data/`.
2. Define the POMP model (rprocess, dmeasure, rinit, covar, partrans).
3. Run a simulation check at MLE starting values (50 trajectories).
4. Run two-phase PMCMC (Phase 1: 5,000 iterations diagonal proposal;
   Phase 2: adaptive MVN proposal, stop at ESS ≥ 200).
5. Assess convergence (trace plots, ESS, posterior densities).
6. Save the posterior chain and run a posterior predictive check.

---

### Data

#### Data / Original_data

Contains the raw RIVM line-list as downloaded:
`COVID-19_casus_landelijk_tm_03102021.csv`

#### Data / Processed_data

Eight CSV files (one per age group per wave) produced by `Data_preprocessing/covid19_age_preprocessing.py`, plus the full weekly all-ages series. These are the direct inputs to all modelling scripts.

| File | Contents |
|---|---|
| `covid19_wave{1-4}_young_NL.csv` | Weekly reported cases, young group (<60), within-wave week index and calendar date |
| `covid19_wave{1-4}_old_NL.csv` | Weekly reported cases, old group (≥60) |
| `covid19_weekly_allages_NL.csv` | Full Feb 2020 – Oct 2021 weekly series, all ages |

#### Data / Degraded_data

Degraded versions of the processed data generated during sensitivity
analysis runs. Produced automatically by `sensitivity_pipeline.R` and saved
here for reproducibility. Organised by wave and degradation type.

---

### Data_preprocessing

| File | Description |
|---|---|
| `covid19_age_preprocessing.py` | Full preprocessing pipeline: loads the raw RIVM line-list, maps ten-year age bands to two groups, excludes 211 records with ambiguous age (0.010% of total), aggregates to ISO calendar weeks, extracts the four wave windows, and writes the eight output CSVs to `Data/Processed_data/`. Requires Python 3.x with `pandas` and `numpy`. |

**To run:**
```bash
python covid19_age_preprocessing.py
```
Input: `Data/Original_data/COVID-19_casus_landelijk_tm_03102021.csv`
Output: all files in `Data/Processed_data/`

---

### Figure_production

Contains standalone scripts for producing thesis figures that are not
generated as a by-product of the sensitivity analysis pipeline.

| File | Description |
|---|---|
| `wave_phase_plots.R` | Phase timeline plots (data only, no model overlay) for all four waves. Saved to `Outcomes/results/phase_plots/`. |
| `wave_phase_ppc_plots.R` | Phase timeline plots with baseline 95% posterior predictive interval overlay. Requires baseline `.rds` files. |
| `fig_seir_inference.py` | Python/matplotlib script producing the two-panel thesis figure: (a) age-structured SEIR model diagram, (b) PMCMC inference pipeline. Outputs `fig_seir_inference.pdf` and `.png`. Requires Python 3.x with `matplotlib`. |

---

### Outcomes

```text
Contains all model outputs. Organised into subfolders mirroring the
analysis structure.
Outcomes/
├── results/
│ ├── wave1/ # sensitivity .rds files, Wave 1
│ ├── wave2/
│ ├── wave3/
│ ├── wave4/
│ ├── identifiability/ # wave1_identifiability_check outputs
│ ├── phase_plots/ # figure_production outputs
│ ├── prior_posterior/ # prior vs posterior figures
│ └── diagnostics/ # trace plots, posterior density plots
└── figures/ # final thesis figures (.pdf and .png)
```

---

## How to reproduce the full analysis

1. **Download raw data** from https://data.rivm.nl/covid-19/ and place in
   `Data/Original_data/`.

2. **Preprocess:**
```bash
   cd Data_preprocessing
   python covid19_age_preprocessing.py
```

3. **Run baseline models** (one per wave, can be run in any order):
```r
   source("Analysis/Wave-by-wave_baseline_model/full_POMP_pipeline_V2_k_fixed_age_groups/covid19_age_wave1.R")
   # repeat for waves 2–4
```

4. **Run sensitivity analysis:**
```r
   setwd("Analysis/Sensitivity_analysis/")
   source("run_wave1_sensitivity.R")
   # repeat for waves 2–4
```

5. **Produce figures:**
```r
   source("sensitivity_plots_wave1.R")
   # repeat for waves 2–4
   source("sensitivity_prior_posterior_plots.R")
```

**Expected runtime:** approximately 4–6 hours per wave for the full
sensitivity analysis (32 PMCMC runs × ~18,000 iterations each, with
Nₚ = 3,000 particles). A single baseline wave run takes approximately
30–60 minutes.





