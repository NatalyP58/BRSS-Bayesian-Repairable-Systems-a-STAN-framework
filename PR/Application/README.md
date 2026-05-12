# Application Module for the Perfect Repair (PR) Model

This folder contains the application workflow for the **Perfect Repair (PR)** model using the sugarcane harvester dataset, as presented in:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The goal of this module is to illustrate how the PR model can be applied to real failure data after preprocessing, adjustment of non-operational periods, and model fitting under frequentist and Bayesian approaches.

---

## Contents of this module

This application script performs the following steps:

1. Load and preprocess the sugarcane harvester data  
2. Select one vehicle for analysis  
3. Build failure-time vectors for the most frequent failure modes  
4. Adjust the time scale by removing off-season periods  
5. Fit the PR model using:
   - frequentist estimation,
   - Bayesian estimation with independent priors,
   - Bayesian estimation with dependent priors  
6. Compare fitted models using AIC, BIC, and DIC  

---

## Main script

### `pr_application_sugarcane.R`

This is the main application script for the PR model.

---

## Data preparation

The script reads the application dataset:

#```r
data <- read.csv("data.csv", sep = ";")

## Off-season adjustment

The function:

`adjust_offseason(x, start_obs, end_obs)`

removes non-operational periods (December–March) from the time scale.

This step:

- eliminates failures occurring during off-season,
- compresses the timeline by removing inactivity gaps,
- ensures that the PR model is fitted only on effective operating time.

The adjusted failure times are stored in:

`times_offseason`

---

## Model functions used

This script depends on the functions defined in:

`functions/`

### Core model functions
- `lambda_PR()` → conditional intensity  
- `Lambda_PR()` → cumulative intensity  
- `log_like_PR()` → log-likelihood  

### Estimation functions
- `frequentist_estimation()` → maximum likelihood estimation  
- `bayesian_estimation_ind()` → Bayesian inference with independent priors  
- `bayesian_estimation_dep()` → Bayesian inference with dependent priors  

### Stan models
- `model.stan.ind` → independent priors  
- `model.stan.dep` → dependent priors  

---

## Parallel estimation

The main estimation routine is:

`run_application_parallel_cluster()`

This function:

- initializes a parallel cluster,
- ensures reproducible random number generation,
- compiles Stan models (if needed),
- runs frequentist and Bayesian estimation,
- returns parameter estimates, fitted objects, and diagnostics.

---

## Output structure

The resulting object:

`results`

contains:

- `freq_est` → frequentist estimates  
- `bayes_est_ind` → Bayesian estimates (independent priors)  
- `bayes_est_dep` → Bayesian estimates (dependent priors)  
- `fit_ind`, `fit_dep` → Stan fit objects  
- `coverage_ind`, `coverage_dep` → (if applicable) coverage indicators  
- success flags and seed information  

Results are saved as:

`results_application01.rds`

---

## Model comparison

### AIC and BIC

The function:

`calc_AIC_BIC(results, times, time_trunc)`

computes:

- AIC and BIC for frequentist estimates,
- AIC and BIC for Bayesian estimates using posterior summaries.

---

### DIC

The function:

`calc_DIC(fit_stan, times, time_trunc)`

computes the Deviance Information Criterion using:

- posterior mean deviance,
- effective number of parameters,
- full posterior distribution.

---

## Required packages

- `tidyverse`
- `dplyr`
- `rstan`
- `coda`
- `parallel`

---

## Notes

- The PR model assumes **perfect repair**, meaning that after each failure the system returns to an as-good-as-new condition.
- The intensity depends only on the time since the last failure.
- The parameterization enforces \( \beta \geq 1 \), ensuring physically meaningful behavior.
- The off-season adjustment is critical in this dataset to avoid bias due to long inactivity periods.
- The implementation is designed for **reproducibility and methodological clarity**, consistent with the manuscript.

---

## Purpose of this folder

This folder demonstrates the full application of the PR model to real-world data.

It complements:

- the **functions module** (core definitions),
- the **simulation module** (controlled experiments),

by providing a **real-data case study**.