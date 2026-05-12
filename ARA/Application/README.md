# Application Module for the Power Law Process (PLP)

This folder contains the application workflow for the **Power Law Process (PLP)** model using the sugarcane harvester dataset, as presented in:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The goal of this module is to illustrate how the PLP model can be applied to real failure data after preprocessing, adjustment of non-operational periods, and model fitting under frequentist and Bayesian approaches.

---

## Contents of this module

This application script performs the following steps:

1. Load and preprocess the sugarcane harvester data  
2. Select one vehicle for analysis  
3. Build failure-time vectors for the most frequent failure modes  
4. Adjust the time scale by removing off-season periods  
5. Fit the PLP model using:
   - frequentist estimation,
   - Bayesian estimation with independent priors,
   - Bayesian estimation with dependent priors  
6. Compare fitted models using AIC, BIC, and DIC  

---

## Main script

### `plp_application_sugarcane.R`

This is the main application script for the PLP model.

---

## Data preparation

The script reads the application dataset:

# ```r
## data <- read.csv("data.csv", sep = ";")

## Functions used in this application

This script relies on several functions defined in external files.  
Below is a concise description of each function to allow standalone use of this module.

### Core PLP functions

- `lambda_PLP(t, param, time_trunc)`  
  Computes the conditional intensity function of the Power Law Process.

- `Lambda_PLP(param, time_trunc)`  
  Computes the cumulative intensity evaluated at the truncation time.

- `log_like_PLP(param, time, time_trunc)`  
  Computes the log-likelihood of the PLP model for multiple systems.

---

### Data generation / preprocessing

- `adjust_offseason(x, start_obs, end_obs)`  
  Removes non-operational periods from the timeline by:
  - eliminating failures occurring during off-season,
  - compressing the time scale to reflect effective operational time.

---

### Estimation methods

- `frequentist_estimation(time, time_trunc)`  
  Computes maximum likelihood estimates using numerical optimization.

- `bayesian_estimation_ind(time, time_trunc, ...)`  
  Performs Bayesian inference using independent priors via Stan.

- `bayesian_estimation_dep(time, time_trunc, ...)`  
  Performs Bayesian inference using correlated priors via Stan (LKJ structure).

---

### Model comparison

- `calc_AIC_BIC(results, times, time_trunc)`  
  Computes AIC and BIC from the PLP log-likelihood evaluated at point estimates.

- `calc_DIC(fit_stan, times, time_trunc)`  
  Computes the Deviance Information Criterion using posterior samples.

---

## Source files

All functions are loaded via:

```r
source_files <- c(
  "plp_model_functions.R",
  "plp_stan_models.R",
  "plp_frequentist_estimation.R",
  "plp_bayesian_estimation_ind.R",
  "plp_bayesian_estimation_dep.R"
)