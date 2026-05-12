# Application Module for the Arithmetic Reduction of Intensity (ARI)

This folder contains the application workflow for the **Arithmetic Reduction of Intensity (ARI)** model using the sugarcane harvester dataset, as presented in:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The goal of this module is to illustrate how the ARI model can be applied to real failure data after preprocessing, adjustment of non-operational periods, and model fitting under frequentist and Bayesian approaches.

---

## Contents of this module

This application script performs the following steps:

1. Load and preprocess the sugarcane harvester data  
2. Select one vehicle for analysis  
3. Build failure-time vectors for the most frequent failure modes  
4. Adjust the time scale by removing off-season periods  
5. Fit the ARI model using:
   - frequentist estimation,
   - Bayesian estimation with independent priors,
   - Bayesian estimation with dependent priors  
6. Compare fitted models using AIC, BIC, and DIC  

---

## Main script

### `ari_application_sugarcane.R`

This is the main application script for the ARI model.

---

## Data preparation

The script reads the application dataset:

```r
data <- read.csv("data.csv", sep = ";")
```
## Functions used in this application

This script relies on several functions defined in external files.  
Below is a concise description of each function to allow standalone use of this module.

### Core ARI functions

- `lambda_ARI(t, param, times, time_trunc)`  
  Computes the conditional intensity function of the ARI model.

- `Lambda_ARI(param, times, time_trunc)`  
  Computes the cumulative intensity evaluated at the truncation time.

- `log_like_ARI(param, time, time_trunc)`  
  Computes the log-likelihood of the ARI model for multiple systems.

### Data preprocessing

- `adjust_offseason(x, start_obs, end_obs)`  
  Removes non-operational periods from the timeline by:
  - eliminating failures occurring during off-season,
  - compressing the time scale to reflect effective operating time.

### Estimation methods

- `frequentist_estimation(time, time_trunc)`  
  Computes maximum likelihood estimates using numerical optimization.

- `bayesian_estimation_ind(time, time_trunc, ...)`  
  Performs Bayesian inference using independent priors via Stan.

- `bayesian_estimation_dep(time, time_trunc, ...)`  
  Performs Bayesian inference using correlated priors via Stan (LKJ structure).

### Model comparison

- `calc_AIC_BIC(results, times, time_trunc)`  
  Computes AIC and BIC from the ARI log-likelihood evaluated at point estimates.

- `calc_DIC(fit_stan, times, time_trunc)`  
  Computes the Deviance Information Criterion using posterior samples.

## Source files

All functions are loaded via:

```r
source_files <- c(
  "ari_model_functions.R",
  "ari_stan_models.R",
  "ari_frequentist_fit.R",
  "ari_bayes_independent_fit.R",
  "ari_bayes_dependent_fit.R"
)
```

## Notes

- The ARI model modifies the failure intensity directly, unlike ARA which modifies the system’s virtual age.
- The parameter \( \rho \in [-1,1] \) governs the repair effect:
  - \( \rho = 0 \): minimal repair (ABAO),
  - \( 0 < \rho < 1 \): imperfect repair,
  - \( \rho = 1 \): optimal repair (intensity is reset after each failure),
  - \( \rho < 0 \): harmful repair (increases failure intensity).
- The application is based on multiple subsystems (failure modes) observed over a common time horizon.
- Off-season periods are explicitly removed to ensure that inference is based on effective operating time.
- Bayesian estimation is implemented using Stan and Hamiltonian Monte Carlo (HMC).
- Model comparison metrics (AIC, BIC, DIC) are used to assess relative model fit.
- Numerical stability of the intensity should be monitored, especially for extreme values of \( \rho \).
