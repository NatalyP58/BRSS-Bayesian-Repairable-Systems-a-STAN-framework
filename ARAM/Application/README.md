# Application Module for the Arithmetic Reduction of Age Modified (ARAM) Model

This folder contains the application workflow for the **Arithmetic Reduction of Age Modified (ARAM)** model using the sugarcane harvester dataset, as presented in:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The goal of this module is to illustrate how the ARAM model can be applied to real failure data after preprocessing, adjustment of non-operational periods, and model fitting under frequentist and Bayesian approaches.

---

## Contents of this module

This application script performs the following steps:

1. Load and preprocess the sugarcane harvester data  
2. Select one vehicle for analysis  
3. Build failure-time vectors for the most frequent failure modes  
4. Adjust the time scale by removing off-season periods  
5. Fit the ARAM model using:
   - frequentist estimation,
   - Bayesian estimation with independent priors,
   - Bayesian estimation with dependent priors  
6. Compare fitted models using AIC, BIC, and DIC  

---

## Main script

### `aram_application_sugarcane.R`

This is the main application script for the ARAM model.

---

## Data preparation

The script reads the application dataset:

```r
data <- read.csv("application/sugarcane_data.csv", sep = ";")
```

---

## Functions used in this application

This script relies on several functions defined in external files.  
Below is a concise description of each function to allow standalone use of this module.

### Core ARAM functions

- `lambda_ARAM(t, param, time_trunc)`  
  Computes the conditional intensity function of the ARAM model.

- `Lambda_ARAM(param, time_trunc)`  
  Computes the cumulative intensity evaluated at the truncation time.

- `log_like_ARAM(param, time, time_trunc)`  
  Computes the log-likelihood of the ARAM model for multiple systems.

---

### Data preprocessing

- `adjust_offseason(x, start_obs, end_obs)`  
  Removes non-operational periods from the timeline by:
  - eliminating failures occurring during off-season,
  - compressing the time scale to reflect effective operating time.

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
  Computes AIC and BIC from the ARAM log-likelihood evaluated at point estimates.

- `calc_DIC(fit_stan, times, time_trunc)`  
  Computes the Deviance Information Criterion using posterior samples.

---

## Source files

All functions are loaded via:

```r
source_files <- c(
  file.path("functions", "aram_model_functions.R"),
  file.path("functions", "aram_stan_models.R"),
  file.path("functions", "aram_frequentist_fit.R"),
  file.path("functions", "aram_bayes_dependent_fit.R"),
  file.path("functions", "aram_bayes_independent_fit.R")
)
```

---

## Notes

- The ARAM (Arithmetic Reduction of Age Modified) model extends the classical ARA framework by allowing both beneficial and detrimental repair effects through a modified virtual-age structure.

- Unlike the classical ARA model, ARAM allows the repair parameter \( \theta \in [-1,1] \), providing greater flexibility for modeling practical repairable-system behavior.

- The model incorporates the transformation function \( h(\theta) \), which modifies the system’s effective age while preserving the interpretation of the repair effect.

- The parameter \( \theta \in [-1,1] \) governs the repair effect:
  - \( \theta = 0 \): minimal repair (ABAO),
  - \( 0 < \theta < 1 \): imperfect beneficial repair,
  - \( \theta = 1 \): perfect repair (AGAN),
  - \( \theta < 0 \): harmful repair.

- The function \( h(\theta) \) is constructed to ensure model stability and depends on the value of the trend parameter \( \beta \).

- The application is based on multiple subsystems (failure modes) observed over a common time horizon.

- Off-season periods are explicitly removed to ensure that inference is based on effective operating time.

- Bayesian estimation is implemented using Stan and Hamiltonian Monte Carlo (HMC).

- Model comparison metrics (AIC, BIC, DIC) are used to assess relative model fit.

- Numerical stability should be monitored when \( \theta \) approaches its boundary values or when the transformed virtual age becomes close to zero.

- The implementation supports both independent and correlated prior structures for Bayesian inference.

---