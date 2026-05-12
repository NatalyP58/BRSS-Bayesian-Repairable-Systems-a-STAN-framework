# Application Module for the Partial Imperfect Repair (PIR) Model

This folder contains the application workflow for the **Partial Imperfect Repair (PIR)** model using the sugarcane harvester dataset, as presented in:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The goal of this module is to illustrate how the PIR model can be applied to real failure data after preprocessing, adjustment of non-operational periods, definition of change points, and model fitting under frequentist and Bayesian approaches.

---

## Contents of this module

This application script performs the following steps:

1. Load and preprocess the sugarcane harvester data  
2. Select one vehicle for analysis  
3. Build failure-time vectors for the most frequent failure modes  
4. Adjust the time scale by removing off-season periods  
5. Define the change-point vector \( \tau \)  
6. Fit the PIR model using:
   - frequentist estimation,
   - Bayesian estimation with independent priors,
   - Bayesian estimation with dependent priors  
7. Compare fitted models using AIC, BIC, and DIC  

---

## Main script

### `pir_application_sugarcane.R`

This is the main application script for the PIR model.

---

## Data preparation

The script reads the application dataset:

```r
data <- read.csv("Application/sugarcane_data.csv", sep = ";")
```

The selected subsystems are renamed for presentation as:

```r
names(times_offseason) <- c("Line Divider", "Diesel Engine", "Elevator", "Transmission")
```

The change-point vector is specified as:

```r
tau <- c(5504.5167, 11633.933)
```

---

## Functions used in this application

This script relies on several functions defined in external files.  
Below is a concise description of each function to allow standalone use of this module.

### Core PIR functions

- `lambda_PIR(t, param, tau, time_trunc)`  
  Computes the conditional intensity function of the PIR model.

- `Lambda_PIR(param, tau, time_trunc)`  
  Computes the cumulative intensity evaluated at the truncation time.

- `log_like_PIR(param, time, tau, time_trunc)`  
  Computes the log-likelihood of the PIR model for multiple systems.

---

### Data preprocessing

- `adjust_offseason(x, start_obs, end_obs)`  
  Removes non-operational periods from the timeline by:
  - eliminating failures occurring during off-season,
  - compressing the time scale to reflect effective operating time.

---

### Estimation methods

- `frequentist_estimation(time, time_trunc, tau)`  
  Computes maximum likelihood estimates using numerical optimization.

- `bayesian_estimation_ind(time, time_trunc, ..., tau)`  
  Performs Bayesian inference using independent priors via Stan.

- `bayesian_estimation_dep(time, time_trunc, ..., tau)`  
  Performs Bayesian inference using correlated priors via Stan (LKJ structure).

---

### Model comparison

- `calc_AIC_BIC(results, times, tau, time_trunc)`  
  Computes AIC and BIC from the PIR log-likelihood evaluated at point estimates.

- `calc_DIC(fit_stan, times, tau, time_trunc)`  
  Computes the Deviance Information Criterion using posterior samples.

---

## Source files

All functions are loaded via:

```r
source_files <- c(
  file.path("functions", "pir_model_functions.R"),
  file.path("functions", "pir_stan_models.R"),
  file.path("functions", "pir_frequentist_fit.R"),
  file.path("functions", "pir_bayes_dependent_fit.R"),
  file.path("functions", "pir_bayes_independent_fit.R")
)
```

---

## Notes

- The PIR (Partial Imperfect Repair) model extends classical repair models by introducing a change-point structure \( \tau \), which partitions the observation window into intervals where the failure dynamics are evaluated.

- The vector \( \tau = (\tau_1, \ldots, \tau_m) \) defines the change points and must be specified as part of the model input.

- In this application, the change-point vector is fixed as:
  \[
  \tau = (5504.5167,\ 11633.933).
  \]

- Unlike ARI or ARAM, the PIR model structures the repair effect through change-point intervals rather than acting directly only on intensity or virtual age.

- The parameter \( \theta \in [-1,1] \) governs the repair effect:
  - \( \theta = 0 \): minimal repair (ABAO),
  - \( 0 < \theta < 1 \): imperfect repair,
  - \( \theta = 1 \): perfect repair (AGAN),
  - \( \theta < 0 \): harmful repair.

- The application is based on multiple subsystems, or failure modes, observed over a common time horizon.

- Off-season periods are explicitly removed to ensure that inference is based on effective operating time.

- Bayesian estimation is implemented using Stan and Hamiltonian Monte Carlo (HMC).

- Model comparison metrics (AIC, BIC, DIC) are used to assess relative model fit.

- Numerical stability should be monitored when \( \theta \) approaches its boundary values or when change points are near the extremes of the observation window.

- The implementation supports both independent and correlated prior structures for Bayesian inference.

---