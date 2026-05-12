# Functions for the Arithmetic Reduction of Age with Memory (ARAM) Model

This folder contains the core R functions used for the **Arithmetic Reduction of Age with Memory (ARAM)** model in the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The files in this folder provide the main components required for simulation, likelihood evaluation, frequentist estimation, and Bayesian inference under the ARAM formulation.

---

## Overview

The ARAM model extends classical imperfect repair models by modifying the **effective age of the system after each repair**, incorporating both the accumulated degradation and a repair-effect parameter.

Unlike the ARA model, ARAM introduces a transformation of the repair parameter to ensure numerical stability and valid model behavior across parameter regions.

The conditional intensity function depends on the adjusted age:

\[
\lambda_{ARAM}(t) = \mu_T \frac{\beta}{T} \left( t - s(\beta)\,h(\theta, \beta)\,t_{i-1} \right)^{\beta - 1},
\]

where:

- \( t_{i-1} \) is the most recent failure time before \( t \),
- \( T \) is the truncation horizon,
- \( \beta > 0 \) is the shape parameter,
- \( \mu_T > 0 \) is the scale parameter,
- \( \theta \in [-1,1] \) controls the repair effect,
- \( h(\theta, \beta) \) is a transformation ensuring valid domain behavior,
- \( s(\beta) \in \{-1,1\} \) is a sign adjustment depending on \( \beta \).

---

## Files included in this folder

### `aram_model_functions.R`
Defines the core mathematical functions for the ARAM model.

**Main functions:**
- `gen_failures_ARAM()`  
  Simulates failure times for \(k\) independent systems using numerical inversion of the cumulative intensity.
- `lambda_ARAM()`  
  Evaluates the ARAM conditional intensity function.
- `Lambda_ARAM()`  
  Evaluates the cumulative intensity at the truncation horizon.
- `log_like_ARAM()`  
  Computes the log-likelihood for multiple independent systems.
- `indicator()`  
  Returns the most recent failure time before a given time \(t\).
- `h()`  
  Transformation of the repair parameter ensuring numerical stability.
- `sign()`  
  Sign correction used in the ARAM formulation.

This file provides the core building blocks used across estimation procedures.

---

### `aram_stan_models.R`
Defines the Stan code used for Bayesian inference under the ARAM model.

**Main Stan model objects:**
- `model.stan.ind`  
  Stan model with **independent priors** on the parameters.
- `model.stan.dep`  
  Stan model with **dependent priors**, using an LKJ-based correlation structure.

These models implement the ARAM likelihood together with the prior specifications described in the manuscript.

---

### `aram_frequentist_fit.R`
Provides the frequentist estimation routine for the ARAM model.

**Main function:**
- `frequentist_estimation()`  
  Computes maximum likelihood estimates via numerical optimization using `optimx`, with constraints:
  - \( \beta > 0 \),
  - \( \mu_T > 0 \),
  - \( \theta \in [-1,1] \).

This function is primarily used for benchmarking against Bayesian estimators.

---

### `aram_bayes_independent_fit.R`
Provides Bayesian estimation under the **independent-prior** specification.

**Main function:**
- `bayesian_estimation_ind()`  
  Runs posterior sampling in Stan, extracts posterior summaries, and optionally computes empirical coverage.

This function expects a compiled Stan model consistent with `model.stan.ind`.

---

### `aram_bayes_dependent_fit.R`
Provides Bayesian estimation under the **dependent-prior** specification.

**Main function:**
- `bayesian_estimation_dep()`  
  Performs posterior sampling using correlated priors via an LKJ structure, and returns posterior summaries and optional coverage.

This function expects a compiled Stan model consistent with `model.stan.dep`.

---

## Required packages

The functions in this folder rely on the following R packages:

- `rstan`
- `coda`
- `optimx`

Additional packages may appear in specific files depending on the workflow.

---

## Recommended workflow

A typical workflow is:

1. Use `aram_model_functions.R` to define intensity, cumulative intensity, likelihood, and simulation functions.
2. Use `aram_stan_models.R` to define the Stan model code.
3. Compile the desired Stan model:
   - `model.stan.ind` for independent priors
   - `model.stan.dep` for dependent priors
4. Run estimation using:
   - `frequentist_estimation()` for maximum likelihood
   - `bayesian_estimation_ind()` for independent priors
   - `bayesian_estimation_dep()` for dependent priors

---

## Notes

- The ARAM model introduces additional flexibility compared to ARA by transforming the repair-effect parameter.
- The transformation \( h(\theta, \beta) \) ensures numerical stability and avoids invalid values in the intensity function.
- The model accommodates both beneficial and harmful repair effects.
- Care must be taken with extreme parameter values, as the transformed intensity can become unstable.
- The Bayesian implementation relies on Hamiltonian Monte Carlo (HMC) via Stan.
- All functions are designed to handle multiple independent systems observed up to a common truncation horizon.
- Naming conventions were standardized in English to ensure consistency across the repository.

---

## Purpose of this folder

This folder contains the reusable functions associated with the ARAM model.

Simulation scripts, application scripts, and datasets should be stored separately from these core functions whenever possible.