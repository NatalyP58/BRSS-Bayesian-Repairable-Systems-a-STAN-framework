# Functions for the Arithmetic Reduction of Intensity (ARI) Model

This folder contains the core R functions used for the **Arithmetic Reduction of Intensity (ARI)** model in the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The files in this folder provide the main components required for simulation, likelihood evaluation, frequentist estimation, and Bayesian inference under the ARI formulation.

---

## Overview

The ARI model defines a class of imperfect repair models in which maintenance actions directly reduce the **failure intensity**, rather than modifying the system’s virtual age.

Under this formulation, the conditional intensity function is defined as

\[
\lambda_{ARI}(t) = \mu_T \frac{\beta}{T} \left( t^{\beta - 1} - \rho \, t_{i-1}^{\beta - 1} \right),
\]

where:
- \( t_{i-1} \) denotes the most recent failure time prior to \( t \),
- \( T \) is the truncation horizon,
- \( \beta > 0 \) is the shape parameter,
- \( \mu_T > 0 \) is the scale parameter,
- \( \rho \in [-1,1] \) is the repair-effect parameter.

The parameter \( \rho \) governs the repair effect:
- \( \rho = 0 \): minimal repair (ABAO),
- \( 0 < \rho < 1 \): effective but imperfect repair,
- \( \rho = 1 \): optimal repair (intensity reset to zero after each event),
- \( \rho < 0 \): harmful repair (increases the failure intensity).

---

## Files included in this folder

### `ari_model_functions.R`
Defines the core mathematical functions for the ARI model.

**Main functions:**
- `gen_failures_ARI()`  
  Simulates failure times for \(k\) independent systems under the ARI model using numerical inversion of the cumulative intensity.
- `lambda_ARI()`  
  Evaluates the ARI conditional intensity function incorporating the repair effect.
- `Lambda_ARI()`  
  Evaluates the cumulative intensity at the truncation horizon, accounting for all inter-failure intervals.
- `log_like_ARI()`  
  Computes the log-likelihood for multiple independent systems.
- `indicator()`  
  Returns the most recent failure time prior to a given time \(t\).

This file provides the fundamental building blocks used by both the frequentist and Bayesian implementations.

---

### `ari_stan_models.R`
Defines the Stan code used for Bayesian inference under the ARI model.

**Main Stan model objects:**
- `model.stan.ind`  
  Stan model with **independent priors**, where:
  - \( \beta > 0 \),
  - \( \mu_T > 0 \),
  - \( \rho \in [-1,1] \) is directly bounded.
  
- `model.stan.dep`  
  Stan model with **dependent priors**, using:
  - an LKJ-based correlation structure,
  - a latent Gaussian representation,
  - transformations ensuring valid parameter support (including \( \rho = \tanh(\cdot) \)).

These models implement the ARI likelihood together with the prior structures described in the manuscript.

---

### `ari_frequentist_fit.R`
Provides the frequentist estimation routine for the ARI model.

**Main function:**
- `frequentist_estimation()`  
  Obtains parameter estimates by numerical maximization of the ARI log-likelihood using `optimx`, with constraints ensuring:
  - \( \beta > 0 \),
  - \( \mu_T > 0 \),
  - \( \rho \in [-1,1] \).

This file is intended for comparison with the Bayesian estimators and for simulation studies.

---

### `ari_bayes_independent_fit.R`
Provides Bayesian estimation under the **independent-prior** specification.

**Main function:**
- `bayesian_estimation_ind()`  
  Runs posterior sampling in Stan, extracts posterior summaries, and optionally evaluates empirical coverage in simulation settings.

This function expects a compiled Stan model consistent with `model.stan.ind`.

---

### `ari_bayes_dependent_fit.R`
Provides Bayesian estimation under the **dependent-prior** specification.

**Main function:**
- `bayesian_estimation_dep()`  
  Runs posterior sampling in Stan under the correlated prior structure, extracts posterior summaries, and optionally evaluates empirical coverage in simulation settings.

This function expects a compiled Stan model consistent with `model.stan.dep`.

---

## Required packages

The functions in this folder rely on the following R packages:

- `rstan`
- `coda`
- `optimx`

---

## Recommended workflow

A typical workflow is:

1. Use `ari_model_functions.R` to define the ARI intensity, cumulative intensity, likelihood, and simulation functions.
2. Use `ari_stan_models.R` to define the Stan model code.
3. Compile the desired Stan model:
   - `model.stan.ind` for independent priors
   - `model.stan.dep` for dependent priors
4. Run estimation using:
   - `frequentist_estimation()` for maximum-likelihood estimation
   - `bayesian_estimation_ind()` for Bayesian estimation with independent priors
   - `bayesian_estimation_dep()` for Bayesian estimation with dependent priors

---

## Notes

- The ARI model modifies the **intensity directly**, unlike ARA which modifies the virtual age.
- The parameter \( \rho \) provides a flexible representation of repair effectiveness.
- Care must be taken to ensure that the intensity remains positive when evaluating the likelihood.
- The Bayesian functions are designed to work with multiple independent systems observed up to a common truncation horizon.
- In simulation studies, the Bayesian routines may also return empirical coverage based on 95% HPD intervals.
- Naming conventions were standardized in English to improve readability and consistency across the repository.

---

## Purpose of this folder

This folder is intended to contain the reusable functions associated with the ARI model.  
Application scripts, simulation scripts, and data files should be stored separately from these core functions whenever possible.