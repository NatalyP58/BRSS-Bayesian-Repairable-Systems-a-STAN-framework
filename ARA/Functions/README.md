# Functions for the Arithmetic Reduction of Age (ARA) Model

This folder contains the core R functions used for the **Arithmetic Reduction of Age (ARA)** model in the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The files in this folder provide the main components required for simulation, likelihood evaluation, frequentist estimation, and Bayesian inference under the ARA formulation.

---

## Overview

The ARA model extends the classical repairable-systems framework by introducing a **repair-effect parameter** that modifies the effective age of the system after each failure.

Under this formulation, the conditional intensity function is defined as

\[
\lambda_{ARA}(t) = \mu_T \frac{\beta}{T} \left(t - \theta t_{i-1}\right)^{\beta - 1},
\]

where:
- \( t_{i-1} \) denotes the most recent failure time before \( t \),
- \( T \) is the truncation horizon,
- \( \beta \geq 1 \) is the shape parameter,
- \( \mu_T > 0 \) is the scale parameter,
- \( \theta \in [-1,1] \) is the repair-effect parameter.

The parameter \( \theta \) controls the effectiveness of repairs:
- \( \theta = 0 \): perfect repair (PR model),
- \( \theta = 1 \): minimal repair,
- intermediate values represent imperfect repair,
- negative values allow modeling over-repair effects.

---

## Files included in this folder

### `ara_model_functions.R`
Defines the core mathematical functions for the ARA model.

**Main functions:**
- `gen_failures_ARA()`  
  Simulates failure times for \(k\) independent systems under the ARA model using numerical inversion of the cumulative intensity.
- `lambda_ARA()`  
  Evaluates the ARA conditional intensity function incorporating the repair effect.
- `Lambda_ARA()`  
  Evaluates the cumulative intensity at the truncation horizon, accounting for all inter-failure intervals.
- `log_like_ARA()`  
  Computes the log-likelihood for multiple independent systems.
- `indicator()`  
  Returns the most recent failure time before a given time \(t\).

This file provides the fundamental building blocks used by both the frequentist and Bayesian implementations.

---

### `ara_stan_models.R`
Defines the Stan code used for Bayesian inference under the ARA model.

**Main Stan model objects:**
- `model.stan.ind`  
  Stan model with **independent priors**, where:
  - \( \beta \geq 1 \) is enforced via a shifted parameterization,
  - \( \mu_T > 0 \),
  - \( \theta \in [-1,1] \) is directly bounded.
  
- `model.stan.dep`  
  Stan model with **dependent priors**, using:
  - an LKJ-based correlation structure,
  - a latent Gaussian representation,
  - transformations ensuring valid parameter support.

These models implement the ARA likelihood together with the prior structures described in the manuscript.

---

### `ara_frequentist_fit.R`
Provides the frequentist estimation routine for the ARA model.

**Main function:**
- `frequentist_estimation()`  
  Obtains parameter estimates by numerical maximization of the ARA log-likelihood using `optimx`, with constraints ensuring:
  - \( \beta \geq 1 \),
  - \( \mu_T > 0 \),
  - \( \theta \in [-1,1] \).

This file is intended for comparison with the Bayesian estimators and for simulation studies.

---

### `ara_bayes_independent_fit.R`
Provides Bayesian estimation under the **independent-prior** specification.

**Main function:**
- `bayesian_estimation_ind()`  
  Runs posterior sampling in Stan, extracts posterior summaries, and optionally evaluates empirical coverage in simulation settings.

This function expects a compiled Stan model consistent with `model.stan.ind`.

---

### `ara_bayes_dependent_fit.R`
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

Additional packages may appear in specific files depending on the workflow.

---

## Recommended workflow

A typical workflow is:

1. Use `ara_model_functions.R` to define the ARA intensity, cumulative intensity, likelihood, and simulation functions.
2. Use `ara_stan_models.R` to define the Stan model code.
3. Compile the desired Stan model:
   - `model.stan.ind` for independent priors
   - `model.stan.dep` for dependent priors
4. Run estimation using:
   - `frequentist_estimation()` for maximum-likelihood estimation
   - `bayesian_estimation_ind()` for Bayesian estimation with independent priors
   - `bayesian_estimation_dep()` for Bayesian estimation with dependent priors

---

## Notes

- The code in this folder is written to support both **methodological illustration** and **reproducible implementation**.
- The ARA model generalizes several classical repair models through the parameter \( \theta \).
- The parameterization enforces \( \beta \geq 1 \), ensuring a non-decreasing failure intensity between failures.
- The Bayesian functions are designed to work with multiple independent systems observed up to a common truncation horizon.
- In simulation studies, the Bayesian routines may also return empirical coverage based on 95% HPD intervals.
- Naming conventions were standardized in English to improve readability and consistency across the repository.

---

## Purpose of this folder

This folder is intended to contain the reusable functions associated with the ARA model.  
Application scripts, simulation scripts, and data files should be stored separately from these core functions whenever possible.