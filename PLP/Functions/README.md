# Functions for the Perfect Repair (PR) Model

This folder contains the core R functions used for the **Perfect Repair (PR)** model in the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The files in this folder provide the main components required for simulation, likelihood evaluation, frequentist estimation, and Bayesian inference under the PR formulation.

---

## Overview

The PR model assumes that after each failure, the system is restored to an **as-good-as-new condition**, meaning that the failure intensity depends only on the time since the last failure.

The implementation follows the parametrization adopted in the manuscript, where the conditional intensity function is defined as

\[
\lambda_{PR}(t) = \mu_T \frac{\beta}{T} (t - t_{i-1})^{\beta - 1},
\]

where \( t_{i-1} \) denotes the most recent failure time before \( t \), \( T \) is the truncation horizon, \( \beta \geq 1 \) is the shape parameter, and \( \mu_T > 0 \) is the scale parameter.

---

## Files included in this folder

### `pr_model_functions.R`
Defines the core mathematical functions for the PR model.

**Main functions:**
- `gen_failures_PR()`  
  Simulates failure times for \(k\) independent systems under the PR model using numerical inversion of the cumulative intensity.
- `lambda_PR()`  
  Evaluates the PR conditional intensity function based on the time since the last failure.
- `Lambda_PR()`  
  Evaluates the cumulative intensity at the truncation horizon, accounting for all inter-failure intervals.
- `log_like_PR()`  
  Computes the log-likelihood for multiple independent systems.
- `indicator()`  
  Returns the most recent failure time before a given time \(t\).

This file provides the fundamental building blocks used by both the frequentist and Bayesian implementations.

---

### `pr_stan_models.R`
Defines the Stan code used for Bayesian inference under the PR model.

**Main Stan model objects:**
- `model.stan.ind`  
  Stan model with **independent priors**, where \( \beta \) is defined via a shifted parameterization to enforce \( \beta \geq 1 \).
- `model.stan.dep`  
  Stan model with **dependent priors**, using an LKJ-based correlation structure and latent variable representation.

These models implement the PR likelihood together with the corresponding prior structures described in the manuscript.

---

### `pr_frequentist_fit.R`
Provides the frequentist estimation routine for the PR model.

**Main function:**
- `frequentist_estimation()`  
  Obtains parameter estimates by numerical maximization of the PR log-likelihood using `optimx`, with constraints ensuring \( \beta \geq 1 \) and \( \mu_T > 0 \).

This file is intended for comparison with the Bayesian estimators and for simulation studies.

---

### `pr_bayes_independent_fit.R`
Provides Bayesian estimation under the **independent-prior** specification.

**Main function:**
- `bayesian_estimation_ind()`  
  Runs posterior sampling in Stan, extracts posterior summaries, and optionally evaluates empirical coverage in simulation settings.

This function expects a compiled Stan model consistent with `model.stan.ind`.

---

### `pr_bayes_dependent_fit.R`
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

1. Use `pr_model_functions.R` to define the PR intensity, cumulative intensity, likelihood, and simulation functions.
2. Use `pr_stan_models.R` to define the Stan model code.
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
- The PR model assumes that the system is fully restored after each failure, which is appropriate in contexts where maintenance actions eliminate accumulated degradation.
- The parameterization enforces \( \beta \geq 1 \), ensuring a non-decreasing failure intensity between failures.
- The Bayesian functions are designed to work with multiple independent systems observed up to a common truncation horizon.
- In simulation studies, the Bayesian routines may also return empirical coverage based on 95% HPD intervals.
- Naming conventions were standardized in English to improve readability and consistency across the repository.

---

## Purpose of this folder

This folder is intended to contain the reusable functions associated with the PR model.  
Application scripts, simulation scripts, and data files should be stored separately from these core functions whenever possible.