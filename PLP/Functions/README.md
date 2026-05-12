# Functions for the Power Law Process (PLP) Model

This folder contains the core R functions used for the Power Law Process (PLP) model in the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The files in this folder provide the main components required for simulation, likelihood evaluation, frequentist estimation, and Bayesian inference under the PLP formulation.

---

# Overview

The Power Law Process (PLP) model is a nonhomogeneous Poisson process commonly used to model repairable systems under minimal repair assumptions. Under this formulation, each repair restores the system to the same stochastic condition it had immediately before failure.

The implementation follows the parametrization adopted in the manuscript, where the conditional intensity function is defined as

\[
\lambda_{PLP}(t) = \mu_T \frac{\beta}{T} t^{\beta - 1},
\]

where \(T\) is the truncation horizon, \(\beta > 0\) is the shape parameter, and \(\mu_T > 0\) is the scale parameter.

The parameter \(\beta\) governs the trend behavior of the failure process:

- \( \beta > 1 \): increasing failure intensity,
- \( \beta = 1 \): homogeneous Poisson process,
- \( \beta < 1 \): decreasing failure intensity.

---

# Files included in this folder

## `plp_model_functions.R`

Defines the core mathematical functions for the PLP model.

### Main functions:

### `gen_failures_PLP()`

Simulates failure times for \(k\) independent systems under the PLP model using inversion of the cumulative intensity function.

### `lambda_PLP()`

Evaluates the PLP conditional intensity function.

### `Lambda_PLP()`

Evaluates the cumulative intensity at the truncation horizon.

### `log_like_PLP()`

Computes the log-likelihood for multiple independent systems.

This file provides the fundamental building blocks used by both the frequentist and Bayesian implementations.

---

## `plp_stan_models.R`

Defines the Stan code used for Bayesian inference under the PLP model.

### Main Stan model objects:

### `model.stan.ind`

Stan model with independent priors for the PLP parameters.

### `model.stan.dep`

Stan model with dependent priors using an LKJ-based correlation structure and latent variable representation.

These models implement the PLP likelihood together with the corresponding prior structures described in the manuscript.

---

## `plp_frequentist_fit.R`

Provides the frequentist estimation routine for the PLP model.

### Main function:

### `frequentist_estimation()`

Obtains parameter estimates by numerical maximization of the PLP log-likelihood using `optimx`.

This file is intended for comparison with the Bayesian estimators and for simulation studies.

---

## `plp_bayes_independent_fit.R`

Provides Bayesian estimation under the independent-prior specification.

### Main function:

### `bayesian_estimation_ind()`

Runs posterior sampling in Stan, extracts posterior summaries, and optionally evaluates empirical coverage in simulation settings.

This function expects a compiled Stan model consistent with `model.stan.ind`.

---

## `plp_bayes_dependent_fit.R`

Provides Bayesian estimation under the dependent-prior specification.

### Main function:

### `bayesian_estimation_dep()`

Runs posterior sampling in Stan under the correlated prior structure, extracts posterior summaries, and optionally evaluates empirical coverage in simulation settings.

This function expects a compiled Stan model consistent with `model.stan.dep`.

---

# Required packages

The functions in this folder rely on the following R packages:

- `rstan`
- `coda`
- `optimx`

Additional packages may appear in specific files depending on the workflow.

---

# Recommended workflow

A typical workflow is:

1. Use `plp_model_functions.R` to define the PLP intensity, cumulative intensity, likelihood, and simulation functions.
2. Use `plp_stan_models.R` to define the Stan model code.
3. Compile the desired Stan model:
   - `model.stan.ind` for independent priors,
   - `model.stan.dep` for dependent priors.
4. Run estimation using:
   - `frequentist_estimation()` for maximum-likelihood estimation,
   - `bayesian_estimation_ind()` for Bayesian estimation with independent priors,
   - `bayesian_estimation_dep()` for Bayesian estimation with dependent priors.

---

# Notes

- The code in this folder is written to support both methodological illustration and reproducible implementation.

- The PLP model corresponds to a minimal repair process, where repairs do not modify the future evolution of the failure intensity.

- The implementation supports both increasing and decreasing failure trends through the parameter \(\beta\).

- The Bayesian functions are designed to work with multiple independent systems observed up to a common truncation horizon.

- In simulation studies, the Bayesian routines may also return empirical coverage based on 95% HPD intervals.

- Naming conventions were standardized in English to improve readability and consistency across the repository.

---

# Purpose of this folder

This folder is intended to contain the reusable functions associated with the PLP model.

Application scripts, simulation scripts, and data files should be stored separately from these core functions whenever possible.