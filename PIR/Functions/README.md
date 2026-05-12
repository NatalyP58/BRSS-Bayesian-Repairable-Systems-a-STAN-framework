# Functions Module for the Partial Imperfect Repair (PIR) Model

This folder contains the core implementation of the **Partial Imperfect Repair (PIR)** model developed in the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

This module provides all fundamental components required to define, simulate, and estimate the PIR model under both frequentist and Bayesian approaches.

---

## Contents of this module

This folder includes the following files:

- `pir_model_functions.R`
- `pir_stan_models.R`
- `pir_frequentist_fit.R`
- `pir_bayes_independent_fit.R`
- `pir_bayes_dependent_fit.R`

---

## 1. Core model functions

File: `pir_model_functions.R`

This file implements the mathematical structure of the PIR model.

### Included functions

- `gen_failures_PIR(k, param, time_trunc, tau)`  
  Simulates failure times for \( k \) independent systems using numerical inversion of the cumulative intensity.

- `lambda_PIR(t, param, tau, time_trunc)`  
  Computes the conditional intensity function of the PIR model.

- `Lambda_PIR(param, tau, time_trunc)`  
  Computes the cumulative intensity evaluated at the truncation horizon.

- `log_like_PIR(param, time, tau, time_trunc)`  
  Computes the log-likelihood for multiple independent systems.

---

## 2. Stan model definitions

File: `pir_stan_models.R`

This file defines the Stan implementations used for Bayesian inference.

### Included models

- `model.stan.ind`  
  PIR model with independent priors:
  - \( \beta \sim \text{Gamma} \)
  - \( \mu_T \sim \text{Gamma} \)
  - \( \theta \sim \text{Uniform}(-1,1) \)

- `model.stan.dep`  
  PIR model with correlated priors using an LKJ structure.

### Key features

- Custom likelihood implemented via `log_like_PIR`
- Explicit handling of change points \( \tau \)
- Parameter transformations for numerical stability

---

## 3. Frequentist estimation

File: `pir_frequentist_fit.R`

### Function

- `frequentist_estimation(time, time_trunc, tau)`

Computes maximum likelihood estimates using numerical optimization (`optimx`).

### Characteristics

- Uses L-BFGS-B method
- Enforces parameter constraints:
  - \( \beta > 0 \)
  - \( \mu_T > 0 \)
  - \( \theta \in [-1,1] \)

---

## 4. Bayesian estimation (independent priors)

File: `pir_bayes_independent_fit.R`

### Function

- `bayesian_estimation_ind(...)`

Performs Bayesian inference using:

- independent priors,
- Stan and Hamiltonian Monte Carlo (HMC).

### Output

- posterior medians (point estimates),
- HPD intervals,
- coverage indicators (for simulation),
- full `stanfit` object.

---

## 5. Bayesian estimation (dependent priors)

File: `pir_bayes_dependent_fit.R`

### Function

- `bayesian_estimation_dep(...)`

Performs Bayesian inference using:

- correlated priors via LKJ,
- latent Gaussian structure,
- Cholesky decomposition.

---

## Model description

The PIR model generalizes imperfect repair processes by introducing:

- a repair-effect parameter \( \theta \),
- a set of change points \( \tau \),
- a flexible structure that modifies both:
  - the effective age of the system,
  - and the scaling of the intensity.

### Interpretation of \( \theta \)

- \( \theta = 0 \): minimal repair (ABAO)  
- \( 0 < \theta < 1 \): imperfect repair  
- \( \theta = 1 \): perfect repair (AGAN-like behavior)  
- \( \theta < 0 \): harmful repair  

---

## Notes

- The PIR model depends explicitly on the vector of change points \( \tau \), which must be provided in all estimation procedures.
- Numerical stability is handled via transformations \( h(\theta, \beta) \) and \( g(\theta) \).
- Simulation is performed using inversion of the cumulative intensity.
- Bayesian estimation relies on Stan and HMC for efficient sampling.
- The module is designed to be fully reusable across simulation and application workflows.

---

## Purpose of this module

This module provides the **core computational layer** of the PIR model, including:

- mathematical formulation,
- simulation tools,
- likelihood evaluation,
- inference methods.

It is designed to ensure **consistency, modularity, and reproducibility** within the overall framework of repairable systems modeling.