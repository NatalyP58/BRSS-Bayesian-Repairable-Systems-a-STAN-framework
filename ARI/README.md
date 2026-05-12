# Repairable Systems Models – ARI Implementation

This repository contains the implementation of the **Arithmetic Reduction of Intensity (ARI)** model developed in the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The goal of this project is to provide a **reproducible and modular framework** for:

- simulation of repairable systems,
- frequentist inference,
- Bayesian inference using Stan,
- application to real-world data.

---

## Repository structure

The project is organized into three main modules:

ARI/  
├── Functions/  
├── Simulation/  
├── Application/  
└── README.md  

---

## 1. Functions

Functions/

This folder contains the **core implementation of the ARI model**.

Includes:

- model definitions (intensity, cumulative intensity, likelihood),
- data simulation functions,
- frequentist estimation,
- Bayesian estimation (independent and dependent priors),
- Stan model specifications.

This is the **core mathematical layer** of the project.

---

## 2. Simulation

Simulation/

This module provides a **controlled experimental framework**.

It allows you to:

- generate synthetic data under known parameters,
- compare estimation methods:
  - frequentist,
  - Bayesian (independent priors),
  - Bayesian (dependent priors),
- evaluate performance using:
  - Bias,
  - Mean Squared Error (MSE),
  - empirical coverage.

This module is used to **validate the methodology**.

---

## 3. Application

Application/

This module applies the ARI model to **real-world data** (sugarcane harvester dataset).

It includes:

- data preprocessing,
- failure-time construction,
- off-season adjustment,
- model estimation,
- model comparison using:
  - AIC,
  - BIC,
  - DIC.

This demonstrates the **practical applicability** of the model.

---

## Workflow overview

A typical workflow is:

1. Define model and estimation functions  
   → Functions/

2. Validate methods through simulation  
   → Simulation/

3. Apply the model to real data  
   → Application/

---

## Statistical model

The ARI model assumes:

- **imperfect repair**, where maintenance actions directly modify the failure intensity,
- the failure intensity evolves over time according to a baseline power law,
- the repair effect depends on the most recent failure time.

The conditional intensity function is given by:

\[
\lambda_{ARI}(t) = \mu_T \frac{\beta}{T} \left( t^{\beta - 1} - \rho \, t_{i-1}^{\beta - 1} \right),
\]

where:

- \( t_{i-1} \) is the most recent failure time prior to \( t \),
- \( \beta > 0 \) controls the trend of failures,
- \( \mu_T > 0 \) is a scale parameter,
- \( \rho \in [-1,1] \) is the repair-effect parameter,
- \( T \) is the truncation horizon.

---

## Estimation methods implemented

The framework includes:

### Frequentist approach
- Maximum likelihood estimation using `optimx`

### Bayesian approaches
- Independent priors  
- Dependent priors (LKJ correlation structure)

Both Bayesian models are implemented in **Stan** using Hamiltonian Monte Carlo (HMC).

---

## Required packages

- rstan
- coda
- optimx
- parallel
- tidyverse
- dplyr

---

## Reproducibility

- Random seeds are controlled in both simulation and application modules
- Parallel execution ensures efficient computation
- Results are stored in `.rds` files for reproducibility

---

## Notes

- The ARI model modifies the **failure intensity directly**, unlike ARA which modifies the system’s virtual age.
- The parameter \( \rho \) governs the repair effect:
  - \( \rho = 0 \): minimal repair (ABAO),
  - \( 0 < \rho < 1 \): imperfect repair,
  - \( \rho = 1 \): optimal repair (intensity reset after each failure),
  - \( \rho < 0 \): harmful repair (increases failure intensity).
- The model allows a flexible representation of repair effects through direct intensity adjustment.
- The off-season adjustment is critical in the application dataset to avoid bias due to long inactivity periods.
- Numerical stability of the intensity should be monitored, especially for extreme values of \( \rho \).
- The implementation is designed for **reproducibility and methodological clarity**, consistent with the manuscript.

---

## Purpose of this repository

This repository provides a **complete and reproducible implementation** of the ARI model, including:

- theoretical formulation,
- simulation-based validation,
- real-data application.

It is intended to support both:

- **methodological research**, and  
- **practical implementation in reliability analysis**.