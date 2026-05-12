# Repairable Systems Models – ARA Implementation

This repository contains the implementation of the **Arithmetic Reduction of Age (ARA)** model developed in the manuscript:

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

ARA/  
├── Functions/  
├── Simulation/  
├── Application/  
└── README.md  

---

## 1. Functions

Functions/

This folder contains the **core implementation of the ARA model**.

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

This module applies the ARA model to **real-world data** (sugarcane harvester dataset).

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

The ARA model assumes:

- a **virtual-age repair mechanism**, meaning that after each failure the effective age of the system is modified according to the repair effect,
- the failure intensity depends on a **modified effective age** of the system,
- the repair effect is controlled by a parameter that adjusts how much past failures influence future behavior.

The intensity function is given by:

\[
\lambda_{ARA}(t) = \mu_T \frac{\beta}{T} \left(t - \theta t_{i-1}\right)^{\beta - 1},
\]

with \(t_0 = 0\) and \(t_{n+1} = T\), where:

- \( \beta \geq 1 \) controls the trend of failures,
- \( \mu_T > 0 \) is a scale parameter,
- \( \theta \in [-1,1] \) is the repair-effect parameter,
- \( t_{i-1} \) is the most recent failure time before \(t\),
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

- The ARA model generalizes classical repair models through the parameter \( \theta \):
  - \( \theta = 0 \): minimal repair (ABAO),
  - \( 0 < \theta < 1 \): efficient imperfect repair,
  - \( \theta = 1 \): perfect repair (AGAN),
  - \( \theta < 0 \): harmful repair.
- The parameterization enforces \( \beta \geq 1 \), ensuring a non-decreasing failure intensity between failures.
- The off-season adjustment is critical in the application dataset to avoid bias due to long inactivity periods.
- The implementation is designed for **reproducibility and methodological clarity**, consistent with the manuscript.

---

## Purpose of this repository

This repository provides a **complete and reproducible implementation** of the ARA model, including:

- theoretical formulation,
- simulation-based validation,
- real-data application.

It is intended to support both:

- **methodological research**, and  
- **practical implementation in reliability analysis**.