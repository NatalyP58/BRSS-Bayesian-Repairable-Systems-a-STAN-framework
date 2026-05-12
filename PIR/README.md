# Repairable Systems Models – PIR Implementation

This repository contains the implementation of the **Partial Imperfect Repair (PIR)** model developed in the manuscript:

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

PIR/  
├── Functions/  
├── Simulation/  
├── Application/  
└── README.md  

---

## 1. Functions

Functions/

This folder contains the **core implementation of the PIR model**.

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

This module applies the PIR model to **real-world data** (sugarcane harvester dataset).

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

The PIR model assumes:

- a **piecewise structure defined by change points** \( \tau = (\tau_1, \ldots, \tau_m) \),
- the failure intensity evolves according to a baseline process (typically PLP),
- the repair effect is incorporated through a parameter \( \theta \),
- the observation window is partitioned into intervals determined by \( \tau \).

The conditional intensity function is given by:

\[
\lambda_{PIR}(t \mid H_{t^-}) = \frac{\mu_T \, \beta}{T + g(\theta)\,\tau_{N(t^-)}} 
\left( t - s(\beta)\,h(\theta)\,\tau_{N(t^-)} \right)^{\beta - 1},
\]

where:

- \( \beta > 0 \) controls the failure trend,
- \( \mu_T > 0 \) is a scale parameter,
- \( \theta \in [-1,1] \) governs the repair effect,
- \( \tau_{N(t^-)} \) is the most recent change point before time \( t \),
- \( h(\theta) \), \( g(\theta) \), and \( s(\beta) \) are transformation functions defining the repair mechanism.

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

- The PIR model generalizes several classical repair models through its flexible structure.
- The repair effect is governed by \( \theta \):
  - \( \theta = 0 \): minimal repair (ABAO),
  - \( 0 < \theta < 1 \): imperfect repair,
  - \( \theta = 1 \): perfect repair (AGAN),
  - \( \theta < 0 \): harmful repair.
- The change-point vector \( \tau \) introduces additional flexibility by structuring the failure process over time.
- The off-season adjustment is critical in the application dataset to avoid bias due to long inactivity periods.
- Numerical stability should be monitored when:
  - \( \theta \) is close to its boundaries,
  - change points are near the extremes of the observation window.
- The implementation is designed for **reproducibility and methodological clarity**, consistent with the manuscript.

---

## Purpose of this repository

This repository provides a **complete and reproducible implementation** of the PIR model, including:

- theoretical formulation,
- simulation-based validation,
- real-data application.

It is intended to support both:

- **methodological research**, and  
- **practical implementation in reliability analysis**.