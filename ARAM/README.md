# Repairable Systems Models – ARAM Implementation

This repository contains the implementation of the **Arithmetic Reduction of Age Modified (ARAM)** model developed in the manuscript:

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

ARAM/  
├── Functions/  
├── Simulation/  
├── Application/  
└── README.md  

---

## 1. Functions

Functions/

This folder contains the **core implementation of the ARAM model**.

Includes:

- model definitions (intensity, cumulative intensity, likelihood),
- virtual-age transformation functions,
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

This module applies the ARAM model to **real-world data** (sugarcane harvester dataset).

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

The ARAM model assumes:

- a repair mechanism based on **modified virtual age**,
- the repair effect acts through a transformation function \( h(\theta) \),
- both beneficial and detrimental repair effects are allowed,
- the failure intensity evolves according to a modified age structure.

The conditional intensity function is given by:

\[
\lambda_{ARAM}(t \mid H_{t^-})
=
\sum_{i=1}^{n+1}
\frac{\mu_T}{T}
\beta
\left(
t -
\operatorname{sign}(\beta-1)
\, h(\theta)\,
t_{i-1}
\right)^{\beta-1}
\mathbb{I}(t_{i-1}<t\le t_i),
\]

where:

- \( \beta > 0 \) controls the trend behavior,
- \( \mu_T > 0 \) is a scale parameter,
- \( \theta \in [-1,1] \) governs the repair effect,
- \( h(\theta) \) modifies the system’s effective age after repair.

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

- The ARAM model extends the classical ARA framework by allowing greater flexibility in the repair effect.

- Unlike classical ARA models restricted to beneficial repairs, ARAM allows:
  - beneficial repair,
  - minimal repair,
  - harmful repair.

- The repair effect is governed by \( \theta \):
  - \( \theta = 0 \): minimal repair (ABAO),
  - \( 0 < \theta < 1 \): imperfect beneficial repair,
  - \( \theta = 1 \): perfect repair (AGAN),
  - \( \theta < 0 \): harmful repair.

- The transformation function \( h(\theta) \) ensures model stability while preserving the interpretation of the repair effect.

- The off-season adjustment is critical in the application dataset to avoid bias due to long inactivity periods.

- Numerical stability should be monitored when:
  - \( \theta \) is close to its boundaries,
  - transformed virtual ages approach zero.

- The implementation is designed for **reproducibility and methodological clarity**, consistent with the manuscript.

---

## Purpose of this repository

This repository provides a **complete and reproducible implementation** of the ARAM model, including:

- theoretical formulation,
- simulation-based validation,
- real-data application.

It is intended to support both:

- **methodological research**, and  
- **practical implementation in reliability analysis**.