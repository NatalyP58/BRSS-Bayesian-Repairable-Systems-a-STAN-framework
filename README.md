# BRSS: Bayesian Repairable Systems – A Stan Framework

This repository contains a modular implementation of repairable-systems models developed for the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The repository provides R and Stan code for simulation studies, frequentist estimation, Bayesian inference, and real-data applications for several repairable-systems models.

---

## Models included

The repository includes implementations for:

- **PLP** – Power Law Process
- **PR** – Perfect Repair model
- **ARA** – Arithmetic Reduction of Age
- **ARI** – Arithmetic Reduction of Intensity
- **ARAM** – Arithmetic Reduction of Age with Memory
- **PIR** – Partial Imperfect Repair model

These models cover different assumptions about repair effects, including minimal repair, perfect repair, imperfect repair, harmful repair, and change-point-based repair dynamics.

---

## Repository structure

Each model folder follows the same modular organization:

```text
Model/
├── Functions/
├── Simulation/
├── Application/
└── README.md
```

---

## Repository modules

### Functions/

Contains the core implementation of each model, including:

- conditional intensity functions,
- cumulative intensity functions,
- log-likelihood functions,
- data generation routines,
- frequentist estimation,
- Bayesian estimation with independent priors,
- Bayesian estimation with dependent priors using LKJ correlation structures,
- Stan model specifications.

### Simulation/

Contains scripts for controlled simulation studies. These scripts:

- generate synthetic repairable-systems data,
- estimate model parameters,
- compare frequentist and Bayesian methods,
- compute bias, mean squared error, and empirical coverage,
- store reproducible outputs in `.rds` files.

### Application/

Contains real-data application scripts using the sugarcane harvester dataset. These scripts include:

- data preprocessing,
- construction of failure-time vectors,
- off-season adjustment,
- model fitting,
- model comparison using AIC, BIC, and DIC.

---

## Bayesian inference

Bayesian estimation is implemented in **Stan** using Hamiltonian Monte Carlo. For each model, two Bayesian specifications are considered:

- independent priors,
- dependent priors based on an LKJ correlation structure.

The dependent-prior formulation allows correlation among transformed model parameters through a Cholesky-factor representation.

---

## Required R packages

Main dependencies include:

- `rstan`
- `coda`
- `optimx`
- `parallel`
- `tidyverse`
- `dplyr`

Additional packages may be used in specific application scripts.

---

## Additional system requirement for Windows

For Windows users, **RTools45** must be installed to compile Stan models with `rstan`.

You can verify the installation in R with:

```r
Sys.which("make")
pkgbuild::check_build_tools(debug = TRUE)
```

---

## Purpose

This repository is intended to provide a reproducible computational framework for repairable-systems modeling, supporting:

- methodological development,
- simulation-based validation,
- Bayesian and frequentist comparison,
- practical reliability analysis using real industrial data.

---

## Notes

Each model folder contains its own `README.md` file with additional details about the corresponding implementation, scripts, workflow, and model-specific assumptions.