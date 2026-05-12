# Simulation Module for the Arithmetic Reduction of Age (ARA)

This folder contains the simulation workflow for the **Arithmetic Reduction of Age (ARA)** model developed in the manuscript:

**A Bayesian Approach to Repairable Systems Models**  
**Author:** Nataly Martinez Riascos  
**Date:** April 15, 2026

The simulation framework is designed to compare three estimation strategies:

- Frequentist estimation
- Bayesian estimation with independent priors
- Bayesian estimation with dependent priors

---

## Contents of this folder

This folder contains two main scripts:

- `ara_simulation_cluster.R`
- `ara_simulation_results.R`

Together, these files define the full simulation pipeline:
1. run simulation experiments,
2. estimate model parameters,
3. store outputs,
4. summarize performance using bias, mean squared error, and empirical coverage.

---

## 1. `ara_simulation_cluster.R`

This script runs the simulation study using parallel cluster execution.

### Purpose

The script repeatedly:
- generates failure-time data under the ARA model,
- fits the ARA model using three estimation methods,
- stores estimates, fitted objects, and coverage indicators,
- saves all simulation outputs in an `.rds` file.

### Main function

#### `run_simulation_parallel_cluster(m, k, param, time_trunc, SEED = NULL)`

Runs the full simulation study in parallel.

### Arguments

- `m`  
  Number of simulation replications.

- `k`  
  Number of independent systems generated in each replication.

- `param`  
  True parameter vector used for data generation.  
  For the ARA model this is:
  - `beta`
  - `mu_T`
  - `theta`

- `time_trunc`  
  Common truncation horizon.

- `SEED`  
  Optional random seed for reproducibility.

---

## Internal workflow

For each replication, the script performs the following steps:

1. **Generate data**
   - `gen_failures_ARA()`

2. **Frequentist estimation**
   - `frequentist_estimation()`

3. **Bayesian estimation with independent priors**
   - `bayesian_estimation_ind()`

4. **Bayesian estimation with dependent priors**
   - `bayesian_estimation_dep()`

---

## Stan-related objects used

The script sources the model and estimation files and uses the following Stan objects:

- `model.stan.ind`
- `model.stan.dep`
- `model.compiled.ind`
- `model.compiled.dep`

---

## Returned object

The function returns a list of length `m`, where each element corresponds to one simulation replication.

Each replication contains:

- `freq_est`  
  Frequentist parameter estimates

- `bayes_est_ind`  
  Bayesian parameter estimates under independent priors

- `bayes_est_dep`  
  Bayesian parameter estimates under dependent priors

- `coverage_ind`  
  Coverage indicators under independent priors

- `coverage_dep`  
  Coverage indicators under dependent priors

- `fit_ind`  
  Stan fit object under independent priors

- `fit_dep`  
  Stan fit object under dependent priors

- `success_freq`  
  Logical flag indicating successful frequentist estimation

- `success_ind`  
  Logical flag indicating successful Bayesian estimation with independent priors

- `success_dep`  
  Logical flag indicating successful Bayesian estimation with dependent priors

- `seed_used`  
  Seed used for that replication

---

## Output file

Simulation outputs are saved as:

```r
results_k101_simulation01.rds