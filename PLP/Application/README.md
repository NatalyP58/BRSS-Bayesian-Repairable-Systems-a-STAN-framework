## Functions used in this application

This script relies on several functions defined in external files.  
Below is a concise description of each function to allow standalone use of this module.

### Core ARAM functions

- `lambda_ARAM(t, param, times, time_trunc)`  
  Computes the conditional intensity function of the ARAM model.

- `Lambda_ARAM(param, times, time_trunc)`  
  Computes the cumulative intensity evaluated at the truncation time.

- `log_like_ARAM(param, time, time_trunc)`  
  Computes the log-likelihood of the ARAM model for multiple systems.

### Data preprocessing

- `adjust_offseason(x, start_obs, end_obs)`  
  Removes non-operational periods from the timeline by:
  - eliminating failures occurring during off-season,
  - compressing the time scale to reflect effective operating time.

### Estimation methods

- `frequentist_estimation(time, time_trunc)`  
  Computes maximum likelihood estimates using numerical optimization.

- `bayesian_estimation_ind(time, time_trunc, ...)`  
  Performs Bayesian inference using independent priors via Stan.

- `bayesian_estimation_dep(time, time_trunc, ...)`  
  Performs Bayesian inference using correlated priors via Stan (LKJ structure).

### Model comparison

- `calc_AIC_BIC(results, times, time_trunc)`  
  Computes AIC and BIC from the ARAM log-likelihood evaluated at point estimates.

- `calc_DIC(fit_stan, times, time_trunc)`  
  Computes the Deviance Information Criterion using posterior samples.

## Source files

All functions are loaded via:

```r
source_files <- c(
  file.path("functions", "aram_model_functions.R"),
  file.path("functions", "aram_stan_models.R"),
  file.path("functions", "aram_frequentist_fit.R"),
  file.path("functions", "aram_bayes_independent_fit.R"),
  file.path("functions", "aram_bayes_dependent_fit.R")
)
```

## Notes

- The ARAM model extends imperfect repair models by modifying the system’s effective age through a transformation, rather than directly modifying intensity (as in ARI).
- The parameter \( \theta \in [-1,1] \) governs the repair effect:
  - \( \theta = 0 \): minimal repair (ABAO),
  - \( 0 < \theta < 1 \): imperfect repair,
  - \( \theta = 1 \): perfect repair (AGAN),
  - \( \theta < 0 \): harmful repair.
- The application is based on multiple subsystems (failure modes) observed over a common time horizon.
- Off-season periods are explicitly removed to ensure that inference is based on effective operating time.
- Bayesian estimation is implemented using Stan and Hamiltonian Monte Carlo (HMC).
- Model comparison metrics (AIC, BIC, DIC) are used to assess relative model fit.
- Numerical stability should be monitored when \( \theta \) is close to its boundary values.