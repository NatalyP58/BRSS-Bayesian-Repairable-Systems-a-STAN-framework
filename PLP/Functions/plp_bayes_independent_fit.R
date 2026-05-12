# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: PLP Bayesian estimation with independent priors
# ------------------------------------------------------------
# This file implements Bayesian estimation for the Partial
# Imperfect Repair (PLP) model using Stan under an
# independent-prior specification.
# ------------------------------------------------------------

# Required packages
library(rstan)
library(coda)

# Stan configuration
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

# ------------------------------------------------------------
# Bayesian estimation for the PLP model with independent priors
# ------------------------------------------------------------
bayesian_estimation_ind <- function(time, time_trunc, beta_r = NULL, mu_T_r = NULL,
                                    theta_r = NULL, tau, compiled_model = NULL,
                                    simulation = TRUE, seed = NULL) {
  
  # If no compiled model is provided, use the global compiled object
  if (is.null(compiled_model)) {
    if (!exists("model.compiled.ind")) {
      stop("model.compiled.ind not found in the global environment")
    }
    compiled_model <- model.compiled.ind
  }
  
  # If simulation = TRUE, the true parameter values must be provided
  if (simulation && (is.null(beta_r) || is.null(mu_T_r) || is.null(theta_r))) {
    stop("If simulation = TRUE, you must provide the true values beta_r, mu_T_r, and theta_r.")
  }
  
  # Flatten failure times and build index vectors for Stan
  time_flat <- unlist(time)
  len_times <- sapply(time, length)
  start_idx <- cumsum(c(1, len_times[-length(len_times)]))
  end_idx   <- cumsum(len_times)
  init_fun <- function() {list(beta = 1.1,mu_T = 1.1) }
  # Data passed to Stan
  stan_data <- list(
    k          = length(time),                             # number of independent systems
    len_times  = array(len_times, dim = length(len_times)),# number of observed failures in each system
    time_flat  = time_flat,                                # stacked vector of all failure times
    start_idx  = array(start_idx, dim = length(start_idx)),# starting index of each system in time_flat
    end_idx    = array(end_idx, dim = length(end_idx)),    # ending index of each system in time_flat
    time_trunc = time_trunc,                               # truncation horizon T
    beta0      = c(0.16, 0.04),                            # Gamma(shape, rate) hyperparameters for beta
    mu0        = c(2, 0.05)                                # Gamma(shape, rate) hyperparameters for mu_T
  )
  
  # Posterior sampling under the independent-prior specification
  fit <- sampling(
    compiled_model,                                           # compiled Stan model
    data    = stan_data,                                      # data list passed to Stan
    chains  = 1,                                              # single MCMC chain
    iter    = 3000,                                           # total number of iterations
    warmup  = 500,                                            # iterations discarded for adaptation
    seed    = if (is.null(seed)) 123L else as.integer(seed),  # fixed seed for reproducibility                                  
    init    = init_fun,                                       # model-specific initial values close to a stable region of the parameter space
    control = list(adapt_delta = 0.95),                       # more conservative HMC adaptation
    refresh = 0                                               # suppress sampling progress output
  )
  
  # Geweke diagnostic: checks convergence by comparing early vs late segments of the chain
  # |z| > 1.96 suggests lack of convergence at the 5% significance level
  zscores <- abs(geweke.diag(As.mcmc.list(fit))[[1]]$z[1:2])
  if (all(!is.na(zscores)) && any(zscores > 1.96)) {
    return(list(success = FALSE, estimate = NULL, coverage = c(0, 0), fit = fit))
  }
  
  # Extract posterior samples
  samples   <- rstan::extract(fit)
  beta_est  <- samples$beta
  mu_T_est  <- samples$mu_T
  theta_est <- samples$theta
  
  # Posterior point estimates
  estimate <- c(median(beta_est), median(mu_T_est))
  
  if (simulation) {
    # 95% HPD intervals used for empirical coverage
    beta_hpd  <- HPDinterval(as.mcmc(matrix(beta_est, ncol = 1)), prob = 0.95)
    mu_T_hpd  <- HPDinterval(as.mcmc(matrix(mu_T_est, ncol = 1)), prob = 0.95)
    
    coverage <- c(
      as.integer(beta_hpd[1] <= beta_r & beta_r <= beta_hpd[2]),
      as.integer(mu_T_hpd[1] <= mu_T_r & mu_T_r <= mu_T_hpd[2]),
      
    )
  } else {
    coverage <- NULL
  }
  
  return(list(success = TRUE, estimate = estimate, coverage = coverage, fit = fit))
}