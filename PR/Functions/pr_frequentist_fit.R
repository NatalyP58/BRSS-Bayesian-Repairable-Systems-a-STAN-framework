# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: PR frequentist fit
# ------------------------------------------------------------
# This file provides the frequentist estimator for the Perfect
# Repair (PR) model using numerical maximization of the
# log-likelihood function.
# ------------------------------------------------------------

# Required packages
library(optimx)
library(pracma)

# ------------------------------------------------------------
# Frequentist estimation for the PR model
# ------------------------------------------------------------
frequentist_estimation <- function(time, time_trunc) {
  
  # Numerical maximization of the PR log-likelihood
  result <- try(
    optimx::optimr(
      par = c(1.1, 1.1),                  # initial values for (beta, mu_T)
      fn = log_like_PR,               # PR log-likelihood to be maximized
      method = "L-BFGS-B",            # bound-constrained optimization method
      time_trunc = time_trunc,        # truncation horizon T
      times = time,                   # observed failure times
      control = list(maximize = TRUE),# direct maximization of the log-likelihood
      lower = c(1, 0),                # enforce beta >= 1 and mu_T > 0
      upper = c(Inf, Inf)             # no finite upper bounds imposed
    )$par,
    silent = TRUE
  )
  
  # Return failure if optimization does not converge properly
  if (inherits(result, "try-error") || any(is.na(result))) {
    return(list(success = FALSE, estimate = NULL))
  }
  
  # Return parameter estimates when optimization is successful
  return(list(success = TRUE, estimate = as.numeric(result)))
}