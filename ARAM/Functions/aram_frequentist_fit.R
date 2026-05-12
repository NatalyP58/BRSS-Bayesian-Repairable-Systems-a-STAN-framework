# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: ARAM frequentist estimation
# ------------------------------------------------------------
# This file provides a frequentist estimator for the Arithmetic
# Reduction of Age with Memory (ARAM) model using numerical
# maximization of the log-likelihood function.
# ------------------------------------------------------------

# Required packages
library(optimx)
#library(pracma)

# ------------------------------------------------------------
# Frequentist estimation for the ARAM model
# ------------------------------------------------------------
frequentist_estimation <- function(time, time_trunc) {
  
  # Numerical maximization of the ARAM log-likelihood
  result <- try(
    optimx::optimr(
      par = c(2, 20, 0.8),            # initial values for (beta, mu_T, theta)
      fn = log_like_ARAM,             # ARAM log-likelihood to be maximized
      method = "L-BFGS-B",            # bound-constrained optimization method
      time_trunc = time_trunc,        # truncation horizon T
      time = time,                    # observed failure times
      control = list(maximize = TRUE),# direct maximization of the log-likelihood
      lower = c(0, 0, -1),            # lower bounds for (beta, mu_T, theta)
      upper = c(Inf, Inf, 1)          # upper bounds for (beta, mu_T, theta)
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

