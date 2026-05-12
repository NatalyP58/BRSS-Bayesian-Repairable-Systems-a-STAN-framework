# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: ARA model functions
# ------------------------------------------------------------
# This file contains the core functions used for the
# Arithmetic Reduction of Age (ARA) model:
#   1. simulation of failure times,
#   2. conditional intensity function,
#   3. cumulative intensity function, and
#   4. log-likelihood evaluation for multiple independent systems.
# ------------------------------------------------------------

# Required package
#library(pracma)

# ------------------------------------------------------------
# Failure-time generation under the ARA model
# ------------------------------------------------------------
gen_failures_ARA <- function(k, param, time_trunc) {
  
  # Interval index of time t within the sequence of failure times
  k_t <- function(t, times, time_trunc) {
    findInterval(t, c(0, times, time_trunc), left.open = TRUE)
  }
  
  # Cumulative intensity up to time t under the ARA model
  lambda_ARA_acum <- function(t, param, times, time_trunc) {
    if (t <= 0) return(0)
    
    beta  <- param[1]
    mu_T  <- param[2]
    theta <- param[3]
    
    kt      <- k_t(t, times, time_trunc)
    times0  <- c(0, times)[1:(kt - 1)]
    times1  <- c(times, time_trunc)[1:(kt - 1)]
    times_k <- ifelse(kt - 1 > 0, c(times, time_trunc)[kt - 1], 0)
    
    if (kt - 1 > 0) {
      lambda1 <- mu_T * (
        (times1 - theta * times0)^beta -
          (times0 - theta * times0)^beta
      ) / time_trunc
      lambda1[is.na(lambda1)] <- 0
    } else {
      lambda1 <- 0
    }
    
    lambda2 <- mu_T * (
      (t - theta * times_k)^beta -
        (times_k - theta * times_k)^beta
    ) / time_trunc
    lambda2[is.na(lambda2)] <- 0
    
    return(sum(lambda1) + lambda2)
  }
  
  lambda_ARA_acum <- Vectorize(lambda_ARA_acum, vectorize.args = "t")
  
  # Objective function used to numerically invert the cumulative intensity
  # It enforces the inversion equation:
  #   Lambda(t + x) - Lambda(t) = -log(1 - u)
  # where u ~ Uniform(0,1).
  objetivo_incremento <- function(x, param, t, times, time_trunc, u) {
    abs(
      lambda_ARA_acum(t + x, param, times, time_trunc) -
        lambda_ARA_acum(t, param, times, time_trunc) +
        log(1 - u)
    )
  }
  
  # Simulate one failure trajectory up to the truncation horizon
  simular_una_trayectoria <- function(param, tau, time_trunc) {
    t_actual <- 0
    times    <- numeric(0)
    
    repeat {
      u        <- runif(1)
      objetivo <- -log(1 - u)
      
      # Remaining cumulative intensity available from the current time to T
      max_Lambda <- lambda_ARA_acum(time_trunc, param, times, time_trunc) -
        lambda_ARA_acum(t_actual, param, times, time_trunc)
      
      # Stop if the simulated jump exceeds the remaining cumulative intensity
      if (objetivo > max_Lambda) break
      
      # Numerically invert the cumulative intensity to obtain the next event time
      resultado <- optimize(
        objetivo_incremento,
        interval   = c(0, time_trunc - t_actual),
        param      = param,
        t          = t_actual,
        times      = times,
        time_trunc = time_trunc,
        u          = u
      )
      
      incremento  <- resultado$minimum
      t_siguiente <- t_actual + incremento
      
      # Stop if the next event falls beyond the truncation horizon
      if (t_siguiente > time_trunc) break
      
      # Store the new failure time and continue the trajectory
      times    <- c(times, t_siguiente)
      t_actual <- t_siguiente
    }
    
    return(times)
  }
  
  # Generate k independent trajectories
  replicate(k, simular_una_trayectoria(param, times, time_trunc), simplify = FALSE)
}

# ------------------------------------------------------------
# Auxiliary function
# ------------------------------------------------------------

# Most recent failure time before t
indicator <- function(vector, t) {
  vals <- vector[vector < t]
  if (length(vals) == 0) return(0)
  max(vals)
}

# ------------------------------------------------------------
# ARA intensity and cumulative intensity
# ------------------------------------------------------------

# ARA conditional intensity function
lambda_ARA <- function(t, param, times, time_trunc) {
  if (t <= 0) return(0)
  
  beta  <- param[1]
  mu_T  <- param[2]
  theta <- param[3]
  
  timess <- indicator(times, t)  # timess = t_{i-1}
  
  mu_T * beta / time_trunc * (t - theta * timess)^(beta - 1)
}

lambda_ARA <- Vectorize(lambda_ARA, vectorize.args = "t")

# ARA cumulative intensity at the truncation horizon
Lambda_ARA <- function(param, times, time_trunc) {
  if (length(times) == 0 || max(times) <= 0) return(0)
  
  beta  <- param[1]
  mu_T  <- param[2]
  theta <- param[3]
  
  times0 <- c(0, times[-length(times)])
  times1 <- times
  timesp <- times[length(times)]
  
  sum(
    mu_T / time_trunc * (
      (times1 - theta * times0)^beta -
        (times0 - theta * times0)^beta
    )
  ) +
    mu_T * (
      (time_trunc - theta * timesp)^beta -
        (timesp - theta * timesp)^beta
    ) / time_trunc
}

# ------------------------------------------------------------
# ARA log-likelihood
# ------------------------------------------------------------

log_like_ARA <- function(param, times, time_trunc) {
  
  log_like_single <- function(times_i) {
    if (length(times_i) == 0) return(0)
    
    log_lambda_sum <- sum(log(lambda_ARA(times_i, param, times_i, time_trunc)))
    lambda_total   <- Lambda_ARA(param, times_i, time_trunc)
    
    log_lambda_sum - lambda_total
  }
  
  # Total log-likelihood for multiple independent systems
  sum(sapply(times, log_like_single))
}