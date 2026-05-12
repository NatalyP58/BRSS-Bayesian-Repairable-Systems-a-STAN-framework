# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: PR model functions
# ------------------------------------------------------------
# This file contains the core functions used for the Perfect
# Repair (PR) model:
#   1. simulation of failure times,
#   2. conditional intensity function,
#   3. cumulative intensity function, and
#   4. log-likelihood evaluation for multiple independent systems.
# ------------------------------------------------------------

# Required package
#library(pracma)

# ------------------------------------------------------------
# Failure-time generation under the PR model
# ------------------------------------------------------------
gen_failures_PR <- function(k, param, time_trunc) {
  
  # Interval index of time t within the sequence of failure times
  k_t <- function(t, times, time_trunc) {
    findInterval(t, c(0, times, time_trunc), left.open = TRUE)
  }
  
  # Cumulative intensity up to time t under the PR model
  lambda_PR_acum <- function(t, param, times, time_trunc) {
    if (t <= 0) return(0)
    
    beta <- param[1]
    mu_T <- param[2]
    
    kt <- k_t(t, times, time_trunc)
    times0 <- c(0, times)[1:(kt - 1)]
    times1 <- c(times, time_trunc)[1:(kt - 1)]
    times_k <- ifelse(kt - 1 > 0, c(times, time_trunc)[kt - 1], 0)
    
    if (kt - 1 > 0) {
      lambda1 <- mu_T * ((times1 - times0)^beta) / time_trunc
      lambda1[is.na(lambda1)] <- 0
    } else {
      lambda1 <- 0
    }
    
    lambda2 <- mu_T * ((t - times_k)^beta) / time_trunc
    lambda2[is.na(lambda2)] <- 0
    
    return(sum(lambda1) + lambda2)
  }
  
  lambda_PR_acum <- Vectorize(lambda_PR_acum, vectorize.args = "t")
  
  # Objective function used to numerically invert the cumulative intensity
  # It enforces the inversion equation:
  #   Lambda(t + x) - Lambda(t) = -log(1 - u)
  # where u ~ Uniform(0,1).
  # The optimizer searches for the increment x such that this condition holds.
  optimize_increment <- function(x, param, t, times, time_trunc, u) {
    abs(
      lambda_PR_acum(t + x, param, times, time_trunc) -
        lambda_PR_acum(t, param, times, time_trunc) +
        log(1 - u)
    )
  }
  # Simulate one failure trajectory up to the truncation horizon
  simulate_one_trajectory <- function(param, time_trunc) {
    t_current <- 0
    times <- numeric(0)
    
    repeat {
      u <- runif(1)
      target <- -log(1 - u)  # exponential waiting-time transform
      
      # Remaining cumulative intensity available from the current time to T
      max_Lambda <- lambda_PR_acum(time_trunc, param, times, time_trunc) -
        lambda_PR_acum(t_current, param, times, time_trunc)
      
      # Stop if the simulated jump exceeds the remaining cumulative intensity
      if (target > max_Lambda) break
      
      # Numerically invert the cumulative intensity to obtain the next event time
      result <- optimize(
        optimize_increment,
        interval   = c(0, time_trunc - t_current),
        param      = param,
        t          = t_current,
        times      = times,
        time_trunc = time_trunc,
        u          = u
      )
      
      increment <- result$minimum
      t_next <- t_current + increment
      
      # Stop if the next event falls beyond the truncation horizon
      if (t_next > time_trunc) break
      
      # Store the new failure time and continue the trajectory
      times <- c(times, t_next)
      t_current <- t_next
    }
    
    return(times)
  }
  
  # Generate k independent trajectories
  replicate(k, simulate_one_trajectory(param, time_trunc), simplify = FALSE)
}

# ------------------------------------------------------------
# Auxiliary function
# ------------------------------------------------------------

# Most recent failure time before t
indicator_PR <- function(vector, t) {
  vals <- vector[vector < t]
  if (length(vals) == 0) return(0)
  max(vals)
}

# ------------------------------------------------------------
# PR intensity and cumulative intensity
# ------------------------------------------------------------

# PR conditional intensity function
lambda_PR <- function(t, param, times, time_trunc) {
  if (t <= 0) return(0)
  
  beta <- param[1]
  mu_T <- param[2]
  
  previous_time <- indicator_PR(times, t)
  
  mu_T * beta / time_trunc * (t - previous_time)^(beta - 1)
}

lambda_PR <- Vectorize(lambda_PR, vectorize.args = "t")

# PR cumulative intensity at the truncation horizon
Lambda_PR <- function(param, times, time_trunc) {
  # Here times = {t_1, ..., t_n}, with t_0 = 0 and t_{n+1} = T
  if (length(times) == 0 || max(times) <= 0) return(0)
  
  beta <- param[1]
  mu_T <- param[2]
  
  times0 <- c(0, times[-length(times)])
  times1 <- times
  times_last <- times[length(times)]
  
  sum(mu_T / time_trunc * ((times1 - times0)^beta)) +
    mu_T * ((time_trunc - times_last)^beta) / time_trunc
}

# ------------------------------------------------------------
# PR log-likelihood
# ------------------------------------------------------------

log_like_PR <- function(param, times, time_trunc) {
  
  log_like_single <- function(times_i) {
    # If no failures are observed, the contribution is zero
    if (length(times_i) == 0) return(0)
    
    log_lambda_sum <- sum(log(lambda_PR(times_i, param, times_i, time_trunc)))
    lambda_total <- Lambda_PR(param, times_i, time_trunc)
    
    log_lambda_sum - lambda_total
  }
  
  # Total log-likelihood for multiple independent systems
  sum(sapply(times, log_like_single))
}