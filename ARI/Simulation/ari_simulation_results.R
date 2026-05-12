# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: ARI simulation results summary
# ------------------------------------------------------------
# This script summarizes simulation results for the Arithmetic
# Reduction of Intensity (ARI) model. It computes bias, mean
# squared error (MSE), and empirical coverage for the
# frequentist estimator and the two Bayesian estimators:
#   1. independent priors,
#   2. dependent priors.
# ------------------------------------------------------------

# ------------------------------------------------------------
# Summary functions
# ------------------------------------------------------------

# Bias for frequentist, independent-prior, and dependent-prior estimates
Bias <- function(est1, est2, est3, real_par) {
  est1 <- sweep(est1, 2, real_par)
  est2 <- sweep(est2, 2, real_par)
  est3 <- sweep(est3, 2, real_par)
  
  return(rbind(
    colMeans(est1, na.rm = TRUE),
    colMeans(est2, na.rm = TRUE),
    colMeans(est3, na.rm = TRUE)
  ))
}

# Mean squared error for frequentist, independent-prior, and dependent-prior estimates
MSE <- function(est1, est2, est3, real_par) {
  est1 <- (sweep(est1, 2, real_par))^2
  est2 <- (sweep(est2, 2, real_par))^2
  est3 <- (sweep(est3, 2, real_par))^2
  
  return(rbind(
    colMeans(est1, na.rm = TRUE),
    colMeans(est2, na.rm = TRUE),
    colMeans(est3, na.rm = TRUE)
  ))
}

# Data frame with bias and MSE summaries
sim_df <- function(sim, param) {
  simul <- data.frame(
    Method = c("frequentist", "independent", "dependent"),
    cbind(
      Bias(sim$freq_est, sim$bayes_est_ind, sim$bayes_est_dep, param),
      MSE(sim$freq_est, sim$bayes_est_ind, sim$bayes_est_dep, param)
    )
  )
  
  colnames(simul) <- c(
    "Method",
    "Bias_beta", "Bias_mu_T", "Bias_rho",
    "MSE_beta", "MSE_mu_T", "MSE_rho"
  )
  
  return(simul)
}

# Data frame with empirical coverage summaries
sim_df_coverage <- function(sim) {
  simul <- data.frame(
    Method = c("independent", "dependent"),
    rbind(sim$coverage_ind_summary, sim$coverage_dep_summary)
  )
  
  colnames(simul) <- c("Method", "HPD_beta", "HPD_mu_T", "HPD_rho")
  return(simul)
}

# ------------------------------------------------------------
# Load simulation results
# ------------------------------------------------------------

results <- readRDS("results_k101_simulation01.rds")

# Frequentist estimates
freq_est <- do.call(rbind, lapply(results, function(res) res$freq_est))

# Bayesian estimates: independent priors
bayes_est_ind <- do.call(rbind, lapply(results, function(res) res$bayes_est_ind))

# Bayesian estimates: dependent priors
bayes_est_dep <- do.call(rbind, lapply(results, function(res) res$bayes_est_dep))

# Coverage: independent priors
coverage_ind <- do.call(rbind, lapply(results, function(res) res$coverage_ind))
coverage_ind_summary <- colMeans(coverage_ind, na.rm = TRUE)

# Coverage: dependent priors
coverage_dep <- do.call(rbind, lapply(results, function(res) res$coverage_dep))
coverage_dep_summary <- colMeans(coverage_dep, na.rm = TRUE)

# Optional: retain fitted Stan objects if needed for later inspection
fit_ind <- lapply(results, function(res) res$fit_ind)
fit_dep <- lapply(results, function(res) res$fit_dep)

# Store results in a single list
sim <- list(
  freq_est             = freq_est,
  bayes_est_ind        = bayes_est_ind,
  bayes_est_dep        = bayes_est_dep,
  coverage_ind         = coverage_ind,
  coverage_dep         = coverage_dep,
  coverage_ind_summary = coverage_ind_summary,
  coverage_dep_summary = coverage_dep_summary
)

# ------------------------------------------------------------
# Compute summaries
# ------------------------------------------------------------

# True parameter values used in the simulation study
param_true <- c(2, 30, 0.5)

# Bias and MSE summary
summary_df <- sim_df(sim, param_true)
summary_df

# Empirical coverage summary
coverage_df <- sim_df_coverage(sim)
coverage_df