# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: ARI application to sugarcane harvester data
# ------------------------------------------------------------
# This script applies the Arithmetic Reduction of Intensity
# (ARI) model to a selected vehicle from the sugarcane
# harvester dataset.
#
# The workflow includes:
#   1. data loading and preprocessing,
#   2. construction of failure-time vectors by subsystem,
#   3. adjustment for off-season periods,
#   4. frequentist and Bayesian estimation,
#   5. model comparison via AIC, BIC, and DIC.
# ------------------------------------------------------------

# Required packages
library(tidyverse)
library(dplyr)

# ------------------------------------------------------------
#  Data loading and preparation-----
# ------------------------------------------------------------

# Read the application dataset
data <- read.csv("application/sugarcane_data.csv", sep = ";")

# Select one vehicle for analysis
vehicle_selected <- "A"  # replace with the desired vehicle

# Keep only observations from the selected vehicle
data_vehicle <- data %>%
  filter(VEICULO == vehicle_selected)

# Convert the failure-start column to POSIXct format
data_vehicle$START_STOPPAGE <- as.POSIXct(
  data_vehicle$INICIO_PARADA,
  format = "%d/%m/%Y %H:%M"
)

# Define the observation window
start_obs <- as.POSIXct("2015-04-01 00:00", format = "%Y-%m-%d %H:%M")
end_obs   <- as.POSIXct("2017-03-31 23:59", format = "%Y-%m-%d %H:%M")

# Compute failure times in hours from the beginning of the observation window
data_vehicle <- data_vehicle %>%
  arrange(START_STOPPAGE) %>%
  mutate(
    time_failure_hours = as.numeric(
      difftime(START_STOPPAGE, start_obs, units = "hours")
    )
  )

# Identify the four most frequent failure modes for the selected vehicle
top4 <- data_vehicle %>%
  group_by(PROBLEMA) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>%
  slice(2:5) %>%
  pull(PROBLEMA)

# Build one failure-time vector per selected failure mode
times <- lapply(top4, function(prob) {
  failure_times <- data_vehicle %>%
    filter(PROBLEMA == prob) %>%
    arrange(time_failure_hours) %>%
    pull(time_failure_hours)
  
  # Remove the first event so that each vector is aligned from the observation origin
  failure_times <- failure_times[failure_times != min(failure_times)]
  return(failure_times)
})

# Name each vector by its corresponding failure mode
names(times) <- top4

# Display failure-time vectors
times

# ------------------------------------------------------------
# Stan and inference setup -----
# ------------------------------------------------------------

# Initial directory setup
current_dir <- getwd()

# Create a safe temporary directory for Stan
safe_temp_dir <- file.path(current_dir, "stan_temp")
if (!dir.exists(safe_temp_dir)) {
  dir.create(safe_temp_dir, recursive = TRUE)
}

# Set critical temporary environment variables
Sys.setenv(TMPDIR = safe_temp_dir)
Sys.setenv(TEMP = safe_temp_dir)
Sys.setenv(TMP = safe_temp_dir)
options(tmpdir = safe_temp_dir)

# Load main packages
library(parallel)
suppressPackageStartupMessages({
  library(rstan)
  library(coda)
})

# Stan configuration
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
options(timeout = 600)

# Load auxiliary scripts
source_files <- c(
  file.path("functions", "ari_model_functions.R"),
  file.path("functions", "ari_stan_models.R"),
  file.path("functions", "ari_frequentist_fit.R"),
  file.path("functions", "ari_bayes_dependent_fit.R"),
  file.path("functions", "ari_bayes_independent_fit.R")
)

for (file in source_files) {
  source(file.path(current_dir, file))
}

# Precompile the main Stan model with independent priors
model.compiled.ind <- tryCatch({
  rstan::stan_model(
    model_code = model.stan.ind,
    auto_write = FALSE,
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error compiling independent model on the main node:", e$message, "\n")
  NULL
})

# Precompile the main Stan model with dependent priors
model.compiled.dep <- tryCatch({
  rstan::stan_model(
    model_code = model.stan.dep,
    auto_write = FALSE,
    verbose = FALSE
  )
}, error = function(e) {
  cat("Error compiling dependent model on the main node:", e$message, "\n")
  NULL
})

# ------------------------------------------------------------
#  Main application function ------
# ------------------------------------------------------------

run_application_parallel_cluster <- function(time, time_trunc, SEED = NULL) {
  
  # Generate a reproducible seed if none is provided
  if (is.null(SEED)) {
    SEED <- as.integer(Sys.time()) %% .Machine$integer.max
    cat("No seed was provided. Generated seed:", SEED, "\n")
  } else {
    cat("Using fixed seed:", SEED, "\n")
  }
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  
  # Use fewer cores for greater stability
  ncores <- min(3, detectCores() - 1)
  cl <- makeCluster(
    ncores,
    type = "PSOCK",
    outfile = file.path(current_dir, "cluster_log1.txt")
  )
  
  # Reproducible RNG streams across workers
  parallel::clusterSetRNGStream(cl, iseed = SEED)
  
  # Export required objects to workers
  clusterExport(
    cl,
    varlist = c(
      "frequentist_estimation",
      "bayesian_estimation_ind",
      "bayesian_estimation_dep",
      "Lambda_ARI",
      "lambda_ARI",
      "log_like_ARI",
      "model.stan.ind",
      "model.stan.dep",
      "current_dir",
      "safe_temp_dir",
      "model.compiled.ind",
      "model.compiled.dep",
      "SEED"
    ),
    envir = environment()
  )
  
  # Worker configuration
  clusterEvalQ(cl, {
    Sys.setenv(TMPDIR = safe_temp_dir)
    Sys.setenv(TEMP = safe_temp_dir)
    Sys.setenv(TMP = safe_temp_dir)
    options(tmpdir = safe_temp_dir, mc.cores = 1, timeout = 600)
    rstan::rstan_options(auto_write = TRUE)
    Sys.setenv(STAN_NUM_THREADS = 1)
    
    suppressPackageStartupMessages({
      library(rstan)
      library(coda)
    })
    
    # Compile independent model if not already available
    if (!is.null(model.compiled.ind)) {
      assign("model.compiled.ind", model.compiled.ind, envir = .GlobalEnv)
    } else {
      assign(
        "model.compiled.ind",
        rstan::stan_model(model_code = model.stan.ind),
        envir = .GlobalEnv
      )
    }
    
    # Compile dependent model if not already available
    if (!is.null(model.compiled.dep)) {
      assign("model.compiled.dep", model.compiled.dep, envir = .GlobalEnv)
    } else {
      assign(
        "model.compiled.dep",
        rstan::stan_model(model_code = model.stan.dep),
        envir = .GlobalEnv
      )
    }
  })
  
  # Frequentist estimation
  freq <- tryCatch({
    frequentist_estimation(time, time_trunc)
  }, error = function(e) {
    cat("Error in frequentist estimation:", e$message, "\n")
    list(success = FALSE)
  })
  
  # Bayesian estimation with independent priors
  bayes_ind <- tryCatch({
    bayesian_estimation_ind(
      time,
      time_trunc,
      beta_r = NULL,
      mu_T_r = NULL,
      rho_r = NULL,
      compiled_model = model.compiled.ind,
      simulation = FALSE,
      seed = SEED
    )
  }, error = function(e) {
    cat("Error in Bayesian estimation with independent priors:", e$message, "\n")
    list(success = FALSE)
  })
  
  # Bayesian estimation with dependent priors
  bayes_dep <- tryCatch({
    bayesian_estimation_dep(
      time,
      time_trunc,
      beta_r = NULL,
      mu_T_r = NULL,
      rho_r = NULL,
      compiled_model = model.compiled.dep,
      simulation = FALSE,
      seed = SEED
    )
  }, error = function(e) {
    cat("Error in Bayesian estimation with dependent priors:", e$message, "\n")
    list(success = FALSE)
  })
  
  stopCluster(cl)
  
  # Final results
  results <- list(
    freq_est      = if (freq$success) freq$estimate else NA,
    bayes_est_ind = if (bayes_ind$success) bayes_ind$estimate else NA,
    fit_ind       = bayes_ind$fit,
    coverage_ind  = if (bayes_ind$success) bayes_ind$coverage else NA,
    bayes_est_dep = if (bayes_dep$success) bayes_dep$estimate else NA,
    fit_dep       = bayes_dep$fit,
    coverage_dep  = if (bayes_dep$success) bayes_dep$coverage else NA,
    success_freq  = freq$success,
    success_ind   = bayes_ind$success,
    success_dep   = bayes_dep$success,
    seed_used     = SEED
  )
  
  return(results)
}

# ------------------------------------------------------------
#   Off-season adjustment----------
# ------------------------------------------------------------

adjust_offseason <- function(x, start_obs, end_obs) {
  # Define off-season intervals (December to March) within the observation window
  years <- seq(
    from = as.numeric(format(start_obs, "%Y")),
    to   = as.numeric(format(end_obs, "%Y")),
    by   = 1
  )
  
  # Build off-season periods for each year
  off_seasons <- lapply(years, function(y) {
    start <- as.POSIXct(paste0(y, "-12-31 00:00"), tz = "UTC")
    end   <- as.POSIXct(paste0(y + 1, "-03-31 00:00"), tz = "UTC")
    c(start, end)
  })
  
  # Keep only off-season intervals that overlap the observation window
  off_seasons <- Filter(function(r) r[1] < end_obs & r[2] > start_obs, off_seasons)
  
  # Convert off-season intervals to hours relative to start_obs
  off_df <- data.frame(
    start = as.numeric(difftime(sapply(off_seasons, `[[`, 1), start_obs, units = "hours")),
    end   = as.numeric(difftime(sapply(off_seasons, `[[`, 2), start_obs, units = "hours"))
  )
  
  # Remove off-season gaps from the time scale
  x_adj <- x
  dur_gap <- 0
  for (i in seq_len(nrow(off_df))) {
    after_gap <- off_df$end[i] - dur_gap
    dur_gap <- off_df$end[i] - off_df$start[i]
    x_adj[x_adj >= after_gap] <- x_adj[x_adj >= after_gap] - dur_gap
  }
  
  # Remove events that fall inside off-season intervals
  for (i in seq_len(nrow(off_df))) {
    x_adj <- x_adj[!(x >= off_df$start[i] & x < off_df$end[i])]
  }
  
  return(x_adj)
}

# Adjust failure times by removing off-season periods
times_offseason <- lapply(times, function(vec) {
  adjust_offseason(vec, start_obs, end_obs)
})

times_offseason

# Truncation horizon in hours
time_trunc <- as.numeric(difftime(end_obs, start_obs, units = "hours"))

# ------------------------------------------------------------
# Main execution--------------
# ------------------------------------------------------------

tryCatch({
  results <- run_application_parallel_cluster(
    time       = times_offseason,
    time_trunc = time_trunc,
    SEED       = 20260418
  )
  
  # Save application results
  saveRDS(results, file = file.path(current_dir, "results_application01.rds"))
  
  # Final cleanup
  unlink(safe_temp_dir, recursive = TRUE)
}, error = function(e) {
  cat("Error in main execution:", e$message, "\n")
}, finally = {
  # Attempt cleanup even if an error occurs
  if (dir.exists(safe_temp_dir)) {
    unlink(safe_temp_dir, recursive = TRUE)
  }
})

# ------------------------------------------------------------
# Information criteria: AIC and BIC------------
# ------------------------------------------------------------

calc_AIC_BIC <- function(results, times, time_trunc) {
  
  # Number of estimated parameters
  k <- length(results$freq_est)
  
  # Total number of observed failures across all systems
  n <- sum(sapply(times, length))
  
  # Frequentist estimates
  loglik_freq <- log_like_ARI(results$freq_est, time = times, time_trunc = time_trunc)
  AIC_freq <- 2 * k - 2 * loglik_freq
  BIC_freq <- k * log(n) - 2 * loglik_freq
  
  # Bayesian estimates with independent priors
  loglik_ind <- log_like_ARI(results$bayes_est_ind, time = times, time_trunc = time_trunc)
  AIC_ind <- 2 * k - 2 * loglik_ind
  BIC_ind <- k * log(n) - 2 * loglik_ind
  
  # Bayesian estimates with dependent priors
  loglik_dep <- log_like_ARI(results$bayes_est_dep, time = times, time_trunc = time_trunc)
  AIC_dep <- 2 * k - 2 * loglik_dep
  BIC_dep <- k * log(n) - 2 * loglik_dep
  
  return(list(
    n        = n,
    k        = k,
    AIC_freq = AIC_freq,
    BIC_freq = BIC_freq,
    AIC_ind  = AIC_ind,
    BIC_ind  = BIC_ind,
    AIC_dep  = AIC_dep,
    BIC_dep  = BIC_dep
  ))
}

ic_vals <- calc_AIC_BIC(results, times = times_offseason, time_trunc = time_trunc)
print(ic_vals)

# ------------------------------------------------------------
# Deviance Information Criterion (DIC)
# ------------------------------------------------------------

calc_DIC <- function(fit_stan, times, time_trunc) {
  
  # Extract posterior draws of the model parameters
  draws <- as.data.frame(fit_stan)
  param_names <- c("beta", "mu_T", "rho")
  draws <- draws[, param_names]
  
  # Deviance at each posterior draw
  dev_samples <- apply(draws, 1, function(par) {
    -2 * log_like_ARI(par, time = times, time_trunc = time_trunc)
  })
  
  # Posterior mean deviance
  D_bar <- mean(dev_samples)
  
  # Deviance at the posterior mean
  par_mean <- colMeans(draws)
  D_hat <- -2 * log_like_ARI(par_mean, time = times, time_trunc = time_trunc)
  
  # Effective number of parameters
  p_D <- D_bar - D_hat
  
  # Deviance Information Criterion
  DIC <- D_bar + p_D
  
  return(list(D_bar = D_bar, D_hat = D_hat, p_D = p_D, DIC = DIC))
}

calc_DIC(results$fit_ind, times = times_offseason, time_trunc = time_trunc)
calc_DIC(results$fit_dep, times = times_offseason, time_trunc = time_trunc)
