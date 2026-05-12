# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: PLP simulation with parallel cluster execution
# ------------------------------------------------------------
# This script runs parallel simulation experiments for the
# Power Law Process (PLP) model, including:
#   1. data generation,
#   2. frequentist estimation,
#   3. Bayesian estimation with independent priors, and
#   4. Bayesian estimation with dependent priors.
# ------------------------------------------------------------

# Initial setup
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

# Load auxiliary scripts using absolute paths
source_files <- c(
  file.path("functions", "plp_model_functions.R"),
  file.path("functions", "plp_stan_models.R"),
  file.path("functions", "plp_frequentist_fit.R"),
  file.path("functions", "plp_bayes_dependent_fit.R"),
  file.path("functions", "plp_bayes_independent_fit.R")
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
# Main parallel simulation function
# ------------------------------------------------------------
run_simulation_parallel_cluster <- function(m, k, param, time_trunc, SEED = NULL) {
  
  # If no seed is provided, generate a reproducible one from system time
  if (is.null(SEED)) {
    SEED <- as.integer(Sys.time()) %% .Machine$integer.max
    cat("No seed was provided. Generated seed:", SEED, "\n")
  } else {
    cat("Using fixed seed:", SEED, "\n")
  }
  
  RNGkind("L'Ecuyer-CMRG")
  set.seed(SEED)
  
  # One seed per simulation replication
  seeds <- SEED + seq_len(m)
  
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
      "gen_failures_PLP",
      "frequentist_estimation",
      "bayesian_estimation_ind",
      "bayesian_estimation_dep",
      "Lambda_PLP",
      "lambda_PLP",
      "log_like_PLP",
      "model.stan.ind",
      "model.stan.dep",
      "current_dir",
      "safe_temp_dir",
      "model.compiled.ind",
      "model.compiled.dep",
      "seeds",
      "SEED"
    ),
    envir = environment()
  )
  
  # Worker configuration
  clusterEvalQ(cl, {
    # Set temporary directory in workers
    Sys.setenv(TMPDIR = safe_temp_dir)
    Sys.setenv(TEMP = safe_temp_dir)
    Sys.setenv(TMP = safe_temp_dir)
    options(tmpdir = safe_temp_dir)
    
    # Load required packages
    suppressPackageStartupMessages({
      library(rstan)
      library(coda)
    })
    
    # Stan configuration for each worker
    options(mc.cores = 1)
    options(timeout = 600)
    rstan::rstan_options(auto_write = TRUE)
    Sys.setenv(STAN_NUM_THREADS = 1)
    
    # Use precompiled independent model when available; otherwise compile locally
    if (!is.null(model.compiled.ind)) {
      assign("model.compiled.ind", model.compiled.ind, envir = .GlobalEnv)
    } else {
      tryCatch({
        assign(
          "model.compiled.ind",
          rstan::stan_model(
            model_code = model.stan.ind,
            auto_write = FALSE,
            verbose = FALSE
          ),
          envir = .GlobalEnv
        )
      }, error = function(e) {
        cat("Error compiling independent model on worker:", e$message, "\n")
        assign("model.compiled.ind", NULL, envir = .GlobalEnv)
      })
    }
    
    # Use precompiled dependent model when available; otherwise compile locally
    if (!is.null(model.compiled.dep)) {
      assign("model.compiled.dep", model.compiled.dep, envir = .GlobalEnv)
    } else {
      tryCatch({
        assign(
          "model.compiled.dep",
          rstan::stan_model(
            model_code = model.stan.dep,
            auto_write = FALSE,
            verbose = FALSE
          ),
          envir = .GlobalEnv
        )
      }, error = function(e) {
        cat("Error compiling dependent model on worker:", e$message, "\n")
        assign("model.compiled.dep", NULL, envir = .GlobalEnv)
      })
    }
    
    # Load core functions inside each worker
    source_files <- c(
      "plp_model_functions.R",
      "plp_bayesian_estimation_ind.R",
      "plp_bayesian_estimation_dep.R"
    )
    
    for (file in source_files) {
      source(file.path(current_dir, file))
    }
  })
  
  # Single simulation run
  single_run <- function(i) {
    set.seed(seeds[i])
    
    cat("Starting simulation", i, "\n")
    
    beta_r <- param[1]
    mu_T_r <- param[2]
    
    # Generate failure-time data
    time <- tryCatch({
      gen_failures_PLP(k, c(beta_r, mu_T_r), time_trunc)
    }, error = function(e) {
      cat("Error generating data:", e$message, "\n")
      return(NULL)
    })
    
    if (is.null(time)) {
      return(list(
        success_freq = FALSE,
        success_ind = FALSE,
        success_dep = FALSE
      ))
    }
    
    # Frequentist estimation
    freq <- tryCatch({
      frequentist_estimation(time, time_trunc)
    }, error = function(e) {
      cat("Error in frequentist estimation:", e$message, "\n")
      list(success = FALSE)
    })
    
    # Bayesian estimation with independent priors
    bayes_ind <- tryCatch({
      if (exists("model.compiled.ind") && !is.null(model.compiled.ind)) {
        bayesian_estimation_ind(
          time,
          time_trunc,
          beta_r,
          mu_T_r,
          compiled_model = model.compiled.ind,
          seed = seeds[i]
        )
      } else {
        stop("Compiled independent model not available")
      }
    }, error = function(e) {
      cat("Error in Bayesian estimation with independent priors:", e$message, "\n")
      list(success = FALSE)
    })
    
    # Bayesian estimation with dependent priors
    bayes_dep <- tryCatch({
      if (exists("model.compiled.dep") && !is.null(model.compiled.dep)) {
        bayesian_estimation_dep(
          time,
          time_trunc,
          beta_r,
          mu_T_r,
          compiled_model = model.compiled.dep,
          seed = seeds[i]
        )
      } else {
        stop("Compiled dependent model not available")
      }
    }, error = function(e) {
      cat("Error in Bayesian estimation with dependent priors:", e$message, "\n")
      list(success = FALSE)
    })
    
    # Collect results
    list(
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
      seed_used     = seeds[i]
    )
  }
  
  # Run simulations in parallel
  results <- tryCatch({
    parLapplyLB(cl, 1:m, single_run)
  }, finally = {
    stopCluster(cl)
  })
  
  # Final summary
  conv_freq <- sum(sapply(results, function(x) x$success_freq))
  conv_ind  <- sum(sapply(results, function(x) x$success_ind))
  conv_dep  <- sum(sapply(results, function(x) x$success_dep))
  
  cat("\n--- FINAL SUMMARY ---\n")
  cat("Completed simulations:", m, "\n")
  cat("Successful frequentist fits:", conv_freq, "\n")
  cat("Successful Bayesian fits (independent priors):", conv_ind, "\n")
  cat("Successful Bayesian fits (dependent priors):", conv_dep, "\n")
  
  return(results)
}

# ------------------------------------------------------------
# Main execution
# ------------------------------------------------------------
tryCatch({
  results <- run_simulation_parallel_cluster(
    m = 1,
    k = 2,
    param = c(2, 30),
    time_trunc = 40,
    SEED = 20260418
  )
  
  # Save simulation results
  saveRDS(results, file = file.path(current_dir, "results_k101_simulation01.rds"))
  
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