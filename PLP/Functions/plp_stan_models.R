# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: PLP Stan models
# ------------------------------------------------------------
# This file defines two Stan implementations for the Power Law
# Process (PLP): one with independent Gamma priors and one with
# correlated priors based on an LKJ structure.
# ------------------------------------------------------------

# Temporary directory used during Stan compilation
Sys.setenv(TMPDIR = "/mnt/nfs/home/natalymr/tmp")
dir.create(Sys.getenv("TMPDIR"), recursive = TRUE, showWarnings = FALSE)

# Stan setup
library(rstan)
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores() - 1)

# ------------------------------------------------------------
# Stan model with independent Gamma priors
# ------------------------------------------------------------
model.stan.ind <- "
functions {
  
  // PLP intensity function
  real lambda_PLP(real t, vector param, real time_trunc) {
    real beta = param[1];
    real mu_T = param[2];
    return mu_T * beta / time_trunc * pow(t, beta - 1);
  }

  // PLP cumulative intensity at the truncation horizon
  real Lambda_PLP(vector param, real time_trunc) {
    real beta = param[1];
    real mu_T = param[2];
    return mu_T / time_trunc * pow(time_trunc, beta);
  }

  // Log-likelihood for k independent systems
  real log_like_PLP(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (i in 1:k) {
      for (j in start_idx[i]:end_idx[i]) {
        log_like += log(lambda_PLP(time_flat[j], param, time_trunc));
      }
    }

    // In the PLP, the cumulative term is the same for each system
    log_like += - k * Lambda_PLP(param, time_trunc);

    return log_like;
  }
}

data {
  int<lower=1> k;                    // number of independent systems
  int<lower=1> len_times[k];         // number of observed failures per system
  vector[sum(len_times)] time_flat;  // stacked failure times
  int<lower=1> start_idx[k];         // start position of each system in time_flat
  int<lower=1> end_idx[k];           // end position of each system in time_flat
  real<lower=0> time_trunc;          // truncation horizon

  vector[2] beta0;                   // Gamma prior hyperparameters for beta
  vector[2] mu0;                     // Gamma prior hyperparameters for mu_T
}

parameters {
  real<lower=0> beta;                // PLP shape parameter
  real<lower=0> mu_T;                // PLP scale parameter
}

transformed parameters {
  vector[2] param = [beta, mu_T]';
}

model {
  // Independent priors
  beta ~ gamma(beta0[1], beta0[2]);
  mu_T ~ gamma(mu0[1], mu0[2]);

  // Likelihood contribution
  target += log_like_PLP(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"

# ------------------------------------------------------------
# Stan model with correlated priors
# ------------------------------------------------------------
model.stan.dep <- "
functions {

  // PLP intensity function
  real lambda_PLP(real t, vector param, real time_trunc) {
    real beta = param[1];
    real mu_T = param[2];
    return mu_T * beta / time_trunc * pow(t, beta - 1);
  }

  // PLP cumulative intensity at the truncation horizon
  real Lambda_PLP(vector param, real time_trunc) {
    real beta = param[1];
    real mu_T = param[2];
    return mu_T / time_trunc * pow(time_trunc, beta);
  }

  // Log-likelihood for k independent systems
  real log_like_PLP(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (i in 1:k) {
      for (j in start_idx[i]:end_idx[i]) {
        log_like += log(lambda_PLP(time_flat[j], param, time_trunc));
      }
    }

    // In the PLP, the cumulative term is shared across systems
    log_like += - k * Lambda_PLP(param, time_trunc);

    return log_like;
  }
}

data {
  int<lower=1> k;                    // number of independent systems
  int<lower=1> len_times[k];         // number of observed failures per system
  vector[sum(len_times)] time_flat;  // stacked failure times
  int<lower=1> start_idx[k];         // start position of each system in time_flat
  int<lower=1> end_idx[k];           // end position of each system in time_flat
  real<lower=0> time_trunc;          // truncation horizon

  real<lower=0> beta0;               // baseline value for beta
  real<lower=0> mu0;                 // baseline value for mu_T
  real<lower=0> sigma_beta0;         // upper bound for beta scale
  real<lower=0> sigma_mu0;           // upper bound for mu_T scale
  real<lower=0> eta;                 // LKJ shape parameter
}

parameters {
  vector[2] Z;                       // latent standard normal vector
  cholesky_factor_corr[2] L;         // Cholesky factor of the prior correlation matrix
  real<lower=0> sigma_beta;          // prior scale for beta
  real<lower=0> sigma_mu_T;          // prior scale for mu_T
}

transformed parameters {
  vector[2] W;
  W = L * Z;                         // correlated latent vector

  // Exponential transformation ensures positivity
  real<lower=0> beta = exp(W[1] * sigma_beta + log(beta0));
  real<lower=0> mu_T = exp(W[2] * sigma_mu_T + log(mu0));

  vector[2] param = [beta, mu_T]';
}

model {
  // Correlated prior specification
  Z ~ normal(0, 1);
  L ~ lkj_corr_cholesky(eta);

  sigma_beta ~ uniform(0, sigma_beta0);
  sigma_mu_T ~ uniform(0, sigma_mu0);

  // Likelihood contribution
  target += log_like_PLP(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"