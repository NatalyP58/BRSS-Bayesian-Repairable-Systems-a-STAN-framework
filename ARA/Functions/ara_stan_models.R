# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: ARA Stan models
# ------------------------------------------------------------
# This file defines the Stan code used for Bayesian inference
# under the Arithmetic Reduction of Age (ARA) model.
#
# Two specifications are included:
#   1. a model with independent priors, and
#   2. a model with dependent priors based on an LKJ structure.
# ------------------------------------------------------------

# Configure a custom temporary directory to avoid disk-space issues
Sys.setenv(TMPDIR = "/mnt/nfs/home/natalymr/tmp")
dir.create(Sys.getenv("TMPDIR"), recursive = TRUE, showWarnings = FALSE)

# Load Stan interface and configure compilation options
library(rstan)
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores())

# ------------------------------------------------------------
# Stan models for the Arithmetic Reduction of Age (ARA) model
# ------------------------------------------------------------

## Independent-prior model ----
# Stan model with independent priors
model.stan.ind <- "
functions {

  // Most recent failure time before t
  // The vector of times must be ordered increasingly
  real indicator(vector times, real t) {
    real max_val = 0;
    for (i in 1:num_elements(times)) {
      if (times[i] >= t) break;
      if (times[i] > max_val)
        max_val = times[i];
    }
    return max_val;
  }

  // ARA conditional intensity function
  real lambda_ARA(real t, vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    real theta = param[3];

    real timess = indicator(times, t);

    return mu_T * beta / time_trunc * pow(t - theta * timess, beta - 1);
  }

  // ARA cumulative intensity at the truncation horizon
  real Lambda_ARA(vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    real theta = param[3];

    int N = num_elements(times);
    vector[N] delta_times_1;
    vector[N] delta_times_2;

    // Terms of the form t_i - theta * t_{i-1}, with t_0 = 0
    delta_times_1 = times - theta * append_row(0, head(times, N - 1));

    // Terms of the form t_{i-1} - theta * t_{i-1}
    delta_times_2 = append_row(0, head(times, N - 1)) -
                    theta * append_row(0, head(times, N - 1));

    // Sum of interval contributions
    real sum_term = sum(pow(delta_times_1, beta)) - sum(pow(delta_times_2, beta));

    // Final contribution from the last failure to T
    real tail_term = pow(time_trunc - theta * times[N], beta) -
                     pow(times[N] - theta * times[N], beta);

    return mu_T / time_trunc * (sum_term + tail_term);
  }

  // Log-likelihood for k independent systems
  real log_like_ARA(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (l in 1:k) {
      int n_l = end_idx[l] - start_idx[l] + 1;
      vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);

      for (i in 1:n_l) {
        log_like += log(lambda_ARA(times_l[i], param, times_l, time_trunc));
      }

      log_like -= Lambda_ARA(param, times_l, time_trunc);
    }

    return log_like;
  }
}

data {
  int<lower=1> k;                    // number of independent systems
  int<lower=1> len_times[k];         // number of observed failures in each system
  vector[sum(len_times)] time_flat;  // stacked failure times across systems
  int<lower=1> start_idx[k];         // starting index of each system in time_flat
  int<lower=1> end_idx[k];           // ending index of each system in time_flat
  real<lower=0> time_trunc;          // truncation horizon

  vector[2] beta0;                   // Gamma prior hyperparameters for beta_raw
  vector[2] mu0;                     // Gamma prior hyperparameters for mu_T
}

parameters {
  real<lower=0> beta_raw;            // shifted parameter so that beta >= 1
  real<lower=0> mu_T;                // scale parameter
  real<lower=-1, upper=1> theta;      // repair effectiveness parameter
}

transformed parameters {
  real beta = beta_raw + 1;
  vector[3] param = [beta, mu_T, theta]';
}

model {
  // Independent priors
  beta_raw ~ gamma(beta0[1], beta0[2]);
  mu_T ~ gamma(mu0[1], mu0[2]);
  theta ~ uniform(-1, 1);

  // Likelihood contribution
  target += log_like_ARA(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"

## Dependent-prior model ----
# Stan model with dependent priors based on an LKJ structure
model.stan.dep <- "
functions {

  // Most recent failure time before t
  // The vector of times must be ordered increasingly
  real indicator(vector times, real t) {
    real max_val = 0;
    for (i in 1:num_elements(times)) {
      if (times[i] >= t) break;
      if (times[i] > max_val)
        max_val = times[i];
    }
    return max_val;
  }

  // ARA conditional intensity function
  real lambda_ARA(real t, vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    real theta = param[3];

    real timess = indicator(times, t);

    return mu_T * beta / time_trunc * pow(t - theta * timess, beta - 1);
  }

  // ARA cumulative intensity at the truncation horizon
  real Lambda_ARA(vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    real theta = param[3];

    int N = num_elements(times);
    vector[N] delta_times_1;
    vector[N] delta_times_2;

    // Terms of the form t_i - theta * t_{i-1}, with t_0 = 0
    delta_times_1 = times - theta * append_row(0, head(times, N - 1));

    // Terms of the form t_{i-1} - theta * t_{i-1}
    delta_times_2 = append_row(0, head(times, N - 1)) -
                    theta * append_row(0, head(times, N - 1));

    // Sum of interval contributions
    real sum_term = sum(pow(delta_times_1, beta)) - sum(pow(delta_times_2, beta));

    // Final contribution from the last failure to T
    real tail_term = pow(time_trunc - theta * times[N], beta) -
                     pow(times[N] - theta * times[N], beta);

    return mu_T / time_trunc * (sum_term + tail_term);
  }

  // Log-likelihood for k independent systems
  real log_like_ARA(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (l in 1:k) {
      int n_l = end_idx[l] - start_idx[l] + 1;
      vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);

      for (i in 1:n_l) {
        log_like += log(lambda_ARA(times_l[i], param, times_l, time_trunc));
      }

      log_like -= Lambda_ARA(param, times_l, time_trunc);
    }

    return log_like;
  }
}

data {
  int<lower=1> k;                    // number of independent systems
  int<lower=1> len_times[k];         // number of observed failures in each system
  vector[sum(len_times)] time_flat;  // stacked failure times across systems
  int<lower=1> start_idx[k];         // starting index of each system in time_flat
  int<lower=1> end_idx[k];           // ending index of each system in time_flat
  real<lower=0> time_trunc;          // truncation horizon

  real<lower=0> beta0;               // baseline value for beta
  real<lower=0> mu0;                 // baseline value for mu_T
  real<lower=0> sigma_beta0;         // upper bound for beta scale
  real<lower=0> sigma_mu0;           // upper bound for mu_T scale
  real<lower=0> eta;                 // LKJ shape parameter
}

parameters {
  vector[3] Z;                       // latent standard normal vector
  cholesky_factor_corr[3] L;         // Cholesky factor of the prior correlation matrix
  real<lower=0> sigma_beta;          // prior scale for beta
  real<lower=0> sigma_mu_T;          // prior scale for mu_T
}

transformed parameters {
  vector[3] W;
  W = L * Z;

  real<lower=1> beta = exp(W[1] * sigma_beta + log(beta0)) + 1;
  real<lower=0> mu_T = exp(W[2] * sigma_mu_T + log(mu0));
  real<lower=-1, upper=1> theta = tanh(W[3]);

  vector[3] param = [beta, mu_T, theta]';
}

model {
  // Dependent prior specification
  Z ~ normal(0, 1);
  L ~ lkj_corr_cholesky(eta);

  sigma_beta ~ uniform(0, sigma_beta0);
  sigma_mu_T ~ uniform(0, sigma_mu0);

  // Likelihood contribution
  target += log_like_ARA(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"