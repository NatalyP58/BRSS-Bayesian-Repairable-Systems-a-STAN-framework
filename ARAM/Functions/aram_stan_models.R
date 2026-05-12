# ------------------------------------------------------------
# Title: A Bayesian Approach to Repairable Systems Models
# Author: Nataly Martinez Riascos
# Date: April 15, 2026
# ------------------------------------------------------------
# File: ARAM Stan models
# ------------------------------------------------------------
# This file defines the Stan code used for Bayesian inference
# under the Arithmetic Reduction of Age with Memory (ARAM)
# model.
#
# Two specifications are included:
#   1. a model with independent priors, and
#   2. a model with dependent priors based on an LKJ structure.
# ------------------------------------------------------------

# Temporary directory used during Stan compilation
Sys.setenv(TMPDIR = "/mnt/nfs/home/natalymr/tmp")
dir.create(Sys.getenv("TMPDIR"), recursive = TRUE, showWarnings = FALSE)

# Stan setup
library(rstan)
rstan_options(auto_write = FALSE)
options(mc.cores = parallel::detectCores() - 1)

# ------------------------------------------------------------
# Stan model with independent priors
# ------------------------------------------------------------
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

  // Transformation used in the ARAM formulation
  real h(real theta, real beta) {
    return ((beta < 1 && theta >= 0) || (beta >= 1 && theta < 0))
           ? log((1 + theta) / (1 - theta)) * 1e4
           : theta;
  }

  // Sign correction used in the ARAM formulation
  int sign(real beta) {
    return (beta >= 1) ? 1 : -1;
  }

  // ARAM conditional intensity function
  real lambda_ARAM(real t, vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    real theta = param[3];

    real timess  = indicator(times, t);
    real h_theta = h(theta, beta);
    real signo   = sign(beta);

    return mu_T * beta / time_trunc * pow(t - signo * h_theta * timess, beta - 1);
  }

  // ARAM cumulative intensity at the truncation horizon
  real Lambda_ARAM(vector param, vector times, real time_trunc) {
    real beta    = param[1];
    real mu_T    = param[2];
    real theta   = param[3];
    real h_theta = h(theta, beta);
    real signo   = sign(beta);

    int N = num_elements(times);
    vector[N] timesPrev;
    vector[N] delta_times_1;
    vector[N] delta_times_2;
    real sum_term;
    real tail_term;

    timesPrev     = append_row(0, head(times, N - 1));
    delta_times_1 = times - signo * h_theta * timesPrev;
    delta_times_2 = timesPrev - signo * h_theta * timesPrev;

    // Sum of interval contributions
    sum_term = sum(pow(delta_times_1, beta) / time_trunc)
               - sum(pow(delta_times_2, beta) / time_trunc);

    // Final contribution from the last failure to T
    tail_term = (pow(time_trunc - signo * h_theta * times[N], beta)
                - pow(times[N] - signo * h_theta * times[N], beta)) / time_trunc;

    return mu_T * (sum_term + tail_term);
  }

  // Log-likelihood for k independent systems
  real log_like_ARAM(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (l in 1:k) {
      int n_l = end_idx[l] - start_idx[l] + 1;
      vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);

      for (i in 1:n_l) {
        log_like += log(lambda_ARAM(times_l[i], param, times_l, time_trunc));
      }

      log_like -= Lambda_ARAM(param, times_l, time_trunc);
    }

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
  real<lower=0> beta;                // ARAM shape parameter
  real<lower=0> mu_T;                // ARAM scale parameter
  real<lower=-1, upper=1> theta;     // repair-effect parameter
}

transformed parameters {
  vector[3] param = [beta, mu_T, theta]';
}

model {
  // Independent priors
  beta  ~ gamma(beta0[1], beta0[2]);
  mu_T  ~ gamma(mu0[1], mu0[2]);
  theta ~ uniform(-1, 1);

  // Likelihood contribution
  target += log_like_ARAM(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"

# ------------------------------------------------------------
# Stan model with dependent priors
# ------------------------------------------------------------
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

  // Transformation used in the ARAM formulation
  real h(real theta, real beta) {
    return ((beta < 1 && theta >= 0) || (beta >= 1 && theta < 0))
           ? log((1 + theta) / (1 - theta)) * 1e4
           : theta;
  }

  // Sign correction used in the ARAM formulation
  int sign(real beta) {
    return (beta >= 1) ? 1 : -1;
  }

  // ARAM conditional intensity function
  real lambda_ARAM(real t, vector param, vector times, real time_trunc) {
    real beta  = param[1];
    real mu_T  = param[2];
    real theta = param[3];

    real timess  = indicator(times, t);
    real h_theta = h(theta, beta);
    real signo   = sign(beta);

    return mu_T * beta / time_trunc * pow(t - signo * h_theta * timess, beta - 1);
  }

  // ARAM cumulative intensity at the truncation horizon
  real Lambda_ARAM(vector param, vector times, real time_trunc) {
    real beta    = param[1];
    real mu_T    = param[2];
    real theta   = param[3];
    real h_theta = h(theta, beta);
    real signo   = sign(beta);

    int N = num_elements(times);
    vector[N] timesPrev;
    vector[N] delta_times_1;
    vector[N] delta_times_2;
    real sum_term;
    real tail_term;

    timesPrev     = append_row(0, head(times, N - 1));
    delta_times_1 = times - signo * h_theta * timesPrev;
    delta_times_2 = timesPrev - signo * h_theta * timesPrev;

    // Sum of interval contributions
    sum_term = sum(pow(delta_times_1, beta) / time_trunc)
               - sum(pow(delta_times_2, beta) / time_trunc);

    // Final contribution from the last failure to T
    tail_term = (pow(time_trunc - signo * h_theta * times[N], beta)
                - pow(times[N] - signo * h_theta * times[N], beta)) / time_trunc;

    return mu_T * (sum_term + tail_term);
  }

  // Log-likelihood for k independent systems
  real log_like_ARAM(vector param, real time_trunc, vector time_flat, int[] start_idx, int[] end_idx, int k) {
    real log_like = 0;

    for (l in 1:k) {
      int n_l = end_idx[l] - start_idx[l] + 1;
      vector[n_l] times_l = segment(time_flat, start_idx[l], n_l);

      for (i in 1:n_l) {
        log_like += log(lambda_ARAM(times_l[i], param, times_l, time_trunc));
      }

      log_like -= Lambda_ARAM(param, times_l, time_trunc);
    }

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
  vector[3] Z;                       // latent standard normal vector
  cholesky_factor_corr[3] L;         // Cholesky factor of the prior correlation matrix
  real<lower=0> sigma_beta;          // prior scale for beta
  real<lower=0> sigma_mu_T;          // prior scale for mu_T
}

transformed parameters {
  vector[3] W;
  W = L * Z;                         // correlated latent vector

  real<lower=0> beta = exp(W[1] * sigma_beta + log(beta0));
  real<lower=0> mu_T = exp(W[2] * sigma_mu_T + log(mu0));
  real<lower=-1, upper=1> theta = tanh(W[3]);

  vector[3] param = [beta, mu_T, theta]';
}

model {
  // Correlated prior specification
  Z ~ normal(0, 1);
  L ~ lkj_corr_cholesky(eta);

  sigma_beta ~ uniform(0, sigma_beta0);
  sigma_mu_T ~ uniform(0, sigma_mu0);

  // Likelihood contribution
  target += log_like_ARAM(param, time_trunc, time_flat, start_idx, end_idx, k);
}
"