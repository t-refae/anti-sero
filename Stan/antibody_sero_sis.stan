// antibody_sero_sis.stan
//   S --lambda(a)--> I --omega--> S
// estimate_omega = 0 with omega_fixed = 0 recovers the SI model

functions {
  matrix sero_trajectory(vector lambda, real omega) {
    int A = num_elements(lambda);
    matrix[A + 1, 2] traj;
    real S = 1.0;
    real I = 0.0;
    traj[1] = [1.0, 0.0];
    for (a in 1:A) {
      real lam = lambda[a];
      real tot = lam + omega;
      if (tot > 1e-10) {                       // guard against 0/0 at old ages
        real Istar = lam / tot;
        I = Istar + (I - Istar) * exp(-tot);
        S = 1 - I;
      }
      traj[a + 1] = [S, I];
    }
    return traj;
  }

  vector log_state_lik(vector z, vector mu_x, vector sig_x,
                       real log_M, real half_log_2pi) {
    int J = num_elements(mu_x);
    vector[J] out;
    for (j in 1:J) {
      out[j] = log_sum_exp(-0.5 * square((z - mu_x[j]) / sig_x[j]) - z)
               - log(sig_x[j]) - half_log_2pi - log_M;
    }
    return out;
  }
}

data {
  int<lower=1> n;
  int<lower=1> age_max;
  array[n] int<lower=0, upper=age_max> ages;

  int<lower=1> n_phi_draws;
  array[n] vector[n_phi_draws] log_phi_draws;

  int<lower=0, upper=1> estimate_omega;
  real<lower=0>         omega_fixed;
  real<lower=0>         omega_prior_sd;
  int<lower=0, upper=1> sigma_by_state;

  int<lower=0, upper=1> foi_piecewise;
  int<lower=0>          n_brackets;
  array[age_max] int<lower=0> bracket_of_age;
}

transformed data {
  int n_sigma = sigma_by_state ? 2 : 1;
  int n_exp   = foi_piecewise  ? 0 : 1;
  int n_pw    = foi_piecewise  ? n_brackets : 0;
  real log_M        = log(n_phi_draws);
  real half_log_2pi = 0.5 * log(2 * pi());
}

parameters {
  ordered[2] mu;                                // (mu_S, mu_I++)
  vector<lower=0>[n_sigma] sigma_phi;
  array[estimate_omega] real<lower=0> omega_raw;
  array[n_exp] real<lower=0> lambda_1;
  array[n_exp] real<lower=0> kappa;
  vector<lower=0>[n_pw] lambda_b;
}

transformed parameters {
  real omega = estimate_omega ? omega_raw[1] : omega_fixed;

  vector[2] mu_x = mu;
  vector[2] sig_x;
  vector[age_max] lambda_long;

  if (sigma_by_state) {
    sig_x = sigma_phi;
  } else {
    sig_x = rep_vector(sigma_phi[1], 2);
  }

  if (foi_piecewise) {
    for (a in 1:age_max){
      lambda_long[a] = lambda_b[bracket_of_age[a]];
    }
  } else {
    for (a in 1:age_max) {
      lambda_long[a] = lambda_1[1] * exp(-kappa[1] * (a - 1));
    }
  }

  matrix[age_max + 1, 2] state_probs = sero_trajectory(lambda_long, omega);
}

model {
  matrix[age_max + 1, 2] log_state;
  for (a in 1:(age_max + 1)) {
    for (j in 1:2) {
      log_state[a, j] = log(fmax(state_probs[a, j], 1e-300));
    }
  }

  for (i in 1:n) {
    target += log_sum_exp(to_vector(log_state[ages[i] + 1])
                          + log_state_lik(log_phi_draws[i], mu_x, sig_x,
                                          log_M, half_log_2pi));
  }

  mu         ~ normal(0, 5);
  sigma_phi  ~ normal(0, 1);
  omega_raw  ~ normal(0, omega_prior_sd);
  lambda_1   ~ normal(0, 1);
  kappa      ~ normal(0, 1);
  lambda_b   ~ normal(0, 1);
}

generated quantities {
  matrix[n, 2] prob_by_group;
  vector[n] log_lik;
  real mean_years_to_reversion = omega > 1e-10 ? inv(omega) : 0;

  real lprior = normal_lpdf(mu | 0, 5)
              + normal_lpdf(sigma_phi | 0, 1)
              + normal_lpdf(omega_raw | 0, omega_prior_sd)
              + normal_lpdf(lambda_1 | 0, 1)
              + normal_lpdf(kappa | 0, 1)
              + normal_lpdf(lambda_b | 0, 1);

  for (i in 1:n) {
    vector[2] lp = to_vector(log(fmax(state_probs[ages[i] + 1], 1e-300)))
                   + log_state_lik(log_phi_draws[i], mu_x, sig_x,
                                   log_M, half_log_2pi);
    real lse = log_sum_exp(lp);
    log_lik[i] = lse;
    prob_by_group[i] = to_row_vector(exp(lp - lse));
  }
}
