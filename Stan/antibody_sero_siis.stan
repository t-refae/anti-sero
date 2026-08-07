// antibody_sero_siis.stan
//   S --lambda(a)--> I++ --psi--> I+ --omega--> S
// estimate_omega = 0 with omega_fixed = 0 recovers the SII model

functions {
  matrix sero_trajectory(vector lambda, real psi, real omega) {
    int A = num_elements(lambda);
    matrix[A + 1, 3] traj;
    matrix[3, 3] G = rep_matrix(0.0, 3, 3);
    G[3, 2] =  psi;    G[2, 2] = -psi;      // I++ -> I+
    G[1, 3] =  omega;  G[3, 3] = -omega;    // I+  -> S
    traj[1] = [1.0, 0.0, 0.0];              // S(0) = 1
    for (a in 1:A) {
      G[1, 1] = -lambda[a];
      G[2, 1] =  lambda[a];                 // S -> I++
      traj[a + 1] = to_row_vector(matrix_exp(G) * to_vector(traj[a]));
    }
    return traj;
  }

  // log of the MC-integrated state-specific likelihood, Eq. (18), for one individual
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
  
  real partial_ll(array[] int idx, int start, int end,
                  array[] vector log_phi_draws, array[] int ages,
                  matrix log_state, vector mu_x, vector sig_x,
                  real log_M, real half_log_2pi) {
    real out = 0;
    for (q in 1:size(idx)) {
      int i = idx[q];
      out += log_sum_exp(to_vector(log_state[ages[i] + 1])
                         + log_state_lik(log_phi_draws[i], mu_x, sig_x,
                                         log_M, half_log_2pi));
    }
    return out;
  }
}

data {
  int<lower=1> n;
  int<lower=1> age_max;
  array[n] int<lower=0, upper=age_max> ages;

  int<lower=1> n_phi_draws;
  array[n] vector[n_phi_draws] log_phi_draws;   // n x n_phi_draws matrix from R

  int<lower=0, upper=1> estimate_omega;
  real<lower=0>         omega_fixed;
  real<lower=0>         omega_prior_sd;
  int<lower=0, upper=1> sigma_by_state;         // 1 = one sigma per serostate

  int<lower=0, upper=1> foi_piecewise;          // 1 = piecewise-constant lambda(a)
  int<lower=0>          n_brackets;
  array[age_max] int<lower=0> bracket_of_age;   // ignored when foi_piecewise = 0
  
  int<lower=1> grainsize;
}

transformed data {
  int n_sigma  = sigma_by_state ? 3 : 1;
  int n_exp    = foi_piecewise  ? 0 : 1;
  int n_pw     = foi_piecewise  ? n_brackets : 0;
  real log_M        = log(n_phi_draws);
  real half_log_2pi = 0.5 * log(2 * pi());
  
  array[n] int ii;
  for (i in 1:n) {
    ii[i] = i;
  }
}

parameters {
  ordered[3] mu;                                // (mu_S, mu_I+, mu_I++)
  vector<lower=0>[n_sigma] sigma_phi;
  real<lower=0> psi;
  array[estimate_omega] real<lower=0> omega_raw;
  array[n_exp] real<lower=0> lambda_1;
  array[n_exp] real<lower=0> kappa;
  vector<lower=0>[n_pw] lambda_b;
}

transformed parameters {
  real omega = estimate_omega ? omega_raw[1] : omega_fixed;

  vector[3] mu_x  = [mu[1], mu[3], mu[2]]';     // indexed (S, I++, I+)
  vector[3] sig_x;
  vector[age_max] lambda_long;

  if (sigma_by_state) {
    sig_x = sigma_phi;
  } else {
    sig_x = rep_vector(sigma_phi[1], 3);
  }

  if (foi_piecewise) {
    for (a in 1:age_max){
      lambda_long[a] = lambda_b[bracket_of_age[a]];
    }
  } else {
    for (a in 1:age_max){
      lambda_long[a] = lambda_1[1] * exp(-kappa[1] * (a - 1));
    }
  }

  matrix[age_max + 1, 3] state_probs = sero_trajectory(lambda_long, psi, omega);
}

model {
  matrix[age_max + 1, 3] log_state;
  
  for (a in 1:(age_max + 1)) {
    for (j in 1:3) {
      log_state[a, j] = log(fmax(state_probs[a, j], 1e-300));
    }
  }

  target += reduce_sum(partial_ll, ii, grainsize,
                       log_phi_draws, ages, log_state, mu_x, sig_x,
                       log_M, half_log_2pi);

  mu         ~ normal(0, 5);
  sigma_phi  ~ normal(0, 1);
  psi        ~ normal(0, 0.5);
  omega_raw  ~ normal(0, omega_prior_sd);
  lambda_1   ~ normal(0, 1);
  kappa      ~ normal(0, 1);
  lambda_b   ~ normal(0, 1);
}

generated quantities {
  matrix[n, 3] prob_by_group;
  vector[n] log_lik;
  real mean_years_to_reversion = omega > 1e-10 ? inv(omega) : 0;

  // log prior, for priorsense power-scaling (item 17)
  real lprior = normal_lpdf(mu | 0, 5)
              + normal_lpdf(sigma_phi | 0, 1)
              + normal_lpdf(psi | 0, 0.5)
              + normal_lpdf(omega_raw | 0, omega_prior_sd)
              + normal_lpdf(lambda_1 | 0, 1)
              + normal_lpdf(kappa | 0, 1)
              + normal_lpdf(lambda_b | 0, 1);

  for (i in 1:n) {
    vector[3] lp = to_vector(log(fmax(state_probs[ages[i] + 1], 1e-300)))
                   + log_state_lik(log_phi_draws[i], mu_x, sig_x,
                                   log_M, half_log_2pi);
    real lse = log_sum_exp(lp);
    log_lik[i] = lse;
    prob_by_group[i] = to_row_vector(exp(lp - lse));
  }
}
