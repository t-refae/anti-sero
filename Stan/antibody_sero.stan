functions {
  array[] real sero_probs_single_group(vector lambda, real psi, int age) {

    if (age == 0) {
      return {1.0, 0.0, 0.0};
    }

    real S = 1;
    real I = 0;
    for (i in 1:age) {
      real fixed_lam = lambda[i];

      real old_S = S;
      S = old_S*exp(-fixed_lam);

      real old_I = I;

      if (psi != fixed_lam) {
          I = exp(-psi)*old_I + (fixed_lam/(psi-fixed_lam))*old_S*(exp(-fixed_lam) - exp(-psi));
      } else {
          I = exp(-psi)*(old_I + psi*old_S);
      }
    }

    real R = 1-(S+I);

    return {S, I, R};
  }
}

data {
  int n;
  array[n] int ages;
  int age_max;

  // antibody posterior draws (log scale)
  int<lower=1> n_phi_draws;
  matrix[n_phi_draws, n] log_phi_draws;
  
  int<lower=0, upper=1> sigma_by_state;
}

transformed data {
  int n_sigma = sigma_by_state ? 3 : 1;
}

parameters {
  ordered[3] mu;
  vector<lower=0>[n_sigma] sigma_phi;
  real<lower=1e-5> psi;
  real<lower=0> lambda_1;
  real<lower=0> kappa;
}

transformed parameters {

  vector[3] mu_x;
  mu_x[1] = mu[1];
  mu_x[2] = mu[3];
  mu_x[3] = mu[2];

  vector[3] sigma_x;
  if(sigma_by_state) {
    sigma_x = sigma_phi;
  } else {
    sigma_x = rep_vector(sigma_phi[1], 3);
  }
  
  vector[age_max] lambda_long;
  for(i in 1:age_max) {
    lambda_long[i] = lambda_1 * exp(-kappa * (i - 1));
  }
}

model {

  matrix[n, 3] log_prob_by_group;

  for(i in 1:n) {
    log_prob_by_group[i] = to_row_vector(log(sero_probs_single_group(lambda_long, psi, ages[i])));
  }

  for(i in 1:n) {
    vector[3] lp_state;

    for(j in 1:3) {
      // MC integrate: E_{phi~posterior}[ lognormal(phi | mu_x[j], sigma) ]
      vector[n_phi_draws] lp_draws;
      for(m in 1:n_phi_draws) {
        // lognormal_lpdf(exp(logphi)|mu,sigma) = normal_lpdf(logphi|mu,sigma) - logphi
        lp_draws[m] = normal_lpdf(log_phi_draws[m, i] | mu_x[j], sigma_x[j]) - log_phi_draws[m, i];
      }
      real log_avg = log_sum_exp(lp_draws) - log(n_phi_draws);
      lp_state[j] = log_prob_by_group[i, j] + log_avg;
    }

    target += log_sum_exp(lp_state);
  }

  // priors
  mu ~ normal(0, 5);
  sigma_phi ~ normal(0, 1);
  lambda_1 ~ normal(0, 1);
  kappa ~ normal(0, 1);
  psi ~ normal(0, 0.5);
}

generated quantities {

  matrix[n, 3] prob_by_group;
  vector[n] log_likelihood;

  {
    matrix[n, 3] log_prob_by_group;

    for(i in 1:n) {
      log_prob_by_group[i] = to_row_vector(log(sero_probs_single_group(lambda_long, psi, ages[i])));
    }

    for(i in 1:n) {
      vector[3] lp_state;

      for(j in 1:3) {
        vector[n_phi_draws] lp_draws;
        for(m in 1:n_phi_draws) {
          lp_draws[m] = normal_lpdf(log_phi_draws[m, i] | mu_x[j], sigma_x[j]) - log_phi_draws[m, i];
        }
        real log_avg = log_sum_exp(lp_draws) - log(n_phi_draws);
        lp_state[j] = log_prob_by_group[i, j] + log_avg;
      }

      real lse = log_sum_exp(lp_state);
      log_likelihood[i] = lse;
      prob_by_group[i] = to_row_vector(exp(lp_state - lse));
    }
  }
}
