functions {
  array[] real sero_probs_single_group(vector lambda, int age) {

    if (age == 0)
      return {1.0, 0.0};

    real S = 1;
    real I = 0;

    for (i in 1:age) {
      real lam = lambda[i];
      real oldS = S;

      S = oldS * exp(-lam);
      I = 1 - S;
    }

    return {S, I};
  }
}

data {
  int n;
  array[n] int ages;
  int age_max;

  int<lower=1> n_phi_draws;
  matrix[n_phi_draws, n] log_phi_draws;
}

parameters {
  ordered[2] mu;
  real<lower=0> sigma;
  real<lower=0> lambda_1;
  real<lower=0> beta;
}

transformed parameters {
  vector[2] mu_x;
  mu_x = mu;

  vector[age_max] lambda_long;
  for(i in 1:age_max)
    lambda_long[i] = lambda_1 * exp(-beta * (i - 1));
}

model {

  matrix[n,2] log_prob_by_group;

  for(i in 1:n)
    log_prob_by_group[i] =
      to_row_vector(log(sero_probs_single_group(lambda_long, ages[i])));

  for(i in 1:n) {

    vector[2] lp_state;

    for(j in 1:2) {

      vector[n_phi_draws] lp_draws;

      for(m in 1:n_phi_draws)
        lp_draws[m] =
          normal_lpdf(log_phi_draws[m,i] | mu_x[j], sigma)
          - log_phi_draws[m,i];

      real log_avg = log_sum_exp(lp_draws) - log(n_phi_draws);

      lp_state[j] = log_prob_by_group[i,j] + log_avg;
    }

    target += log_sum_exp(lp_state);
  }

  mu ~ normal(0,5);
  sigma ~ normal(0,1);
  lambda_1 ~ normal(0,1);
  beta ~ normal(0,1);
}

generated quantities {

  matrix[n,2] prob_by_group;

  for(i in 1:n) {

    vector[2] lp_state;

    for(j in 1:2) {

      vector[n_phi_draws] lp_draws;

      for(m in 1:n_phi_draws)
        lp_draws[m] =
          normal_lpdf(log_phi_draws[m,i] | mu_x[j], sigma)
          - log_phi_draws[m,i];

      real log_avg = log_sum_exp(lp_draws) - log(n_phi_draws);

      lp_state[j] = log(sero_probs_single_group(lambda_long, ages[i])[j]) + log_avg;
    }

    real lse = log_sum_exp(lp_state);
    prob_by_group[i] = to_row_vector(exp(lp_state - lse));
  }
}
