functions {
  matrix sero_trajectory(vector lambda) {
    int A = num_elements(lambda);
    matrix[A + 1, 2] traj;
    real S = 1;
    traj[1] = [1.0, 0.0];
    for (a in 1:A) {
      S = S * exp(-lambda[a]);
      traj[a + 1] = [S, 1 - S];
    }
    return traj;
  }

  real partial_ll(array[] int idx, int start, int end,
                  array[] vector z, array[] vector neg_z,
                  array[] int ages, matrix log_traj,
                  vector mu_x, vector inv_sigma, vector log_norm) {
    real out = 0;
    int J = num_elements(mu_x);
    for (q in 1:size(idx)) {
      int i = idx[q];
      vector[J] lp;
      for (j in 1:J)
        lp[j] = log_traj[ages[i] + 1, j] + log_norm[j]
                + log_sum_exp(-0.5 * square((z[i] - mu_x[j]) * inv_sigma[j])
                              + neg_z[i]);
      out += log_sum_exp(lp);
    }
    return out;
  }
}

data {
  int n;
  array[n] int ages;
  int age_max;

  int<lower=1> n_phi_draws;
  matrix[n_phi_draws, n] log_phi_draws;
  
  int<lower=0, upper=1> sigma_by_state;
  int<lower=1> grainsize;
}

transformed data {
  int n_sigma = sigma_by_state ? 2 : 1;
  real log_M        = log(n_phi_draws);
  real half_log_2pi = 0.5 * log(2 * pi());

  array[n] vector[n_phi_draws] z;
  array[n] vector[n_phi_draws] nz;
  array[n] int ii;
  for (i in 1:n) {
    z[i]  = col(log_phi_draws, i);
    nz[i] = -z[i];
    ii[i] = i;
  }
}

parameters {
  ordered[2] mu;
  vector<lower=0>[n_sigma] sigma_phi;
  real<lower=0> lambda_1;
  real<lower=0> kappa;
}

transformed parameters {
  vector[2] mu_x;
  mu_x = mu;
  
  vector[2] sigma_x;
  if (sigma_by_state) {
    sigma_x = sigma_phi;
  } else {
    sigma_x = rep_vector(sigma_phi[1], 2);
  }

  vector[age_max] lambda_long;
  for(i in 1:age_max)
    lambda_long[i] = lambda_1 * exp(-kappa * (i - 1));
}

model {
  matrix[age_max + 1, 2] log_traj = log(sero_trajectory(lambda_long));
  vector[2] inv_sigma;
  vector[2] log_norm;
  for (j in 1:2) {
    inv_sigma[j] = inv(sigma_x[j]);
    log_norm[j]  = -log(sigma_x[j]) - half_log_2pi - log_M;
  }

  target += reduce_sum(partial_ll, ii, grainsize,
                       z, nz, ages, log_traj, mu_x, inv_sigma, log_norm);

  mu ~ normal(0,5);
  sigma_phi ~ normal(0,1);
  lambda_1 ~ normal(0,1);
  kappa ~ normal(0,1);
}

generated quantities {
  matrix[n,2] prob_by_group;
  vector[n] log_likelihood;
  {
    matrix[age_max + 1, 2] log_traj = log(sero_trajectory(lambda_long));

    vector[2] inv_sigma;
    vector[2] log_norm;
    for (j in 1:2) {
      inv_sigma[j] = inv(sigma_x[j]);
      log_norm[j]  = -log(sigma_x[j]) - half_log_2pi - log_M;
    }

    for (i in 1:n) {
      vector[2] lp_state;
      for (j in 1:2) {
        lp_state[j] = log_traj[ages[i] + 1, j] + log_norm[j]
                      + log_sum_exp(-0.5 * square((z[i] - mu_x[j]) * inv_sigma[j])
                                    + nz[i]);
      }
      real lse = log_sum_exp(lp_state);
      log_likelihood[i] = lse;
      prob_by_group[i] = to_row_vector(exp(lp_state - lse));
    }
  }
}
