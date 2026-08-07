functions {
  matrix sero_trajectory(vector lambda, real psi) {
    int A = num_elements(lambda);
    matrix[A + 1, 3] traj;
    real S = 1;
    real I = 0;
    traj[1] = [1.0, 0.0, 0.0];
    for (a in 1:A) {
      real fixed_lam = lambda[a];
      real old_S = S;
      real old_I = I;
      S = old_S * exp(-fixed_lam);
      if (psi != fixed_lam) {
        I = exp(-psi)*old_I + (fixed_lam/(psi-fixed_lam))*old_S*(exp(-fixed_lam) - exp(-psi));
      } else {
        I = exp(-psi)*(old_I + psi*old_S);
      }
      traj[a + 1] = [S, I, 1 - (S + I)];
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

  // antibody posterior draws (log scale)
  int<lower=1> n_phi_draws;
  matrix[n_phi_draws, n] log_phi_draws;
  
  int<lower=0, upper=1> sigma_by_state;
  int<lower=1> grainsize;
}

transformed data {
  int n_sigma = sigma_by_state ? 3 : 1;
  
  
  real log_M        = log(n_phi_draws);
  real half_log_2pi = 0.5 * log(2 * pi());

  array[n] vector[n_phi_draws] z;
  array[n] vector[n_phi_draws] nz;
  for (i in 1:n) {
    z[i]  = col(log_phi_draws, i);
    nz[i] = -z[i];
  }
  
  array[n] int ii;
  for (i in 1:n) {
    ii[i] = i;
  }
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
  matrix[age_max + 1, 3] log_traj = log(sero_trajectory(lambda_long, psi));

  vector[3] inv_sigma;
  vector[3] log_norm;
  for (j in 1:3) {
    inv_sigma[j] = inv(sigma_x[j]);
    log_norm[j]  = -log(sigma_x[j]) - half_log_2pi - log_M;
  }

  target += reduce_sum(partial_ll, ii, grainsize,
                       z, nz, ages, log_traj, mu_x, inv_sigma, log_norm);

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
    matrix[age_max + 1, 3] log_traj = log(sero_trajectory(lambda_long, psi));

    vector[3] inv_sigma;
    vector[3] log_norm;
    for (j in 1:3) {
      inv_sigma[j] = inv(sigma_x[j]);
      log_norm[j]  = -log(sigma_x[j]) - half_log_2pi - log_M;
    }

    for (i in 1:n) {
      vector[3] lp_state;
      for (j in 1:3) {
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
