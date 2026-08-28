functions {
  real partial_binom(array[] int idx, int start, int end,
                     vector log_phi, array[,] int z,
                     int n_replicates, vector log_d, real k0, real k1) {
    real out = 0;
    for (q in 1:size(idx)) {
      int i = idx[q];
      out += binomial_lpmf(z[i] | n_replicates, Phi(k0 + k1 * (log_phi[i] - log_d)));
    }
    return out;
  }
}

data {
  int<lower=1> n_individuals;                       // number of serosurveyed individuals
  int<lower=1> n_dilutions;                         // number of dilutions
  int<lower=1> n_replicates;                        // number of replicates
  vector<lower=1>[n_dilutions] d;                   // dilution factor for each observation
  array[n_individuals, n_dilutions] int<lower=0, upper=n_replicates> z; // survival count at each dilution

  int<lower=0, upper=1> estimate_k0;                // 0: hold k0 at k0_fixed
  real k0_fixed;                                    // rungs: -2, -1, 0, 1, 2
  real<lower=0> k0_prior_sd;                        // used only when estimate_k0 = 1
  real<lower=0> phi_prior_scale;                    // 500 in antibody_mech.stan
  real<lower=0> k1_prior_scale;                     //   1 in antibody_mech.stan
  int<lower=1> grainsize;
}

transformed data {
  vector[n_dilutions] log_d = log(d);
  array[n_individuals] int ii;
  for (i in 1:n_individuals) ii[i] = i;
}

parameters {
  array[estimate_k0] real k0_free;                  // length 0 when k0 is fixed
  real<lower=0> k1;                                 // slope
  vector<lower=0>[n_individuals] phi;
}

transformed parameters {
  real k0 = estimate_k0 ? k0_free[1] : k0_fixed;
  vector[n_individuals] log_phi = log(phi);
}

model {
  target += reduce_sum(partial_binom, ii, grainsize,
                       log_phi, z, n_replicates, log_d, k0, k1);

  k1  ~ cauchy(0, k1_prior_scale);
  phi ~ cauchy(0, phi_prior_scale);

  if (estimate_k0) {
    k0_free[1] ~ normal(0, k0_prior_sd);
  }
}

generated quantities {
  
  // for priorsense
  real lprior_phi = cauchy_lpdf(phi | 0, phi_prior_scale);
  real lprior_k1  = cauchy_lpdf(k1  | 0, k1_prior_scale);
  real lprior_k0  = estimate_k0 ? normal_lpdf(k0_free[1] | 0, k0_prior_sd) : 0;
  real lprior     = lprior_phi + lprior_k1 + lprior_k0;
  
  
  
  real mean_log_phi = mean(log_phi);

  vector[n_individuals] log_LD50;
  vector[n_individuals] log_lik;
  for (i in 1:n_individuals) {
    log_LD50[i] = log_phi[i] - k0 / k1;
    log_lik[i]  = binomial_lpmf(z[i] | n_replicates,
                                Phi(k0 + k1 * (log_phi[i] - log_d)));
  }
}
