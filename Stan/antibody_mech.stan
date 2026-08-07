functions {
  real partial_binom(array[] int idx, int start, int end,
                     vector log_phi, array[,] int z,
                     int n_replicates, vector log_d, real k1) {
    real out = 0;
    for (q in 1:size(idx)) {
      int i = idx[q];
      out += binomial_lpmf(z[i] | n_replicates, Phi(k1 * (log_phi[i] - log_d)));
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
  // int k0;
  int<lower=1> grainsize;
}

transformed data {
  vector[n_dilutions] log_d = log(d);
  array[n_individuals] int ii;
  for (i in 1:n_individuals) ii[i] = i;
}

parameters { 
  // real k0;          // intercept
  real<lower=0> k1; // slope 

  // phi on the original scale, constrained to (1, 2^11) = (1, 2048)
  // vector<lower=1, upper=2048>[n_individuals] phi;
  
  vector<lower=0>[n_individuals] phi;
} 

transformed parameters {
  vector[n_individuals] log_phi;
  log_phi = log(phi);
}

model {
  target += reduce_sum(partial_binom, ii, grainsize,
                       log_phi, z, n_replicates, log_d, k1);

  k1  ~ cauchy(0, 1);
  phi ~ cauchy(0, 500);
}

generated quantities {
  vector[n_individuals] log_LD50;
  for (i in 1:n_individuals) {
    log_LD50[i] = log_phi[i];
  }
}
