// antibody_mech.stan with k0 made estimable (fixed or free)
//
//   estimate_k0 = 0 : hold k0 at the constant k0_fixed
//   estimate_k0 = 1 : estimate k0 freely, with prior normal(0, k0_prior_sd)
//
//   transform_prior = 1 : use the exactly transformed prior on phi,
//                         cauchy(0, phi_prior_scale * exp(-k0 / k1)).
//                         Because the half-Cauchy is a scale family and
//                         phi' = phi * exp(-k0/k1), this makes the model an
//                         exact reparametrisation of the k0 = 0 model
//
//   transform_prior = 0 : leave the prior at cauchy(0, phi_prior_scale),
//                         i.e. what happens if k0 is simply freed without
//                         adjusting anything else
//
// With estimate_k0 = 0, k0_fixed = 0, transform_prior = 1 and
// phi_prior_scale = 500 this file is identical to antibody_mech.stan.

data {
  int<lower=1> n_individuals;                       // number of serosurveyed individuals
  int<lower=1> n_dilutions;                         // number of dilutions
  int<lower=1> n_replicates;                        // number of replicates
  vector<lower=1>[n_dilutions] d;                   // dilution factor for each observation
  array[n_individuals, n_dilutions] int<lower=0, upper=n_replicates> z; // survival count at each dilution

  int<lower=0, upper=1> estimate_k0;
  int<lower=0, upper=1> transform_prior;
  real k0_fixed;
  real<lower=0> k0_prior_sd;
  real<lower=0> phi_prior_scale;                    // 500 in antibody_mech.stan
}
parameters {
  array[estimate_k0] real k0_free;                  // length 0 when k0 is fixed
  real<lower=0> k1;                                 // slope
  vector<lower=0>[n_individuals] phi;
}
transformed parameters {
  real k0 = estimate_k0 ? k0_free[1] : k0_fixed;
  vector[n_individuals] log_phi;
  log_phi = log(phi);
}
model {
  {
    array[n_individuals] vector[n_dilutions] eta;
    for (i in 1:n_individuals) {
      eta[i] = k0 + k1 * (log_phi[i] - log(d));
    }

    for (i in 1:n_individuals) {
      z[i] ~ binomial(n_replicates, Phi(eta[i]));
    }
  }

  k1 ~ cauchy(0, 1);
  phi ~ cauchy(0, phi_prior_scale * (transform_prior ? exp(-k0 / k1) : 1.0));

  if (estimate_k0){
    k0_free[1] ~ normal(0, k0_prior_sd);
  }
}
generated quantities {
  vector[n_individuals] log_LD50;
  vector[n_individuals] log_lik;                    // leave-one-serum-out
  for (i in 1:n_individuals) {
    log_LD50[i] = log_phi[i] - k0 / k1;             // invariant to k0 by construction
    log_lik[i]  = binomial_lpmf(z[i] | n_replicates,
                                Phi(k0 + k1 * (log_phi[i] - log(d))));
  }
}
