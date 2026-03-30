data {
  int<lower=1> n_individuals;                       // number of serosurveyed individuals
  int<lower=1> n_dilutions;                         // number of dilutions
  int<lower=1> n_replicates;                        // number of replicates
  vector<lower=1>[n_dilutions] d;                   // dilution factor for each observation 
  array[n_individuals, n_dilutions] int<lower=0, upper=n_replicates> z; // survival count at each dilution
  // int k0;
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
  {
    array[n_individuals] vector[n_dilutions] eta;
    for (i in 1:n_individuals) {
      eta[i] = k1 * (log_phi[i] - log(d));
    }
    
    for (i in 1:n_individuals) {
      z[i] ~ binomial(n_replicates, Phi(eta[i]));
    }
  }

  // k0 ~ cauchy(0,10);
  
  // k1 ~ cauchy(0,10); 
  
  k1 ~ cauchy(0,1);

  phi ~ cauchy(0, 500);
} 

generated quantities {
  vector[n_individuals] LD50;
  for (i in 1:n_individuals) {
    LD50[i] = log_phi[i];
  }
}
