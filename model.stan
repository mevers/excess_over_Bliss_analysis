data {
  int<lower=1> N;                   // Number of observations
  int<lower=1> J;                   // Number of drug+dose combinations
  int<lower=1> K;                   // Number of stages
  int<lower=1,upper=J> drug[N];     // Drug+dose of measurement
  int<lower=1,upper=K> stage[N];    // Stage of measurement
  vector[N] eob;                    // Excess over Bliss values
}

parameters {
  matrix[J, K] mu;                  // Mean eob per drug+dose
  real<lower=0,upper=100> sigma;
  real mu_eob;
  real<lower=0,upper=100> sigma_eob;
}

model {

  // Pooled CI and priors for hyperparameters
  to_vector(mu) ~ normal(mu_eob, sigma_eob);
  sigma ~ cauchy(0, 2.5);

  mu_eob ~ normal(0, 10);
  sigma_eob ~ cauchy(0, 2.5);

  // Likelihood
  for (i in 1:N)
    eob[i] ~ normal(mu[drug[i], stage[i]], sigma);
}
