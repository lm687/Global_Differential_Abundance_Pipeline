// Based on https://mc-stan.org/docs/2_18/stan-users-guide/latent-dirichlet-allocation.html
data {
  int<lower=2> K;               // num signature communities
  int<lower=2> V;               // num signatures/mutation categories
  int<lower=1> M;               // num samples
  int<lower=1> N;               // total exposure count
  int<lower=1,upper=V> w[N];    // mutation n
  int<lower=1,upper=M> doc[N];  // sample ID for mutation n
  vector<lower=0>[K] alpha;     // signature community prior
  vector<lower=0>[V] beta;      // signature/mutation category prior
}
parameters {
  simplex[K] theta[M];   // signature community distribution for sample m
  simplex[V] phi[K];     // signature/mutation category distribution for signature community k
}
model {
  for (m in 1:M)
    theta[m] ~ dirichlet(alpha);  // prior
  for (k in 1:K)
    phi[k] ~ dirichlet(beta);     // prior
  for (n in 1:N) {
    real gamma[K];
    for (k in 1:K)
      gamma[k] = log(theta[doc[n], k]) + log(phi[k, w[n]]);
    target += log_sum_exp(gamma);  // likelihood;
  }
}

