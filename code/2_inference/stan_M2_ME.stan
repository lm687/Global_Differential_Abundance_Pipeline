
// Stan model for the multinomial with mixed effects

data {
  int<lower=0> n; // number of samples. \in N^1
  int<lower=0> d; // number of signatures. \in N^1
  int w[2*n,d]; // number of mutations attributed to each signature. \in N^{Ns,Nk}
  int p; // number of covariates
  matrix[p, 2*n] x; // covariate matrix
  matrix[n, 2*n] Z; // matrix for random effects
}

parameters {
  matrix[p,d-1] beta; // coefficients for fixed effects
  vector[n] u; // coefficients for random effects
  real<lower=0> sigma_u; // sd for random effect coefficients
}



transformed parameters {
  matrix[2*n,d-1] thetaprime;
  simplex[d] theta[2*n];
  for(j in 1:(d-1)){
    thetaprime[,j] = x'*beta[,j] + Z'*u;
  }
  for(l in 1:(2*n)){
    theta[l] = softmax(append_row(to_vector(thetaprime[l,]), 0));
  }
}

model {
  
  // prior for standard deviation of random effects
  sigma_u ~ gamma(5, 5);
  
  // prior for random effects
  to_vector(u) ~ normal(0, sqrt(sigma_u));
  
  
  for(d_it in 1:(d-1)){
    beta[,d_it] ~ uniform(-5, 5);
  }
  
  for (l in 1:(2*n) ) {
    w[l,] ~ multinomial(theta[l]);
  }
}
