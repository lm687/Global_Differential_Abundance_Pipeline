
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
  matrix[n, d-1] ularge; // matrix of coefficients for random effects
  cov_matrix[d-1] sigma_RE; // variance for random effect coefficients
}

transformed parameters {
  matrix[2*n,d-1] thetaprime;
  simplex[d] theta[2*n];

  thetaprime = x'*beta + Z'*ularge;

  for(l in 1:(2*n)){
    theta[l] = softmax(append_row(to_vector(thetaprime[l,]), 0));
  }
  
}

model {
  
  // prior for covariance matrix of random effects
  // sigma_RE ~ inv_wishart(d**2/2 - 1/2*d -1 , diag_matrix(rep_vector(0.1, d-1));
  sigma_RE ~ inv_wishart(d-1, diag_matrix(rep_vector(1.0,d-1)));

  // sigma_RE ~ inv_wishart(8, diag_matrix(rep_vector(0.1, d-1));
  
  // random effects are drawn from a multi-variate normal
  for (i in 1:(n) ) {
    ularge[i,] ~ multi_normal(rep_vector(0, d-1), sigma_RE);
  }

  // prior for beta: uniform
  for(d_it in 1:(d-1)){
    beta[,d_it] ~ uniform(-5, 5);
  }
  
  // multinomial draws from the probabilities
  for (l in 1:(2*n) ) {
    w[l,] ~ multinomial(theta[l]);
  }
}

