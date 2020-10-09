
// Stan model for the multinomial with mixed effects
// IN PROGRESS

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
  vector[n] ularge; // matrix of coefficients for random effects
  cov_matrix[d-1] sigma_RE; // variance for random effect coefficients
}

transformed parameters {
  // matrix[n,d-1] u_trans = rep_matrix(u, d-1);
  matrix[2*n,d-1] thetaprime;
  simplex[d] theta[2*n];

  for(j in 1:(d-1)){
    thetaprime[,j] = x'*beta[,j] + Z'*ularge;
  }

  for(l in 1:(2*n)){
    theta[l] = softmax(append_row(to_vector(thetaprime[l,]), 0));
  }
}

model {
  
  // prior for variance of random effects
  identity <- diag_matrix(rep_vector(1.0,dim)); 
  sigma_RE ~ inv_wishart(5, identity);
  
  // prior for random effects
  ularge ~ multi_normal(0, sigma_RE);

  for(d_it in 1:(d-1)){
    beta[,d_it] ~ uniform(-5, 5);
  }
  
  for (l in 1:(2*n) ) {
    w[l,] ~ multinomial(theta[l]);
  }
}

