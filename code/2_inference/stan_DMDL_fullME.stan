
// Stan model for the Dirichlet-multinomial with mixed effects
// Equivalent to stan_DMSL_fullME.stan, but this time with two lambdas; one for each group

functions {
  real dirichlet_multinomial_lpmf(int[] y, vector alpha) {
    real alpha_plus = sum(alpha);
    return lgamma(alpha_plus) + sum(lgamma(alpha + to_vector(y)))
    - lgamma(alpha_plus+sum(y)) - sum(lgamma(alpha));
  }
}

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
  real lambda1; // first parameter of overdispersion
  real lambda2; // second parameter of overdispersion
}

transformed parameters {
  matrix[2*n,d-1] thetaprime;
  simplex[d] theta[2*n];
  matrix<lower=0>[2*n,d] alpha; 
  real[2*n] loglambda_vec;
  
  // repeat the first overdispersion parameter n times, then the second overdispersion parameter
  loglambda_vec = append_col(rep_array(n, loglambda1),rep_array(n, loglambda2));
  
  thetaprime = x'*beta + Z'*ularge;
  
  for(l in 1:(2*n)){
    theta[l] = softmax(append_row(to_vector(thetaprime[l,]), 0));
  }
  
  for(l in 1:(2*n)){
    alpha[l,] = to_row_vector(softmax(to_vector(theta[l,])))* loglambda_vec[l];
  }
  
}

model {
    
  //  prior for overdispersion parameters
  loglambda1 ~ gamma(30, 5);
  loglambda2 ~ gamma(30, 5);
  
  
  // prior for covariance matrix of random effects
  sigma_RE ~ inv_wishart(d-1, diag_matrix(rep_vector(1.0,d-1)));
  
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
    w[l,] ~ dirichlet_multinomial(to_vector(alpha[l,]));
  }
}

