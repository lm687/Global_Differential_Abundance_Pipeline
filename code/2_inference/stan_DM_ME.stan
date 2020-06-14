// Stan model for the dirichlet-multinomial
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
  vector[n] u; // coefficients for random effects
  real<lower=0> sigma_u[d-1]; // sd for random effect coefficients
  real<lower=0> overdispersion_scalar[2]; // for the dirichlet
}

transformed parameters {
  matrix[2*n,d] alpha_mean;
  matrix<lower=0>[2*n,d] alpha; 
  vector<lower=0>[2*n] overdispersion_scalars = append_row(rep_vector(overdispersion_scalar[1], n),rep_vector(overdispersion_scalar[2], n));
  
  alpha_mean = append_col( (x'*beta + Z'*rep_matrix(u, d-1)), rep_vector(0, 2*n));
  for(l in 1:(2*n)){
    alpha[l,] = to_row_vector(softmax(to_vector(alpha_mean[l,])))* overdispersion_scalars[l];
  }
}

model {
  
  for(j in 1:(d-1)){
    sigma_u ~ gamma(5, 5);
  }
  
    overdispersion_scalar ~ normal(1,1);
  

  for(d_it in 1:(d-1)){
    beta[,d_it] ~ uniform(-5, 5);
  }
  
// prior for random effects
  to_vector(u) ~ normal(0, sqrt(sigma_u));
  
  for (l in 1:(2*n) ) {
    w[l,] ~ dirichlet_multinomial(to_vector(alpha[l,]));
  }
}


