
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
  // int<lower=0> m[2*n]; // number of mutations per sample \in N^{Ns}
  int<lower=0> d; // number of signatures. \in N^1
  int w[2*n,d]; // number of mutations attributed to each signature. \in N^{Ns,Nk}
  int p; // number of covariates
  matrix[p, 2*n] x; // covariate matrix
  matrix[n, 2*n] Z; // matrix for random effects
}

parameters {
  matrix[p,d] beta; // coefficients for fixed effects
  vector[n] u; // coefficients for random effects
  real<lower=0> sigma_u; // sd for random effect coefficients
}



transformed parameters {
  vector[d] prior_alpha;
  matrix<lower=0>[2*n,d] alpha;
  
    alpha = exp(x'*beta + Z'*u);

}

model {

  for(d_it in 1:d){
    beta[,d_it] ~ normal(0, 1);
  }

  for(i in 1:(2*n)){
    for(j in 1:d){
      alpha[i,j] ~ gamma(5,5);  //gamma(prior_alpha); // prior for the alphas
    }
  }
  
  // prior for random effects
   u ~ normal(0, sigma_u^2);

  for (i2 in 1:(2*n) ) {
    w[i2,] ~ dirichlet_multinomial(to_vector(alpha[i2,]));
  }
}
