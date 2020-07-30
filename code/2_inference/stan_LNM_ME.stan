
data{
  int<lower=1> d; // number of compositions
  int<lower=1> n; // number of samples
  int w[2*n,d]; // this is the matrix of observations, which are counts (or count-like)
  matrix[2, 2*n] x; // predictor
  matrix[n,2*n] Z; // random effects id
}


parameters{
  cov_matrix[d-1] Sigma; // covariance matrix for the normal
  vector[d-1] v[2*n]; // a realisation of the multivariate normal
  real<lower=0> var_u; // variance for random effects
  matrix[2,d-1] beta; // FE coefficients
  vector[n] u; // RE coeffcicients
}

transformed parameters{
  matrix[2*n,d-1] mu; // parameter for the normal. Unconstrained as it's a logR of two parts
  simplex[d] theta[2*n]; // this is the parameter in the simplex that is fed into the multinomial
  vector[d] full_v[2*n]; // transformed version of v, with an extra 0 for the last column

  for(i in 1:2*n){
    for(j in 1:(d-1)){
      full_v[i,j] = v[i,j];
    }
    full_v[i,d] = 0;
  }

  for(l in 1:2*n){
    theta[l,] = softmax(full_v[l,]);
  }

  mu = x'*beta + Z'*rep_matrix(u, d-1);

}

model {
  // prior for variance of random effects
  var_u ~ gamma(5, 5);
  
  // prior for random effects
  u ~ normal(0, sqrt(var_u));

  for(d_it in 1:(d-1)){
    beta[,d_it] ~ uniform(-5, 5);
  }

  for(i in 1:2*n){
    v[i,] ~ multi_normal(mu[i,], Sigma);
  }
  
  for(i in 1:2*n){
    w[i,] ~ multinomial(theta[i,]);
  }
  
}
