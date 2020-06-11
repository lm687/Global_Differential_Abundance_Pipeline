
data{
  int<lower=1> d; // number of compositions
  int<lower=1> n; // number of samples
  int W[2*n,d]; // this is the matrix of observations, which are counts (or count-like)
  matrix[2, 2*n] x; // predictor
  matrix[n,2*n] Z; // random effects id
}


parameters{
  cov_matrix[d-1] Sigma; // covariance matrix for the normal
  vector[d-1] v[2*n]; // a realisation of the multivariate normal
  real<lower=0> sigmasqrt_uk; // variance for random effects
  matrix[2,d-1] beta; // FE coefficients
  matrix[n,d-1] u; // RE coeffcicients
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

  for(i in 1:2*n){
    theta[i,] = softmax(full_v[i,]);
  }

  mu = x'*beta + Z'*u;

}

model {
  // prior for random effects
  for(j in 1:(d-1)){
    u[,j] ~ normal(0, sigmasqrt_uk);
  }
  
  for(i in 1:2*n){
    v[i,] ~ multi_normal(mu[i,], Sigma);
  }
  
  for(i in 1:2*n){
    W[i,] ~ multinomial(theta[i,]);
  }
  
}
