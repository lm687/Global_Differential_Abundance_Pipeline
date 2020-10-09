simulate_from_DM = function(beta_coefs, RE_coefs, lambda, coefficient_overdispersion){
  alphabar = softmax(cbind(sapply(1:(length(beta_coefs)/2),
                                  function(some_dummy_idx){
                                    give_z_matrix(length(RE_coefs) * 2) %*% RE_coefs}) +
                             give_x_matrix(length(RE_coefs) * 2) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  alpha_mat = alphabar*rep(exp(lambda), each=n)*coefficient_overdispersion
  return(t(apply(alpha_mat, 1, MCMCpack::rdirichlet, n=1)))
}

simulate_from_M_RE = function(beta_coefs, RE_coefs){
  RE_coefs = matrix(RE_coefs, ncol=length(beta_coefs)/2)
  n = length(RE_coefs)/(length(beta_coefs)/ 2)
  theta = softmax(cbind( (give_z_matrix(2*n) %*% RE_coefs) +
                           (give_x_matrix(2*n)) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  return(theta)
}

simulate_from_DM_RE = function(beta_coefs, RE_coefs, lambda, coefficient_overdispersion){
  RE_coefs = matrix(RE_coefs, ncol=length(beta_coefs)/2)
  n = length(RE_coefs)/(length(beta_coefs)/ 2)
  alphabar = softmax(cbind( (give_z_matrix(2*n) %*% RE_coefs) +
                              (give_x_matrix(2*n)) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  alpha_mat = alphabar*rep(exp(lambda), each=n)*coefficient_overdispersion
  return(t(apply(alpha_mat, 1, MCMCpack::rdirichlet, n=1)))
}

simulate_from_DM_RE_altpar = function(beta_coefs, RE_coefs, lambda){
  RE_coefs = matrix(RE_coefs, ncol=length(beta_coefs)/2)
  n = length(RE_coefs)/(length(beta_coefs)/ 2)
  alphabar = softmax(cbind( (give_z_matrix(2*n) %*% RE_coefs) +
                              (give_x_matrix(2*n)) %*% matrix(beta_coefs, nrow=2) , 0))
  ## here there is a single lambda
  alpha_mat = alphabar*rep(1/exp(lambda), each=n)
  return(t(apply(alpha_mat, 1, MCMCpack::rdirichlet, n=1)))
}