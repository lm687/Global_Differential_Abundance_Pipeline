#-------------------------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)
library(scales) ## for alpha transparency in plots
source("helper_functions.R")
set.seed(1245)
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#

softmax = function(x){
  if(is.null(dim(x))){
    ## vector
    sum_x = sum(exp(x))
    exp(x)/sum_x
  }else{
    ## matrix
    sum_x = rowSums(exp(x))
    sweep(exp(x), 1, sum_x, '/')
  }
}

reflect_matrix = function(m){
  m[nrow(m):1,]
}

compile("FE_multinomial.cpp", "-std=gnu++17")
# compile("sea.cpp", "-std=gnu++17")
dyn.load(dynlib("FE_multinomial"))
# set.seed(123)

give_data_and_parameters = function(n, d){
  beta = rbind(rgamma(d-1, 1, 1),rgamma(d-1, 6, 2))
  # beta = rbind(rgamma(d-1, 1, 1), 0)
  # beta[2,d-1] = 0 ## this is set to zero as it is redundant. It is not redundant!
  x = cbind(1, sample(c(0,1), n, replace = TRUE))
  # x = cbind(1, sample(c(0), n, replace = TRUE))
  # x = x[order(x[,2]),]
  true_theta = softmax(cbind(x %*% beta + matrix(rnorm(n*(d-1), mean = 0, sd=0.1), ncol=(d-1)), 0))
  # image(t(true_theta)) 
  # pheatmap::pheatmap(true_theta)
  data <- list(Y = t(apply(true_theta, 1, function(p) rmultinom(1, 300, p))), n=n, d=d, x=x)
  # dev.off(); image(t(data$Y))
  
  return(list(data=data, true_theta=true_theta, beta=beta, x=x))
  
}

inference_random_start = function(bool_identical_initial_beta, d, data, x, true_beta){
  
  if(bool_identical_initial_beta){
    parameters <- list( beta = array(rep(runif(1, min = -4, max = 4), 2*(d-1)), dim = c(2,d-1)) )
  }else{
    parameters <- list(
      # beta = array(rbind(rep(1,d-1), rep(NA,d-1)), dim = c(2,d-1))
      # beta = array(rep(runif(1, min = -4, max = 4), 2*(d-1)), dim = c(2,d-1))
      beta = array(c(rep(runif(1, min = -4, max = 4), (d-1)),
                     rep(runif(1, min = -4, max = 4), (d-1))),
                   dim = c(2,d-1))
    )
  }
  
  obj <- MakeADFun(data, parameters, DLL="FE_multinomial")
  obj$hessian <- TRUE
  opt <- do.call("optim", obj)
  opt
  opt$hessian ## <-- FD hessian from optim
  obj$he()    ## <-- Analytical hessian
  sdreport(obj)
  
  # betas_estimate = rbind(opt$par[grep('beta_intersect', names(opt$par))],
  #                        c(opt$par[grep('beta_slope', names(opt$par))]))
  betas_estimate = matrix(opt$par, nrow=2)
  # betas_estimate
  # parameters$beta
  # beta
  
  return(list(x=x, betas_estimate = matrix(opt$par, nrow=2),
              true_beta = true_beta,
              initial_beta = parameters$beta))
}

replicate_inference = function(n=100, d=2, nreplicas = 30){
  simulated_data = give_data_and_parameters(n = n, d = d)
  inference = c(replicate(nreplicas, inference_random_start(F, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE),
                replicate(nreplicas, inference_random_start(T, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE))
  return(inference)
}

create_plot = function(n=100, d=2){
  simulated_data = give_data_and_parameters(n = n, d = d)
  
  nreplicas = 30
  inference = c(replicate(nreplicas, inference_random_start(F, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE),
                replicate(nreplicas, inference_random_start(T, d = d, data = simulated_data$data, x=simulated_data$x, true_beta=simulated_data$beta), simplify = FALSE))
  
  lims = lapply(1:2, function(beta_idx) sapply(list(c(sapply(inference, function(i) i$betas_estimate)[beta_idx,],
                       sapply(inference, function(i) i$initial_beta)[beta_idx,],
                       sapply(inference, function(i) i$true_beta)[beta_idx,])), function(j) c(min(j), max(j))))
  plot(t(sapply(inference, function(i) i$betas_estimate)), col=alpha(factor(rep(c(1,2), each=nreplicas)), 0.2),
       pch='E', xlab='Beta intercept', ylab='Beta slope', xlim=lims[[1]], ylim=lims[[2]])
  points(t(sapply(inference, function(i) i$initial_beta)), col=alpha(factor(rep(c(1,2), each=nreplicas)), 0.2), pch=1, cex=0.5)
  points(t(sapply(inference, function(i) i$true_beta)), col='blue', pch='T')
  
  return(inference)
}

#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#

pdf("../../../../CDA_in_Cancer/results/ProjectDA/TMB/bias_beta_FE_multinomial.pdf",
    width = 10)
par(mfrow=c(2,3))
results_inference = replicate(6, create_plot())
dev.off()
## two examples of sets of coefficients whihc should give equally good fits

## I only select the first one from each dataset run since the estimates are the same for any initial value
plot(t(sapply(1:ncol(results_inference), function(dataset_idx) c(results_inference[2,dataset_idx][[1]]$betas_estimate[1,],
                                                               results_inference[2,dataset_idx][[1]]$true_beta[1,]))),
     xlab = "True intercept", ylab='Estimated intercept', main='There is a bias in underestimating\n the intercept')
abline(coef = c(0,1), col='blue', lty='dashed')
#-------------------------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------------------------#
## Now width d>2
results_inference_multivariate = replicate(50, replicate_inference(n = 100, d = 4))
results_inference_multivariate

par(mfrow=c(1,2))
plot(t(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) c(results_inference_multivariate[2,dataset_idx][[1]]$betas_estimate[1,1],
                                                                              results_inference_multivariate[2,dataset_idx][[1]]$true_beta[1,1]))),
     xlab = "True intercept", ylab='Estimated intercept', main='Intercept is also biased\n in the multivariate case\n (one logR)')
abline(coef = c(0,1), col='blue', lty='dashed')

plot(t(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) c(results_inference_multivariate[2,dataset_idx][[1]]$betas_estimate[2,1],
                                                                              results_inference_multivariate[2,dataset_idx][[1]]$true_beta[2,1]))),
     xlab = "True intercept", ylab='Estimated intercept', main='Slope is not well estimated\n in the multivariate case\n (one logR)')
abline(coef = c(0,1), col='blue', lty='dashed')

## pairs plots
# for intercept
intercepts_multivariate = t(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) c(results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[1,],
                                                                         results_inference_multivariate[1,dataset_idx][[1]]$true_beta[1,])))
colnames(intercepts_multivariate) = apply(expand.grid(paste0('beta intercept ', 1:(ncol(intercepts_multivariate)/2)), c(' estimate', ' true')), 1, paste0, collapse="")
pairs(intercepts_multivariate)
plot_pairs_with_identity(intercepts_multivariate, pch=19)

## Compute the bias
#dev.off()
par(mfrow=c(2,3))
# intercept
apply(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) results_inference_multivariate[1,dataset_idx][[1]]$true_beta[1,] - results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[1,]),
      1, hist)
# slope
apply(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) results_inference_multivariate[1,dataset_idx][[1]]$true_beta[2,] - results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[2,]),
      1, hist)

bias_intercept = apply(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) results_inference_multivariate[1,dataset_idx][[1]]$true_beta[1,] - results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[1,]),
                   1, mean)
bias_slope = apply(sapply(1:ncol(results_inference_multivariate), function(dataset_idx) results_inference_multivariate[1,dataset_idx][[1]]$true_beta[2,] - results_inference_multivariate[1,dataset_idx][[1]]$betas_estimate[2,]),
                       1, mean)
bias_intercept
bias_slope
#-------------------------------------------------------------------------------------------------------------------#


# #-------------------------------------------------------------------------------------------------------------------#
# par(mfrow=c(1,3))
# image(t(softmax(cbind(results_inference[1,5][[1]]$x %*% results_inference[1,5][[1]]$true_beta, 0))), main='True')
# image(t(softmax(cbind(results_inference[1,5][[1]]$x %*% results_inference[1,5][[1]]$betas_estimate, 0))), main='Estimate #1')
# image(t(softmax(cbind(results_inference[2,5][[1]]$x %*% results_inference[2,5][[1]]$betas_estimate, 0))), main='Estimate #2')
# #-------------------------------------------------------------------------------------------------------------------#
# 
# #-------------------------------------------------------------------------------------------------------------------#
# par(mfrow=c(1,3))
# image(t(reflect_matrix(true_theta)), main='True draws')
# image(t(reflect_matrix(softmax(cbind(x %*% betas_estimate, 0)))), main='Recovered')
# image(t(reflect_matrix(obj$simulate()$theta)), main='Recovered 2')
# 
# 
# par(mfrow=c(1,2))
# image(t(data$Y), main='True')
# image(apply(softmax(cbind(x %*% betas_estimate, 0)), 1, function(p) rmultinom(1, 300, p)), main='Recovered')
# 
# idx_first_group = which(x[,2] == 0)
# idx_second_group = which(x[,2] == 1)
# head(x[idx_first_group,] %*% matrix(opt$par, nrow=2))
# ## first group
# head(matrix(x[idx_first_group,1]) %*% betas_estimate[1,])
# head(matrix(x[idx_first_group,1]) %*% beta[1,])
# softmax(cbind(head(matrix(x[idx_first_group,1]) %*% betas_estimate[1,]), 0))
# softmax(cbind(head(matrix(x[idx_first_group,1]) %*% beta[1,]), 0))
# ## second group
# head((x[idx_second_group,]) %*% betas_estimate)
# head((x[idx_second_group,]) %*% beta)
# softmax(cbind(head((x[idx_second_group,]) %*% betas_estimate), 0))
# softmax(cbind(head((x[idx_second_group,]) %*% beta), 0))
# 
# z = rnorm(2*(d-1))
# z %*% opt$hessian %*% matrix(z)
# sdreport(obj)
# opt$hessian
# 
# # closure = function(v) v/sum(v)
# # 
# # theta_est = closure(exp(c( opt$par, 0)))
# # 
# # plot(true_theta, theta_est)
# #-------------------------------------------------------------------------------------------------------------------#
