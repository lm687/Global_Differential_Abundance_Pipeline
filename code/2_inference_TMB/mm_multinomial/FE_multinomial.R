rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
library(TMB)

set.seed(1245)

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

create_plot = function(){
  n = 100
  d = 2
  beta = rbind(rgamma(d-1, 1, 1),rgamma(d-1, 6, 2))
  # beta = rbind(rgamma(d-1, 1, 1), 0)
  # beta[2,d-1] = 0 ## this is set to zero as it is redundant. It is not redundant!
  x = cbind(1, sample(c(0,1), n, replace = TRUE))
  # x = cbind(1, sample(c(0), n, replace = TRUE))
  # x = x[order(x[,2]),]
  true_theta = softmax(cbind(x %*% beta + matrix(rnorm(n*(d-1), mean = 0.5, sd=0.01), ncol=(d-1)), 0))
  # image(t(true_theta)) 
  # pheatmap::pheatmap(true_theta)
  data <- list(Y = t(apply(true_theta, 1, function(p) rmultinom(1, 300, p))), n=n, d=d, x=x)
  # dev.off(); image(t(data$Y))
  
  inference_random_start = function(bool_identical_initial_beta){
    
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
      #beta_intersect=rep(0,d-1),
      #beta_slope= c(rep(0,d-1)))
      #beta_slope= c(rep(0,d-2)))
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
    betas_estimate
    parameters$beta
    beta
  
    return(list(x=x, betas_estimate = matrix(opt$par, nrow=2),
                true_beta = beta,
                initial_beta = parameters$beta))
  }
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-
  
  # c(opt$par[grep('beta_slope', names(opt$par))], 0))
  ## there is a problem with the slope betas
  
  #plot(as.vector(beta), betas_estimate)
  
  nreplicas = 30
  inference = c(replicate(nreplicas, inference_random_start(F), simplify = FALSE),
                replicate(nreplicas, inference_random_start(T), simplify = FALSE))
  
  lims = lapply(1:2, function(beta_idx) sapply(list(c(sapply(inference, function(i) i$betas_estimate)[beta_idx,],
                       sapply(inference, function(i) i$initial_beta)[beta_idx,],
                       sapply(inference, function(i) i$true_beta)[beta_idx,])), function(j) c(min(j), max(j))))
  plot(t(sapply(inference, function(i) i$betas_estimate)), col=alpha(factor(rep(c(1,2), each=nreplicas)), 0.2),
       pch='E', xlab='Beta intercept', ylab='Beta slope', xlim=lims[[1]], ylim=lims[[2]])
  points(t(sapply(inference, function(i) i$initial_beta)), col=alpha(factor(rep(c(1,2), each=nreplicas)), 0.2), pch=1, cex=0.5)
  points(t(sapply(inference, function(i) i$true_beta)), col='blue', pch='T')
  
  return(inference)

}


par(mfrow=c(2,6))
results_inference = replicate(12, create_plot())
## two examples of sets of coefficients whihc should give equally good fits

results_inference[2,4][[1]]$betas_estimate
results_inference[2,4][[1]]$true_beta

results_inference[2,5][[1]]$betas_estimate

par(mfrow=c(1,3))
image(t(softmax(cbind(results_inference[1,5][[1]]$x %*% results_inference[1,5][[1]]$true_beta, 0))), main='True')
image(t(softmax(cbind(results_inference[1,5][[1]]$x %*% results_inference[1,5][[1]]$betas_estimate, 0))), main='Estimate #1')
image(t(softmax(cbind(results_inference[2,5][[1]]$x %*% results_inference[2,5][[1]]$betas_estimate, 0))), main='Estimate #2')

#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-

par(mfrow=c(1,3))
image(t(reflect_matrix(true_theta)), main='True draws')
image(t(reflect_matrix(softmax(cbind(x %*% betas_estimate, 0)))), main='Recovered')
image(t(reflect_matrix(obj$simulate()$theta)), main='Recovered 2')


par(mfrow=c(1,2))
image(t(data$Y), main='True')
image(apply(softmax(cbind(x %*% betas_estimate, 0)), 1, function(p) rmultinom(1, 300, p)), main='Recovered')

idx_first_group = which(x[,2] == 0)
idx_second_group = which(x[,2] == 1)
head(x[idx_first_group,] %*% matrix(opt$par, nrow=2))
## first group
head(matrix(x[idx_first_group,1]) %*% betas_estimate[1,])
head(matrix(x[idx_first_group,1]) %*% beta[1,])
softmax(cbind(head(matrix(x[idx_first_group,1]) %*% betas_estimate[1,]), 0))
softmax(cbind(head(matrix(x[idx_first_group,1]) %*% beta[1,]), 0))
## second group
head((x[idx_second_group,]) %*% betas_estimate)
head((x[idx_second_group,]) %*% beta)
softmax(cbind(head((x[idx_second_group,]) %*% betas_estimate), 0))
softmax(cbind(head((x[idx_second_group,]) %*% beta), 0))

z = rnorm(2*(d-1))
z %*% opt$hessian %*% matrix(z)
sdreport(obj)
opt$hessian

# closure = function(v) v/sum(v)
# 
# theta_est = closure(exp(c( opt$par, 0)))
# 
# plot(true_theta, theta_est)
