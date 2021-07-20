## are there bias in standard multinomial regression from nnet?
library(nnet)
library(reshape2)
library(ggplot2)
source("~/Documents/PhD/GlobalDA/code/2_inference_TMB/helper_TMB.R")

d <- 5
n <- 200
lambda_nmuts <- 150

beta = rbind(runif(d-1), runif(d-1))
x = cbind(rep(1, n*2), rep(c(0,1), each=n))

give_est = function(){
  probs = softmax(cbind(x %*% beta, 0))
  Y = t(apply(probs, 1, rmultinom, size = rpois(n = 1, lambda = lambda_nmuts), n = 1))
  image(Y)
  
  ## nnet uses first category as baseline
  res_multinom = nnet::multinom(formula = Y[,ncol(Y):1]~x[,2])
  
  return(list(est_beta=coef(res_multinom)[(d-1):1,],
       true_beta=t(beta)))
}

ests = replicate(n = 100, give_est())

bias = sapply(ests['est_beta',], function(i) as.vector(i) - as.vector(t(beta)))
bias = melt(bias)
bias$type_beta = rep(c('Intercept', 'Slope'), (d-1)*n/2)
bias$idx_param = rep(rep(1:(d-1), each=2), n/2)
bias$param = paste0(bias$type_beta, bias$idx_param)

ggplot(bias, aes(x=param, group=param, y=value))+geom_boxplot()+facet_wrap(.~type_beta, scales = "free_x")


ggplot(melt(sapply(ests['est_beta',], as.vector)), aes(x=Var1, group=Var1, y=value))+geom_violin()+
  geom_point(data = cbind.data.frame(Var1=1:(2*(d-1)), value=as.vector(t(beta))), aes(x=Var1, y=value))
