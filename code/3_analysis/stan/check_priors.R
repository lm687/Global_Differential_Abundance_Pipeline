## looking at posteriors and pairs plots of the stan runs for simulated datasets
## trying to find if something should re re-parametrised
setwd("~/Documents/PhD/GlobalDA/code/3_analysis/")
source("../2_inference_TMB/helper_TMB.R")
load_stan = function(stanfile){
  load(stanfile)
  return(list(W=W, X=X, Z=Z, d=d, fit_stan=fit_stan))
}
## this is a DM model
stan1 = load_stan("../../data/assessing_models_simulation/inference_results/20200625_n10_nlambda100_d3_beta_intensity0_Nits2000_ModelDMROO.RData")

post = as.matrix(stan1$fit_stan)

## check if priors are not well suited
## priors
## 1. beta
# for(d_it in 1:(d-1)){
#   beta[,d_it] ~ uniform(-5, 5);
# }
par(mfrow=c(1, ncol(python_like_select_colnames(post, 'beta'))))
sapply(1:ncol(python_like_select_colnames(post, 'beta')), function(idx){
  lims = c(min(min(python_like_select_colnames(post, 'beta')[,idx]),
               -5),
           max(max(python_like_select_colnames(post, 'beta')[,idx]),
               5))
  plot(density(python_like_select_colnames(post, 'beta')[,idx]), xlim=c(lims[1]-0.2, lims[2]+0.2))
  abline(v=-5, col='blue'); abline(v=5, col='blue')
})

## 2. var RE
# var_u ~ gamma(5, 5);
par(mfrow=c(1,1))
python_like_select_colnames(post, 'sigma_u') ## although it's called sigma_u, this is in fact the variance.
## the confusing name was solved in a later version of the code
plot(density(rgamma(1e3, shape = 5, rate = 5)))
points(density(python_like_select_colnames(post, 'sigma_u')), pch=16, cex=0.4, col='blue')
## in blue, posterior. in black, prior


## Looking at the pairs plot to identify why it runs to slowly
par(mfrow=c(1,1))
plot(post[,"lp__"], type='l')

pdf("~/Desktop/pairs_example_DM.png", width=20, heigh=20)
pairs(post)
dev.off()
# python_like_select_colnames(post, 'sigma_u')

