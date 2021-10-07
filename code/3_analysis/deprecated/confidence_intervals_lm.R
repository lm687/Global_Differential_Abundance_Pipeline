##' confidence intervals of a simple linear regression
##' in 95% of cases, it should be the case that the confidence
##' interval contains the true value

source("/Users/morril01/Documents/PhD/GlobalDA/code/2_inference_TMB/helper_TMB.R")

## y = a + bx, and we only look at b (we set a_true=0 in all simulations)
give_run = function(i){
  a = 0 
  b = runif(1)
  ## create data
  N = 200
  x = runif(N)
  y = cbind(1, x) %*% c(a,b) + rnorm(n = N, mean = 0, sd = 0.02)
  res_lm = lm(y~x)
  list(true=c(a,b), estimate=  summary(res_lm)$coefficients[,'Estimate'],
       stderr=summary(res_lm)$coefficients[,'Std. Error'])
}
nreps = 3000
replicates = replicate(nreps, give_run())


params_in_CI_res = apply(replicates, 2, function(j) give_params_in_CI(vec_est = j[[2]][2], vec_stderr = j[[3]][2], vec_true = j[[1]][2]) )
table(params_in_CI_res)/length(params_in_CI_res)
## that is what I expect