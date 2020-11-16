## How can it be that the multinomial does give results that include the observed values when the DM equivalent doesn't?!
## (which to me makes no sense)

rm(list = ls())
setwd("~/Documents/PhD/GlobalDA/code/3_analysis/")
source("../2_inference_TMB/helper_TMB.R")
source("../3_analysis/helper/helper_simulate_from_model.R")

coefficient_overdispersion = 1000 ## value used in TMB for better convergence

M = readRDS("../../data/robjects_cache/tmb_results/fullRE_M_Kidney-RCC.clearcell_nucleotidesubstitution1.RDS")
DM = readRDS("../../data/robjects_cache/tmb_results/fullRE_DM_Kidney-RCC.clearcell_nucleotidesubstitution1.RDS")

dataset_TMB = list(M=M, DM=DM)

#-------------------------------------------------------------------------------------------------#
folder_roo = "../../data/roo/"
count_objects = sapply(list.files(folder_roo, full.names = TRUE), readRDS)
names(count_objects) = gsub("_ROO.RDS", "", list.files(folder_roo))
#-------------------------------------------------------------------------------------------------#

count_object = count_objects$`Kidney-RCC.clearcell_nucleotidesubstitution1`

simulate_theta_and_accuracy = function(dataset, model){

  beta_coefs = python_like_select_name(dataset$par.fixed, 'beta')
  RE_coefs = dataset$par.random
  lambda = python_like_select_name(dataset$par.fixed, 'log_lambda')
  
  if(model %in% c('DM', 'modified_DM')){
    sim_thetas = replicate(1e3, simulate_from_DM_RE(beta_coefs, RE_coefs, lambda, coefficient_overdispersion))
  }else if(model == 'DM_altpar'){
    sim_thetas = replicate(1e3, simulate_from_DM_RE_altpar(beta_coefs, RE_coefs, lambda))
  }else if(model == 'M'){
    sim_thetas = replicate(1e3, simulate_from_M_RE(beta_coefs, RE_coefs))
  }else{
    stop('Indicate a correct <model>')
  }
  
  ## from observed
  matrices = slot(count_object, 'count_matrices_active')
  if(sum(sapply(matrices, length)) == 0){
    matrices = slot(count_object, 'count_matrices_all')
  }
  
  matrices = do.call('rbind', lapply(matrices, round))
  vec_observed = as.vector(matrices)
  
  mut_toll = rowSums(matrices)
  ## simulate with multinomial from these thetas
  multinom_draws0 = lapply(1:dim(sim_thetas)[3], function(simulation_idx){
    sapply(1:nrow(sim_thetas[,,simulation_idx]), function(i) rmultinom( n = 1, size = mut_toll[i], 
                                                                 prob = sim_thetas[i,,simulation_idx]) )
  })
  multinom_draws = sapply(multinom_draws0, function(i) as.vector(t(i)))
  
  ## check if real values fall in the confidence interval, for each of the thetas
  
  # multinom_draws[[1]]
  # multinom_draws[[2]]
  in_conf_int = sapply(1:dim(multinom_draws)[1], function(i){
    confint_bounds = quantile(multinom_draws[i,], probs = c(0.025, 0.975))
    (vec_observed[i] >= confint_bounds[1]) & (vec_observed[i] <= confint_bounds[2])
  })
  
  ## plotting with subsample
  subset_pars = 1:length(vec_observed) #sample(1:length(vec_observed), size = 10)
  ccccc=reshape2::melt(multinom_draws[subset_pars,sample(1:1000, size = 100)])

  ## the ones that worked
  ccccc_success=reshape2::melt(multinom_draws[subset_pars,sample(1:1000, size = 100)][in_conf_int,])

  ## the ones that didn't work
  ccccc_failed=reshape2::melt(multinom_draws[subset_pars,sample(1:1000, size = 100)][!in_conf_int,])

  return(list(sim_thetas=sim_thetas, draws=multinom_draws, in_conf_int=in_conf_int, melt_all=ccccc,
              melt_success=ccccc_success, melt_failed=ccccc_failed, vec_observed_subset=vec_observed[subset_pars],
              subset_pars_idx=subset_pars, multinom_draws0=multinom_draws0))

}
res_draws = list(simulate_theta_and_accuracy(dataset_TMB$M, 'M'),
                 simulate_theta_and_accuracy(dataset_TMB$DM, 'DM'))
names(res_draws) = c('M', 'DM')

## What model does best?
mean(res_draws$M$in_conf_int)
mean(res_draws$DM$in_conf_int)

par(mfrow=c(1,2))
plot(res_draws$M$melt_all$Var1, res_draws$M$melt_all$value)
points(1:length(res_draws$M$subset_pars_idx), res_draws$M$vec_observed_subset, col='red')
plot(res_draws$DM$melt_all$Var1, res_draws$DM$melt_all$value)
points(1:length(res_draws$DM$subset_pars_idx), res_draws$DM$vec_observed_subset, col='red')

## the ones that worked
plot(res_draws$M$melt_success$Var1, res_draws$M$melt_success$value)
points(1:sum(res_draws$M$in_conf_int), res_draws$M$vec_observed_subset[res_draws$M$in_conf_int], col='red')

## the ones that didn't work
plot(res_draws$M$melt_failed$Var1, res_draws$M$melt_failed$value)
points(1:sum(!(res_draws$M$in_conf_int)), res_draws$M$vec_observed_subset[!(res_draws$M$in_conf_int)], col='red')

plot(res_draws$M$melt_failed$Var1, log(res_draws$M$melt_failed$value))
points(1:sum(!(res_draws$M$in_conf_int)), log(res_draws$M$vec_observed_subset[!(res_draws$M$in_conf_int)]), col='red')

length(res_draws$M$melt_all$Var1)
length(res_draws$DM$melt_all$Var1)

## sort them
require(dplyr)
ordered_pars = order(res_draws$M$melt_all %>% group_by(Var1) %>% dplyr::summarize(a=mean(value)) %>% dplyr::select(a))
par(mfrow=c(1,2))
sorted_M = do.call('rbind', lapply(ordered_pars, function(sorted_idx) res_draws$M$melt_all[res_draws$M$melt_all$Var1 == sorted_idx,]))
plot(rep(1:length(ordered_pars), each=100), sorted_M$value)
points(1:length(res_draws$M$subset_pars_idx), res_draws$M$vec_observed_subset[ordered_pars], col='red')
## colour by success/failure. it's not related to the absolute value of the parameter
plot(rep(1:length(ordered_pars), each=100), sorted_M$value, col=as.factor(rep(res_draws$M$in_conf_int[ordered_pars], each=100)), cex=0.2)

exp(python_like_select_name(DM$par.fixed, 'log_lambda'))
plot(python_like_select_name(M$par.fixed, 'beta'),
     python_like_select_name(DM$par.fixed, 'beta'))

plot(python_like_select_name(M$par.fixed, 'beta'),
     python_like_select_name(DM$par.fixed, 'beta'))

par(mfrow=c(1,1))
plot(c(python_like_select_name(M$par.fixed, 'beta'), M$par.random),
     c(python_like_select_name(DM$par.fixed, 'beta'), DM$par.random))
abline(coef = c(0,1), lty='dashed')

## Simulating using M with DM's overdispersion
give_mod_DM = function(new_lambda=NULL){
  modified_DM = DM
  modified_DM$par.fixed[1:(length(modified_DM$par.fixed)-2)] = M$par.fixed
  if(!is.null(new_lambda)){
    modified_DM$par.fixed[(length(modified_DM$par.fixed)-1):length(modified_DM$par.fixed)] = rep(new_lambda, 2) #M$par.fixed
  }else{
    ## leave the lambdas that you get from the DM run
  }
  modified_DM$par.random = M$par.random
  modM_draws = simulate_theta_and_accuracy(modified_DM, model = 'DM')
  return(modM_draws)
}

modM_draws = give_mod_DM()
modM_concentrated = give_mod_DM(new_lambda = 200)

## What model does best?
mean(res_draws$M$in_conf_int)
mean(res_draws$DM$in_conf_int)

## DM and modified DM don't fail in the same simulated samples
table(res_draws$DM$in_conf_int,
      modM_draws$in_conf_int)

## Does changing lambda help?
## What model does best?
mean(res_draws$M$in_conf_int)
mean(res_draws$DM$in_conf_int)
mean(modM_draws$in_conf_int) ## that does not make any sense as we have the same mean as in M, plus overdispersion!
mean(modM_concentrated$in_conf_int) ## that does not make any sense as we have the same mean as in M, plus overdispersion!

frac_inconfint_lambdas = sapply(seq(-10, 10, length.out = 20), function(new_lambda_it){
  .draws = give_mod_DM(new_lambda_it)
  return(mean(.draws$in_conf_int))
})

require(ggplot2)
ggplot(cbind.data.frame(x=1:length(frac_inconfint_lambdas), frac_inconfint_lambdas=frac_inconfint_lambdas),
       aes(x=x, y=frac_inconfint_lambdas))+
  geom_point()+geom_line()+ labs(x='Overdispersion coefficient (higher <=> more concentrated)',
                                 y = 'Fraction of samples with observed theta in simulated range')

## The more overdispersed the better

## Seeing exactly what type of data you get under each lambda value
which(!res_draws$DM$in_conf_int & res_draws$M$in_conf_int)[1]

idx_sample = 3

modM_draws$vec_observed_subset[idx_sample]
res_draws$M$vec_observed_subset[idx_sample]
res_draws$DM$vec_observed_subset[idx_sample]

par(mfrow=c(1,3))
hist(res_draws$M$draws[idx_sample,], main='M')
abline(v=modM_draws$vec_observed_subset[idx_sample], col='red')

hist(res_draws$DM$draws[idx_sample,], main='DM')
abline(v=modM_draws$vec_observed_subset[idx_sample], col='red')

hist(modM_draws$draws[idx_sample,], main='Modified DM')
abline(v=modM_draws$vec_observed_subset[idx_sample], col='red')


in_conf_int_fun = function(multinom_draws, vec_observed){
  confint_bounds = quantile(multinom_draws, probs = c(0.025, 0.975))
  (vec_observed >= confint_bounds[1]) & (vec_observed <= confint_bounds[2])
}
in_conf_int_fun(modM_draws$draws[idx_sample,], vec_observed = modM_draws$vec_observed_subset[idx_sample])
in_conf_int_fun(res_draws$M$draws[idx_sample,], vec_observed = modM_draws$vec_observed_subset[idx_sample])
in_conf_int_fun(res_draws$DM$draws[idx_sample,], vec_observed = modM_draws$vec_observed_subset[idx_sample])

dim(modM_draws$sim_thetas[idx_sample,,])
dim(res_draws$M$sim_thetas[idx_sample,,])

## there should be NO VARIATION WHATSOEVER IN THE THETAS IN THE MULTINOMIAL
res_draws$M$sim_thetas[idx_sample,1,]
res_draws$M$sim_thetas[idx_sample,1,]

## and there isn't
## it should be the case that the y axis (multinomial) shows no variation
par(mfrow=c(3,dim(res_draws$M$sim_thetas)[2]))
sapply(1:dim(res_draws$M$sim_thetas)[2], function(j){
  plot(modM_draws$sim_thetas[idx_sample,j,],
       res_draws$M$sim_thetas[idx_sample,j,],
  xlab='Modified M theta', ylab='M thetas')
abline(coef = c(0,1), lty='dashed')
})
## and now the same for draws
sapply(1:dim(res_draws$M$sim_thetas)[2], function(j){
  plot(sapply(1:dim(res_draws$M$sim_thetas)[3],
              function(i) modM_draws$multinom_draws0[[i]][j,idx_sample]),
       sapply(1:dim(res_draws$M$sim_thetas)[3],
              function(i) res_draws$M$multinom_draws0[[i]][j,idx_sample]),
       xlab='Modified M draws', ylab='M draws')
  abline(coef = c(0,1), lty='dashed')
})
sapply(1:dim(res_draws$M$sim_thetas)[2], function(j){
  plot(density(sapply(1:dim(res_draws$M$sim_thetas)[3],
              function(i) modM_draws$multinom_draws0[[i]][j,idx_sample])),
       main='In blue, M draws.\nIn black, modified M draws')
  lines(density(sapply(1:dim(res_draws$M$sim_thetas)[3],
              function(i) res_draws$M$multinom_draws0[[i]][j,idx_sample])), col='blue')
  abline(v=res_draws$M$vec_observed_subset[idx_sample], lty='dashed')
  abline(coef = c(0,1), lty='dashed')
})

## The variability in theta is (only moderately?) transmitted to the draws





