#-------------------------------------------------------------------------------------------------#
rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
# setwd("~/Documents/PhD/GlobalDA/code/2_inference_TMB/mm_multinomial/")
library(TMB)
library(scales)
library(uuid)
library(ROCR)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(xtable)
library(mvtnorm)
source("../1_create_ROO/roo_functions.R")
source("../2_inference/helper/helper_DA_stan.R") ## for normalise_rw
source("mm_multinomial/helper_functions.R")
source("helper_TMB.R")
uuid = uuid::UUIDgenerate()
re_run_test = FALSE ## use cache
# set.seed(1245)
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
TMB::compile("mm_multinomial/ME_LNM.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/ME_LNM"))
TMB::compile("mm_multinomial/ME_multinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/ME_multinomial"))
TMB::compile("mm_multinomial/ME_dirichletmultinomial.cpp", "-std=gnu++17")
dyn.load(dynlib("mm_multinomial/ME_dirichletmultinomial"))
#-------------------------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------------------------#
datasets_files = list.files("../../data/assessing_models_simulation/datasets/", full.names = TRUE)


# a = wrapper_run_TMB(ct = datasets_files[[1]], typedata = "simulation", model = "M", simulation = TRUE)

if(re_run_test){
  runs_M = lapply(datasets_files, function(i)  try(wrapper_run_TMB(i, model = "M", typedata = "simulation", simulation = TRUE)))
  names(runs_M) = gsub(".RDS", "", basename(datasets_files))
  saveRDS(object = runs_M, file = paste0("../../data/robjects_cache/tmb_results_simulations/simulation_runs_M_", uuid, ".RDS"))
}else{
  runs_M = readRDS("../../data/robjects_cache/tmb_results_simulations/simulation_runs_M_17b8ae57-29f2-4ffe-951b-2bad65ec6387.RDS")
}
if(re_run_test){
  runs_DM = lapply(datasets_files, function(i)  try(wrapper_run_TMB(i, model = "DM", typedata = "simulation", simulation = TRUE)))
  names(runs_DM) = gsub(".RDS", "", basename(datasets_files))
  saveRDS(object = runs_DM, file = paste0("../../data/robjects_cache/tmb_results_simulations/simulation_runs_DM_", uuid, ".RDS"))
}else{
  runs_DM = readRDS("../../data/robjects_cache/tmb_results_simulations/simulation_runs_DM_17b8ae57-29f2-4ffe-951b-2bad65ec6387.RDS")
}
if(re_run_test){
  runs_LNM = lapply(datasets_files, function(i)  try( wrapper_run_TMB(i, model = "LNM", typedata = "simulation", simulation = TRUE)))
  names(runs_LNM) = gsub(".RDS", "", basename(datasets_files))
  saveRDS(object = runs_LNM, file = paste0("../../data/robjects_cache/tmb_results_simulations/simulation_runs_LNM_", uuid, ".RDS"))
}else{
  runs_LNM = readRDS("../../data/robjects_cache/tmb_results_simulations/simulation_runs_LNM_f4c2696a-252f-461f-840c-34008c0394be.RDS")
}

# t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE){
#   if( equal.variance==FALSE ) 
#   {
#     se <- sqrt( (s1^2/n1) + (s2^2/n2) )
#     # welch-satterthwaite df
#     df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
#   } else
#   {
#     # pooled standard deviation, scaled by the sample sizes
#     se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
#     df <- n1+n2-2
#   }      
#   t <- (m1-m2-m0)/se 
#   dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
#   names(dat) <- c("Difference of means", "Std Error", "t",      "p-value")
#   return(dat) 
# }

# t.test2(m1 = colMeans(props[[1]]), m2 = colMeans(props[[2]]), s1 = )

# Compositional::hotel2T2(x1 = matrix(runif(20*4), ncol=4), x2 = matrix(runif(10*4), ncol=4))

wrapper_run_ttest = function(i){
  x = readRDS(i)
  x = x[[1]]@count_matrices_all
  props = sapply(x, normalise_rw, simplify = FALSE)
  return(Compositional::hotel2T2(x1 = compositions::ilr(props[[1]]), x2 = compositions::ilr(props[[2]]))$pvalue)
}
runs_ttest_proportions = lapply(datasets_files, function(i)  try(wrapper_run_ttest(i)))
pvals_ttest = as.numeric(unlist(runs_ttest_proportions))

## assess if there was good convergence
runs_M[[10]]
sapply(runs_M, typeof)

# which_betas = which(names(runs_M[[1]]$par.fixed) == "beta")
# R = runs_M[[1]]$par.fixed[which_betas]
# 
# wald_TMB_wrapper(runs_M[[1]])

pvals_M = sapply(runs_M, wald_TMB_wrapper)
pvals_DM = sapply(runs_DM, wald_TMB_wrapper)
pvals_LNM = sapply(runs_LNM, wald_TMB_wrapper)

pdf("../../results/TMB/pvals_simulation.pdf", height = 3, width = 10)
par(mfrow=c(1,3))
plot(sort(na.omit(pvals_M)), type = "h", main="M")
plot(sort(na.omit(pvals_DM)), type = "h", main='DM')
plot(sort(na.omit(pvals_LNM)), type = "h", main="LNM")
dev.off()

## adjust p-values
## TO DO
pvals_M_adj = pvals_M
pvals_DM_adj = pvals_DM
pvals_LNM_adj = pvals_LNM
pvals_ttest_adj = pvals_ttest

sapply(list(M=pvals_M_adj <= 0.05, DM=pvals_DM_adj <= 0.05, LNM=pvals_LNM_adj <= 0.05, ttest=pvals_ttest_adj <= 0.05), table)

VennDiagram::venn.diagram(as.vector(table(pvals_M <= 0.05, pvals_DM <= 0.05, pvals_LNM <= 0.05)))


## which simulations were actually differentially abundant?
datasets = lapply(datasets_files, readRDS)
DA_bool = ( sapply(datasets, function(i) i$beta_gamma_shape) > 0 )
sapply(datasets, function(i) i$beta[1,])
table(DA_bool) ## most are differentially abundant

## Remove outlier datasets
plot(density(sapply(datasets, function(i) max(i$beta[2,]))))
remove_idx_datasets = which(sapply(datasets, function(i) max(i$beta[2,])) > 5)
datasets = datasets[-remove_idx_datasets]
pvals_M_adj = pvals_M_adj[-remove_idx_datasets]
pvals_DM_adj = pvals_DM_adj[-remove_idx_datasets]
pvals_LNM_adj = pvals_LNM_adj[-remove_idx_datasets]
pvals_ttest_adj = pvals_ttest_adj[-remove_idx_datasets]
runs_M = runs_M[-remove_idx_datasets]
runs_DM = runs_DM[-remove_idx_datasets]
runs_LNM = runs_LNM[-remove_idx_datasets]
DA_bool = DA_bool[-remove_idx_datasets]

summarise_DA_detection = function(true, predicted){
  ## remove NAs
  which_na = which(is.na(predicted))
  if(length(which_na) > 0){ ## some NA
    true = true[-which_na]
    predicted = predicted[-which_na]
  }
  
  FPs = sum(!true & predicted)/sum(predicted)
  TPs = sum(true & predicted)/sum(predicted)
  TNs = sum(!true & !predicted)/sum(!predicted)
  FNs = sum(true & !predicted)/sum(!predicted)
  total_pos = sum(true | predicted)
  Power = TPs/total_pos
  Sensitivity = TPs / (TPs + FNs)
  Specificity = TNs / (TNs + FPs)
  pred <- prediction(as.numeric(true), as.numeric(predicted))
  AUC = try(performance(pred, "auc")@y.values[[1]])
  return(c(FP=FPs, TP=TPs, Power=Power, AUC=AUC, Specificity=Specificity, Sensitivity=Sensitivity))
}

## they are all pretty horrendous in terms of FP
summarise_DA_detection(true = DA_bool, predicted = pvals_M_adj <= 0.05)
summarise_DA_detection(true = DA_bool, predicted = pvals_DM_adj <= 0.05)
summarise_DA_detection(true = DA_bool, predicted = pvals_LNM_adj <= 0.05)
summarise_DA_detection(true = DA_bool, predicted = pvals_ttest_adj <= 0.05)

all_summary = do.call('rbind', list(M=summarise_DA_detection(true = DA_bool, predicted = pvals_M_adj <= 0.05),
                                    DM=summarise_DA_detection(true = DA_bool, predicted = pvals_DM_adj <= 0.05),
                                    LNM=summarise_DA_detection(true = DA_bool, predicted = pvals_LNM_adj <= 0.05),
                                    ttest=summarise_DA_detection(true = DA_bool, predicted = pvals_ttest_adj <= 0.05)))
xtable::xtable(all_summary, digits=c(0,4,4,4,4,4,4))

table(DA_bool, pvals_DM_adj <= 0.05)

## here there should be a negative correlation (higher effect size <=> lower p-value) (which there is to some extend)
par(mfrow=c(1,3))
plot( sapply(datasets, function(i) i$beta_gamma_shape), pvals_M_adj)
plot( sapply(datasets, function(i) i$beta_gamma_shape), pvals_DM_adj)
plot( sapply(datasets, function(i) i$beta_gamma_shape), pvals_LNM_adj)

## compare estimated betas to actual betas
par(mfrow=c(1,1))
plot((unlist(sapply(datasets, function(i) i$beta[2,]))),
     (unlist(sapply(runs_M, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta"))))))
abline(coef = c(0, 1))

for(j in which(sapply(runs_DM, typeof) == "character")){
  runs_DM[[j]] = list(par.fixed=c(beta=rep(NA, 2*(datasets[[j]]$d-1)))) ## *2 for slope and intercept
}

for(j in which(sapply(runs_LNM, typeof) == "character")){
  runs_LNM[[j]] = list(par.fixed=c(beta=rep(NA, 2*(datasets[[j]]$d-1)))) ## *2 for slope and intercept
}

## here there is only the intercept
df_beta_recovery = cbind.data.frame(beta_true = unlist(sapply(datasets, function(i) i$beta[2,])),
                                    idx = rep(1:length(datasets) , unlist(sapply(datasets, function(i) i$d))-1),
                                    d =  rep(unlist(sapply(datasets, function(i) i$d)), unlist(sapply(datasets, function(i) i$d))-1),
                                    n =  rep(unlist(sapply(datasets, function(i) i$n)), unlist(sapply(datasets, function(i) i$d))-1),
                                    beta_gamma_shape =  rep(unlist(sapply(datasets, function(i) i$beta_gamma_shape)), unlist(sapply(datasets, function(i) i$d))-1),
                                    beta_est_M = unlist(sapply(runs_M, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta")))),
                                    beta_est_DM = unlist(sapply(runs_DM, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta")))),
                                    beta_est_LNM = unlist(sapply(runs_LNM, function(i) select_slope_2(python_like_select_name(i$par.fixed, "beta")))),
                                    pvals_M_adj=rep(pvals_M_adj, unlist(sapply(datasets, function(i) i$d))-1),
                                    pvals_DM_adj=rep(pvals_DM_adj, unlist(sapply(datasets, function(i) i$d))-1),
                                    pvals_LNM_adj=rep(pvals_LNM_adj, unlist(sapply(datasets, function(i) i$d))-1))
df_beta_recovery$bool_zero_true_beta = factor((df_beta_recovery$beta_true == 0), levels = c(TRUE, FALSE))

pairs(df_beta_recovery[,c('pvals_M_adj', 'pvals_DM_adj', 'pvals_LNM_adj')], main='Pairs plot of p-values')
pairs(df_beta_recovery[,c('beta_est_M', 'beta_est_DM', 'beta_est_LNM')], main='Pairs plot of betas')

table(na.omit(df_beta_recovery$pvals_M_adj <= 0.05))
table(na.omit(df_beta_recovery$pvals_DM_adj <= 0.05))
table(na.omit(df_beta_recovery$pvals_LNM_adj <= 0.05))

#----------------------------------------------------------------------------------------------------------#
## zero slopes are not well detected
ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=beta_gamma_shape))+geom_point()+
  geom_abline(intercept = 0, slope = 1)

ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=beta_gamma_shape))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~idx,scales = "free")

ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_DM), col=beta_gamma_shape))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~idx,scales = "free")

ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_LNM), col=beta_gamma_shape))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~idx,scales = "free")
#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#
## plotting separately the differentially abundant (right) and non-differentially abundant (left) datasets
ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=n))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")

## Multinomial
p <- ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=(pvals_M_adj<0.05)))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
gp <- ggplotGrob(p); facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]; gp$widths[facet.columns] <- gp$widths[facet.columns] * c(1,4); grid::grid.draw(gp)

## Dirichlet-Multinomial
## there are some very extreme values of beta
plot(density(log(abs(na.omit(df_beta_recovery$beta_est_DM))))) ## all betas
plot(density(na.omit(sapply(runs_DM, function(i) max(python_like_select_name(i$par.fixed, "beta"))))))
df_beta_recovery[df_beta_recovery$idx %in% which(sapply(runs_DM, function(i) max(abs(python_like_select_name(i$par.fixed, "beta")))) > 10),'beta_est_DM'] = NA
p <- ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_DM), col=(pvals_DM_adj)<0.05))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
gp <- ggplotGrob(p); facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]; gp$widths[facet.columns] <- gp$widths[facet.columns] * c(1,4); grid::grid.draw(gp)

## LNM
df_beta_recovery[df_beta_recovery$idx %in% which(sapply(runs_LNM, function(i) max(abs(python_like_select_name(i$par.fixed, "beta")))) > 20),'beta_est_LNM'] = NA
p <- ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_LNM), col=(pvals_LNM_adj)<0.05))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+facet_wrap(.~bool_zero_true_beta, scales = "free_x")
gp <- ggplotGrob(p); facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]; gp$widths[facet.columns] <- gp$widths[facet.columns] * c(1,4); grid::grid.draw(gp)
#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#
## colour code per patient
ggplot(df_beta_recovery,
       aes(x=(beta_true), y=(beta_est_M), col=factor(idx)))+geom_point()+
  geom_abline(intercept = 0, slope = 1)+guides(color=FALSE)
#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#
## there are some simulated betas which are extreme and I have removed those

## a clear problem is some datasets with no zero or near-zero differential abundance having a very small p-value
df_beta_recovery[df_beta_recovery$beta_true == 0,][1:10,]
## e.g. idx=6
df_beta_recovery[df_beta_recovery$idx == 6,]


simulate_theta_from_TMB_M = function(tmb_object, dataset_object){
  ## there is a column of zeros because the intercepts have been absorbed by the random effects intercept
  theta_sim = softmax(
    cbind(t(dataset_object$X_sim) %*% t(matrix(python_like_select_name(tmb_object$par.fixed, 'beta'), ncol = 2)) +
    replicate(n = dataset_object$d-1, tmb_object$par.random %*% dataset_object$Z_sim, simplify = TRUE), 0)
    )
  list(theta_sim[1:dataset_object$n,], theta_sim[(dataset_object$n+1):(2*dataset_object$n),])
}

simulate_theta_from_TMB_DM = function(tmb_object, dataset_object){
  ## there is a column of zeros because the intercepts have been absorbed by the random effects intercept
  alphabar_sim = softmax(
    cbind(t(dataset_object$X_sim) %*% t(matrix(python_like_select_name(tmb_object$par.fixed, 'beta'), ncol = 2)) +
            replicate(n = dataset_object$d-1, tmb_object$par.random %*% dataset_object$Z_sim, simplify = TRUE), 0)
  )
  theta_sim = t(apply(alphabar_sim, 1, function(i) MCMCpack::rdirichlet(1, i*exp(tmb_object$par.fixed['log_lambda']))))
  list(theta_sim[1:dataset_object$n,], theta_sim[(dataset_object$n+1):(2*dataset_object$n),])
}

simulate_theta_from_TMB_LNM = function(tmb_object, dataset_object){
  ## there is a column of zeros because the intercepts have been absorbed by the random effects intercept
  mu_sim = t(dataset_object$X_sim) %*% t(matrix(python_like_select_name(tmb_object$par.fixed, 'beta'), ncol = 2)) +
            replicate(n = dataset_object$d-1, tmb_object$par.random %*% dataset_object$Z_sim, simplify = TRUE)
  
  cov_matrix = matrix(NA, nrow=(dataset_object$d-1), ncol=(dataset_object$d-1))
  fill_cov_matrix = matrix(python_like_select_name(tmb_object$par.fixed, 'cov_par'))
  ## fill in the order that TMB's UNSTRUCTURED_CORR_t saves the covariances
  cov_matrix[unlist(sapply(2:nrow(cov_matrix), function(rw) seq(from = rw,length.out = (rw-1), by = nrow(cov_matrix) )))] = fill_cov_matrix
  cov_matrix[unlist(sapply(2:nrow(cov_matrix), function(cl) seq(from = (((cl-1)*nrow(cov_matrix))+1),length.out = (cl-1), by = 1 )))] = fill_cov_matrix
  diag(cov_matrix) = 1 ## because this is how it was simulated
  theta_sim2 = softmax(cbind(t(apply(mu_sim, 1, function(mu_it) mvtnorm::rmvnorm(n = 1, mean = mu_it, sigma = cov_matrix))),
                     0))
  list(theta_sim2[1:dataset_object$n,], theta_sim2[(dataset_object$n+1):(2*dataset_object$n),])
}

par(mfrow=c(4,2))
sapply(datasets[[6]][[1]]@count_matrices_all, function(i) image(t(i), main='true'))
sapply(simulate_theta_from_TMB_M(runs_M[[6]], datasets[[6]]) , function(i) image(t(i), main='sim M est'))
sapply(simulate_theta_from_TMB_DM(runs_DM[[6]], datasets[[6]]) , function(i) image(t(i), main='sim DM est'))
sapply(simulate_theta_from_TMB_LNM(runs_LNM[[6]], datasets[[6]]) , function(i) image(t(i), main='sim LNM est'))

simulate_theta_from_TMB_DM(runs_M[[6]], datasets[[6]])
runs_DM[[6]]$par.random

## another problem is clearly differential abundance runs having a very large p-value
#----------------------------------------------------------------------------------------------------------#

#----------------------------------------------------------------------------------------------------------#
# Compare to stan runs
folder_stan = "../../data/assessing_models_simulation/inference_results/"
stanfiles = list.files(folder_stan, full.names = TRUE)
load_stan = function(stanfile){
  load(stanfile)
  return(list(W=W, X=X, Z=Z, d=d, fit_stan=fit_stan))
}
# names(stan_simA) = gsub(".RData", "", basename(stanfiles))
# rhats = sapply(stan_simA, function(i) rstan::Rhat(as.matrix(i$fit_stan)))

stan_simA_M = gsub(".RData", "", basename(stanfiles[grepl("ModelMROO", stanfiles)]))
stan_simA_DM = stanfiles[grepl("ModelDMROO", stanfiles)]
stan_simA_DM = gsub(".RData", "", basename(stan_simA_DM[!grepl("Nits1000_ModelDMROO", stan_simA_DM)]))
stan_simA_LNM = gsub(".RData", "", basename(stanfiles[grepl("ModelLNMROO", stanfiles)]))

# min(rhats)
## supossing this is well computed, none of these runs converged

mod_names_stan_M = stan_simA_M %>% gsub("_nlambda", "_", .) %>% gsub("_d", "_", .) %>% gsub("_n", "_", .) %>% 
  gsub("beta_intensity", "", .) %>% 
  gsub("_Nits1000_ModelMROO", "", .) 
mod_names_stan_DM = stan_simA_DM %>% gsub("_nlambda", "_", .) %>% gsub("_d", "_", .) %>% gsub("_n", "_", .) %>% 
  gsub("beta_intensity", "", .) %>% 
  gsub("_Nits2000_ModelDMROO", "", .) 
mod_names_stan_LNM = stan_simA_LNM %>% gsub("_nlambda", "_", .) %>% gsub("_d", "_", .) %>% gsub("_n", "_", .) %>% 
  gsub("beta_intensity", "", .) %>% 
  gsub("_Nits2000_ModelLNMROO", "", .) 
mod_names_M = gsub("_dataset", "", names(runs_M))
mod_names_DM = gsub("_dataset", "", names(runs_DM))
mod_names_LNM = gsub("_dataset", "", names(runs_LNM))

# dim(python_like_select_colnames(as.matrix(stan_simA[[1]]$fit_stan), "beta\\[2,"))

stan_simA_M = stan_simA_M[match(mod_names_M, mod_names_stan_M)]
stan_simA_DM = stan_simA_DM[match(mod_names_DM, mod_names_stan_DM)]
stan_simA_LNM = stan_simA_LNM[match(mod_names_LNM, mod_names_stan_LNM)]

sapply(list(stan_simA_M, stan_simA_DM, stan_simA_LNM), names)

idx = (1:length(stan_simA_M))[51]

give_plot_stan_vs_TMB_density = function(stan_list, TMB_list, nameout){
pdf(paste0("../../results/assessing_models/", nameout), width = 10, height = 3)
sapply(1:length(stan_list), function(idx){
  if(!is.na(stan_list[[idx]])){
    .xx = load_stan(paste0(folder_stan, stan_list[[idx]], ".RData"))
    post_betaslope = python_like_select_colnames(as.matrix(.xx$fit_stan), "beta\\[2,")
    ML_betaslope = select_slope_2(python_like_select_name(runs_M[[idx]]$par.fixed, "beta"))
    
    par(mfrow=c(1,ncol(post_betaslope)))
    for(i in 1:ncol(post_betaslope)){
      plot(density(post_betaslope[,i]))
      abline(v=ML_betaslope[i])
      # abline(v=ML_betaslope_2[i], col='red')
    }
  }
})
dev.off()
}

give_plot_stan_vs_TMB_density(stan_list = stan_simA_M, TMB_list = runs_M, nameout = "simA_M_stan_TMB_betaslope.pdf")
give_plot_stan_vs_TMB_density(stan_list = stan_simA_DM, TMB_list = runs_DM, nameout = "simA_DM_stan_TMB_betaslope.pdf")
give_plot_stan_vs_TMB_density(stan_list = stan_simA_LNM, TMB_list = runs_LNM, nameout = "simA_LNM_stan_TMB_betaslope.pdf")


## scatterplot for betas for mean of posterior and ML from TMB

means_stan = function(stanfiles_list){
  lapply(stanfiles_list, function(stanfile){
  if(!is.na(stanfile)){
    stanfile = paste0(folder_stan, stanfile, '.RData')
    cat(stanfile, '\n')
    load(stanfile)
    betas = (python_like_select_colnames(as.matrix(fit_stan), 'beta\\[2'))
    if(is.null(dim(betas))){ ## one column
      beta_means=mean(betas)
    }else{
      beta_means=colMeans(betas)
    }
    return(beta_means)
  }else{return(NA)}
  })
}

means_stan_M = means_stan(stan_simA_M)
means_stan_DM = means_stan(stan_simA_DM)
means_stan_LNM = means_stan(stan_simA_LNM)
## add as many NAs as necessary for runs for which we don't have data
for(j in which(sapply(means_stan_M, typeof) == "logical")){ ## they are NA
  means_stan_M[[j]] = rep(NA,  (datasets[[j]]$d-1))
}
for(j in which(sapply(means_stan_DM, typeof) == "logical")){ ## they are NA
  means_stan_DM[[j]] = rep(NA,  (datasets[[j]]$d-1))
}
for(j in which(sapply(means_stan_LNM, typeof) == "logical")){ ## they are NA
  means_stan_LNM[[j]] = rep(NA,  (datasets[[j]]$d-1))
}

df_beta_recovery = cbind(df_beta_recovery,
                         means_stan_M = as.vector(unlist(means_stan_M)),
                         means_stan_DM = as.vector(unlist(means_stan_DM)),
                         means_stan_LNM = as.vector(unlist(means_stan_LNM)))

ggplot(df_beta_recovery, aes(x=beta_est_M, means_stan_M))+geom_point()+geom_abline(coef = c(0,1))
ggsave("../../results/assessing_models/simA_M_scatter.pdf")
ggplot(df_beta_recovery, aes(x=beta_est_DM, means_stan_DM, shape=factor(idx)))+geom_point()+geom_abline(coef = c(0,1))+
  guides(shape=FALSE)
ggsave("../../results/assessing_models/simA_DM_scatter.pdf")
ggplot(df_beta_recovery, aes(x=beta_est_LNM, means_stan_LNM))+geom_point()+geom_abline(coef = c(0,1))
ggsave("../../results/assessing_models/simA_LNM_scatter.pdf")

## In the DM there are still
